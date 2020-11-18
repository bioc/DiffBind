#setwd("20151110_BirgitN_RG_FusionChIP/delivered/")
#setwd("/Users/stark01/SynologyDrive/Work/RELA")

# library(DiffBind)
# library("profileplyr")
# require("Rsamtools")
# library(rtracklayer)
# library(BiocParallel)
# library(grid)



pv.plotProfile <- function(pv, mask, sites, maxSites=1000, 
                           scores, annotate, labels=DBA_ID,
                           doPlot=TRUE, returnVal="profileplyr",
                           ...) {
  
  if (!requireNamespace("profileplyr",quietly=TRUE)) {
    stop("Package profileplyr not installed",call.=FALSE)
  }
  
  ## Which samples
  if(missing(mask)) {
    mask <- 1:length(pv$peaks)
  } else if(is(mask,"logical")) {
    mask <- which(mask)
  }
  samples <- pv$class["bamRead",mask]
  
  ## Sample labels
  sampnames <- pv$class[labels,mask]
  names(samples) <- sampnames
  
  ## Which sites
  if(missing(sites)) {
    sites <- NULL
  }
  
  groups <- FALSE
  # sites should end up a GRangesList
  if(is.null(sites)) { # all sites
    # Get all sites from binding matrix
  }  else if (is(sites,"vector")) { # specified sites
    # Get specified sites from binding matrix
  } else if(is(sites, "GRanges")) { # Supplied sites
    # Use supplied sites
  } else if(is(sites,"GRangesList")) { # Supplied groups of sites
    groups <- TRUE
    groupnames <- names(sites)
  } else if(is(sites,"DBA")) { # Report DBA object
    
    # check if report-based
    if(is.null(sites$resultObject)) {
      stop("sites is not a report-based DBA object.")
    } else if(sites$resultObject==FALSE) {
      stop("sites is not a report-based DBA object.")
    }
    
    # Get site groups from report DBA
    groups <- TRUE
    
    # Get each peakset, form GRangesList
    peaks <- sites$peaks
    for(i in 1:length(peaks)) {
      peaks[[i]] <- GRanges(peaks[[i]])
      names(mcols(peaks[[i]])) <- "score"
    }
    peaks <- GRangesList(peaks)
    
    # Groupnames from report DBA $ class
    names(peaks) <- groupnames <- pv.groupNames(sites)
    sites <- peaks
  }
  
  if(is(sites,"GRanges")) {
    sites <- GRangesList(sites)
  }
  
  # Limit sites
  for(i in 1:length(sites)) {
    sites[[i]] <- sites[[i]][1:min(maxSites,length(sites[[i]])),]
  }
  
  # Scores
  
  # save sites in bedfiles
  bedfiles <- NULL
  for(i in 1:length(sites)) {
    bedfile <-tempfile(as.character(Sys.getpid()))
    rtracklayer::export.bed(sites[[i]],con=bedfile)
    bedfiles <- c(bedfiles,bedfile)
  }
  
  ## Generate profiles
  profiles <- pv.profiles(pv, samples, bedfiles) #, 
  #style, nOfWindows, bin_size, distanceAround)
  
  # Add group labels
  if(!groups) {
    groupnames <- names(sites)[1]
    if(is.null(groupnames)) {
      rowData(profiles)$sgGroup <- "Sites"
    } else {
      rowData(profiles)$sgGroup <- groupnames
    }
  } else {
    fns <- strsplit(bedfiles,"/")
    grlabels <- as.character(rowData(profiles)$sgGroup)
    for(i in 1:length(fns)) {
      fn <- fns[[i]]
      grp <- fn[length(fn)]
      grlabels[grlabels==grp] <- groupnames[i]
    }
    rowData(profiles)$"Binding Sites" <- factor(grlabels, ordered=TRUE)
    profiles@params$rowGroupsInUse <- "Binding Sites"
  }

  # delete bedfiles
  unlink(bedfiles)
  
  
  # Annotation
  profiles <- pv.annotate(pv, profiles, annotate)
  
  # Set up Groups
  if(groups) {
    profiles <- pv.groupProfiles(profiles)
  }
  
  # Generate plot
  if(doPlot) {
    if(returnVal=="HeatmapList") {
      return_ht_list <- TRUE
    } else {
      return_ht_list <- FALSE
    }
    profilehm <- pv.profileHeatmap(profiles, 
                                   sample_names=sampnames,
                                   group_names=groupnames,
                                   return_ht_list,
                                   ...)
  }
  
  #return
  if(returnVal == "profileplyr") {
    return(profiles)
  }
  
  if(returnVal == "EnrichedHeatmap") {
    return(profilehm)
  }
  
  if(returnVal == "HeatmapLKist") {
    return(profilehm)
  }
  
}

pv.profiles <- function(pv, samples, sites, 
                        style="point", nOfWindows=100, bin_size=25, 
                        distanceAround=1000) {
  
  # Check bam files: exist, PE/SE, BAI
  
  if(is.null(pv$config$singleEnd)) {
    bfile <- pv.BamFile(samples[1], bIndex=TRUE)
    pv$config$singleEnd	<- !suppressMessages(
      Rsamtools::testPairedEndBam(bfile))
  }
  paired <- !pv$config$singleEnd
  
  # Get MulticoreParam
  
  if(pv$config$RunParallel) {
    if(is.null(pv$config$cores)) {
      cores <- BiocParallel::multicoreWorkers()
    } else {
      cores <- pv$config$cores
    }
    param <- BiocParallel::MulticoreParam(workers <- cores)
  } else {
    param <-  BiocParallel::SerialParam()
  }
  
  # Call profileplyr
  message("Generating profiles...")
  profiles <- NULL
  res <- tryCatch(
    suppressMessages(
      profiles <- bplapply(samples,
                           profileplyr::BamBigwig_to_chipProfile,
                           sites, 
                           format="bam", paired=paired,
                           style=style,nOfWindows=nOfWindows,
                           bin_size=bin_size, 
                           distanceAround=distanceAround,
                           BPPARAM=param)),
    error=function(x){stop("profileplyr error: ",x)}
  )
  
  if(length(profiles) > 1) {
    profileobjs <- profileplyr::as_profileplyr(profiles[[1]])
    for(i in 2:length(profiles) ) {
      profileobjs <- c(profileobjs, profileplyr::as_profileplyr(profiles[[i]]))
    }
    proplyrObject <- profileobjs
  } else {
    proplyrObject <- profiles
  }
  
  if(!is(proplyrObject,"profileplyr")) {
    proplyrObject <- profileplyr::as_profileplyr(proplyrObject)
  }
  
  proplyrObject <- profileplyr::orderBy(proplyrObject,"score")
  
  return(proplyrObject)
}


pv.annotate <- function(pv, profiles, annotate) {
  return(profiles)
}

pv.groupColors <- c(crukBlue,crukCyan,crukMagenta, crukGrey)

pv.profileHeatmap <- function(profiles, sample_names, group_names,
                              return_ht_list=FALSE, ...){
  if(is.null(group_names)) {
    include_group_annotation <- FALSE
    group_anno_color <- NULL
  } else {
    include_group_annotation <- TRUE
    group_anno_color <- pv.groupColors[1:length(group_names)]
  }
  
  hm <- profileplyr::generateEnrichedHeatmap(profiles,
                                             matrices_pos_line=FALSE,
                                             decreasing=TRUE,
                                             sample_names=sample_names,
                                             include_group_annotation =
                                               include_group_annotation,
                                             group_anno_color =  
                                               group_anno_color,
                                             group_anno_column_names_gp =
                                               grid::gpar(col="white"),
                                             return_ht_list=return_ht_list,
                                             ...)
  return(hm)
}

pv.groupNames <- function(pv) {
  if(length(unique(pv$class[DBA_ID,])) == ncol(pv$class)) {
    return(pv$class[DBA_ID,])
  }
  
  labels <- NULL
  for(i in c(DBA_ID,DBA_TISSUE, DBA_FACTOR, 
             DBA_CONDITION, DBA_TREATMENT)) {
    if(length(unique(pv$class[i,])) > 1) {
      labels <- paste(labels,pv$class[i,],sep="_")
    }
  }
  
  labels <- sub('.', '',labels)
  
  return(labels)
}

pv.groupProfiles <- function(profiles) {
  return(profiles)
}

