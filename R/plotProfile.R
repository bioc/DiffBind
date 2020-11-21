#setwd("20151110_BirgitN_RG_FusionChIP/delivered/")
#setwd("/Users/stark01/SynologyDrive/Work/RELA")

# library(DiffBind)
# library("profileplyr")
# require("Rsamtools")
# library(rtracklayer)
# library(BiocParallel)
# library(grid)


pv.style          <- "point"
pv.nOfWindows     <- 100
pv.bin_size       <- 20
pv.distanceAround <- 1500

pv.plotProfile <- function(pv, mask, sites, maxSites=1000, 
                           scores, annotate=TRUE, labels=DBA_ID,
                           doPlot=TRUE, returnVal="profileplyr",
                           ...) {
  
  if (!requireNamespace("profileplyr",quietly=TRUE)) {
    stop("Package profileplyr not installed",call.=FALSE)
  }
  
  if(is(pv,"profileplyr")) {
    profiles <- pv
    sampnames <- names(assays(profiles))
    gnames <- unique(rowData(profiles)$"Binding Sites")
    if(is.null(gnames)) {
      gnames <- rowData(profiles)$sgGroup
      groups <- FALSE
    } else {
      groups <- TRUE
    }
    groupnames <- as.character(gnames)
  } else {
    
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
      # User supplied sites
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
    rownames(profileplyr::sampleData(profiles)) <-
      names(assays(profiles)) <- sampnames
    
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
      rowData(profiles)$"Binding Sites" <- factor(grlabels, levels=groupnames)
      profiles@params$rowGroupsInUse <- "Binding Sites"
    }
    
    # delete bedfiles
    unlink(bedfiles)
  }
  
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
    doAnnotation <- annotate != FALSE
    profilehm <- pv.profileHeatmap(profiles, 
                                   sample_names=sampnames,
                                   group_names=groupnames,
                                   annotate=doAnnotation,
                                   return_ht_list,
                                   ...)
  }
  
  #return
  gc()
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
                        style=pv.style, nOfWindows=pv.nOfWindows, 
                        bin_size=pv.bin_size, 
                        distanceAround=pv.distanceAround) {
  
  # check config
  if(!is.null(pv$config$pp.style)) {
    style <- pv$config$pp.style
  }
  if(!is.null(pv$config$pp.nOfWindows)) {
    nOfWindows <- pv$config$pp.nOfWindows
  }
  if(!is.null(pv$config$pp.bin_size)) {
    bin_size <- pv$config$pp.bin_size
  }
  if(!is.null(pv$config$pp.distanceAround)) {
    distanceAround <- pv$config$pp.distanceAround
  }
  
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
  
  if(!is.null(rowData(profiles)$Features)) {
   return(profiles)
  }
  
  if(is(annotate,"logical")){
    if(!annotate) {
      return(profiles)
    } else {
      annotate <- pv.genomes(pv$class["bamRead",1],pv$chrmap)
      if(annotate=="BSgenome.Hsapiens.UCSC.hg19") {
        annotate <- "hg19"
      } else if(annotate=="BSgenome.Hsapiens.UCSC.hg38") {
        annotate <- "hg38"
      } else if(annotate=="BSgenome.Mmusculus.UCSC.mm9") {
        annotate <- "mm9"
      } else if(annotate=="BSgenome.Mmusculus.UCSC.mm10") {
        annotate <- "mm10"
      } else {
        message("Must specify transcriptome")
        return(profiles)
      }
      message("Annotate using: ",annotate)
    }
  }
  
  message("Annotating...")
  suppressMessages(
    profiles <- profileplyr::annotateRanges(profiles, TxDb=annotate, 
                                            verbose=FALSE)
  )
  profiles <- pv.setAnno(profiles)
  
  return(profiles)
  
}

pv.setAnno <- function(dataset) {
  Features <-  as.character(rowData(dataset)$annotation_short)
  Features[Features=="3p UTR"] <- "Gene Body"
  Features[Features=="5p UTR"] <- "Gene Body"
  Features[Features=="Downstream"] <- "Gene Body"
  Features[Features=="Exon"] <- "Gene Body"
  Features[Features=="Intron"] <- "Gene Body"
  Features[Features=="Distal Intergenic"] <- "Intergenic"
  Features <- factor(Features,levels=c("Promoter","Gene Body","Intergenic"))
  rowData(dataset)$Features <- Features
  
  Distances <- rowData(dataset)$distanceToTSS
  Distances[Distances < -100000] <- -1000000
  Distances[Distances > 100000]  <- 1000000
  rowData(dataset)$"TSS Distance" <- Distances
  
  return(dataset)
}

pv.groupColors <- c(crukBlue,crukCyan,crukMagenta, crukGrey, 
                    "green", "forestgreen")

pv.profileHeatmap <- function(profiles, sample_names, group_names,
                              annotate=FALSE,
                              return_ht_list=FALSE, ...){
  message("Plotting...")
  
  if(is.null(group_names)) {
    include_group_annotation <- FALSE
    group_anno_color <- NULL
  } else {
    include_group_annotation <- TRUE
    group_anno_color <- pv.groupColors[1:length(group_names)]
  }
  
  if(is.null(rowData(profiles)$Features)) {
    annotate <- FALSE
  }
  if(annotate) {
    extra_annotation_columns <- c("Features","TSS Distance" )
    extra_anno_color = list(c(crukBlue,crukCyan,crukMagenta),
                            colorRampPalette(c("red", grDevices::rgb(.99,.99,.99), 
                                               "green"))(n =9))
  } else {
    extra_annotation_columns <- NULL
    extra_anno_color <- NULL
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
                                             extra_annotation_columns = 
                                               extra_annotation_columns,
                                             extra_anno_color =
                                               extra_anno_color,
                                             return_ht_list=return_ht_list,
                                             raster_device="png",
                                             raster_quality = 10,
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

pv.mergeProfiles <- function(profiles, bMean=TRUE) {
  
  profile <- profiles[[1]]
  
  merged <- assay(profile)
  for(i in 2:length(profiles)) {
    merged <- merged + assay(profiles[[i]])
  }
  
  if(bMean) {
    merged <- merged / length(profiles)
  }
  
  assay(profile) <- merged
  
  return(profile)
  
}

