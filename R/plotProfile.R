pv.style          <- "point"
pv.nOfWindows     <- 100
pv.bin_size       <- 20
pv.distanceAround <- 1500
pv.distanceUp     <- 1000
pv.distanceDown   <- 1000


pv.plotProfile <- function(pv, mask, sites, maxSites=1000, labels,
                           scores="Score", absScores=TRUE, annotate=TRUE, 
                           normalize=TRUE,merge=DBA_REPLICATE,
                           doPlot=TRUE, returnVal="profileplyr",
                           ...) {
  
  if (!suppressWarnings(requireNamespace("profileplyr",quietly=TRUE))) {
    stop("Package profileplyr not installed",call.=FALSE)
  }
  
  if(missing(sites)) {
    sites <- NULL
  }
  
  if(!is(pv,"profileplyr")) { # Need to compute profiles
    
    if(missing(mask)) {
      mask <- NULL
    }
    
    if(missing(labels)) {
      labels <- NULL
    }
    
    sitelabels <- NULL
    if(is(labels,"list")) {
      if(!is.null(labels$sites)) {
        sitelabels <- labels$sites
      }
      if(!is.null(labels$samples)) {
        labels <- labels$samples
      } else {
        labels <- NULL
      }
    }
    
    groups <- FALSE
    condition1 <- condition2 <- NULL
    
    # Check for default contrast 
    if(is.null(sites)) { 
      if(!is.null(pv$contrasts[[1]])) {
        if(pv$config$AnalysisMethod == DBA_EDGER) {
          if(!is.null(pv$contrasts[[1]]$edgeR)) {
            sites <- 1
          }
        } else {
          if(!is.null(pv$contrasts[[1]]$DESeq2)) {
            sites <- 1
          }
        }
      }
    }
    
    # If sites is unspecified
    if(is.null(sites)) { 
      # Get all sites from binding matrix and scramble
      sites <- pv.peakMatrix_toGR(pv,pv$binding)
      sites <- sites[sample(1:length(sites)),]
      scores <- NULL
    }  else if(is(sites,"numeric")) {
      # Generate Report-based
      con <- sites
      sites <- dba.report(pv, contrast=con,
                          bDB=TRUE,bGain=TRUE,bLoss=TRUE,bAll=FALSE)
      #sites$peaks[[2]]$Score <- abs(sites$peaks[[2]]$Score)
      
      if(is.null(mask)) {
        if(!is.null(pv$contrasts[[con]]$group1)) {
          condition1 <- which(pv$contrasts[[con]]$group1)
          condition2 <- which(pv$contrasts[[con]]$group2)
          name1 <- pv$contrasts[[con]]$name1
          name2 <- pv$contrasts[[con]]$name2
          mask <- list(name1=condition1,name2=condition2)
          names(mask) <- c(name1,name2)
        }
      }
    }
    
    ## Which samples
    samplegroups <- NULL
    if(is.null(mask)) {
      mask <- 1:length(pv$peaks)
    } else if(is(mask,"logical")) {
      mask <- which(mask)
    } else if(is(mask,"list")) {
      samplegroups <- lapply(mask,function(x){if(is(x,"logical")){which(x)}else{x}})
      mask <- unlist(samplegroups)
    }
    
    samples <- pv$class["bamRead",mask]
    
    # Normalization factors
    normfacs <- rep(1, length(mask))
    if(is(normalize,"logical")) {
      if(normalize != FALSE) {
        normfacs <- dba.normalize(pv,bRetrieve=TRUE)$norm.factors[mask]
      }
    } else {
      if(length(normalize) == length(mask)) {
        normfacs <- normalize
      } else {
        stop("normalize must include one factor for each sample.")
      }
    }
    
    ## Sample labels
    if(is.null(labels)) {
      sampnames <- pv$class[DBA_ID,mask]
    } else if(is(labels[1],"character")) {
      sampnames  <- pv$class[DBA_ID,mask]
    } else if(length(labels) == 1) {
      sampnames <- pv$class[labels,mask]
    } else {
      sampnames <-  apply(pv$class[labels,mask], 2, paste,collapse="_")
    }
    names(samples) <- sampnames
    
    ## Get sites
    if (is(sites,"vector")) { # specified sites
      # Get specified sites from binding matrix
      sites  <- pv.peakMatrix_toGR(pv,pv$binding[sites,])
      scores <- NULL
    } else if(is(sites, "GRanges")) { # Supplied sites
      # User supplied sites
      # just pass them through
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
        #names(mcols(peaks[[i]])) <- "score"
      }
      peaks <- GRangesList(peaks)
      
      # Groupnames from report DBA $ class
      names(peaks) <- groupnames <- pv.groupNames(sites)
      sites <- peaks
    }
    
    if(is(sites,"GRanges")) {
      sites <- GRangesList(sites)
    }
    
    if(is(sites,"GRangesList")){ 
      if(!is.null(sitelabels)) {
        if(length(sitelabels) != length(sites)) {
          warning("Should be same number of site group labels as site groups",
                  .call=FALSE)
        } else {
          names(sites) <- groupnames <- sitelabels
        }
      }
    }
    
    # Limit sites
    for(i in 1:length(sites)) {
      sites[[i]] <- sites[[i]][1:min(maxSites,length(sites[[i]])),]
    }
    
    # save sites in bedfiles
    bedfiles <- NULL
    for(i in 1:length(sites)) {
      bedfile <- tempfile(as.character(Sys.getpid()))
      rtracklayer::export.bed(sites[[i]],con=bedfile)
      bedfiles <- c(bedfiles,bedfile)
    }
    
    # Generate mergelist from attributes if required
    if(!is.null(merge)) {
      merge     <- pv.mergelist(pv, mask, merge, labels)
      sampnames <- names(merge)
    } 
    
    ## Generate profiles
    profiles <- pv.profiles(pv, samples, bedfiles,
                            mergelist = merge, normfacs = normfacs,
                            ...) 
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
      rowData(profiles)$"Binding Sites" <- factor(grlabels, 
                                                  levels=unique(grlabels),
                                                  labels=groupnames)
      profiles@params$rowGroupsInUse <- "Binding Sites"
    }
    
    # Add samplegroup labels
    if(!is.null(samplegroups)) {
      # convert to merged samples 
      if(!is.null(merge)) {
        if(!is.null(merge)) {
          conditions <- NULL
          sgroups <- samplegroups
          sampnum <- 1
          for(sgroup in 1:length(sgroups)) {
            sgroups[[sgroup]] <- sampnum:(sampnum+length(sgroups[[sgroup]])-1)
            sampnum <- sampnum + length(sgroups[[sgroup]])
          }
          for(tomerge in merge) {
            gnum <- which(unlist(lapply(sgroups,function(x){
              if(tomerge[1] %in% x){TRUE}else{FALSE}})))
            conditions <- c(conditions,names(sgroups)[gnum])
          }
          gnames <- NULL
          for(gname in names(samplegroups)) {
            gnames <- c(gnames,rep(gname, sum(conditions %in% gname)))
          }
          conditions <- gnames
        }
      }
      else {
        conditions <- NULL
        for(sgroup in 1:length(samplegroups)) {
          sname <- names(samplegroups)[sgroup]
          conditions <- c(conditions,rep(sname,length(samplegroups[[sgroup]])))
        }
      }
      # add samplegroup labels
      metadata(profiles) <- c(metadata(profiles),
                              list("Sample Group"=factor(conditions, 
                                                         levels=unique(conditions))))
    }
    
    # delete bedfiles
    unlink(bedfiles)
    
  } else {  # Profiles already defined -- just plotting
    profiles <- pv
    sampnames <- names(assays(profiles))
    gnames <- unique(rowData(profiles)$"Binding Sites")
    if(is.null(gnames)) {
      gnames <- unique(rowData(profiles)$sgGroup)
      groups <- FALSE
    } else {
      groups <- TRUE
    }
    groupnames <- as.character(gnames)
    
    # Normalization factors
    normfacs <- rep(1, length(sampnames))
    if(!is(normalize,"logical")) {
      if(length(normalize) == length(sampnames)) {
        normfacs <- normalize
        if(length(normfacs) != length(assays(profiles))) {
          warning("Wrong number of normalization factors, skipping.",call.=FALSE)
        } else {
          for(i in 1:length(normfacs)) {
            assay(profiles,i) <-
              assay(profiles,i) / normfacs[i]
          }
        }
      } else {
        stop("normalize must include one factor for each sample.")
      }
    }
  }
  
  
  # Scores
  profiles <- pv.siteScore(profiles, sites, scores, absScores)
  
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
                                   samples_names=sampnames,
                                   group_names=groupnames,
                                   annotate=doAnnotation,
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

pv.profiles <- function(pv, samples, sites, mergelist=NULL, normfacs=NULL,
                        style=pv.style, nOfWindows=pv.nOfWindows, 
                        bin_size=pv.bin_size, 
                        distanceAround=pv.distanceAround, 
                        distanceUp=pv.distanceUp, distanceDown=pv.distanceDown,
                        ...) {
  
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
  if(!is.null(pv$config$pp.distanceUp)) {
    distanceUp <- pv$config$pp.distanceUp
  }
  if(!is.null(pv$config$pp.distanceDown)) {
    distanceDown <- pv$config$pp.distanceDown
  }
  
  args <- pv.sepProfilingArgs(list(...), remove=FALSE)
  if(length(args) > 0) {
    if(!is.null(args$nOfWindows)) {
      nOfWindows <- args$nOfWindows
    }
    if(!is.null(args$bin_size)) {
      bin_size <- args$bin_size
    }
    if(!is.null(args$distanceAround)) {
      distanceAround <- args$distanceAround
    }
    if(!is.null(args$distanceUp)) {
      distanceDown <- args$distanceDown
    }
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
    param <- BiocParallel::MulticoreParam(workers=cores)
  } else {
    param <-  BiocParallel::SerialParam()
  }
  
  # Call profileplyr
  message("Generating profiles...")
  profiles <- NULL
  if(is.null(profiles)) { # NULL check for debugging
    res <- tryCatch(
      suppressMessages(
        profiles <- bplapply(samples,
                             profileplyr::BamBigwig_to_chipProfile,
                             sites, 
                             format="bam", paired=paired,
                             style=style,nOfWindows=nOfWindows,
                             bin_size=bin_size, 
                             distanceAround=distanceAround,
                             distanceUp=distanceUp, distanceDown=distanceDown,
                             BPPARAM=param)),
      error=function(x){stop("profileplyr error: ",x)}
    )
  }
  
  # Normalize
  if(!is.null(normfacs)) {
    if(length(normfacs) != length(profiles)) {
      warning("Wrong number of normalization factors, skipping.",call.=FALSE)
    } else {
      for(i in 1:length(normfacs)) {
        assay(profiles[[i]],1) <-
          assay(profiles[[i]],1) / normfacs[i]
      }
    }
  }
  
  if(!is.null(mergelist)) {
    profiles <- pv.mergeProfiles(profiles, mergelist)
  }
  
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
  
  return(proplyrObject)
}

pv.siteScore <- function(profiles, sites, scores, absScores) {
  
  if(!is.null(sites)) { # Add metadata
    psites <- mcols(profiles)$name
    if(is(sites,"GRangesList")) {
      names(sites) <- NULL
    }
    sites  <- unlist(sites)
    ssites <- names(sites)
    
    matches <- match(psites, ssites)
    if(sum(is.na(matches)) > 0) {
      mcols(profiles) <- cbind(mcols(profiles), mcols(sites))      
    } else {
      mcols(profiles) <- cbind(mcols(profiles),
                               mcols(sites[match(psites, ssites),]))
    }
  }
  
  profiles@params$mcolToOrderBy <- NULL
  if(!is.null(scores)) {
    if(!scores %in% names(mcols(profiles))) {
      if(!is.null(profiles@params$mcolToOrderBy)) {
        message(scores," not a valid score column, using mean signal.")
      }
    } else {
      profiles@params$mcolToOrderBy <- "score"
      whichscore <- match(scores, names(mcols(profiles)))
      mcols(profiles)$score <- mcols(profiles)[whichscore][,1]   
      if(absScores) {
        mcols(profiles)$score <- abs(mcols(profiles)$score)
      } 
    }
  } 
  
  profiles <- profileplyr::orderBy(profiles,profiles@params$mcolToOrderBy)
  
  return(profiles)
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
  
  meta <- metadata(profiles)
  
  message("Annotating...")
  suppressMessages(
    profiles <- profileplyr::annotateRanges(profiles, TxDb=annotate, 
                                            verbose=FALSE)
  )
  profiles <- pv.setAnno(profiles)
  
  metadata(profiles) <- c(meta,list(annotate=annotate))
  
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

pv.conditionColors <- list(c("white",crukMagenta), 
                           c("white",crukBlue), 
                           c("white",crukGrey), 
                           c("white",crukCyan))

pv.profileHeatmap <- function(profiles, samples_names, group_names,
                              annotate=FALSE,
                              ret_ht_list=FALSE, ...){
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
  
  arguments <- list(profiles)
  
  args <- list(...)
  args <- pv.sepProfilingArgs(args, remove=TRUE)
  addarg <- NULL
  
  if(!"matrices_color" %in% names(args)) {
    conditions <- metadata(profiles)$"Sample Group"
    if(!is.null(conditions)) {
      addarg <- pv.addArg(addarg, "matrices_color",
                          pv.conditionColors[conditions])
    }
  }
  
  addarg <- pv.addArg(addarg, "matrices_pos_line",FALSE, args)
  addarg <- pv.addArg(addarg, "decreasing",TRUE, args)
  addarg <- pv.addArg(addarg, "sample_names",samples_names, args)
  addarg <- pv.addArg(addarg, "include_group_annotation",
                      include_group_annotation, args)
  addarg <- pv.addArg(addarg, "group_anno_color",
                      group_anno_color, args)
  addarg <- pv.addArg(addarg, "group_anno_column_names_gp",
                      grid::gpar(col="white"), args)
  addarg <- pv.addArg(addarg, "extra_annotation_columns",
                      extra_annotation_columns, args)
  addarg <- pv.addArg(addarg, "extra_anno_color", extra_anno_color, args)
  addarg <- pv.addArg(addarg, "return_ht_list",ret_ht_list, args)
  addarg <- pv.addArg(addarg, "raster_device","png", args)
  addarg <- pv.addArg(addarg, "raster_quality","10", args)
  
  
  if(!is.null(addarg)){
    arguments <- c(arguments, addarg)
  }  
  
  hm <- do.call(profileplyr::generateEnrichedHeatmap, c(arguments,args))
  
  return(hm)
}


pv.addArg <- function(addarg, param, val, args=NULL) {
  
  if(!is.null(args)) {
    if(param %in% names(args)) {
      return(addarg)
    }
  }
  
  if(is.null(addarg)) {
    addarg <- list(x=val)
    names(addarg) <- param
  } else {
    addarg <- pv.listadd(addarg, val)
    names(addarg)[length(addarg)] <- param
  }
  return(addarg)
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

pv.mergelist <- function(pv, samples, merge, labels=NULL){
  
  if(!is(merge,"list")) { # Make list from excluded attributes
    attributes <- c(DBA_TISSUE,DBA_FACTOR,DBA_CONDITION,DBA_TREATMENT,DBA_REPLICATE)
    atts <- pv$class[attributes[! attributes %in% merge],samples]
    allmerge <- apply(atts,2,paste,collapse="")
    tomerge  <- unique(allmerge)
    if(is.null(labels)) {
      labels <- -merge
    }
    
    merge <- NULL
    for(val in tomerge) {
      merge <- pv.listadd(merge,which(allmerge %in% val))
    }
  }
  
  # Attach labels to mergelist
  
  names(merge) <- pv.makeMergeLabel(pv, merge, samples,labels)
  
  # check uniqueness of labels
  names(merge) <- make.unique(as.character(names(merge)))
  
  # return list of merge groups with label names
  return(merge)
}

pv.makeMergeLabel <- function(pv, merge, samples, labels=NULL) {
  
  if(is(labels[1],"character")) {
    if(length(labels) != length(merge)) {
      warning("Should be same number of sample labels as samples.",call.=FALSE)
    } else {
      return(labels)
    }
  }
  
  if(is.null(names(merge))) {
    
    atts <- c(DBA_TISSUE,DBA_FACTOR,DBA_CONDITION,DBA_TREATMENT,DBA_REPLICATE)
    
    #samples <- samples[unlist(merge)]
    
    if(is.null(labels)) {
      labels <- atts
    }
    
    if(labels[1] < 0) {
      labels <-  atts[which(!atts %in% abs(labels))]
    }
    
    # Strip common attributes
    if(length(labels) > 1) {
      uatts <- apply(pv$class[labels,samples],1,unique)
      for(i in length(labels):1) {
        if(length(uatts[[i]])==1) {
          labels <- labels[-i]
        }
      }
    }
    
    res <- lapply(merge, pv.doMakeMergeLabel, pv$class, samples, labels)
    
    return(res)
  }
}

pv.doMakeMergeLabel <- function(tomerge, meta, samples, labels) {
  res <- NULL
  for(label in labels) {
    res <- c(res, unique(meta[label,samples[tomerge]]))
  }
  if(length(res) > 1) {
    res <- paste(res, collapse="_")
  }
  return(res)
}


pv.mergeProfiles <- function(profiles, mergelist, bMean=TRUE) {
  
  totalsamps <- length(profiles)
  mergelist <- mergelist[order(unlist(lapply(mergelist,min)))]
  mergesamps <- sort(unlist(mergelist))
  if(length(mergesamps) != totalsamps) {
    toadd <- which(!(1:totalsamps %in% mergesamps))
  } else toadd <- NULL
  
  resultlist <- NULL
  
  for(tomerge in mergelist) {
    
    if(length(tomerge) > 1) {
      
      adding <- TRUE
      while(adding){
        if(length(toadd) > 0) {
          if(toadd[1] < tomerge[1]) {
            resultlist <- pv.listadd(resultlist, profiles[[toadd[1]]])
            toadd <- toadd[-1]
          } else {
            adding <- FALSE
          }
        } else adding <- FALSE
      }
      
      profile <- profiles[[tomerge[1]]]
      
      merged <- assay(profile)
      for(i in 2:length(tomerge)) {
        merged <- merged + assay(profiles[[tomerge[i]]])
      }
      
      if(bMean) {
        merged <- merged / length(tomerge)
      }
      
      assay(profile) <- merged
      
    } else {
      profile <- profiles[[tomerge[1]]]
    }
    
    resultlist <- pv.listadd(resultlist, profile)
  }
  
  
  while(length(toadd) > 0) {
    resultlist <- pv.listadd(resultlist, profiles[[toadd[1]]])
    toadd <- toadd[-1]
  }
  
  return(resultlist)
  
}

pv.ProfilingArgs <- c("style","nOfWindows","bin_size",
                      "distanceAround","distanceUp","distanceDown")

pv.sepProfilingArgs <- function(arglist, remove=FALSE) {
  profiling <- which(names(arglist) %in% pv.ProfilingArgs)
  if(length(profiling) > 0) {
    proargs <- arglist[profiling]
    plotargs <- arglist[!profiling]
  } else {
    proargs <- NULL
    plotargs <- arglist
  }

  if(remove) {
    return(plotargs)
  } else {
    return(proargs)
  }
}

