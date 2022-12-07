#####################################
## pv_model.R -- make pv object    ##
## 20 October 2009                 ##
## 3 February 2011 -- packaged     ##
## Rory Stark                      ##
## Cancer Research UK              ##
#####################################
PV_DEBUG <- FALSE

## pv.peakset -- add a peakset to the model
pv.peakset <- function(pv = NULL,peaks, sampID, tissue, factor,condition, treatment, replicate,
                       control, peak.caller, peak.format, reads = 0, consensus =FALSE, 
                       readBam, controlBam, scoreCol = NULL, bLowerScoreBetter = NULL,
                       bRemoveM = TRUE, bRemoveRandom = TRUE,
                       minOverlap = 2,bFast = FALSE,bMakeMasks = TRUE,
                       skipLines = 1, filter = NULL, counts = NULL, spikein) {
  zeroVal <- -1
  bLog <- FALSE
  
  if (missing(peaks)) {
    peaks <- NULL
  }
  
  if(length(peaks) == 0) {
    peaks <- NULL
  }
  
  # if (!is.null(pv$peaks) && length(peaks) == 0) {
  #   peaks <- 1:length(pv$peaks)
  # }
  
  if (missing(counts))
    counts <- NULL
  if (!is.null(counts)) {
    res <- pv.peaksetCounts(
      pv = pv,peaks = peaks,counts = counts,
      sampID = sampID,tissue = tissue,factor = factor,
      condition = condition, treatment = treatment,replicate = replicate
    )
    return(res)
  }
  
  if (missing(peak.format))
    peak.format <- NULL
  if (missing(scoreCol))
    scoreCol <- NULL
  if (missing(bLowerScoreBetter))
    bLowerScoreBetter <- NULL
  if (missing(filter))
    filter <- NULL
  
  bConsensus <- FALSE
  if (is.numeric(consensus)) {
    ## Add a set of consensus peaksets
    bConsensus <- TRUE
    pv <-
      pv.consensusSets(pv,peaks = peaks,minOverlap = minOverlap,
                       attributes = consensus, tissue,factor,condition,treatment,
                       replicate,control,peak.caller, readBam, controlBam,
                       spikein)
    
  } else {
    ## add a specific consensus peakset
    if (is.vector(peaks) && length(peaks) > 1) {
      # consensus
      bConsensus <- TRUE
      pv <-
        pv.consensus(pv,peaks,minOverlap = minOverlap,bFast = bFast)
      if (!is.null(minOverlap)) {
        nset <- length(pv$peaks)
        if (!missing(sampID)) {
          pv$class[PV_ID,nset] <- sampID
          colnames(pv$class)[nset] <- sampID
        }
      }
      
      if (!missing(tissue))
        pv$class[PV_TISSUE,nset] <- tissue
      if (!missing(factor))
        pv$class[PV_FACTOR,nset] <- factor
      if (!missing(condition))
        pv$class[PV_CONDITION,nset] <- condition
      if (!missing(treatment))
        pv$class[PV_TREATMENT,nset] <- treatment
      if (!missing(replicate))
        pv$class[PV_REPLICATE,nset] <- replicate
      if (!missing(control))
        pv$class[PV_CONTROL,nset] <- control
      if (!missing(peak.caller))
        pv$class[PV_CALLER,nset] <- peak.caller
      if (!missing(readBam))
        pv$class[PV_BAMREADS,nset] <- readBam
      if (!missing(controlBam))
        pv$class[PV_BAMCONTROL,nset] <- controlBam
      if (!missing(spikein))
        pv$class[PV_SPIKEIN,nset] <- spikein
    }
  }
  if (bConsensus) {
    if (bMakeMasks) {
      pv$masks <- pv.mask(pv)
    }
    return(pv)
  }
  
  if (missing(tissue))
    tissue <- ''
  if (missing(factor))
    factor <- ''
  if (missing(condition))
    condition <- ''
  if (missing(treatment))
    treatment <- ''
  if (missing(replicate))
    replicate <- ''
  if (missing(control))
    control <- ''
  if (missing(peak.caller))
    peak.caller <- ''
  if (missing(readBam)) {
    readBam <- NA
    #warning("No bam file specified.",call.=FALSE)
  }
  if (length(readBam) == 0) {
    readBam <- NA
    #warning("No bam file specified.",call.=FALSE)
  }
  if (missing(controlBam))
    controlBam <- NA
  if (length(controlBam) == 0)
    controlBam <- NA
  if (missing(spikein))
    spikein <- NA
  if (length(spikein) == 0)
    spikein <- NA
  
  if (!is.null(peaks) && length(peaks) <= 1) {
    if (is.na(peaks) || length(peaks)==0) {
      peaks <- NULL
    } else {
      if (is.character(peaks)) {
        if (peaks == "" || peaks == " ")
          peaks <- NULL
      }
    }
  }
  if (is.null(peaks)) {
    peaks <- matrix(0,0,4)
  } else {
    if (is.character(peaks)) {
      # Read in peaks from a file
      if(file.info(peaks)$size > 0) {
        pcaller <- strtrim(peak.caller,6)
        if (is.null(peak.format)) {
          peak.format <- pcaller
        }
        if (is.null(scoreCol)) {
          scoreCol <- pv.defaultScoreCol(peak.format)
        }
        peaks <- pv.readPeaks(peaks,peak.format,skipLines)
      } else {
        peaks <- matrix(0,0,4)
        scoreSave <- scoreCol <- 0
      }
    } else {
      peaks <- pv.peaksort(peaks)
      if (is.null(scoreCol))
        scoreCol <- pv.defaultScoreCol(peak.format)
      if (is.null(bLowerScoreBetter))
        bLowerScoreBetter <- FALSE
    }
    
    scoreSave <- scoreCol
    if ( (nrow(peaks) > 0) & 
         ( (ncol(peaks) < scoreSave) | (ncol(peaks) == 3))){
      peaks <- cbind(peaks[,1:3],1)
      colnames(peaks)[ncol(peaks)] <- "score"
      scoreCol <- 0
    }
    
    if (is.null(bLowerScoreBetter)) {
      if (peak.caller == "report") {
        bLowerScoreBetter <- TRUE
      } else {
        bLowerScoreBetter <- FALSE
      }
    }
    
    if (is.null(filter) && peak.caller == "bayes") {
      filter <- 0.5
    }
    
    if (scoreCol > 0) {
      if (!missing(filter)) {
        if (!is.null(filter)) {
          if (bLowerScoreBetter) {
            tokeep <- peaks[,scoreCol] <= filter
          } else {
            tokeep <- peaks[,scoreCol] >= filter
          }
          peaks <- peaks[tokeep,]
        }
      }
      peaks[,scoreCol] <-
        pv.normalizeScores(peaks,scoreCol,zeroVal = zeroVal,bLog = bLog)
      if (bLowerScoreBetter) {
        peaks[,scoreCol] <- 1 - peaks[,scoreCol]
      }
      peaks <- peaks[,c(1:3,scoreCol)]
    }
    
    if (bRemoveM) {
      idx <- peaks[,1] != "chrM"
      peaks <- peaks[idx,]
      if (sum(!idx) > 0) {
        peaks[,1] <- as.factor(as.character(peaks[,1]))
      }
    }
    
    if (bRemoveRandom) {
      for (i in c(1:22,"X","Y")) {
        ch <- sprintf("chr%s_random",i)
        idx <- peaks[,1] != ch
        peaks <- peaks[idx,]
        if (sum(!idx) > 0) {
          peaks[,1] <- as.factor(as.character(peaks[,1]))
        }
      }
    }
    
    newchrs <- as.character(peaks[,1])
    pv$chrmap  <- sort(unique(c(pv$chrmap,newchrs)))
    #         peaks[,1] <- factor(peaks[,1],pv$chrmap)
    peaks[,1] <- as.character(peaks[,1])
  }
  
  colnames(peaks)[1:4] <- c("Chr","Start","End","Score")
  
  pv$peaks <- pv.listadd(pv$peaks,peaks)
  
  if (missing(sampID)) {
    if (is.null(pv)) {
      sampID <- 1
    } else if (is.null(pv$peaks)) {
      sampID <- 1
    } else {
      sampID <- length(pv$peaks)
    }
  }
  clascol <-
    cbind(
      NULL,c(
        sampID,tissue,factor,condition,consensus,peak.caller,control,
        reads,replicate,readBam,controlBam,treatment,spikein
      )
    )
  colnames(clascol) <- sampID
  pv$class <- cbind(pv$class,clascol)
  rownames(pv$class) <-
    c(
      "ID","Tissue","Factor","Condition", "Consensus",
      "Peak caller","Control","Reads","Replicate","bamRead",
      "bamControl","Treatment", "Spikein"
    )
  pv$merged  <- NULL
  pv$binding <- NULL
  if (bMakeMasks) {
    pv$masks <- pv.mask(pv)
  }
  return(pv)
}

## pv.vectors -- build the binding expression vectors and do clustering/PCA
pv.vectors <- function(pv,mask,minOverlap = 2,attributes,bAllSame = FALSE,
                       merge = TRUE) {
  if (missing(attributes)) {
    if (is.null(pv$attributes)) {
      attributes <- PV_ID
    } else {
      attributes <- pv$attributes
    }
  }
  
  if(missing(mask) && merge) {
    mask <- 1:length(pv$peaks)
  }
  
  if(is.null(pv$config$mergeOverlap)) {
    maxGap <- as.integer(-1)
  } else {
    maxGap <- as.integer(-pv$config$mergeOverlap)
  }
  
  called <- SN <- allcalled <- NULL
  if (!missing(mask)) {
    if (is.logical(mask)) {
      mask <- which(mask)
    }
    peaks <- NULL
    for (i in mask) {
      #if(nrow(pv$peaks[[i]]) > 0) {
      peaks <- pv.listadd(peaks,pv$peaks[[i]])
      #}
    }
    class      <- pv$class[,mask]
    chrmap     <- pv$chrmap
    config     <- pv$config
    samples    <- pv$samples
    if(!is.null(pv$called)) {
      if(length(mask) <= ncol(pv$called)) {
        called <- pv$called[,mask]
        if(!is.null(pv$allcalled)) {
          allcalled <- pv$allcalled[,mask]
        } 
      }
    }
    if(!is.null(pv$SN)) {
      if(length(mask)==length(SN)) {
        SN <- pv$SN[mask]
      }
    }
    score      <- pv$score
    summits    <- pv$summits
    minCount   <- pv$minCount
    #contrasts  <- pv$contrasts
    blacklist  <- pv$blacklist
    greylist   <- pv$greylist
    peaks.blacklisted <- pv$peaks.blacklisted
    resultObject <- pv$resultObject
    #annotation <- pv$annotation
    
    pv <- NULL
    pv$peaks      <- peaks
    pv$class      <- class
    pv$chrmap     <- chrmap
    pv$config     <- config
    pv$samples    <- samples
    pv$called     <- called
    pv$allcalled  <- allcalled
    #pv$contrasts  <- contrasts
    pv$score      <- score
    pv$SN         <- SN
    pv$summits    <- summits
    pv$minCount   <- pv$minCount
    pv$blacklist  <- blacklist
    pv$greylist   <- greylist
    pv$peaks.blacklisted <- peaks.blacklisted
    pv$resultObject <- resultObject
    #pv$annotation <- annotation
  }
  
  if (is.vector(pv$class)) {
    pv$class <- matrix(pv$class,length(pv$class),1)
    colnames(pv$class) <- pv$class[1,]
  }
  
  peaks <- pv$peaks
  
  numvecs <- length(peaks)
  ncols <- numvecs + 3
  
  if (minOverlap > 0 && minOverlap < 1) {
    minOverlap <- ceiling(numvecs * minOverlap)
  }
  
  npeaks <- 0
  defval <- -1
  
  if (!bAllSame) {
    if (sum(vapply(peaks,nrow,1)) > 0) {
      if (merge) {
        allp <-lapply(peaks,
                      function(x) {
                        y <- x[,1:3]
                        colnames(y) <- c("chr","start","end")
                        y})
        allpeaks <- NULL
        for(el in allp) { # check for empty peaksets
          if(nrow(el)>0) { 
            allpeaks <- pv.listadd(allpeaks,el)
          }
        }
        allpeaks <- bind_rows(allpeaks)
      } else {
        allpeaks <- data.frame(pv$merged)
        allpeaks[,1] <- pv$chrmap[allpeaks[,1]]
      }
    } else {
      allpeaks <- matrix(0,0,4)
    }
    allnames <- NULL
    if (nrow(allpeaks) > 0) {
      res  <- pv.merge(allpeaks,peaks,pv$class, maxgap=maxGap)
      pv$totalMerged <- do.nrow(res$merged)
      rownames(res$merged) <- 1:do.nrow(res$merged)
      allnames <- res$chrmap
      pv$called  <- res$included
      pv$allcalled <- NULL
      pv$merged <- res$merged[,1:3]
      if ((ncol(res$merged) > 4) && (minOverlap > 1)) {
        olaps <- pv.overlaps(pv,minOverlap)
        pv$binding <- res$merged[olaps,]
        pv$called  <- pv$called[olaps,]
        pv$allcalled <- res$included
      }  else {
        pv$binding <- res$merged
      }
    } else {
      pv$merged  <- matrix(0,0,3 + length(pv$peaks))
      pv$overlaps <- NULL
      colnames(pv$merged) <- colnames(allpeaks)
      pv$binding   <- pv$merge
      pv$totalMerged <- 0
      pv$called <- NULL
    }
  } else {
    ## ALL SAME
    result <- matrix(0,nrow(pv$peaks[[1]]),length(pv$peaks) + 3)
    if (is.character(pv$peaks[[1]][1,1])) {
      result[,1] <- match(pv$peaks[[1]][,1],pv$chrmap)
    }
    result[,2] <- Numeric(pv$peaks[[1]][,2])
    result[,3] <- Numeric(pv$peaks[[1]][,3])
    for (i in 1:numvecs) {
      result[,i + 3] <- pv$peaks[[i]][,4]
    }
    colnames(result) <-
      c("CHR","START","END",pv$class[PV_ID,1:numvecs])
    pv$binding <- result
    pv$merged <- pv$binding[,1:3]
    pv$totalMerged <- do.nrow(pv$binding)
    pv$called <- called
    pv$allcalled <- allcalled
    allnames <- pv$chrmap
  }
  
  pv$attributes <- attributes
  pv$minOverlap <- minOverlap
  
  if (is.null(allnames)) {
    allnames <- pv$chrmap[pv$binding[,1]]
  }
  
  pv$binding <- pv.check1(pv$binding)
  if (do.nrow(pv$binding) > 0) {
    vnames <- allnames[pv$binding[,1]]
  }
  if (!is.null(allnames)) {
    newmap <- sort(unique(allnames))
  } else {
    newmap <- NULL
  }
  if (do.nrow(pv$binding) > 0) {
    pv$binding[,1] <- match(vnames,newmap)
    if (is.unsorted(unique(pv$binding[,1]))) {
      pv$binding <- pv.peaksort(pv$binding)
    }
    rownames(pv$binding) <- 1:do.nrow(pv$binding)
  }
  
  pv$merged <- pv.check1(pv$merged)
  pv$called <- pv.check1(pv$called)
  pv$chrmap <- newmap
  
  pv$hc <- NULL
  pv$pc <- NULL
  pv$masks <- pv.mask(pv)
  pv$config <- as.list(pv$config)
  pv.gc()
  return(pv)
}

## pv.model -- build model, e.g. from sample sheet
pv.model <- function(model,mask,minOverlap=2,
                     samplesheet='sampleSheet.csv',config=data.frame(RunParallel=FALSE),
                     caller="raw",format, scorecol, bLowerBetter, skipLines=0,bAddCallerConsensus=TRUE,
                     bRemoveM=TRUE, bRemoveRandom=TRUE,filter,
                     attributes, dir) {
  
  if(missing(format))       format       <- NULL
  if(missing(scorecol))     scorecol     <- NULL
  if(missing(bLowerBetter)) bLowerBetter <- NULL
  if(missing(filter))       filter       <- NULL   
  
  if(!missing(model)){
    ChIPQCobj <- model$ChIPQCobj
  } else ChIPQCobj <- NULL
  
  if(!missing(model)) {
    if(missing(attributes)) {
      if(is.null(model$attributes)) {   
        attributes <- PV_ID
      } else {
        attributes <- model$attributes
      }
    }
    config <- as.list(model$config)
    allsame <- pv.allSame(model)
    model <- pv.vectors(model,mask=mask,minOverlap=minOverlap,
                        attributes=attributes,bAllSame=allsame)
    model$config <- as.list(config)
    model$ChIPQCobj <- ChIPQCobj
    model$class[DBA_REPLICATE,is.na(model$class[DBA_REPLICATE,])] <- ""
    if(!missing(mask)) {
      if(length(model$config$fragmentSize) > 1) {
        model$config$fragmentSize <- model$config$fragmentSize[mask]
      }
    }
    return(model)
  }
  
  if(missing(attributes)) {   
    attributes <- PV_ID
  }
  
  if(is.character(samplesheet)) {
    if(!missing(dir)) {
      samplesheet <- file.path(dir,samplesheet)
    }
    ext <- file_ext(samplesheet)
    if (ext %in% c("xls","xlsx")) {
      if (requireNamespace("XLConnect",quietly=TRUE)) {
        samples <- XLConnect::readWorksheetFromFile(samplesheet,sheet=1)
      } else {
        stop("Package XLConnect is needed to read Excel-format sample sheets.",call.=FALSE)
      }
    } else {
      samples <- read.table(samplesheet,sep=',',stringsAsFactors=FALSE,header=TRUE, 
                            comment.char="")
    }
    samples <- stripSpaces(samples)
  } else samples <- samplesheet
  
  if(is.null(samples$SampleID)){
    samples$SampleID <- 1:nrow(samples)
  }
  if(is.null(samples$Tissue)){
    samples$Tissue <- ""
  } 
  if(is.null(samples$Factor)){
    samples$Factor <- ""
  }
  if(is.null(samples$Condition)){
    samples$Condition <- ""
  }
  if(is.null(samples$Treatment)){
    samples$Treatment <- ""
  }
  if(is.null(samples$Replicate)){
    samples$Replicate <- ""
  }
  
  if(sum(is.na(samples$SampleID)))  samples$SampleID[is.na(samples$SampleID)]=""
  if(sum(is.na(samples$Tissue)))    samples$Tissue[is.na(samples$Tissue)]=""
  if(sum(is.na(samples$Factor)))    samples$Factor[is.na(samples$Factor)]=""
  if(sum(is.na(samples$Condition))) samples$Condition[is.na(samples$Condition)]=""
  if(sum(is.na(samples$Treatment))) samples$Treatment[is.na(samples$Treatment)]=""
  if(sum(is.na(samples$Replicate))) samples$Replicate[is.na(samples$Replicate)]=""
  
  # prepend working directory to file paths
  if(!missing(dir)) {
    if(!is.null(samples$Peaks)) {
      samples$Peaks <- sapply(samples$Peaks,
                              function(x){file.path(dir,x)})
    }
    if(!is.null(samples$bamReads)) {
      samples$bamReads <- sapply(samples$bamReads,
                                 function(x){file.path(dir,x)})
    }
    if(!is.null(samples$bamControl)) {
      samples$bamControl <- sapply(samples$bamControl,
                                   function(x){file.path(dir,x)})
    }
    if(!is.null(samples$Spikein)) {
      samples$Spikein <- sapply(samples$Spikein,
                                function(x){file.path(dir,x)})
    }
  }
  
  model <- NULL
  if(is.character(config)) {
    if(!is.null(config)) {
      config  <- read.table(config,colClasses='character',sep=',',header=TRUE)
      x <- config$DataType
      if(!is.null(x)) {
        if(x=="DBA_DATA_FRAME")	{
          config$DataType <- DBA_DATA_FRAME
        } else if(x=="DBA_DATA_RANGEDDATA"){
          config$DataType <- DBA_DATA_RANGEDDATA            
        } else {
          config$DataType <- DBA_DATA_GRANGES
        } 
      }
      x <- config$RunParallel
      if(!is.null(x)) {
        if(x=="FALSE") {
          config$RunParallel=FALSE
        } else {
          config$RunParallel=TRUE
        }
      }
    }
    config <- as.list(config)
  }
  if(is.null(config$parallelPackage)){
    config$parallelPackage=DBA_PARALLEL_MULTICORE
  } else if (config$parallelPackage == "DBA_PARALLEL_MULTICORE") {
    config$parallelPackage=DBA_PARALLEL_MULTICORE
  } else if (config$parallelPackage == "DBA_PARALLEL_RLSF") {
    config$parallelPackage=DBA_PARALLEL_RLSF  
  }
  
  if(is.null(config$AnalysisMethod)){
    config$AnalysisMethod <- DBA_DESEQ2
  } else if(is.character(config$AnalysisMethod[1])){
    x <- strsplit(config$AnalysisMethod[1],',')
    if(length(x[[1]])==1) {
      config$AnalysisMethod <- pv.getMethod(config$AnalysisMethod)
    }	 else if (length(x[[1]])==2) {
      #config$AnalysisMethod <- c(pv.getMethod(x[[1]][1]),pv.getMethod(x[[1]][2]))	
      config$AnalysisMethod <- pv.getMethod(x[[1]][1])
    }
  }
  
  model$config <- as.list(config)
  curcontrol=1
  for(i in 1:nrow(samples)) {
    if(is.null(samples$PeakCaller[i])) {
      peakcaller  <- caller
    } else if(is.na(samples$PeakCaller[i])) {
      peakcaller  <- caller
    } else {
      peakcaller <- as.character(samples$PeakCaller[i])
    }
    if(is.null(samples$PeakFormat[i])) {
      peakformat  <- format
    } else if(is.na(samples$PeakFormat[i])) {
      peakformat  <- format
    } else {
      peakformat <- as.character(samples$PeakFormat[i])
    } 
    if(is.null(samples$ScoreCol[i])) {
      peakscores  <- scorecol
    } else if(is.na(samples$ScoreCol[i])) {
      peakscores  <- scorecol
    } else {
      if(is.factor(samples$ScoreCol[i])) {
        peakscores <- as.integer(as.character(samples$ScoreCol[i]))	
      } else {
        peakscores <- as.integer(samples$ScoreCol[i])
      }
    }
    if(is.null(samples$LowerBetter[i])) {
      peaksLowerBetter  <- bLowerBetter
    } else if(is.na(samples$LowerBetter[i])) {
      peaksLowerBetter  <- bLowerBetter
    } else {
      peaksLowerBetter <- as.logical(samples$LowerBetter[i])
    }
    if(is.null(samples$Filter[i])) {
      peakfilter  <- filter
    } else if(is.na(samples$Filter[i])) {
      peakfilter  <- filter
    } else {
      if(is.factor(samples$Filter[i])) {
        peakfilter <- as.integer(as.character(samples$Filter[i]))	
      } else {
        peakfilter <- as.integer(samples$Filter[i])
      }
    }
    
    controlid  <- pv.controlID(samples,i,model$class,curcontrol)
    if(is.numeric(controlid)) {
      curcontrol <- controlid+1
      controlid <- sprintf("Control%d",controlid)
    }
    counts <- samples$Counts[i]
    if(!is.null(counts)) {
      if(is.na(counts)[1]) {
        counts <- NULL
      } else if (counts == "") {
        counts <- NULL
      }
    }
    if(!is.null(counts)) {
      peakcaller <- 'counts'
    }
    
    if(is.null(samples$Peaks[i])) {
      peakid <- NULL
    } else {
      peakid <- as.character(samples$Peaks[i])
    }
    
    message(as.character(samples$SampleID[i]),' ',
            as.character(samples$Tissue[i]),' ',
            as.character(samples$Factor[i]),' ',
            as.character(samples$Condition[i]),' ',
            as.character(samples$Treatment[i]),' ',
            as.integer(samples$Replicate[i]),' ',peakcaller)
    
    model <- pv.peakset(model,
                        peaks       = peakid,
                        sampID      = as.character(samples$SampleID[i]),
                        tissue      = as.character(samples$Tissue[i]),
                        factor      = as.character(samples$Factor[i]),
                        condition   = as.character(samples$Condition[i]),
                        treatment   = as.character(samples$Treatment[i]),
                        consensus   = FALSE,
                        peak.caller = peakcaller,
                        peak.format = peakformat,
                        scoreCol    = peakscores,
                        bLowerScoreBetter = peaksLowerBetter,
                        control     = controlid,
                        reads       = NA,
                        replicate   = as.integer(samples$Replicate[i]),
                        readBam     = as.character(samples$bamReads[i]),
                        controlBam  = as.character(samples$bamControl[i]),
                        spikein     = as.character(samples$Spikein[i]),
                        filter      = peakfilter,
                        counts      = counts,
                        bRemoveM=bRemoveM, bRemoveRandom=bRemoveRandom,
                        skipLines=skipLines)
  }
  
  model$samples <- samples
  
  if(bAddCallerConsensus){
    model <- pv.add_consensus(model)
  }
  
  model <- pv.vectors(model,mask=mask,minOverlap=minOverlap,
                      attributes=attributes,
                      bAllSame <- (peakcaller=="counts")) 
  
  model$config <- as.list(model$config)
  model$ChIPQCobj <- ChIPQCobj
  model$class[DBA_REPLICATE,is.na(model$class[DBA_REPLICATE,])]=""
  return(model)
}
