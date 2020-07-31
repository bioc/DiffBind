#####################################
## pv_core.R -- Peak Vectorization ##
## 20 October 2009                 ##
## 3 February 2011 -- packaged     ##
## Rory Stark                      ##
## Cancer Research UK              ##
#####################################

############################
## Top Level Entry Points ##
############################
## pv.peakset       -- add a peakset w/ attributes to the model
## pv.vectors       -- build the binding expression vectors and do clustering/PCA analysis
## pv.list          -- list peaksets w/ attributes in model
## pv.consensus     -- add a consensus peakset based on peaksets already in model

## pv.mask          -- create a mask to define a subset of peaksets in a model
## pv.whichSites    -- create a mask of sites belonging to specific peakset(s)

## pv.plotClust     -- hierarchical cluster plot
## pv.sort          -- sort binding sites (e.g. for heatmap)

## pv.overlap       -- generate overlapping/unique peaksets

## pv.occupancy     -- generate occupancy statistics for peaksets in a model
## pv.plotScatter   -- scatter plots of contrasts


################
## Constants  ##
################
PV_GROUP      <- 0
PV_ID         <- 1
PV_TISSUE     <- 2
PV_FACTOR     <- 3
PV_CONDITION  <- 4
PV_CONSENSUS  <- 5
PV_CALLER     <- 6
PV_CONTROL    <- 7
PV_READS      <- 8
PV_REPLICATE  <- 9
PV_BAMREADS   <- 10
PV_BAMCONTROL <- 11
PV_TREATMENT  <- 12
PV_INTERVALS  <- 13
PV_SN_RATIO   <- 14

PV_DEBUG <- FALSE

###########################################

## pv.peakset -- add a peakset to the model
pv.peakset <- function(pv = NULL,peaks, sampID, tissue, factor,condition, treatment, replicate,
                       control, peak.caller, peak.format, reads = 0, consensus =FALSE, 
                       readBam, controlBam, scoreCol = NULL, bLowerScoreBetter = NULL,
                       bRemoveM = TRUE, bRemoveRandom = TRUE,
                       minOverlap = 2,bFast = FALSE,bMakeMasks = TRUE,
                       skipLines = 1, filter = NULL, counts = NULL) {
  zeroVal <- -1
  bLog <- FALSE
  
  if (missing(peaks)) {
    peaks <- NULL
  }
  
  if (!is.null(pv$peaks) && length(peaks) == 0) {
    peaks <- 1:length(pv$peaks)
  }
  
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
      pv.consensusSets(
        pv,peaks = peaks,minOverlap = minOverlap,attributes = consensus,
        tissue,factor,condition,treatment,replicate,control,peak.caller,
        readBam, controlBam
      )
    
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
  if (missing(readBam))
    readBam <- NA
  if (length(readBam) == 0)
    readBam <- NA
  if (missing(controlBam))
    controlBam <- NA
  if (length(controlBam) == 0)
    controlBam <- NA
  
  if (!is.null(peaks) && length(peaks) <= 1) {
    if (is.na(peaks)) {
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
        scoreCol <- 0
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
        pv.normalize(peaks,scoreCol,zeroVal = zeroVal,bLog = bLog)
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
        reads,replicate,readBam,controlBam,treatment
      )
    )
  colnames(clascol) <- sampID
  pv$class <- cbind(pv$class,clascol)
  rownames(pv$class) <-
    c(
      "ID","Tissue","Factor","Condition", "Consensus",
      "Peak caller","Control","Reads","Replicate","bamRead",
      "bamControl","Treatment"
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
  
  called <- SN <- NULL
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
      if(length(mask)<=ncol(pv$called)) {
        called   <- pv$called[,mask]
      }
    }
    if(!is.null(pv$SN)) {
      if(length(mask)==length(SN)) {
        SN       <- pv$SN[mask]
      }
    }
    score      <- pv$score
    summits    <- pv$summits
    #contrasts  <- pv$contrasts
    blacklist  <- pv$blacklist
    greylist   <- pv$greylist
    peaks.blacklisted <- pv$peaks.blacklisted
    #annotation <- pv$annotation
    
    pv <- NULL
    pv$peaks      <- peaks
    pv$class      <- class
    pv$chrmap     <- chrmap
    pv$config     <- config
    pv$samples    <- samples
    pv$called     <- called
    #pv$contrasts  <- contrasts
    pv$score      <- score
    pv$SN         <- SN
    pv$summits    <- summits
    pv$blacklist  <- blacklist
    pv$greylist   <- greylist
    pv$peaks.blacklisted <- peaks.blacklisted
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
      res  <- pv.merge(allpeaks,peaks,pv$class)
      pv$called <- res$included
      pv$totalMerged <- nrow(res$merged)
      rownames(res$merged) <- 1:nrow(res$merged)
      allnames <- res$chrmap
      pv$merged <- res$merged[,1:3]
      if ((ncol(res$merged) > 4) && (minOverlap > 1)) {
        pv$binding <- res$merged[pv.overlaps(pv,minOverlap),]
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
    result[,2] <- pv$peaks[[1]][,2]
    result[,3] <- pv$peaks[[1]][,3]
    for (i in 1:numvecs) {
      result[,i + 3] <- pv$peaks[[i]][,4]
    }
    colnames(result) <-
      c("CHR","START","END",pv$class[PV_ID,1:numvecs])
    pv$binding <- result
    pv$merged <- pv$binding[,1:3]
    pv$totalMerged <- nrow(pv$binding)
    pv$called <- called
    allnames <- pv$chrmap
  }
  
  pv$attributes <- attributes
  pv$minOverlap <- minOverlap
  
  if (is.null(allnames)) {
    allnames <- pv$chrmap[pv$binding[,1]]
  }
  
  pv$binding <- pv.check1(pv$binding)
  if (nrow(pv$binding) > 0) {
    vnames <- allnames[pv$binding[,1]]
  }
  if (!is.null(allnames)) {
    newmap <- sort(unique(allnames))
  } else {
    newmap <- NULL
  }
  if (nrow(pv$binding) > 0) {
    pv$binding[,1] <- match(vnames,newmap)
    if (is.unsorted(unique(pv$binding[,1]))) {
      pv$binding <- pv.peaksort(pv$binding)
    }
    rownames(pv$binding) <- 1:nrow(pv$binding)
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

## pv.list -- list attributes of samples in model
pv.deflist <- c(
  PV_ID,PV_TISSUE,PV_FACTOR,PV_CONDITION,PV_TREATMENT,
  PV_REPLICATE,PV_CALLER,PV_INTERVALS,PV_READS,PV_SN_RATIO
)

pv.list <- function(pv,mask,bContrasts = FALSE, bDesign=FALSE,
                    attributes = pv.deflist, th = 0.05) {
  
  if (!missing(mask)) {
    if (!is.logical(mask)) {
      tmp  <- rep (FALSE,length(pv$peaks))
      tmp[mask] <- TRUE
      mask <- tmp
    }
  }
  
  if (bContrasts) {
    return(pv.listContrasts(pv,th = th))
  }
  
  if (bDesign) {
    if(is.null(pv$design)) {
      warning("No design present.",call.=FALSE)
      return(NULL)
    } else {
      return(pv$design)
    }
  }
  
  if (missing(attributes)) {
    attributes <- pv.deflist
  }
  
  if (missing(mask)) {
    mask <- rep(TRUE,ncol(pv$class))
  }
  
  if (!is.logical(mask)) {
    newm <- rep(FALSE,length(pv$peaks))
    for (ps in mask) {
      newm[ps] <- TRUE
    }
    mask <- newm
  }
  
  if (PV_INTERVALS %in% attributes) {
    attributes <- attributes[-which(attributes %in% PV_INTERVALS)]
    bIntervals <- TRUE
  } else
    bIntervals <- FALSE
  
  if (PV_SN_RATIO %in% attributes) {
    attributes <- attributes[-which(attributes %in% PV_SN_RATIO)]
    bSN <- TRUE
  } else
    bSN <- FALSE
  
  
  res <- t(pv$class[attributes,mask])
  colnames(res) <- sapply(attributes,pv.attname,pv)
  rownames(res) <- which(mask)
  
  if (bIntervals) {
    intervals <- NULL
    for (i in 1:length(mask)) {
      if (mask[i]) {
        intervals <- c(intervals,nrow(pv$peaks[[i]]))
      }
    }
    res <- cbind(res,intervals)
    colnames(res)[ncol(res)] <- 'Intervals'
  }
  
  if(PV_READS %in% attributes) {
    movecol <- which(attributes %in% PV_READS)
    res <- cbind(res,res[,movecol])
    res <- res[,-movecol]
    colnames(res)[ncol(res)] <- "Reads"
  }
  
  if (bSN) {
    if (!is.null(pv$SN)) {
      res <- cbind(res,pv$SN)
      colnames(res)[ncol(res)] <- 'FRiP'
    }
  }
  
  j <- ncol(res)
  if (nrow(res) > 1) {
    for (i in j:1) {
      x <- unique(res[,i])
      if (colnames(res)[i] == 'Peak caller') {
        if (all.equal(attributes,pv.deflist) == TRUE) {
          if (length(x) == 1) {
            res <- res[,-i]
          }
        }
      } else if (length(x) == 1) {
        if (is.na(x)) {
          res <- res[,-i]
        } else if (x[1] == "") {
          res <- res[,-i]
        }
      }
    }
  }
  
  return(data.frame(res))
}

## pv.consensus -- add a consensus peakset based on peaksets already in model
pv.consensus <- function(pv,sampvec,minOverlap = 2,
                         bFast = FALSE,sampID) {
  if (missing(sampvec)) {
    sampvec <- 1:length(pv$peaks)
  }
  if (is.null(sampvec)) {
    sampvec <- 1:length(pv$peaks)
  }
  if (is(sampvec,"logical")) {
    sampvec <- which(sampvec)
  }
  
  tmp <- NULL
  if (bFast & (max(sampvec) <= ncol(pv$class)))  {
    if ((length(sampvec) < length(pv$peaks)) ||
        (pv$totalMerged != nrow(pv$binding))) {
      pv <- pv.vectors(pv,sampvec,minOverlap = 1)
      sampvec <- 1:length(pv$peaks)
    } else {
      pv <- pv.check(pv)
    }
    sites <- pv.whichSites(pv,sampvec,bUseAllvecs = FALSE)
    tmp$binding <- pv$binding[sites,c(1:3,(sampvec + 3))]
  } else {
    peaklist  <- NULL
    classlist <- NULL
    chrs <- NULL
    for (samp in sampvec) {
      peaklist  <- pv.listadd(peaklist,pv$peaks[[samp]])
      chrs <- sort(unique(c(chrs,unique(pv$peaks[[samp]][,1]))))
    }
    tmp$peaks  <- peaklist
    tmp$class  <- pv$class[, sampvec]
    tmp$chrmap <- chrs
    tmp <- pv.vectors(tmp,minOverlap = minOverlap)
  }
  
  if (minOverlap > 0 && minOverlap < 1) {
    minOverlap <- ceiling(length(tmp$peaks) * minOverlap)
  }
  
  goodvecs <-
    apply(tmp$binding[,4:ncol(tmp$binding)],1,pv.minOverlap,minOverlap)
  
  tmp$binding  <- tmp$binding[goodvecs,]
  
  mean.density <-
    apply(tmp$binding[,4:ncol(tmp$binding)],1,pv.domean)
  tmp$binding <- cbind(tmp$binding[,1:3],mean.density)
  
  #kludge to get peakset in correct format
  tmpf <- tempfile(as.character(Sys.getpid())) #tmpdir='.')
  pv.do_peaks2bed(tmp$binding, tmp$chrmap,tmpf)
  tmp$binding <- pv.readbed(tmpf)
  unlink(tmpf)
  
  if (length(unique(pv$class[PV_REPLICATE, sampvec])) == 1) {
    replicate <- unique(pv$class[PV_REPLICATE, sampvec])
  } else {
    replicate <- pv.catstr(pv$class[PV_REPLICATE, sampvec])
  }
  if (missing(sampID)) {
    sampID <- pv.catstr(pv$class[PV_ID, sampvec])
  }
  pv <- pv.peakset(
    pv, peaks = tmp$binding,
    sampID       =  sampID,
    tissue       =  pv.catstr(pv$class[PV_TISSUE, sampvec]),
    factor       =  pv.catstr(pv$class[PV_FACTOR, sampvec]),
    condition    =  pv.catstr(pv$class[PV_CONDITION, sampvec]),
    treatment    =  pv.catstr(pv$class[PV_TREATMENT, sampvec]),
    peak.caller  =  pv.catstr(pv$class[PV_CALLER, sampvec]),
    control      =  pv.catstr(pv$class[PV_CONTROL, sampvec]),
    reads        =  mean(as.numeric(pv$class[PV_READS, sampvec])),
    replicate    =  replicate,
    consensus    =  TRUE,
    readBam      =  pv.getoneorNA(pv$class[PV_BAMREADS, sampvec]),
    controlBam   =  pv.getoneorNA(pv$class[PV_BAMCONTROL, sampvec]),
    scoreCol     =  0
  )
  
  pv$binding    <- NULL
  pv$hc         <- NULL
  pv$pc         <- NULL
  
  return(pv)
}

pv.consensusSets <- function(pv,peaks = NULL,minOverlap,attributes,
                             tissue,factor,condition,treatment,replicate,control,peak.caller,
                             readBam, controlBam)	{
  if (is.character(peaks)) {
    stop(
      "\"peaks\" parameter can not be a filename when \"consensus\" specifies attributes",
      call. =FALSE
    )
  } 
  
  if (is.null(peaks)) {
    peaks <- rep(TRUE,ncol(pv$class))
  }
  
  include <- FALSE
  exclude <- FALSE
  if (sum(attributes < 0))
    exclude <- TRUE
  if (sum(attributes > 0))
    include <- TRUE
  
  if (include & exclude) {
    stop(
      'Consensus attributes must be all inclusive (positive) or all exclusive (negative)',
      call. = FALSE
    )
  }
  
  if (exclude) {
    atts <- NULL
    if (!(-PV_TISSUE    %in% attributes))
      atts <- c(atts,PV_TISSUE)
    if (!(-PV_FACTOR    %in% attributes))
      atts <- c(atts,PV_FACTOR)
    if (!(-PV_CONDITION %in% attributes))
      atts <- c(atts,PV_CONDITION)
    if (!(-PV_TREATMENT %in% attributes))
      atts <- c(atts,PV_TREATMENT)
    if (!(-PV_REPLICATE %in% attributes))
      atts <- c(atts,PV_REPLICATE)
    if (!(-PV_CALLER    %in% attributes))
      atts <- c(atts,PV_CALLER)
    attributes <- atts
  }
  
  numatts <- length(attributes)
  class <- pv$class[attributes,]
  if (is.vector(class)) {
    class <- matrix(class,1,length(class))
  }
  sampids <- pv$class[PV_ID,]
  
  specs <- unique(class[,peaks],MARGIN = 2)
  if (is.vector(specs)) {
    specs <- matrix(specs,1,length(specs))
  }
  
  if (ncol(specs) == ncol(class)) {
    warning(
      'All peaksets unique for specified attributes; no consensus peaksets added.',
      call. = FALSE
    )
    return(pv)
  }
  for (i in 1:ncol(specs)) {
    cand <- class %in% specs[,i]
    if (is.vector(cand)) {
      cand <- matrix(cand,numatts,ncol(class))
    }
    samples <-
      apply(cand,MARGIN = 2,function(x) {
        sum(x) == numatts
      }) & peaks
    diffatts <-
      apply(class,MARGIN = 1,function(x) {
        length(unique(x)) > 1
      })
    if (sum(samples) > 1) {
      message('Add consensus: ',paste(specs[diffatts,i],collapse = " "))
      if (length(unique(sampids[samples])) == 1) {
        sampid <- sampids[samples][1]
      } else
        sampid <- paste(specs[diffatts,i],collapse = ":")
      pv <-
        pv.consensus(pv,samples,sampID = sampid,minOverlap = minOverlap)
      sampnum <- ncol(pv$class)
      if (pv$class[PV_ID,sampnum] == "")
        pv$class[PV_ID,sampnum] <- "ALL"
      if (!missing(tissue))
        pv$class[PV_TISSUE,sampnum]     <- tissue
      if (!missing(factor))
        pv$class[PV_FACTOR,sampnum]     <- factor
      if (!missing(condition))
        pv$class[PV_CONDITION,sampnum]  <- condition
      if (!missing(treatment))
        pv$class[PV_TREATMENT,sampnum]  <- treatment
      if (!missing(replicate))
        pv$class[PV_REPLICATE,sampnum]  <- replicate
      if (!missing(control))
        pv$class[PV_CONTROL,sampnum]    <- control
      if (!missing(peak.caller))
        pv$class[PV_CALLER,sampnum]     <- peak.caller
      if (!missing(readBam))
        pv$class[PV_BAMREADS,sampnum]   <- readBam
      if (!missing(controlBam))
        pv$class[PV_BAMCONTROL,sampnum] <- controlBam
    }
  }
  
  return(pv)
}


## pv.mask -- create a mask to define a subset of peaksets in a model
pv.mask <-
  function(pv,attribute,value,combine = 'or',mask,merge = 'or',bApply = FALSE) {
    numsamps <- ncol(pv$class)
    
    if (missing(mask)) {
      if ((merge == 'or') | (merge == 'and')) {
        mask <- rep(FALSE, numsamps)
      } else {
        mask <- rep(TRUE, numsamps)
      }
    }
    
    if (missing(attribute)) {
      masks <- NULL
      for (att in c(PV_TISSUE,PV_FACTOR,PV_CONDITION,PV_TREATMENT,PV_CALLER)) {
        vals <- unique(pv$class[att,])
        for (v in vals) {
          res <- list(x = pv.mask(pv,att,v))
          names(res) <- v
          masks <- c(masks,res)
        }
      }
      res <- list(x = pv.mask(pv,PV_CONSENSUS,TRUE))
      if (sum(res[[1]])) {
        names(res) <- "Consensus"
        masks <- c(masks,res)
      }
      reps <- unique(pv$class[PV_REPLICATE,])
      reps <- reps[!is.na(reps)]
      if (length(reps > 1)) {
        for (rep in reps) {
          res <- list(x = pv.mask(pv,PV_REPLICATE,rep))
          names(res) <- sprintf("Replicate.%s",rep)
          masks <- c(masks,res)
        }
      }
      
      masks$All  <- rep(TRUE,ncol(pv$class))
      masks$None <- rep(FALSE,ncol(pv$class))
      
      return(masks)
    }
    
    curmask <- NULL
    for (v in value) {
      newmask <- pv$class[attribute,] == v
      if (is.null(curmask))
        curmask <- newmask
      if ((combine == 'or') | (combine == 'nor')) {
        curmask <- curmask | newmask
      }
      if ((combine == 'and') | (combine == 'nand')) {
        curmask <- curmask & newmask
      }
    }
    
    if ((combine == 'nor') | (combine == 'nand')) {
      curmask <- !curmask
    }
    
    if (merge == 'or') {
      mask <- mask | curmask
    }
    
    if (merge == 'and') {
      mask <- mask & curmask
    }
    
    if (merge == 'nor') {
      mask <- !(mask | curmask)
    }
    
    if (merge == 'nand') {
      mask <- !(mask & curmask)
    }
    if (bApply) {
      pv <- dba(pv, mask)
      return(pv)
    } else {
      return(mask)
    }
  }

## pv.whichSites -- return index vector of sites belonging to a specific peakset
pv.whichSites <-
  function(pv,pnum,combine = "or",minVal = -1,bUseAllvecs = FALSE) {
    if (bUseAllvecs) {
      stop("Internal error: bUseAllvecs=FALSE in pv.whichsites, please report.")
    } else {
      vecs <- pv$binding
    }
    if (length(pnum) == 1) {
      res <- vecs[,pnum + 3] > minVal
    } else {
      res <- vecs[,pnum[1] + 3] > minVal
      for (p in 2:length(pnum)) {
        newvec <- vecs[,pnum[p] + 3] > minVal
        if (combine == 'or') {
          res <- res | newvec
        }
        if (combine == 'and') {
          res <- res & newvec
        }
        if (combine == 'nor') {
          res <- !(res | newvec)
        }
        if (combine == 'nand') {
          res <- !(res & newvec)
        }
      }
    }
    return(res)
  }

## pv.plotClust  -- hierarchical cluster plot
pv.plotClust <-
  function(pv,mask,numSites,sites,attributes = pv$attributes,distMeth = "pearson") {
    if (missing(mask)) {
      mask <- rep(TRUE,ncol(pv$class))
    }
    if (missing(sites)) {
      sites <- rep(TRUE,nrow(pv$binding))
    }
    if (missing(numSites)) {
      numSites <- length(sites)
    }
    pv$binding <- pv$binding[sites,mask][1:numSites,]
    pv$class   <- pv$class[,mask]
    pv         <-
      pv.analysis(pv,attributes,bPCA = FALSE,distMeth = distMeth)
    plot(pv$hc)
    return(pv$hc)
  }


## pv.overlap -- generate overlapping/unique peaksets
pv.overlap <- function(pv,mask,bFast = FALSE,minVal = 0) {
  if (!missing(mask)) {
    if (!is.logical(mask)) {
      peaksets <- mask
    }	else {
      peaksets <- which(mask)
    }
    if (length(peaksets) <= 4) {
      A <- peaksets[1]
      B <- peaksets[2]
      if (length(peaksets) >= 3) {
        C <- peaksets[3]
      }
      if (length(peaksets) == 4) {
        D <- peaksets[4]
      }
    } else {
      warning('Too many peaksets in mask.',call.=FALSE)
      return(NULL)
    }
  } else {
    stop('Must specify mask for peaksets to overlap.',call. = FALSE)
  }
  
  if (pv$totalMerged != nrow(pv$binding)) {
    pv <- pv.vectors(pv,mask = peaksets,minOverlap=1)
    peaksets <- 1:length(peaksets)
    A <- 1; B <- 2; C <- 3; D <- 4
  }
  if (is.null(pv$called)) {
    stop("Called masks not present; re-run dba.count",call.=FALSE)
  }
  maskA <- pv$called[,A]
  maskB <- pv$called[,B]
  if (length(peaksets) >= 3)
    maskC <- pv$called[,C]
  if (length(peaksets) == 4)
    maskD <- pv$called[,D]
  
  pv$binding <- data.frame(pv$binding)
  pv$binding[,1] <- pv$chrmap[pv$binding[,1]]
  if (length(peaksets) < 4) {
    if (length(peaksets) < 3) {
      res <- pv.contrast2(pv$binding,A,B,v1 = maskA,v2 = maskB)
    } else {
      res <-
        pv.contrast3(pv$binding,A,B,C,v1 = maskA,v2 = maskB,v3 = maskC)
    }
  } else {
    res <-
      pv.contrast4(
        pv$binding,A,B,C,D,v1 = maskA,v2 = maskB,v3 = maskC,v4 = maskD
      )
  }
  
  return(res)
}


## pv.sort  - sort binding sites (e.g. for heatmap)
pv.sort <- function(pv,fun = sd,mask,...) {
  if (missing(mask)) {
    mask <- rep(TRUE,ncol(pv$class))
  }
  
  scores <- apply(pv$binding[,c(FALSE,FALSE,FALSE,mask)],1,fun,...)
  ranked <- order(scores,decreasing = TRUE)
  
  pv$binding   <- pv$binding[ranked,]
  
  return(pv)
}

pv.overlapRate <- function(pv,mask = mask) {
  if (!missing(mask)) {
    if (!is.null(mask)) {
      pv <- dba(pv,mask)
    }
  }
  if (is.null(pv$called)) {
    stop("No peak overlap information available",call.=FALSE)
  }
  sums <- apply(pv$called,1,sum)
  res <- sapply(1:length(pv$peaks),function(x)
    sum(sums >= x))
  return(res)
}


## pv.occupancy-- generate co-occupancy stats from peaksets in a model
pv.occupancy <-
  function(pv,mask,sites,byAttribute,Sort = 'inall',CorMethod = "pearson",
           labelAtts = pv$attributes,bPlot = FALSE,minVal = 0,bCorOnly = FALSE,
           bNonZeroCors = FALSE,chrmask) {
    pv <- pv.check(pv)
    
    vecs <- pv$binding
    
    if (missing(mask)) {
      mask <- rep(TRUE,ncol(pv$class))
    } else if (is.null(mask)) {
      mask <- rep(TRUE,ncol(pv$class))
    } else {
      if (!is.logical(mask)) {
        tmp  <- rep (FALSE,length(pv$peaks))
        tmp[mask] <- TRUE
        mask <- tmp
      }
    }
    if (missing(sites)) {
      sites <- 1:nrow(vecs)
    }
    
    res <- NULL
    if (missing(byAttribute)) {
      if (length(sites) < nrow(vecs)) {
        pv$binding <- vecs[sites,]
      }
      if (!missing(chrmask)) {
        chrindex <- match(chrmask,pv$chrmap)
        vecindex <- pv$binding[,1] == chrmask
        pv$binding <- pv$binding[vecindex,]
      }
      res <-
        pv.pairs(
          pv,mask=mask,CorMethod=CorMethod,bPlot=bPlot,minVal=minVal,bCorOnly =
            bCorOnly,bNonZeroCors=bNonZeroCors
        )
      if (!is.null(nrow(res))) {
        if (!missing(labelAtts)) {
          res <- pv.overlapToLabels(pv,res,labelAtts)
        }
      }
    } else {
      ## by attribute
      vals <- unique(pv$class[byAttribute,mask])
      for (i in 1:length(vals)) {
        comps <- which(pv$class[byAttribute,] %in% vals[i])
        vmask <- rep(FALSE,length(mask))
        vmask[comps] <- TRUE
        if (sum(vmask) > 1) {
          res <- rbind(
            res,pv.occupancy(
              pv,mask=vmask,sites=sites,
              Sort=Sort,CorMethod=CorMethod,minVal =
                minVal,bCorOnly=bCorOnly
            )
          )
        }
      }
    }
    
    if (!is.null(nrow(res))) {
      if (Sort == 'cor') {
        res <- res[pv.orderfacs(res[,6],decreasing=TRUE),]
      } else if (Sort == 'percent') {
        res <- res[pv.orderfacs(res[,7],decreasing=TRUE),]
      } else {
        res <- res[pv.orderfacs(res[,5],decreasing=TRUE),]
      }
    }
    
    return(res)
  }


pv.isConsensus <- function(DBA) {
  if(is.null(DBA$class)) {
    return(FALSE)
  }
  if(sum(DBA$class[DBA_CALLER, ]=="counts") != length(DBA$peaks)) {
    return(FALSE)
  } 
  if(is.null(DBA$peaks)) {
    return(FALSE)
  }   
  if (sum(sapply(DBA$peaks,function(x)nrow(x)!=nrow(DBA$peaks[[1]])))) {
    return(FALSE)
  }
  return(TRUE)
}
