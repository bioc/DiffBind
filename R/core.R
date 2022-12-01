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
PV_SPIKEIN    <- 13

PV_DEBUG <- FALSE

###########################################



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
  
  
  res <- matrix(pv$class[attributes,mask], sum(mask), length(attributes),
                byrow=TRUE)
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
    cols <- colnames(res)
    res <- cbind(res,res[,movecol])
    res <- res[,-movecol]
    if(is.null(dim(res))) {
      cols <- cols[-movecol]
      cols <- c(cols, "Reads")
      res <- matrix(res,1,length(res))
      colnames(res) <- cols
    }
    colnames(res)[ncol(res)] <- "Reads"
  }
  
  if (bSN) {
    if (!is.null(pv$SN)) {
      res <- cbind(res,pv$SN[mask])
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
  
  res <- data.frame(res)
  
  if("Reads" %in% names(res)) {
    res$Reads <- as.numeric(res$Reads)
  }
  
  if("FRiP" %in% names(res)) {
    res$FRiP <- as.numeric(res$FRiP)
  }
  
  if("Intervals" %in% names(res)) {
    res$Intervals <- as.numeric(res$Intervals)
  }
  
  return(res)
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
        (pv$totalMerged != do.nrow(pv$binding))) {
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
  
  if (do.nrow(tmp$binding) == 0) {
    message(">> zero peaks in consensus, skipping.")
    return(pv)
  }
  
  if (minOverlap > 0 && minOverlap < 1) {
    minOverlap <- ceiling(length(tmp$peaks) * minOverlap)
  }
  
  if(do.nrow(tmp$binding) == 1) {
    pscores <- matrix(tmp$binding[,4:ncol(tmp$binding)],1)
  } else {
    pscores <- tmp$binding[,4:ncol(tmp$binding)]
  }
  
  goodvecs <-
    apply(pscores,1,pv.minOverlap,minOverlap)
  
  tmp$binding  <- matrix(tmp$binding[goodvecs,],sum(goodvecs))
  pscores      <- matrix(pscores[goodvecs,],sum(goodvecs))
  
  mean.density <- apply(pscores,1,pv.domean)
  if(do.nrow(tmp$binding) == 1) {
    tmp$binding <- matrix(c(tmp$binding[1,1:3],mean.density),1)
  } else {
    tmp$binding  <- cbind(tmp$binding[,1:3],mean.density)
  }
  
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
    spikein      =  pv.getoneorNA(pv$class[PV_SPIKEIN, sampvec]),
    scoreCol     =  0
  )
  
  pv$binding    <- NULL
  pv$hc         <- NULL
  pv$pc         <- NULL
  
  return(pv)
}

pv.consensusSets <- function(pv,peaks = NULL,minOverlap,attributes,
                             tissue,factor,condition,treatment,replicate,control,peak.caller,
                             readBam, controlBam, spikein)	{
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
      if (!missing(spikein))
        pv$class[PV_SPIKEIN,sampnum]    <- spikein
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
  
  if (pv$totalMerged != do.nrow(pv$binding)) {
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
  if(!is.null(pv$allcalled)) {
    pv$called <- pv$allcalled
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
  
  if(length(unique(sapply(DBA$peaks,nrow))) > 1) {
    return(FALSE)
  }
  
  return(TRUE)
}
