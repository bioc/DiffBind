PV_BLACKLIST_HG19   <- 200
PV_BLACKLIST_HG38   <- 201
PV_BLACKLIST_GRCH37 <- 202
PV_BLACKLIST_GRCH38 <- 203
PV_BLACKLIST_MM9    <- 204
PV_BLACKLIST_MM10   <- 205
PV_BLACKLIST_CE10   <- 206
PV_BLACKLIST_CE11   <- 207
PV_BLACKLIST_DM3    <- 208
PV_BLACKLIST_DM6    <- 209

pv.BlackGreyList <- function (DBA, blacklist, greylist,
                              Retrieve, cores) {
  
  allsame <- pv.allSame(DBA)
  
  if(!missing(Retrieve)) {
    if(!missing(blacklist) || !missing(greylist)) 
      stop('blacklist and greylist parameters must be missing to retrieve blacklisted peaks.',
           call.=FALSE)
    
    if(Retrieve==DBA_BLACKLIST){
      if(is.null(DBA$blacklist)) {
        stop("No blacklist.",call.=FALSE)
      }
      return(DBA$blacklist)
    } else if(Retrieve==DBA_GREYLIST) {
      if(is.null(DBA$greylist)) {
        stop("No greylist",call.=FALSE)
      }
      return(DBA$greylist)
    } else if(Retrieve==DBA_BLACKLISTED_PEAKS) {
      if(is.null(DBA$peaks.blacklisted)) {
        stop("No blacklisted peaks.",call.=FALSE)
      }
      return(DBA$peaks.blacklisted)
    } else {
      stop("Invalid value for Retrieve parameter.",call.=FALSE)
    }
  }
  
  if(allsame) {
    originalPeaks    <- nrow(DBA$binding)
  } else {
    originalPeaks     <- sum(unlist(lapply(DBA$peaks,nrow)))
    originalMerged    <- nrow(DBA$merged)
    originalConsensus <- nrow(DBA$binding)
  }
  
  if(!missing(blacklist)) {
    message("Applying blacklist...")
    DBA <- pv.blacklist(DBA, blacklist=blacklist, allsame=allsame)
  }
  
  if(!missing(greylist)) {
    DBA <- pv.greylist(DBA, greylist=greylist, allsame=allsame, cores=cores)
  }
  
  
  if(allsame) {
    endPeaks <- nrow(DBA$peaks[[1]])
  } else {
    endPeaks     <- sum(unlist(lapply(DBA$peaks,nrow)))
  }
  
  if(endPeaks < originalPeaks) {
    
    if(allsame) {
      DBA <- pv.removeBlacklistedPeaks(DBA)
      DBA <- pv.reNormalize(DBA)
    } else {
      DBA <- dba(DBA)
    } 
    
    endMerged    <- nrow(DBA$merged)
    endConsensus <- nrow(DBA$binding)
    
    if(!missing(blacklist) && !missing(greylist)) {
      if(allsame) {
        DBA$peaks.blacklisted <- GRangesList(lapply(DBA$peaks.blacklisted,
                                                    function(x){
                                                      x[,c("cReads","Reads","Score")]
                                                    }))
        msgstring <- sprintf("Removed %d (of %d) consensus peaks.",
                             originalPeaks-endPeaks,originalPeaks)
        DBA$SN <- pv.Signal2Noise(DBA)
      } else {
        mergedRemoved    <- originalMerged - endMerged
        consensusRemoved <- originalConsensus - endConsensus
        msgstring <- sprintf("Removed: %d merged (of %d) and %d (of %d) consensus.",
                             mergedRemoved, originalMerged, 
                             consensusRemoved, originalConsensus)
      }
      message(msgstring)
    }
  } else {
    message("No intervals removed.")
  }
  
  return(DBA)
  
}

pv.blacklist <- function(pv, blacklist, allsame=FALSE) {
  
  if(missing(blacklist)) {
    pv$config$blacklist <- NULL
    return(pv)
  }
  
  if(!is(blacklist,"GRanges")) {
    
    ce10.blacklist   <- dm3.blacklist    <- NULL
    grch37.blacklist <- grch38.blacklist <- NULL
    hg19.blacklist   <- hg38.blacklist   <- NULL
    mm10.blacklist   <- mm9.blacklist    <- NULL
    ce11.blacklist   <- dm6.blacklist    <- NULL
    
    if(blacklist==PV_BLACKLIST_HG19) {
      load(system.file("data/hg19.blacklist.RData",package="GreyListChIP"),
           envir = environment())
      blacklist <- hg19.blacklist
    } else if(blacklist==PV_BLACKLIST_HG38) {
      load(system.file("data/hg38.blacklist.RData",package="GreyListChIP"),
           envir = environment())
      blacklist <- hg38.blacklist
    } else if(blacklist==PV_BLACKLIST_GRCH37) {
      load(system.file("data/grch37.blacklist.RData",package="GreyListChIP"),
           envir = environment())
      blacklist <- grch37.blacklist
    } else if(blacklist==PV_BLACKLIST_GRCH38) {
      load(system.file("data/grch38.blacklist.RData",package="GreyListChIP"),
           envir = environment())
      blacklist <- grch38.blacklist
    } else if(blacklist==PV_BLACKLIST_MM9) {
      load(system.file("data/mm9.blacklist.RData",package="GreyListChIP"),
           envir = environment())
      blacklist <- mm9.blacklist
    } else if(blacklist==PV_BLACKLIST_MM10) {
      load(system.file("data/mm10.blacklist.RData",package="GreyListChIP"),
           envir = environment())
      blacklist <- mm10.blacklist
    } else if(blacklist==PV_BLACKLIST_CE10) {
      load(system.file("data/ce10.blacklist.RData",package="GreyListChIP"),
           envir = environment())
      blacklist <- ce10.blacklist
    } else if(blacklist==PV_BLACKLIST_CE11) {
      load(system.file("data/ce11.blacklist.RData",package="GreyListChIP"),
           envir = environment())
      blacklist <- ce11.blacklist
    } else if(blacklist==PV_BLACKLIST_DM3) {
      load(system.file("data/dm3.blacklist.RData",package="GreyListChIP"),
           envir = environment())
      blacklist <- dm3.blacklist
    } else if(blacklist==PV_BLACKLIST_DM6) {
      load(system.file("data/dm6.blacklist.RData",package="GreyListChIP"),
           envir = environment())
      blacklist <- dm6.blacklist
    }
  }
  
  # check that at least one chr matches
  snames <- as.character(unique(seqnames(blacklist)))
  if(sum(pv$chrmap %in% snames)==0) {
    warning('Blacklist does not overlap any peak chromosomes!',call.=FALSE)
  }
  
  # Apply blacklist to peaksets
  
  if(allsame) {
    totalPeaks <- nrow(pv$peaks[[1]])
  } else {
    totalPeaks <- sum(unlist(lapply(pv$peaks,nrow)))
  }
  
  pv <- pv.applyBlacklist(pv, blacklist)
  
  if(allsame) {
    totalRemoved   <- totalPeaks - nrow(pv$peaks[[1]])
  } else {
    totalRemoved   <- totalPeaks - sum(unlist(lapply(pv$peaks,nrow)))
  }
  msgstring <- sprintf("Removed: %d of %d intervals.",
                       totalRemoved, totalPeaks)
  message(msgstring)
  
  pv$blacklist <- blacklist
  
  return(pv)
  
}

pv.applyBlacklist <- function(pv, blacklist) {
  
  blacklisted <- GRangesList(lapply(lapply(pv$peaks,pv.doBlacklist,
                                           blacklist,bReturnFiltered=TRUE),
                                    GRanges))
  names(blacklisted) <- pv$class[PV_ID,]
  if(is(pv$peaks.blacklisted,"GRangesList")) {
    for(i in 1:length(blacklisted)) {
      pv$peaks.blacklisted[[i]] <- 
        sort(c(pv$peaks.blacklisted[[i]], blacklisted[[i]]))
    }
  } else {
    pv$peaks.blacklisted <- blacklisted
  }
  
  pv$peaks <- lapply(pv$peaks,pv.doBlacklist,
                     blacklist,bReturnFiltered=FALSE)
  return(pv)
}

pv.doBlacklist <- function(peakset, blacklist, bReturnFiltered=FALSE){
  peakset <- GRanges(peakset)
  if(!bReturnFiltered) {
    peakset <- peakset[!peakset %over% blacklist]
  } else {
    peakset <- peakset[peakset %over% blacklist]
  }
  peakset <- data.frame(peakset)
  if(nrow(peakset)==0) {
    return(NULL)
  }
  peakset <- peakset[,-c(4,5)]
  peakset[,1] <- as.character(peakset[,1])
  return(peakset)
}

pv.greylist <- function(pv, greylist, allsame=FALSE, 
                        cores=pv$config$cores) {
  
  if(is(greylist,"list")) {
    greylist <- greylist$master
  }
  
  if(!is(greylist,"GRanges")) {
    controls <- pv$class[PV_BAMCONTROL,]
    if(sum(is.na(controls))==length(controls)) {
      stop("No control reads specified.",call.=FALSE)
    } else {
      if(length(unique(controls))==1 && controls[1]=="") {
        stop("No control reads specified.",call.=FALSE)        
      }
    }
    whichcontrols <- !duplicated(controls)
    whichcontrols <- whichcontrols & !is.na(controls)
    whichcontrols <- whichcontrols & controls != ""
    controls      <- controls[whichcontrols]
    controlnames  <- pv$class[PV_CONTROL,whichcontrols]
    if(length(unique(controlnames))==1 && is.na(unique(controlnames))) {
      controlnames <- controls
    } else {
      if(length(unique(controls))==1 && controls[1]=="") {
        controlnames <- controls      
      }
    }
    
    
    if(is(greylist,"BSgenome")){
      ktype <- seqinfo(greylist)
    } else if(is(greylist,"Seqinfo")) {
      ktype <- greylist
    } else {
      
      dba.ktypes <- NULL
      load(system.file("extra/ktypes.rda", package="DiffBind"),envir = environment())
      
      if(greylist==PV_BLACKLIST_HG19) {
        ktype <- dba.ktypes$BSgenome.Hsapiens.UCSC.hg19
      } else if(greylist==PV_BLACKLIST_HG38) {
        ktype <- dba.ktypes$BSgenome.Hsapiens.UCSC.hg38
      }  else if(greylist==PV_BLACKLIST_GRCH38) {
        ktype <- dba.ktypes$BSgenome.Hsapiens.NCBI.GRCh38
      } else if(greylist==PV_BLACKLIST_MM9) {
        ktype <- dba.ktypes$BSgenome.Mmusculus.UCSC.mm9
      } else if(greylist==PV_BLACKLIST_MM10) {
        ktype <- dba.ktypes$BSgenome.Mmusculus.UCSC.mm10
      } else if(greylist==PV_BLACKLIST_CE10) {
        ktype <- dba.ktypes$BSgenome.Celegans.UCSC.ce10
      } else if(greylist==PV_BLACKLIST_CE11) {
        ktype <- dba.ktypes$BSgenome.Celegans.UCSC.ce11
      } else if(greylist==PV_BLACKLIST_DM3) {
        ktype <- dba.ktypes$BSgenome.Dmelanogaster.UCSC.dm3
      } else if(greylist==PV_BLACKLIST_DM6) {
        ktype <- dba.ktypes$BSgenome.Dmelanogaster.UCSC.dm6
      } else {
        stop("Unknown genome ID.",call.=FALSE)
      }
      rm(dba.ktypes)
    }
    
    ## GreyListChIP for each control    
    controllist <- NULL
    for (i in 1:length(controls)) {
      message(sprintf("Building greylist: %s",controlnames[i]))
      greylist <- pv.makeGreylist(ktype,controls[i],cores,pval=pv$config$greylist.pval)
      controllist <- pv.listadd(controllist, greylist)
    }
    ## Merge
    if(length(controllist)==1) {
      greylist <- controllist[[1]]
      controllist <- NULL
    } else {
      names(controllist) <- controlnames
      greylist <- pv.mergeGreyLists(controllist)
    }
    
    gc(verbose=FALSE)
    
  } else {
    controls <- controllist <- NULL 
  }
  
  message(sprintf("Master greylist: %d ranges, %d bases",
                  length(greylist),sum(width(greylist))))   
  
  # check that at least one chr matches
  snames <- as.character(unique(seqnames(greylist)))
  if(sum(pv$chrmap %in% snames)==0) {
    warning('Greylist does not overlap any peak chromosomes!',call.=FALSE)
  }
  
  # Apply greylist to peaksets
  
  if(allsame) {
    totalPeaks <- nrow(pv$peaks[[1]])
  } else {
    totalPeaks <- sum(unlist(lapply(pv$peaks,nrow)))
  }
  
  pv <- pv.applyBlacklist(pv, greylist)
  
  if(allsame) {
    totalRemoved   <- totalPeaks - nrow(pv$peaks[[1]])
  } else {
    totalRemoved   <- totalPeaks - sum(unlist(lapply(pv$peaks,nrow)))
  }
  msgstring <- sprintf("Removed: %d of %d intervals.",
                       totalRemoved, totalPeaks)
  message(msgstring)
  
  if(length(controls)<=1 || is.null(controllist)) {
    pv$greylist <- greylist
  } else {
    pv$greylist  <- list(master=greylist,controls=controllist)
    pv$greylist$controls <- GRangesList(pv$greylist$controls)
  }
  
  return(pv)
}

pv.makeGreylist <- function(ktype,bamfile,parallel,pval=.999){
  gl <- new("GreyList",karyotype=ktype)
  gl <- GreyListChIP::countReads(gl, bamfile)
  usecores <- 1
  if(!is.null(parallel)) {
    if(parallel != FALSE) {
      usecores <- parallel
    }
  }
  if(is.null(pval)) {
    pval=.999
  }
  gl <- GreyListChIP::calcThreshold(gl,p=pval,cores=usecores)
  gl <- GreyListChIP::makeGreyList(gl)
  
  res <- gl@regions
  rm(gl)
  return(res)
}

pv.mergeGreyLists <- function(controllist) {
  
  greylist <- controllist[[1]]
  
  message(sprintf("%s: %d ranges, %d bases", names(controllist)[1],
                  length(greylist),sum(width(greylist))))    
  
  for(i in 2:length(controllist)) {
    gl <- controllist[[i]]
    message(sprintf("%s: %d ranges, %d bases",
                    names(controllist)[i],length(gl), sum(width(gl))))
    greylist <- union(greylist,gl)
  }
  
  return(greylist)
}

pv.removeBlacklistedPeaks <- function(DBA) {
  bl <- data.frame(DBA$peaks.blacklisted[[1]])[,1:3]
  bl[,1] <- match(bl[,1],DBA$chrmap)
  toremove    <- which(GRanges(data.frame(DBA$merged)) %over% GRanges(bl))
  DBA$merged  <- DBA$merged[-toremove,]
  DBA$binding <- DBA$binding[-toremove,]  
  
  return(DBA)
}





