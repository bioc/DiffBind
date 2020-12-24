PV_BLACKLIST_HG19   <- "BSgenome.Hsapiens.UCSC.hg19"
PV_BLACKLIST_HG38   <- "BSgenome.Hsapiens.UCSC.hg38"
PV_BLACKLIST_GRCH37 <- "BSgenome.Hsapiens.NCBI.GRCh37"
PV_BLACKLIST_GRCH38 <- "BSgenome.Hsapiens.NCBI.GRCh38"
PV_BLACKLIST_MM9    <- "BSgenome.Mmusculus.UCSC.mm9"
PV_BLACKLIST_MM10   <- "BSgenome.Mmusculus.UCSC.mm10"
PV_BLACKLIST_CE10   <- "BSgenome.Celegans.UCSC.ce10"
PV_BLACKLIST_CE11   <- "BSgenome.Celegans.UCSC.ce11"
PV_BLACKLIST_DM3    <- "BSgenome.Dmelanogaster.UCSC.dm3"
PV_BLACKLIST_DM6    <- "BSgenome.Dmelanogaster.UCSC.dm6"

pv.BlackGreyList <- function (DBA, blacklist, greylist,
                              Retrieve, cores) {
  
  isConsensus <- pv.isConsensus(DBA) 
  
  if(!missing(Retrieve)) {
    
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
  
  if(isConsensus) {
    originalPeaks    <- nrow(DBA$binding)
  } else {
    originalPeaks     <- sum(unlist(lapply(DBA$peaks,nrow)))
    originalMerged    <- nrow(DBA$merged)
    originalConsensus <- nrow(DBA$binding)
  }
  
  doblacklist <- dogreylist <- getgenome <- FALSE
  genome <- NULL
  if(is(blacklist,"logical")) {
    if(blacklist == TRUE) {
      doblacklist <- TRUE
      getgenome <- TRUE
    } 
  } else {
    genome <- blacklist
    doblacklist <- TRUE
  }
  
  havegreylist <- FALSE
  if(is(greylist,"logical")) {
    if(greylist == TRUE) {
      dogreylist <- TRUE
      getgenome <- TRUE
    } 
  } else if (is(greylist,"list") ||
             is(greylist,"GRanges") || is(greylist,"GRangesList")) {
    havegreylist <- TRUE
    dogreylist   <- TRUE
    if(is.null(genome)) {
      genome <- FALSE
    }
  } else {
    if(is.null(genome)) {
      genome <- greylist
    }
    dogreylist <- TRUE
  }
  
  if(getgenome) {
    if(doblacklist) {
      blacklist <- genome <- pv.genomes(DBA$class["bamRead",],DBA$chrmap)
    } 
    if(dogreylist) {
      if(!havegreylist) {
        if(is.null(genome)) {
          greylist <- genome <- pv.genomes(DBA$class["bamControl",],DBA$chrmap)
        } else {
          greylist <- genome
        }
      }
    }
  } 
  
  if(doblacklist || dogreylist) {
    if(is.null(genome)) {
      message("No genome detected.")
      return(DBA)
    } else {
      if(is(genome,"character")) {
        message("Genome detected: ",strsplit(genome,"BSgenome.")[[1]][2])
      }
    }
  }
  
  if(doblacklist) {
    message("Applying blacklist...")
    DBA <- pv.blacklist(DBA, blacklist=blacklist, isConsensus=isConsensus)
  }
  
  if(dogreylist) {
    DBA <- pv.greylist(DBA, greylist=greylist, isConsensus=isConsensus, cores=cores)
  }
  
  if(isConsensus) {
    endPeaks <- nrow(DBA$peaks[[1]])
  } else {
    endPeaks     <- sum(unlist(lapply(DBA$peaks,nrow)))
  }
  
  if(endPeaks < originalPeaks) {
    
    if(isConsensus) {
      DBA <- pv.removeBlacklistedPeaks(DBA)
      DBA <- pv.reNormalize(DBA)
    } else {
      DBA <- dba(DBA)
    } 
    
    endMerged    <- nrow(DBA$merged)
    endConsensus <- nrow(DBA$binding)
    
    if(!missing(blacklist) && !missing(greylist)) {
      if(isConsensus) {
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

pv.blacklist <- function(pv, blacklist, isConsensus=FALSE) {
  
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
    } else {
      message("No blacklist found for ",blacklist)
      pv$blacklist <- NULL
      return(pv)
    }
  }
  
  # check that at least one chr matches
  snames <- as.character(unique(seqnames(blacklist)))
  if(sum(pv$chrmap %in% snames)==0) {
    warning('Blacklist does not overlap any peak chromosomes!',call.=FALSE)
  }
  
  # Apply blacklist to peaksets
  
  if(isConsensus) {
    totalPeaks <- nrow(pv$peaks[[1]])
  } else {
    totalPeaks <- sum(unlist(lapply(pv$peaks,nrow)))
  }
  
  pv <- pv.applyBlacklist(pv, blacklist)
  
  if(isConsensus) {
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
    suppressWarnings(peakset <- peakset[!peakset %over% blacklist])
  } else {
    suppressWarnings(peakset <- peakset[peakset %over% blacklist])
  }
  peakset <- data.frame(peakset)
  if(nrow(peakset)==0) {
    return(NULL)
  }
  peakset <- peakset[,-c(4,5)]
  peakset[,1] <- as.character(peakset[,1])
  return(peakset)
}

pv.greylist <- function(pv, greylist, isConsensus=FALSE, 
                        cores=pv$config$cores) {
  
  if(is(greylist,"list")) {
    greylist <- greylist$master
  }
  
  if(!is(greylist,"GRanges")) {
    controls <- pv$class[PV_BAMCONTROL,]
    if(pv.noControls(controls)) {
      message("No control reads specified, unable to generate greylist.")
      return(pv)
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
      load(system.file("extra/ktypes.rda", package="DiffBind"),
           envir = environment())
      
      ktype <- match(greylist,names(dba.ktypes))
      if(is.na(ktype)) {
        message("No known BSgenome: ",greylist)
        pv$greylist <- NULL
        return(pv)
      }
      ktype <- dba.ktypes[[ktype]]
      
      rm(dba.ktypes)
    }
    
    ## GreyListChIP for each control    
    controllist <- pv.makeGreylists(pv,ktype,controls,cores,
                                    pval=pv$config$greylist.pval)
    
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
  
  if(isConsensus) {
    totalPeaks <- nrow(pv$peaks[[1]])
  } else {
    totalPeaks <- sum(unlist(lapply(pv$peaks,nrow)))
  }
  
  pv <- pv.applyBlacklist(pv, greylist)
  
  if(isConsensus) {
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

pv.makeGreylists <- function(pv,ktype,bamfiles,parallel,pval=.999){
  
  usecores <- 1
  if(!is.null(parallel)) {
    if(parallel != FALSE) {
      usecores <- parallel
    }
  }
  if(usecores > 1) {
    param <- BiocParallel::MulticoreParam(workers <- usecores)
  } else {
    param <-  BiocParallel::SerialParam()
  }
  
  message("Counting control reads for greylist...")
  res <- tryCatch(
    suppressMessages(
      gllist <- bplapply(bamfiles,
                         pv.countGreylist,
                         pv, ktype, 
                         BPPARAM=param)),
    error=function(x){stop("GreyListChIP error: ",x)}
  )
  
  controllist <- NULL
  for (i in 1:length(gllist)) {
    message(sprintf("Building greylist: %s",bamfiles[i]))
    greylist <- pv.makeGreylist(pv,ktype,gllist[[i]],usecores,
                                pval=pv$config$greylist.pval)
    controllist <- pv.listadd(controllist, greylist)
  }
  return(controllist)
}

pv.countGreylist <- function(bamfile,pv,ktype) {
  gl <- new("GreyList",karyotype=ktype[pv$chrmap,])
  gl <- GreyListChIP::countReads(gl, bamfile)
  return(gl)
}

pv.makeGreylist <- function(pv,ktype,gl,usecores,pval=.999){
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

pv.noControls <- function(controls) {
  if(sum(is.na(controls))==length(controls)) {
    return(TRUE)
  } else {
    if(length(unique(controls))==1 && controls[1]=="") {
      return(TRUE)        
    }
  }
  return(FALSE)
}

pv.genome <- function(bamfile, chrmap=NULL, ref=NULL, dba.ktypes=NULL) {
  
  header <- Rsamtools::scanBamHeader(bamfile)[[1]]$targets
  
  if(is.null(chrmap)) {
    chrmap <- names(header)
  }
  chrs <- match(chrmap,names(header))
  if(is.na(chrs[1])) {
    message('No matching chromosomes found for file: ',bamfile)
    return(NULL)
  }
  chrnames <- names(header)[chrs]
  
  if(is.null(dba.ktypes)) {
    load(system.file("extra/ktypes.rda", package="DiffBind"),
         envir = environment())
    remove.ktypes <- TRUE
  } else {
    remove.ktypes <- FALSE    
  }
  
  if(!is.null(ref)) {
    dba.ktypes <- list(dba.ktypes[[match(ref,names(dba.ktypes))]])
    names(dba.ktypes) <- ref
  }
  
  res <- NULL
  for(knum in 1:length(dba.ktypes)) {
    matches <- mismatches <- 0
    ktype <- dba.ktypes[[knum]]
    matchChrs <- match(chrnames, names(ktype)) 
    if(!is.na(matchChrs[1])) {
      for(chr in chrs) {
        bamlen <- header[chr]
        krec   <- ktype[names(header)[chr],]
        klen <- GenomeInfoDb::seqlengths(krec)
        if(!is.na(klen)) {
          if(bamlen == klen ) {
            matches <- matches + 1
          } else {
            mismatches <- mismatches + 1
          }
        }
      }
    }
    if(matches > 0 && mismatches == 0) {
      res <- names(dba.ktypes)[knum]
    }
  }
  
  if(remove.ktypes) {
    dba.ktypes <- NULL
  }
  
  return(res)
}

pv.genomes <- function(bamfiles, chrmap=NULL) {
  dba.ktypes <- NULL
  
  load(system.file("extra/ktypes.rda", package="DiffBind"),
       envir = environment())
  
  ref <- NULL
  for(bamfile in bamfiles) {
    newref <- pv.genome(bamfile, chrmap, ref, dba.ktypes)
    if(is.null(ref)) {
      ref <- newref
    } else if(is.null(newref)) {
      return(NULL)      
    }  else if (ref != newref) {
      message("bam files match different reference genomes:",ref," ", newref)
      return(NULL)
    }
  }
  
  dba.ktypes <- NULL
  return(ref)
}




