#####################################
## pv_counts.R -- count-dependant  ##
## 20 October 2009                 ##
## 3 February 2011 -- packaged     ##
## Rory Stark                      ##
## Cancer Research UK              ##
#####################################
PV_DEBUG <- FALSE


## pv.counts -- add peaksets with scores based on read counts
PV_RES_RPKM             <- 1
PV_RES_RPKM_FOLD        <- 2
PV_RES_READS            <- 3
PV_RES_READS_FOLD       <- 4
PV_RES_READS_MINUS      <- 5
PV_RES_RPKM_MINUS       <- 20
PV_SCORE_RPKM           <- PV_RES_RPKM
PV_SCORE_RPKM_FOLD      <- PV_RES_RPKM_FOLD
PV_SCORE_RPKM_MINUS     <- PV_RES_RPKM_MINUS
PV_SCORE_READS          <- PV_RES_READS
PV_SCORE_READS_FOLD     <- PV_RES_READS_FOLD
PV_SCORE_READS_MINUS    <- PV_RES_READS_MINUS
PV_SCORE_CONTROL_READS  <- 18
PV_SCORE_TMM_MINUS_FULL       <- 6
PV_SCORE_TMM_MINUS_EFFECTIVE  <- 7
PV_SCORE_TMM_READS_FULL       <- 8
PV_SCORE_TMM_READS_EFFECTIVE  <- 9
PV_SCORE_TMM_MINUS_FULL_CPM       <- 10
PV_SCORE_TMM_MINUS_EFFECTIVE_CPM  <- 11
PV_SCORE_TMM_READS_FULL_CPM       <- 12
PV_SCORE_TMM_READS_EFFECTIVE_CPM  <- 13
PV_SCORE_READS_FULL               <- 14
PV_SCORE_READS_EFFECTIVE          <- 15
PV_SCORE_READS_MINUS_FULL         <- 16
PV_SCORE_READS_MINUS_EFFECTIVE    <- 17
PV_SCORE_SUMMIT                   <- 101
PV_SCORE_SUMMIT_ADJ               <- 102
PV_SCORE_SUMMIT_POS               <- 103
PV_SCORE_NORMALIZED               <- 104

PV_READS_DEFAULT   <- 0
PV_READS_BAM       <- 3
PV_READS_BED       <- 1

PV_DEFAULT_FILTER <- 1

pv.counts <- function(pv,peaks,minOverlap=2,defaultScore=PV_SCORE_NORMALIZED,
                      bLog=TRUE,insertLength=0,
                      bOnlyCounts=TRUE,bCalledMasks=TRUE,
                      minMaxval= PV_DEFAULT_FILTER,filterFun=max,
                      bParallel=FALSE,bUseLast=FALSE,bWithoutDupes=FALSE, 
                      bScaleControl=FALSE, bSignal2Noise=TRUE,
                      bLowMem=FALSE, readFormat=PV_READS_DEFAULT, 
                      summits, minMappingQuality=0,
                      maxGap=-1, bRecentered=FALSE, minCount=0,
                      bSubControl=FALSE) {
  
  pv <- pv.check(pv)
  
  if(sum(is.na(pv$class[PV_BAMREADS,]))) {
    stop("Can't count: some peaksets are not associated with a .bam file.",call.=FALSE)
  }
  
  pv$class[PV_BAMCONTROL,pv$class[PV_BAMCONTROL,]==""]=NA
  
  if(minOverlap >0 && minOverlap <1) {
    minOverlap <- ceiling(length(pv$peaks) * minOverlap)	
  }
  
  bRecenter <- FALSE
  saveLowMem <- bLowMem
  if(is.logical(summits) && summits==FALSE) {
    summits <- -1
  } 
  if(summits != -1) {
    if(bLowMem==TRUE) {
      bLowMem <- bRecentered
      if(bRecentered == TRUE) {
        #message("Summit info will be lost when using summarizeOverlaps.")
        summits <- -1
      } else {
        message("Computing summits...")
      }
    }
    if(is.logical(summits) && summits==TRUE) {
      summits=0
    }
    if (summits>0) {
      bRecenter=TRUE
    } 
  }
  
  bed <- NULL
  called <- NULL
  if(!missing(peaks)) { # peaks provided
    if(is.vector(peaks)) {
      if(is.character(peaks)){ # peaks in file
        tmp <- pv.peakset(NULL,peaks)
        pv$chrmap <- tmp$chrmap
        peaks <- tmp$peaks[[1]]
        if(is.character(peaks[1,1])){
          peaks[,1] <- factor(peaks[,1],pv$chrmap)
        }
      } else { # peaks based on mask
        saveCorPlot <- pv$config$bCorPlot
        pv$config$bCorPlot <- FALSE
        tmp <- dba(pv,mask=peaks,minOverlap=minOverlap)
        pv$config$bCorPlot <- saveCorPlot 
        pv$chrmap <- tmp$chrmap
        bed <- tmp$binding[,1:3]
        called <- tmp$called[pv.overlaps(tmp,minOverlap),]
      }
    } else { # peaks provided
      pv$chrmap <- unique(as.character(peaks[,1]))
      rownames(peaks) <- 1:nrow(peaks)
      if(is.character(peaks[1,1])){
        peaks[,1] <- factor(peaks[,1],pv$chrmap)
      }
      if(nrow(peaks)==nrow(pv$called)) {
        called <- pv$called
      }
      peaks <- pv.peaksort(peaks)
    }
    if(is.null(bed)) {
      colnames(peaks)[1:3] <- c("CHR","START","END")
      bed <- mergePeaks(peaks[,1:3],maxGap)
    }
  } else { # no peaks provided, derive consensus
    overlaps <- pv.overlaps(pv,minOverlap)
    if(minOverlap == pv$minOverlap) {
      bed <- pv$binding[,1:3]
    } else {
      bed <- pv$merged[overlaps,]
    }
    called <- pv$called[overlaps,]
  }
  
  bed <- pv.check1(bed)
  called <- pv.check1(called)
  
  if(nrow(bed)==0) {
    stop("Zero peaks to count!",call.=FALSE)
  }
  
  bed <- as.data.frame(pv.peaksort(bed,pv$chrmap))
  bed[,1] <- pv$chrmap[bed[,1]]
  
  numChips <- ncol(pv$class)
  chips  <- unique(pv$class[PV_BAMREADS,])
  chips  <- unique(chips[!is.na(chips)])
  inputs <- pv$class[PV_BAMCONTROL,]
  inputs <- unique(inputs[!is.na(inputs)])
  todo   <- unique(c(chips,inputs))
  
  if(!pv.checkExists(todo)) {
    stop('Some read files could not be accessed. See warnings for details.',call.=FALSE)
  }
  
  if(length(insertLength)==1) {
    insertLength <- rep(insertLength,length(todo))
  }
  if(length(insertLength)<length(todo)) {
    warning('Fewer fragment sizes than libraries -- using mean fragment size for missing values',call.=FALSE)
    insertLength <- c(insertLength,rep(mean(insertLength),length(todo)-length(insertLength)))
  }
  if(length(insertLength)>length(todo)) {
    warning('More fragment sizes than libraries',call.=FALSE)
  }
  
  todorecs <- NULL
  for(i in 1:length(todo)) {
    newrec =NULL
    newrec$bamfile <- todo[i]
    newrec$insert <- insertLength[1]
    todorecs <- pv.listadd(todorecs,newrec)
  }
  
  yieldSize <- 5000000
  mode      <- "IntersectionNotEmpty"
  singleEnd <- TRUE
  
  scanbamparam <- NULL
  addfuns <- NULL
  if(bLowMem){
    
    requireNamespace("Rsamtools",quietly=TRUE)
    
    addfuns <- c("BamFileList","summarizeOverlaps","ScanBamParam","scanBamFlag","countBam","SummarizedExperiment")   
    if (insertLength[1] !=0) {
      #warning("fragmentSize ignored when bUseSummarizeOverlaps is TRUE in dba.count",call.=FALSE)
    }
    bAllBam <- TRUE
    # for(st in todo) {
    #   if(substr(st,nchar(st)-3,nchar(st)) != ".bam")	{
    #     bAllBam <- FALSE
    #     warning(st,": not a .bam",call.=FALSE)	
    #   } else if(file.access(paste(st,".bai",sep=""))==-1) {
    #     bAllBam <- FALSE
    #     warning(st,": no associated .bam.bai index",call.=FALSE)	
    #   }
    # }
    # if(!bAllBam) {
    #   stop('All files must be BAM (.bam) with associated .bam.bai index when UseSummarizeOverlaps is TRUE in dba.count',call.=FALSE)	
    # }
    if(!is.null(pv$config$yieldSize)) {
      yieldSize <- pv$config$yieldSize	
    }
    if(!is.null(pv$config$intersectMode)) {
      mode <- pv$config$intersectMode	
    }
    
    if(is.null(pv$config$singleEnd)) {
      bfile <- pv.BamFile(todo[1], bIndex=TRUE)
      pv$config$singleEnd	<- !suppressMessages(
        Rsamtools::testPairedEndBam(bfile))
    } 
    singleEnd <-  pv$config$singleEnd
    if(singleEnd) {
      message("Reads will be counted as Single-end.")
    } else {
      message("Reads will be counted as Paired-end.")        
    }
    
    if(!is.null(pv$config$fragments)) {
      fragments <- pv$config$fragments   
    } else fragments <- FALSE
    
    scanbamparam <- pv$config$scanbamparam 	
  }
  
  if(!bUseLast) {
    pv <- dba.parallel(pv)
    if((pv$config$parallelPackage>0) && bParallel) {   	     
      params  <- dba.parallel.params(pv$config,c("pv.do_getCounts","pv.getCounts",
                                                 "pv.bamReads","pv.BAMstats",
                                                 "fdebug",addfuns))            
      results <- dba.parallel.lapply(pv$config,params,todorecs,
                                     pv.do_getCounts,bed,bWithoutDupes=bWithoutDupes,
                                     bLowMem,yieldSize,mode,singleEnd,
                                     scanbamparam,readFormat,
                                     summits,fragments,minMappingQuality,minCount)
    } else {
      results <- NULL
      for(job in todorecs) {
        message('Sample: ',job,' ')
        results <- pv.listadd(results,pv.do_getCounts(job,bed,bWithoutDupes=bWithoutDupes,
                                                      bLowMem,yieldSize,mode,singleEnd,
                                                      scanbamparam,readFormat,
                                                      summits,fragments,
                                                      minMappingQuality,minCount))
      }	
    }
    if(PV_DEBUG){
      #save(results,file='dba_last_result.RData')
    }
  } else {
    if(PV_DEBUG) {
      load('dba_last_result.RData')
    } else {
      warning("Can't load last result: debug off")
    }
  }
  
  pv.gc()
  
#  if(minMaxval>0) {
    redoScore <- defaultScore
    if(bSubControl) {
      defaultScore <- PV_SCORE_RPKM_MINUS
    } else {
      defaultScore <- PV_SCORE_RPKM
    }
#  } else redoScore <- 0
  if(defaultScore == redoScore) redoScore <- 0
  
  errors <- vapply(results,function(x) if(is.list(x)) return(FALSE) else return(TRUE),TRUE)
  if(sum(errors)) {
    errors <- which(errors)
    for(err in errors) {
      if(is(results[[err]],"try-error")) {
        warning(strsplit(results[[err]][1],'\n')[[1]][2],call.=FALSE)   
      } else {
        warning(results[[err]],call.=FALSE)
      }
    }
    stop("Error processing one or more read files. Check warnings().",call.=FALSE)
  }
  
  allchips <- unique(pv$class[c(PV_BAMREADS,PV_BAMCONTROL),])
  numAdded <- 0
  for(chipnum in 1:numChips) {
    if (pv.nodup(pv,chipnum)) {
      jnum <- which(todo %in% pv$class[PV_BAMREADS,chipnum])
      cond <- results[[jnum]]
      if(length(cond$counts)==0){
        warning('ERROR IN PROCESSING ',todo[jnum],call.=FALSE)
      }
      if(length(cond$libsize)==0){
        warning('ERROR IN PROCESSING ',todo[jnum],call.=FALSE)
      }         
      if(!is.na(pv$class[PV_BAMCONTROL,chipnum])) {
        cnum <- which(todo %in% pv$class[PV_BAMCONTROL,chipnum])
        cont <- results[[cnum]]
        if(length(cont$counts)==0){
          warning('ERROR IN PROCESSING ',todo[cnum],call.=FALSE)
        }
        if(bScaleControl==TRUE) {
          if(cond$libsize>0) {
            scale <- cond$libsize / cont$libsize
            if(scale > 1) scale <- 1
            if(scale != 0) {
              cont$counts <- ceiling(cont$counts * scale)
            }	
          }   	        
        }
      } else {
        cont <- NULL
        cont$counts <- rep(minCount,length(cond$counts))	
        cont$rpkm   <- rep(minCount,length(cond$rpkm))   
      }
      
      cond$counts[cond$counts<minCount] <- minCount
      rpkm_fold   <- cond$rpkm   / cont$rpkm
      reads_fold  <- cond$counts / cont$counts
      rpkm_minus  <- cond$rpkm   - cont$rpkm
      rpkm_minus[rpkm_minus<minCount] <- minCount
      reads_minus <- cond$counts - cont$counts
      reads_minus[reads_minus<minCount] <- minCount
      
      if(bLog) {
        rpkm_fold  <- log2(rpkm_fold)
        reads_fold <- log2(reads_fold)
      }
      if(defaultScore == PV_RES_RPKM) {
        scores <- cond$rpkm
      } else if (defaultScore == PV_RES_RPKM_FOLD) {
        scores <- rpkm_fold
      } else if (defaultScore == PV_RES_RPKM_MINUS) {
        scores <- rpkm_minus
      } else if (defaultScore == PV_RES_READS) {
        scores <- cond$counts    
      } else if (defaultScore == PV_RES_READS_FOLD) {
        scores <- reads_fold
      } else if (defaultScore == PV_RES_READS_MINUS) {
        scores <- reads_minus
      }
      
      if (summits != -1) {
        res <- cbind(bed,scores,cond$rpkm,cond$counts,cont$rpkm,cont$counts,cond$summits,cond$heights)
        colnames(res) <- c("Chr","Start","End","Score","RPKM","Reads","cRPKM","cReads","Summits","Heights")
      } else {
        res <- cbind(bed,scores,cond$rpkm,cond$counts,cont$rpkm,cont$counts)
        colnames(res) <- c("Chr","Start","End","Score","RPKM","Reads","cRPKM","cReads")
      }
      pv <- pv.peakset(pv,
                       peaks       = res,
                       sampID      = pv$class[PV_ID,chipnum],
                       tissue      = pv$class[PV_TISSUE,chipnum],
                       factor      = pv$class[PV_FACTOR,chipnum],
                       condition   = pv$class[PV_CONDITION,chipnum],
                       treatment   = pv$class[PV_TREATMENT,chipnum],
                       consensus   = TRUE,
                       peak.caller = 'counts',
                       control     = pv$class[PV_CONTROL,chipnum],
                       reads       = cond$libsize, #pv$class[PV_READS,chipnum],
                       replicate   = pv$class[PV_REPLICATE,chipnum],
                       readBam     = pv$class[PV_BAMREADS,chipnum],
                       controlBam  = pv$class[PV_BAMCONTROL,chipnum],
                       spikein     = pv$class[PV_SPIKEIN,chipnum],
                       scoreCol    = 0,
                       bRemoveM = FALSE, bRemoveRandom=FALSE,bMakeMasks=FALSE)
      numAdded <- numAdded + 1
    }                  
  }
  
  pv.gc()
  
  if(redoScore > 0) {
    pv$score <- redoScore
  } else {
    pv$score <- defaultScore
  }
  
  pv$maxFilter <- minMaxval
  pv$filterFun <- filterFun
  pv$bSubControl <- bSubControl
  
  if(bOnlyCounts) {
    numpeaks <- length(pv$peaks)
    if(bRecenter) {
      res <- pv.Recenter(pv,summits,(numpeaks-numAdded+1):numpeaks,called)
      if(redoScore>0) {
        defaultScore <- redoScore
        redoScore <- 0
      }
      res <- pv.counts(res,peaks=res$merged,defaultScore=defaultScore,bLog=bLog,insertLength=insertLength,
                       bOnlyCounts=TRUE,bCalledMasks=TRUE,minMaxval=minMaxval,filterFun=filterFun,
                       bParallel=bParallel,bWithoutDupes=bWithoutDupes,bScaleControl=bScaleControl,
                       bSignal2Noise=bSignal2Noise,bLowMem=saveLowMem,readFormat=readFormat,summits=0,
                       bRecentered=TRUE,minMappingQuality=minMappingQuality)
      pv.gc()
      return(res)
    } else {
      savecalled <- pv$called
      if(ncol(savecalled) == numAdded) {
        pv$called <- NULL
      }
      res <- pv.vectors(pv,(numpeaks-numAdded+1):numpeaks,minOverlap=1,bAllSame=TRUE)
      if(is.null(res$called)) {
        if(nrow(savecalled)==nrow(res$peaks[[length(res$peaks)]])) {
          res$called <- savecalled
        } else if(!is.null(called)) {
          if (nrow(called)==nrow(res$peaks[[length(res$peaks)]])) {
            res$called <- called
          }
        }
      }
      if(bRecentered) {
        called <- pv$called
      } 
    }   
    if(!missing(minMaxval)) {
      data <- pv.check1(res$binding[,4:ncol(res$binding)])
      maxs <- apply(data,1,filterFun)
      tokeep <- maxs>=minMaxval
      if(sum(tokeep)<length(tokeep)) {
        if(sum(tokeep)>1) {
          res$binding <- pv.check1(res$binding[tokeep,])
          rownames(res$binding) <- 1:sum(tokeep)
          for(i in 1:length(res$peaks)) {
            res$peaks[[i]] <- res$peaks[[i]][tokeep,]
            rownames(res$peaks[[i]]) <- 1:sum(tokeep)
          } 
          res <- pv.vectors(res,minOverlap=1,bAllSame=pv.allSame(res))
        } else {
          stop('No sites have activity greater than minMaxval',call.=FALSE)
        }
      }
    } else {
      minMaxval <- PV_DEFAULT_FILTER
    }
    if(redoScore > 0) {
      res <- pv.setScore(res,redoScore,minMaxval=minMaxval,filterFun=filterFun,
                         bSignal2Noise=bSignal2Noise)
    }
    if(is.null(res$called)) {
      res$called <- called
    }
  } else {
    if(redoScore > 0) {
      res <- pv.setScore(res,redoScore,minMaxval=minMaxval,filterFun=filterFun,
                         bSignal2Noise=bSignal2Noise)	
    } 
    res <- pv.vectors(pv,bAllSame=pv.allSame(pv))   
  }
  
  if(bSignal2Noise) {
    res$SN <- pv.Signal2Noise(res)
  }
  
  res$minCount <- minCount
  
  pv.gc()
  return(res)	
}

pv.nodup <- function(pv,chipnum) {
  
  
  if(is.null(pv$class[PV_BAMREADS,chipnum])){
    return(FALSE)
  }
  
  if(is.na(pv$class[PV_BAMREADS,chipnum])){
    return(FALSE)
  }   
  
  if(chipnum == 1) {
    return(TRUE)
  }
  
  chips <- pv$class[PV_BAMREADS,1:(chipnum-1)] == pv$class[PV_BAMREADS,chipnum]
  conts <- pv$class[PV_BAMCONTROL,1:(chipnum-1)] == pv$class[PV_BAMCONTROL,chipnum]
  
  conts[is.na(conts)] <- TRUE
  
  if(sum(chips&conts)>0) {
    return(FALSE)
  } else {
    return(TRUE)
  }
  
}

pv.checkExists <- function(filelist){
  res <- file.access(filelist,mode=4)
  for(i in 1:length(filelist)) {
    if(res[i]==-1) {
      warning(filelist[i]," not accessible",call.=FALSE)	
    }	
  }
  return(sum(res)==0)
}

pv.do_getCounts <- function(countrec,intervals,bWithoutDupes=FALSE,
                            bLowMem=FALSE,yieldSize,mode,singleEnd,scanbamparam,
                            fileType=0,summits,fragments,minMappingQuality=0,
                            minCount=0) {
  
  res <- pv.getCounts(bamfile=countrec$bamfile,intervals=intervals,insertLength=countrec$insert,
                      bWithoutDupes=bWithoutDupes,
                      bLowMem=bLowMem,yieldSize=yieldSize,mode=mode,singleEnd=singleEnd,
                      scanbamparam=scanbamparam,
                      fileType=fileType,summits=summits,fragments=fragments,
                      minMappingQuality=minMappingQuality,minCount=minCount)
  pv.gc()
  return(res)
  
}
pv.getCounts <- function(bamfile,intervals,insertLength=0,bWithoutDupes=FALSE,
                         bLowMem=FALSE,yieldSize,mode,singleEnd,scanbamparam,
                         fileType=0,summits=-1,fragments,minMappingQuality=0,
                         minCount=0) {
  
  bufferSize <- 1e6
  fdebug(sprintf('pv.getCounts: ENTER %s',bamfile))
  
  if(bLowMem) {
    if(minMappingQuality>0) {
      #warning('minMappingQuality ignored for summarizeOverlaps, set in ScanBamParam.',call.=FALSE)
    }
    res <- pv.getCountsLowMem(bamfile,intervals,bWithoutDupes,mode,yieldSize,
                              singleEnd,fragments,
                              scanbamparam,minCount=minCount,minQC=minMappingQuality)
    return(res)
  }
  
  fdebug("Starting croi_count_reads...")
  result <- cpp_count_reads(bamfile,insertLength,fileType,bufferSize,
                            intervals,bWithoutDupes,summits,minMappingQuality,
                            minVal=minCount)
  fdebug("Done croi_count_reads...")
  fdebug(sprintf("Counted %d reads...",result$libsize))
  return(result)
}

pv.filterRate <- function(pv,vFilter,filterFun=max) {
  if(!is.numeric(vFilter)) {
    stop('Filter value must be a numeric vector to retrieve filter rate',call.=FALSE)	
  }
  maxs <- apply(pv$binding[,4:ncol(pv$binding)],1,filterFun)
  res <- NULL
  for(filter in vFilter) {
    tokeep <- maxs >= filter
    res <- c(res,sum(tokeep))	
  }
  return(res)
}

pv.getCountsLowMem <- function(bamfile,intervals,bWithoutDups=FALSE,
                               mode="IntersectionNotEmpty",yieldSize=5000000,
                               singleEnd=TRUE,fragments=FALSE,params=NULL,
                               minCount=0, minQC=0) {
  
  intervals <- pv.peaks2DataType(intervals,DBA_DATA_GRANGES)
  
  bf  <- pv.BamFile(bamfile)
  bfl <- Rsamtools::BamFileList(bf,yieldSize=yieldSize)
  
  if(is.null(params)) {
    if(bWithoutDups==FALSE) {
      Dups <- NA
    } else {
      Dups <- FALSE   
    }
    params  <- 
      Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isDuplicate=Dups),
                              mapqFilter=minQC)
  }
  
  counts  <- assay(
    GenomicAlignments::summarizeOverlaps(features=intervals,reads=bfl,
                                         ignore.strand=TRUE,singleEnd=singleEnd,
                                         fragments=fragments,param=params))
  counts[counts<minCount] <- minCount
  libsize <-
    Rsamtools::countBam(bfl,param=params)$records
  if(singleEnd == FALSE) {
    libsize <- libsize/2
  }
  rpkm    <- (counts/(width(intervals)/1000))/(libsize/1e+06)
  
  return(list(counts=counts,rpkm=rpkm,libsize=libsize))
}

pv.BamFile <- function(bamfile,bIndex=TRUE) {
  bai <- NULL
  trybai <- paste(substr(bamfile, 1,nchar(bamfile)-4),".bai", sep="")
  if(file.access(trybai,mode=4) > -1) {
    bai <- trybai
  } else {
    trybai <- paste(bamfile,".bai",sep="") 
    if(file.access(trybai,mode=4) > -1) {
      bai <- trybai
    } 
  }
  if(is.null(bai)) {
    if(!bIndex) {
      stop(bamfile, " has no associated .bai.", call.=FALSE)
    } else {
      message("Indexing ",bamfile)
      bai <- Rsamtools::indexBam(bamfile)
    }
  }
  bamfileobj <- Rsamtools::BamFile(bamfile, index=bai)
  
  return(bamfileobj)
}

pv.Recenter <- function(pv,summits,peakrange,called=NULL) {
  peaklist <- pv$peaks[peakrange]
  if(is.null(peaklist[[1]]$Summits)) {
    stop('Summits not available; re-run dba.count with summits=TRUE',call.=FALSE)   
  }
  
  message('Re-centering peaks...')
  
  positions <- sapply(peaklist,function(x)x$Summits)
  heights   <- sapply(peaklist,function(x) pmax.int(1,x$Heights))
  
  if(!is.null(called)) {
    called <- split(called,rep(1:ncol(called),each=nrow(called)))
    heights <- heights * sapply(called,function(x)x)
  }
  
  centers <- sapply(1:nrow(positions),
                    function(x)round(weighted.mean(positions[x,],heights[x,])))
  starts  <- centers-summits
  ends    <- centers+summits
  
  bed <- peaklist[[1]][,1:3]
  bed[,2] <- starts
  bed[,3] <- ends
  
  pv$peaks <- peaklist
  pv$class <- pv$class[,peakrange]
  if(!is.null(called)) {
    peaklist <- lapply(called,function(x)cbind(bed,x))
    res <- pv.merge(bed,peaklist,pv$class)
    pv$merged <- res$merged[,1:3]
    pv$called <- res$merged[,4:ncol(res$merged)]
    pv$chrmap <- res$chrmap
  } else {
    res <- pv.merge(bed)
    pv$merged <- res$merged
    pv$chrmap <- res$chrmap
  }
  pv$merged <- data.frame(pv$merged)
  pv$merged[,1] <- pv$chrmap[pv$merged[,1]]
  pv$binding <- pv$merged
  return(pv)
}

makepv <- function(peakset,called) {
  peaklist <- lapply(called,function(x)cbind(peakset,x))
}

pv.controlID <- function(samples,i,class, curnum){
  makeID <- FALSE
  if(is.null(samples$ControlID[i])) {
    makeID <- TRUE
  } else if(is.na(samples$ControlID[i])) {
    makeID <- TRUE
  } else {
    return(as.character(samples$ControlID[i]))
  }
  newid <- NULL
  if(makeID) {
    if(!is.null(samples$bamControl[i])) {
      if(!is.na(samples$bamControl[i])) {
        if(!samples$bamControl[i]=="") {
          if(i==1) {
            newid <- 1
          } else {
            res <- samples$bamReads %in% samples$bamControl[i]
            if(sum(res)) {
              return(samples$sampID[which(res)[1]])
            }
            res <- samples$bamControl[1:(i-1)] %in% samples$bamControl[i]
            if(sum(res)) {
              newid <- class[PV_CONTROL,which(res)[1]]
            } else {
              return(curnum)
            }
          }
        }
      }
    }  
  }
  
  if(!is.null(newid)) {
    res <- newid
  } else {
    res <- ""
  }
  return(res)
}

pv.resetCounts <- function(pv,counts,minCount=0) {
  
  if(!is(counts,"data.frame")) {
    stop("New counts must be passed as a data.frame.",call.=FALSE)
  }
  if(sum(pv$class[PV_CALLER,]=="counts") != ncol(pv$class)) {
    stop("All peaks must have counts to replace.")
  }
  if(nrow(pv$binding) != nrow(counts)) {
    stop("All samples must have same number of peaks as binding matrix.",call.=FALSE)
  }
  if(ncol(pv$binding) != ncol(counts)) {
    stop("All samples must have same number of samples as binding matrix.",call.=FALSE)
  }
  if(sum(pv$class[PV_ID,] == colnames(counts[,4:ncol(counts)])) != 
     ncol(pv$class)) {
    stop("All samples must have same IDs, and be in same order, as binding matrix.",call.=FALSE)
  }
  
  for(sample in 4:ncol(counts)) {
    snum <- sample - 3
    pv$peaks[[snum]]$Reads <- sapply(round(counts[,sample]),
                                     function(x){max(minCount,x)})
  }
  if(!is.null(pv$score)) {
    scoreVal <- pv$score
    pv$score <- NULL
    pv <- dba.count(pv,peaks=NULL,score=scoreVal)
  }
  return(pv)
} 

pv.get_reads <- function(pv,peaksets,bSubControl=FALSE,numReads){
  
  if(is.null(bSubControl)) {
    bSubControl <- FALSE
  }
  if(missing(peaksets)) {
    peaksets <- rep(TRUE,length(pv$peaks))
  }
  reads <- NULL
  if(!is.null(pv$peaks_alt)) {
    peaklist <- pv$peaks_alt
  } else {
    peaklist <- pv$peaks
  }
  
  if(is.logical(peaksets)) {
    peaksets <- which(peaksets)
  }
  
  if(!missing(numReads)) {
    numReads <- min(length(peaklist[[1]]$Reads),numReads)
    for(peakset in peaksets) {
      reads <- cbind(reads,as.integer(peaklist[[peakset]]$Reads[1:numReads]))
      if(bSubControl) {
        reads[,ncol(reads)] <- reads[,ncol(reads)] - 
          as.integer(peaklist[[peakset]]$cReads[1:numReads])
      }
    }
  } else {
    for(peakset in peaksets) {
      reads <- cbind(reads,as.integer(peaklist[[peakset]]$Reads))
      if(bSubControl) {
        reads[,ncol(reads)] <- reads[,ncol(reads)] - 
          as.integer(peaklist[[peakset]]$cReads)
      }
    }
  }
  
  if(!is.null(pv$minCount)) {
    reads[reads<pv$minCount] <- pv$minCount
  } else {
    reads[reads<0] <- 0    
  }
  
  rownames(reads) <- 1:nrow(reads)
  
  return(reads)
}

pv.setScore <- function(pv,score,bLog=FALSE,minMaxval=0,rescore=FALSE,
                        filterFun=max,bSignal2Noise=TRUE) {
  
  doscore <- TRUE
  
  if(rescore == TRUE) {
    if(!is.null(pv$score)) {
      if(pv$score == score) {
        if(!is.null(pv$maxFilter)) {
          if(pv$maxFilter == minMaxval) {
            return(pv)	
          }	
        } 
        doscore <- FALSE
      }	
    }
  }
  
  if(doscore) {
    
    pv <- pv.doSetScore(pv, score, noSub=!rescore)
    
    # if(rescore && minMaxval > 0) {
    #   pv <- pv.doFilter(pv, minMaxval, filterFun, bSignal2Noise)
    #   pv$score <- NULL
    #   pv <- pv.doSetScore(pv, score)
    # }
    
  }
  
  pv$score <- score
  pv$maxFilter <- minMaxval
  pv$filterFun <- filterFun
  # pv$bSubControl <- bSubControl
  
  return(pv)
  
}

pv.doSetScore <- function(pv,score,bLog=FALSE,rescore=TRUE,
                          bSignal2Noise=TRUE, noSub=FALSE) {  
  
  if ((score >= DBA_SCORE_TMM_MINUS_FULL) && (score <= DBA_SCORE_TMM_READS_EFFECTIVE_CPM) ) {
    bCPM=FALSE
    if(score == DBA_SCORE_TMM_MINUS_FULL || score == DBA_SCORE_TMM_MINUS_FULL_CPM) {
      bMinus   <- TRUE
      bFullLib <- TRUE
      if(score == DBA_SCORE_TMM_MINUS_FULL_CPM) {
        bCPM=TRUE
      }
    }
    if(score == DBA_SCORE_TMM_MINUS_EFFECTIVE || score == DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM) {
      bMinus   <- TRUE
      bFullLib <- FALSE
      if(score == DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM) {
        bCPM=TRUE
      }
    }
    if(score == DBA_SCORE_TMM_READS_FULL || score == DBA_SCORE_TMM_READS_FULL_CPM) {
      bMinus   <- FALSE
      bFullLib <- TRUE	
      if(score == DBA_SCORE_TMM_READS_FULL_CPM) {
        bCPM=TRUE
      }
    }
    if(score == DBA_SCORE_TMM_READS_EFFECTIVE || score == DBA_SCORE_TMM_READS_EFFECTIVE_CPM) {
      bMinus   <- FALSE
      bFullLib <- FALSE	
      if(score == DBA_SCORE_TMM_READS_EFFECTIVE_CPM) {
        bCPM=TRUE
      }
    }
    
    pv$binding[,4:ncol(pv$binding)] <- pv.normTMM(pv,bMinus=bMinus,bFullLib=bFullLib,bCPM=bCPM)
    
    for(i in 1:length(pv$peaks)) {
      colnum <- 3+i
      pv$peaks[[i]]$Score <- pv$binding[,colnum]
    }
    
  } else {
    
    if(score %in% c(PV_SCORE_READS_FULL,PV_SCORE_READS_MINUS_FULL, 
                    PV_SCORE_READS_EFFECTIVE, PV_SCORE_READS_MINUS_EFFECTIVE)) {
      
      pv$binding[,4:ncol(pv$binding)] <- pv.normLibsize(pv,score)
      
      for(i in 1:length(pv$peaks)) {
        colnum <- 3+i
        pv$peaks[[i]]$Score <- pv$binding[,colnum]
      }
      
    } else if(score == PV_SCORE_NORMALIZED) {
      
      pv$binding[,4:ncol(pv$binding)] <- 
        pv.countsMA(pv, method=DBA_DESEQ2, contrast=NULL, 
                    bNormalized=TRUE, bCountsOnly=TRUE, filter=0,
                    noSub=noSub)
      
      for(i in 1:length(pv$peaks)) {
        colnum <- 3+i
        pv$peaks[[i]]$Score <- pv$binding[,colnum]
      }
      
    } else {
      
      for(i in 1:length(pv$peaks)) {
        colnum <- 3+i
        if(score == PV_SCORE_RPKM) {
          pv$binding[,colnum] <- pv$peaks[[i]]$RPKM	
        }   		
        if(score == PV_SCORE_RPKM) {
          pv$binding[,colnum] <- pv$peaks[[i]]$RPKM	
        } else if(score == PV_SCORE_RPKM_FOLD) {
          numer <- pv$peaks[[i]]$RPKM
          numer[numer==0] <- 1
          denom <- pv$peaks[[i]]$cRPKM
          denom[denom==0] <- 1
          pv$binding[,colnum] <- numer/denom
          if(bLog) {
            pv$binding[,colnum] <- log2(pv$binding[,colnum])	
          }
        } else if(score == PV_SCORE_RPKM_MINUS) {
          pv$binding[,colnum] <- pv$peaks[[i]]$RPKM-pv$peaks[[i]]$cRPKM	
        } else if(score == PV_SCORE_READS) {
          pv$binding[,colnum] <- pv$peaks[[i]]$Reads	
        }  else if(score == PV_SCORE_CONTROL_READS) {
          pv$binding[,colnum] <- pv$peaks[[i]]$cReads	
        } else if(score == PV_SCORE_READS_FOLD) {
          numer <- pv$peaks[[i]]$Reads
          numer[numer==0] <- 1
          denom <- pv$peaks[[i]]$cReads
          denom[denom==0] <- 1
          pv$binding[,colnum] <- numer/denom
          if(bLog) {
            pv$binding[,colnum] <- log2(pv$binding[,colnum])	
          }
        } else if(score == PV_SCORE_READS_MINUS) {
          pv$binding[,colnum] <- pv$peaks[[i]]$Reads-pv$peaks[[i]]$cReads	
        } else if(score == PV_SCORE_SUMMIT || score == PV_SCORE_SUMMIT_ADJ) {
          if(is.null(pv$peaks[[i]]$Heights)) {
            warning('DBA_SCORE_SUMMIT not available; re-run dba.count with summits=TRUE',call.=FALSE)   
          } else {
            pv$binding[,colnum] <- pv$peaks[[i]]$Heights
            if (score == PV_SCORE_SUMMIT_ADJ) {
              pv$binding[,colnum] <- pv$binding[,colnum] * pv.normFactors(pv)[i]   
            }
          }
        } else if(score == PV_SCORE_SUMMIT_POS) {
          if(is.null(pv$peaks[[i]]$Summits)) {
            warning('DBA_SCORE_SUMMIT_POS not available; re-run dba.count with summits=TRUE',call.=FALSE)   
          } else {
            pv$binding[,colnum] <- pv$peaks[[i]]$Summits
          }
        }
        checkscores <- pv$binding[,colnum]
        if(is.null(pv$minCount)) {
          minRead <- 0
        } else {
          minRead <- pv$minCount
        }
        pv$binding[checkscores<minRead,colnum] <- minRead
        pv$peaks[[i]]$Score <- pv$binding[,colnum]
      }
    }
  }
  pv$score <- score
  
  if(bSignal2Noise) {
    pv$SN <- pv.Signal2Noise(pv)
  }
  
  return(pv)
}

pv.doFilter <- function(pv,minMaxval, filterFun, bSignal2Noise=TRUE) {
  
  if(!missing(minMaxval)) {
    maxs <- apply(pv$binding[,4:ncol(pv$binding)],1,filterFun)
    maxs[is.na(maxs)] <- 0
    tokeep <- maxs>=minMaxval
    if(sum(tokeep)>1 && (sum(tokeep) < length(tokeep)) ) {
      pv$binding <- pv$binding[tokeep,]
      rownames(pv$binding) <- 1:sum(tokeep)
      for(i in 1:length(pv$peaks)) {
        pv$peaks[[i]] <- pv$peaks[[i]][tokeep,]
        rownames(pv$peaks[[i]]) <- 1:sum(tokeep)
      }
      
      pv <- pv.vectors(pv,minOverlap=1,bAllSame=TRUE)
      
      if(!is.null(pv$contrasts)) {
        for(i in 1:length(pv$contrasts)) {
          pv$contrasts[[i]]$edgeR=NULL
          pv$contrasts[[i]]$DESeq=NULL
        }
      }
    } else {
      if(sum(tokeep)<2) {
        stop('No sites have activity greater than filter value',call.=FALSE)
      }
    }
    pv$maxFilter <- minMaxval
  }
  
  if(bSignal2Noise) {
    pv$SN <- pv.Signal2Noise(pv)
  }
  
  return(pv)
}

pv.countsMA <- function(pv, method, contrast, 
                        bNormalized=FALSE, bCountsOnly=FALSE,
                        filter=0, filterFun=max, noSub=FALSE) {
  
  #message("pv.countsMA")
  if(is.null(contrast)) {
    contrast<-NULL
    contrast$group1 <- rep(FALSE,ncol(pv$class))
    contrast$group1[1] <- TRUE
    contrast$group2 <- !contrast$group1
  }
  
  if(bNormalized) {
    
    design <- pv$design
    if(is.null(design)) {
      pv$design <- "~1"
    }
    
    savescore <- pv$score
    pv$score <- DBA_SCORE_READS
    
    nofilter <- FALSE
    if(is.null(filter)) {
      filter <- pv$maxFilter
      nofilter <- TRUE
    }
    
    if(noSub) {
      dosub <- FALSE
    } else {
      dosub <- is.null(pv$greylist)
    }
    
    if(method == DBA_DESEQ2) {
      
      if(is.null(pv$norm$DESeq2)) {
        pv <- dba.normalize(pv, method=DBA_DESEQ2, background=FALSE,
                            bSubControl=dosub)
      }
      
      if(is.null(pv$DESeq2$DEdata)) {
        
        if(nofilter)  {
          filter <- pv$norm$DESeq2$filter.val
          filterFun <- pv$norm$DESeq2$filter.fun
        }
        
        dedata <- pv.DEinitDESeq2(pv,bSubControl=pv$norm$DESeq2$bSubControl,
                                  bFullLibrarySize=pv$norm$DESeq2$lib.sizes,
                                  filter=filter,filterFun=filterFun) 
      } else {
        dedata <- pv$DESeq2$DEdata
      }
      counts <- DESeq2::counts(dedata, normalized = bNormalized)
      
    } else if(method == DBA_EDGER) {
      
      if(is.null(pv$norm$edgeR)) {
        pv <- dba.normalize(pv, method=DBA_EDGER, background=FALSE,
                            bSubControl=dosub)
      }
      
      if(nofilter)  {
        filter <- pv$norm$edgeR$filter.val
        filterFun <- pv$norm$edgeR$filter.fun
      }
      
      if(is.null(pv$edgeR$DEdata)) {
        dedata <- pv.DEinitedgeR(pv,bSubControl=pv$norm$edgeR$bSubControl,
                                 filter=filter, filterFun=filterFun)
      } else {
        dedata <-pv$norm$edgeR 
      }
      counts <- pv.edgeRCounts(pv,method,bNormalized)
    }
    counts <- cbind(counts[,contrast$group1], counts[,contrast$group2])
    pv$design <- design
    pv$score <- savescore
  } else {
    counts <- pv.DEinit(pv,contrast$group1, contrast$group2, bRawCounts=TRUE)
  }
  
  if(bCountsOnly) {
    return(counts)
  }
  
  num1 <- sum(contrast$group1)
  num2 <- sum(contrast$group2)
  samps <- num1 + num2
  
  conc <- log2(apply(counts,1,mean))
  
  if(num1>1) {
    con1 <- log2(apply(counts[,1:num1],1,mean))
  } else {
    con1 <- log2(counts[,1])
  }
  
  if(num2>1) {
    con2 <- log2(apply(counts[,(num1+1):ncol(counts)],1,mean))
  } else {
    con2 <- log2(counts[,num1+1])
  }
  
  con1[con1<0] <- 0
  con2[con2<0] <- 0
  conc[conc<0] <- 0
  
  fold <- con1 - con2
  
  res <- NULL
  res$Conc <- conc
  res$Fold <- con1 - con2
  res$Con1 <- con1
  res$Con2 <- con2
  
  return(res)
}

