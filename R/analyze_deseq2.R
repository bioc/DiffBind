pv.DESeq2 <- function(pv,group1,group2,label1="Group 1",label2="Group 2",
                      bSubControl=TRUE,bFullLibrarySize=FALSE,bTagwise=TRUE,
                      bGLM=TRUE,blockList=NULL,filter=0, filterFun=max){
  
  if (!requireNamespace("DESeq2",quietly=TRUE)) {
    stop("Package DESeq2 not installed",call.=FALSE)
  }
  res <- NULL
  
  targets <- NULL
  if(!is.null(blockList)) {
    targets <- pv.blockFactors(pv,group1,group2,label1,label2,blockList)
    if(is.null(targets)){
      return(res)	
    }
  } 
  res$DEdata <- pv.DEinit(pv,group1,group2,label1,label2,method='DESeq2',
                          bSubControl=bSubControl,
                          bFullLibrarySize=bFullLibrarySize,targets=targets,
                          filter=filter, filterFun=filterFun)
  res$counts <- DESeq2::counts(res$DEdata)
  
  if (!is.null(blockList)) {
    message('DESeq2 multi-factor analysis')
    attr <-  blockList[[1]]$attribute
    if(attr=='Replicate') {
      DESeq2::design(res$DEdata) <- formula(~Replicate + group)
    } else if(attr=='Tissue') {
      DESeq2::design(res$DEdata) <- formula(~Tissue + group)    
    } else if(attr=='Factor') {
      DESeq2::design(res$DEdata) <- formula(~Factor + group)
    } else if(attr=='Condition') {
      DESeq2::design(res$DEdata) <- formula(~Condition + group)
    } else if(attr=='Caller') {
      DESeq2::design(res$DEdata) <- formula(~Caller + group)
    } else if(attr=='Treatment') {
      DESeq2::design(res$DEdata) <- formula(~Treatment + group)
    } else if(attr=='Block') {
      DESeq2::design(res$DEdata) <- formula(~Block + group)
    } else {
      warning('Unsupported blocking attribute: ',attr,call.=FALSE)
      return(NULL)  
    }
    fdebug(sprintf('pv.DESeq blocking analysis: %d db (%s/%s)',
                   sum(res$de$padj<0.05),label1,label2))
  } 
  
  # if(!is.null(pv$norm$DESeq2)) {
  #   message("DESeq2: Using size factors from dba.normalize()")
  #   DESeq2::sizeFactors(res$DEdata) <- pv$norm$DESeq2$norm.facs
  # } else {
  #   message("DESeq2: NO size factors from dba.normalize()")
  if(!bFullLibrarySize) {
    res$DEdata <- DESeq2::estimateSizeFactors(res$DEdata)
  }
  # }
  
  res$facs   <- DESeq2::sizeFactors(res$DEdata)
  res$DEdata <- DESeq2::estimateDispersions(res$DEdata,fitType='local')
  res$DEdata <- DESeq2::nbinomWaldTest(res$DEdata)
  res$de     <- DESeq2::results(res$DEdata)
  
  res$de$pvalue[is.na(res$de$pvalue)]=1
  res$de$padj[is.na(res$de$padj)]=1
  
  res$de <- res$de[order(res$de$padj),]
  res$de <- cbind(as.numeric(rownames(res$de)),res$de$pvalue,res$de$padj)
  colnames(res$de) <- c("id","pval","padj")
  rownames(res$de) <- res$de[,1]
  res$de <- data.frame(res$de)
  
  res$bSubControl      <- bSubControl
  res$bFullLibrarySize <- bFullLibrarySize
  
  return(res)
  
}

pv.DESeq2design <- function(pv,
                            bSubControl=TRUE,bFullLibrarySize=TRUE,
                            existing=NULL,
                            filter=0, filterFun=max){
  
  if (!requireNamespace("DESeq2",quietly=TRUE)) {
    stop("Package DESeq2 not installed",call.=FALSE)
  }
  res <- NULL
  res$design <- pv$design
  
  res$DEdata <- pv.DEinitDESeq2(pv,
                                bSubControl=bSubControl,
                                bFullLibrarySize=bFullLibrarySize,
                                filter=filter, filterFun=filterFun)
  
  res$facs   <- DESeq2::sizeFactors(res$DEdata)
  
  #Can we avoid re fitting model?
  if(!is.null(existing$DEdata)) {
    if(!identical(existing$facs, res$facs)) {
      existing <- NULL
    } else if(bSubControl != existing$bSubControl) {
      existing <- NULL
    } else if(!is.null(pv$config$DESeq2$fitType)) {
      if(pv$config$DESeq2$fitType != existing$fitType) {
        existing <- NULL
      }
    }
  }
  
  if(is.null(existing$DEdata)) {
    if(!is.null(pv$config$DESeq2$fitType)) {
      fittype <- pv$config$DESeq2$fitType
    } else {
      fittype <- 'local'
    }
    res$DEdata <- DESeq2::estimateDispersions(res$DEdata,fitType=fittype)
    res$DEdata <- DESeq2::nbinomWaldTest(res$DEdata)
    res$bSubControl      <- bSubControl
    res$bFullLibrarySize <- bFullLibrarySize
    res$fitType          <- fittype
  } else {
    res <- existing
  }
  
  res$names <- DESeq2::resultsNames(res$DEdata)
  
  return(res)
}

pv.DEinitDESeq2 <- function(pv,
                            bSubControl=FALSE,bFullLibrarySize=FALSE,
                            bRawCounts=FALSE,
                            filter=0,filterFun=max,
                            numReads) {
  
  if (requireNamespace("DESeq2",quietly=TRUE)) {
    DESeq2 <- TRUE 
  } else {
    stop("Package DESeq2 not installed",call.=FALSE)
  }    
  
  srcmask <- pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts")
  if(sum(srcmask)==0) {
    stop('Error: no samples present with read counts',call.=FALSE)
  }
  counts  <- pv.get_reads(pv, srcmask,bSubControl=bSubControl,numReads=numReads)
  
  if(is.null(filter)) {
    filter <- 0
  }
  if(filter > 0){
    scores <- apply(counts,1,filterFun)
    keep   <- scores > filter
    counts <- counts[keep,]
    rownames(counts) <- which(keep)
  } else {
    rownames(counts) <- as.character(1:nrow(counts))
  }
  
  colnames(counts) <- pv$class[PV_ID,]
  
  if(bRawCounts) {
    return(counts)	
  }
  
  if(is(bFullLibrarySize,"logical")) {
    if(bFullLibrarySize) {
      libsize <- as.numeric(pv$class[PV_READS,])
    } else {
      libsize <- colSums(counts)
    }
  } else {
    libsize <- bFullLibrarySize
    bFullLibrarySize <- TRUE
  }
  
  meta <- pv.getMeta(pv)
  if(!pv.isDesign(pv$design)) {
    designstr <- "~ 1"
  } else {
    designstr <- pv$design
  }
  res <- suppressMessages(
    DESeq2::DESeqDataSetFromMatrix(counts,meta,formula(designstr))
  )
  
  if(!is.null(pv$norm$DESeq2)) {
    # message("DESeq2: Using size factors from dba.normalize()")
    if(!is.null(pv$norm$DESeq2$norm.facs)) {
      DESeq2::sizeFactors(res) <- pv$norm$DESeq2$norm.facs
    }
    if(pv$norm$DESeq2$norm.method == PV_NORM_OFFSETS ||
       pv$norm$DESeq2$norm.method == PV_NORM_OFFSETS_ADJUST) {
      if(!is.null(pv$norm$offsets$offsets)) {
        # message("DESeq2: use offsets")
        offsets <- assay(pv$norm$offsets$offsets, "offsets")
        if(pv$norm$offsets$offset.method == PV_OFFSETS_LOESS ||
           pv$norm$DESeq2$norm.method == PV_NORM_OFFSETS_ADJUST) {
          offsets <- pv.offsetsAdjust(pv, offsets, res)
        }
        offsets <- offsets[1:nrow(res),]
        DESeq2::normalizationFactors(res) <- offsets
      } else {
        stop('Internal error: no offsets available.',call.=FALSE)
      }
    }
  } else {
    # message("DESeq2: NO size factors from dba.normalize()")
    if(!bFullLibrarySize) {
      res <- DESeq2::estimateSizeFactors(res)
    } else {
      DESeq2::sizeFactors(res) <- libsize/min(libsize)
    }
  }
  
  return(res)
}

pv.allDESeq2 <- function(pv,block,bSubControl=FALSE,bFullLibrarySize=FALSE,
                         bTagwise=TRUE,bGLM=FALSE,bParallel=FALSE,
                         filter=0,filterFun=max) {
  
  if (!requireNamespace("DESeq2",quietly=TRUE)) {
    stop("Package DESeq2 not installed",call.=FALSE)
  }
  
  if(is.null(pv$contrasts)) {
    if(missing(block)) {
      pv <- pv.contrast(pv)
    } else {
      pv <- pv.contrast(pv,block=block)
    }
  }
  
  if(bParallel) {
    pv <- dba.parallel(pv)
    jobs <- NULL
    blocks <- NULL
  }
  
  fdebug('pv.allDESeq2: for each contrast')
  reslist <- NULL
  
  if(bParallel && (pv$config$parallelPackage > 0)) {   
    params <- dba.parallel.params(pv$config,c('pv.DESeq2_parallel','pv.DESeq2'))
    
    reslist  <- dba.parallel.lapply(pv$config,params,1:length(pv$contrasts),
                                    pv.DESeq2_parallel,pv,NULL, 
                                    bSubControl,bFullLibrarySize,
                                    bTagwise=bTagwise,bGLM=bGLM,
                                    filter=filter,filterFun=filterFun)
    blist <- NULL
    for(i in 1:length(pv$contrasts)) {
      if(!is.null(pv$contrasts[[i]]$blocklist)) {
        blist <- c(blist,i)
      }
    }
    if(length(blist > 0)) {
      bres <-  dba.parallel.lapply(pv$config,params,blist,
                                   pv.DESeq2_parallel,pv,TRUE, 
                                   bSubControl,bFullLibrarySize,
                                   bTagwise=bTagwise,bGLM=bGLM,
                                   filter=filter,filterFun=filterFun)
      for(i in 1:length(blist)) {
        fdebug(sprintf('pv.allDESeq2: contrast %d gets bres %d (%d db)',
                       blist[i],i,sum(bres[[i]]$de$padj<0.05)))
        reslist[[blist[i]]]$block <- bres[[i]]
      }    
    }     
  } else { # SERIAL
    for(i in 1:length(pv$contrast)) { 	
      if(!is.null(pv$contrasts[[i]]$contrast)){
        res <- pv.DESeq2Contrast(pv,pv$contrasts[[i]])
      } else {
        res <- pv.DESeq2(pv,pv$contrasts[[i]]$group1,pv$contrasts[[i]]$group2,
                         pv$contrasts[[i]]$name1,pv$contrasts[[i]]$name2,
                         bSubControl=bSubControl,
                         bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise,
                         bGLM=bGLM, filter=filter,filterFun=filterFun)
        
        if(!is.null(pv$contrasts[[i]]$blocklist)) {
          res$block <- pv.DESeq2(pv,pv$contrasts[[i]]$group1,pv$contrasts[[i]]$group2,
                                 pv$contrasts[[i]]$name1,pv$contrasts[[i]]$name2,
                                 bSubControl=bSubControl,
                                 bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise,
                                 blockList=pv$contrasts[[i]]$blocklist,
                                 filter=filter,filterFun=filterFun)   
        }
      }
      reslist <- pv.listadd(reslist,res)   
    }
  }  
  
  return(reslist)
}

pv.DESeq2_parallel <- function(contrast,pv,blockList,bSubControl,
                               bFullLibrarySize,bTagwise=TRUE,bGLM=FALSE,
                               filter=0,filterFun=max) {
  crec <- pv$contrasts[[contrast]]
  if(!is.null(blockList)) {
    blockList <- crec$blocklist
  }
  if(!is.null(crec$contrast)){
    res <- pv.DESeq2Contrast(pv,crec)
  } else {
    res <- pv.DESeq2(pv,crec$group1,crec$group2,crec$name1,crec$name2,
                     bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,
                     bTagwise=bTagwise,bGLM=bGLM,blockList=blockList,
                     filter=filter,filterFun=filterFun)
    pv.gc()
  }
  return(res)
}

pv.DESeq2Contrast <- function(pv, contrast, shrink=TRUE) {
  
  if(is.null(pv$DESeq2)) {
    stop("Can not test contrast: model has not been run",call.=FALSE)
  }
  
  res <- NULL
  if(contrast$contrastType!="byname") {
    res$de     <- DESeq2::results(pv$DESeq2$DEdata, contrast=contrast$contrast)
  } else {
    res$de     <- DESeq2::results(pv$DESeq2$DEdata, name=contrast$contrast)
  }
  
  if(shrink) {
    if(contrast$contrastType=="byname") {
      res$de <- suppressMessages(DESeq2::lfcShrink(pv$DESeq2$DEdata,
                                                   coef=contrast$contrast,
                                                   res=res$de,type="apeglm"))
    } else  if(contrast$contrastType=="results1") {
      res$de <- suppressMessages(DESeq2::lfcShrink(pv$DESeq2$DEdata,
                                                   coef=contrast$contrast[[1]],
                                                   res=res$de,type="apeglm"))
    } else {
      res$de <- suppressMessages(DESeq2::lfcShrink(pv$DESeq2$DEdata,
                                                   contrast=contrast$contrast,
                                                   res=res$de,type="ashr"))
    }
  }
  
  res$de$pvalue[is.na(res$de$pvalue)]=1
  res$de$padj[is.na(res$de$padj)]=1
  res$de <- res$de[order(res$de$padj),]
  res$de <- cbind(as.numeric(rownames(res$de)),
                  res$de$pvalue,res$de$padj,
                  res$de$log2FoldChange)
  colnames(res$de) <- c("id","pval","padj","fold")
  
  rownames(res$de) <- res$de[,1]
  res$de <- data.frame(res$de)
  
  return(res)
  
}

pv.stripDESeq2 <- function(drec) {
  if(!is.null(drec)) {
    drec$counts     <- NULL
    drec$DEdata     <- NULL
    drec$fullFit    <- NULL
    drec$reducedFit <- NULL
    
    if(!is.null(drec$block)) {
      drec$block <- pv.stripDESeq2(drec$block)	
    }	
  }
  return(drec)
}
