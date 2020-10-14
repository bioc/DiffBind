pv.DBA <- function(pv,method='edgeR',bTagwise=TRUE,
                   minMembers=3,bParallel=FALSE, block) {
  
  if(bParallel) {
    setParallel <- TRUE
    bParallel <- FALSE
  } else {
    setParallel <- FALSE	
  }
  
  if(is.null(pv$contrasts)) {
    message("Adding contrasts(s)...")
    if(missing(block)) {
      pv <- pv.contrast(pv,minMembers=minMembers)
    } else {
      pv <- pv.contrast(pv,minMembers=minMembers,block=block)
    }
  }
  
  if(is.null(pv$contrasts)) {
    stop('Unable to perform analysis: no contrasts specified.')	
  }
  
  noreps <- FALSE
  for(contrast in pv$contrasts) {
    if(!is.null(contrast$group1)) {
      if(sum(contrast$group1)<2) {
        noreps <- TRUE
      }
      if(sum(contrast$group2)<2) {
        noreps <- TRUE
      }
    }
  }
  if(noreps) {
    warning("Some groups have no replicates. Results may be unreliable.",call.=FALSE)	
  }
  
  if(bParallel) {
    pv <- dba.parallel(pv)
    jobs <- NULL
    numjobs <- 0
  }
  
  
  results <- NULL
  
  if(pv.isDesign(pv$design)) { ## establish design-based analyses
    
    pv <- pv.designAnalysis (pv, method, tagwise=bTagwise, setParallel)
    
  }
  
  if('edgeRGLM' %in% method) {
    
    if(bParallel && (pv$config$parallelPackage > 0)) {
      numjobs <- numjobs + 1
      params <- dba.parallel.params(pv$config,c('pv.allDEedgeR','pv.DEedgeR','pv.contrast','pv.listadd'))
      fdebug('submit job: pv.all')
      jobs <- pv.listadd(jobs,dba.parallel.addjob(pv$config,params,
                                                  pv.allDEedgeR,pv,
                                                  bFullLibrarySize=(pv$norm$edgeR$lib.method!=DBA_LIBSIZE_PEAKREADS),
                                                  bParallel=TRUE,
                                                  bSubControl=pv$norm$edgeR$bSubControl,
                                                  bTagwise=bTagwise,bGLM=TRUE,
                                                  filter=pv$norm$edgeR$filter.val,filterFun=pv$norm$edgeR$filter.fun))
    } else {
      results <- pv.listadd(results, pv.allDEedgeR(pv,block=block,
                                                   bSubControl=pv$norm$edgeR$bSubControl,
                                                   bFullLibrarySize=(pv$norm$edgeR$lib.method!=DBA_LIBSIZE_PEAKREADS),
                                                   bParallel=setParallel,bTagwise=bTagwise,bGLM=TRUE,
                                                   filter=pv$norm$edgeR$filter.val,filterFun=pv$norm$edgeR$filter.fun))
    }
  }
  
  if('DESeq2' %in% method) {
    if (!requireNamespace("DESeq2",quietly=TRUE)) {
      stop("Package DESeq2 not installed",call.=FALSE)
    }

    if(bParallel && (pv$config$parallelPackage > 0)) {
      numjobs <- numjobs + 1
      params <- dba.parallel.params(pv$config,c('pv.DESeq2'))
      jobs <- pv.listadd(jobs,dba.parallel.addjob(pv$config,params,pv.allDESeq2,pv,
                                                  bSubControl=pv$norm$DESeq2$bSubControl,
                                                  bFullLibrarySize=(pv$norm$DESeq2$lib.method!=DBA_LIBSIZE_PEAKREADS),
                                                  bTagwise=bTagwise,bGLM=FALSE,
                                                  bParallel=TRUE,
                                                  filter=pv$norm$DESeq2$filter.val,filterFun=pv$norm$DESeq2$filter.fun))
    } else {
      results <- pv.listadd(results,pv.allDESeq2(pv,bSubControl=pv$norm$DESeq2$bSubControl,
                                                 bFullLibrarySize=(pv$norm$DESeq2$lib.method!=DBA_LIBSIZE_PEAKREADS),
                                                 bTagwise=bTagwise,bGLM=FALSE,bParallel=setParallel,
                                                 filter=pv$norm$DESeq2$filter.val,filterFun=pv$norm$DESeq2$filter.fun))
    }
  }
  
  
  if(bParallel && (pv$config$parallelPackage > 0)) {
    results <- dba.parallel.wait4jobs(pv$config,jobs)
  }
  
  jnum <- 1
  if('edgeRGLM' %in% method) {
    edger <- results[[jnum]]
    for(i in 1:length(edger)){
      pv$contrasts[[i]]$edgeR <- edger[[i]]
    }
    jnum <- jnum+1
  }
  
  if( ('DESeq2' %in% method) || ('DESeq2GLM' %in% method) ) {
    deseq2 <- results[[jnum]]
    for(i in 1:length(deseq2)){
      pv$contrasts[[i]]$DESeq2 <- deseq2[[i]]
    }
    jnum <- jnum+1
  } 
  
  # pv$filter    <- filter
  # pv$filterFun <- filterFun
  
  fdebug(sprintf('Exit pv.DBA: %f',pv$contrasts[[1]]$edgeR$counts[7,1]))
  return(pv)
}

pv.DEinit <- function(pv,mask1,mask2,group1=1,group2=2,method='edgeR',
                      bSubControl=FALSE,bFullLibrarySize=FALSE,removeComps=0,
                      bRawCounts=FALSE,targets=NULL,
                      filter=0,filterFun=max) {
  
  fdebug('enter pv.DEinit')
  
  edgeR  <- FALSE
  DESeq2 <- FALSE
  if(method == 'edgeR') {
    edgeR <- TRUE   
  } else if (method == 'DESeq2') {
    if (requireNamespace("DESeq2",quietly=TRUE)) {
      DESeq2 <- TRUE 
    } else {
      stop("Package DESeq2 not installed",call.=FALSE)
    }    
  } else {
    warning('Invalid method: ',method,call.=FALSE)
    return(NULL)
  }
  
  srcmask <- pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts")
  
  fdebug(sprintf('pv.DEinit: %s %2.0f %s %2.0f srcmask %2.0f',
                 group1,sum(mask1),group2,sum(mask2),sum(srcmask)))
  g1 <- which(mask1 & srcmask)
  g2 <- which(mask2 & srcmask)
  
  s1 <- pv.get_reads(pv,g1,bSubControl=bSubControl)
  s2 <- pv.get_reads(pv,g2,bSubControl=bSubControl)
  
  counts <- cbind(s1,s2)
  
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
  
  colnames(counts) <- c(pv$class[PV_ID,mask1],pv$class[PV_ID,mask2])
  
  if(bRawCounts) {
    return(counts)	
  }
  
  groups <- factor(c(rep(group1,length(g1)),rep(group2,length(g2))))
  
  
  if(bFullLibrarySize) {
    libsize <- as.numeric(pv$class[PV_READS,c(g1,g2)])
  } else {
    libsize <- colSums(counts)
  }
  
  if(!bRawCounts) {
    
    if(edgeR) {
      
      normfacs <- NULL
      
      res <- edgeR::DGEList(counts,lib.size=libsize,norm.factors=normfacs,
                            group=groups,genes=as.character(1:nrow(counts)))
      rownames(res$counts) <- 1:nrow(res$counts)
      fdebug(sprintf('DGEList counts: %f',res$counts[7,1]))
    }
    
    if(DESeq2) {
      colnames(counts) <- NULL
      if(is.null(targets)) {     	
        res <- DESeq2::DESeqDataSetFromMatrix(counts,data.frame(groups),formula(~ groups))
      } else {
        res <- DESeq2::DESeqDataSetFromMatrix(counts,data.frame(targets),formula(~ group))     		
      }
      
      # if(!is.null(pv$norm$DESeq2)) {
      #   message("DESeq2: Using size factors from dba.normalize()")
      #   DESeq2::sizeFactors(res) <- pv$norm$DESeq2$norm.facs
      # } else {
      #   message("DESeq2: NO size factors from dba.normalize()")
      
      if(bFullLibrarySize) {
        DESeq2::sizeFactors(res) <- libsize/min(libsize)
      }  else {
        res <- DESeq2::estimateSizeFactors(res)
      }
      # }
      
    }
  }                
  return(res)
}

pv.blockFactors <- function(pv,group1,group2,label1,label2,blockList) {
  samples <- group1 | group2
  targets <- NULL
  for(block in blockList) {
    samps <- block$samples & samples &  
      (pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts"))
    IDs   <- pv$class[PV_ID,samps]
    groups <- NULL
    wsamps <- which(samps)
    for(sw in wsamps) {
      if(group1[sw]) {
        groups <- c(groups,label1)
      } else {
        groups <- c(groups,label2)
      }   
    }
    block <- rep(block$label, sum(samps))
    block <- cbind(groups,block)
    rownames(block) <- IDs
    targets <- rbind(targets,block)
  }
  
  snames <- c(colnames(pv$class[,group1]),colnames(pv$class[,group2]))
  if(length(unique(snames))!=sum(group1|group2)){
    warning('Error: all samples must have unique IDs for blocking analysis',call.=FALSE)
    return(NULL)	
  }
  tnames <- rownames(targets)
  #if(length(snames)!=length(tnames)){
  #   warning('Error: all samples must be matched for blocking analysis')
  #   return(res)	
  #}  	 
  
  newt <- targets
  for(i in 1:nrow(targets)) {
    idx <- match(tnames[i],snames)
    newt[idx,] <- targets[i,]	
  }
  targets <- newt
  rownames(targets) <- snames
  
  colnames(targets) <- c("group",blockList[[1]]$attribute)
  
  targets <- data.frame(targets)
  targets[[1]] <- factor(targets[[1]])
  targets[[2]] <- factor(targets[[2]])
  
  return(targets)
  
}

pv.normFactors <- function(pv,bMinus=FALSE,bFullLib=TRUE) {
  if(length(pv$peaks)<2) {
    warning('Unable to TMM normalize -- not enough peaksets',call.=FALSE)
    return(pv)   
  }
  g1     <- rep(FALSE,length(pv$peaks))
  g1[1]  <- TRUE
  
  savenames <- pv$class[PV_ID,]
  pv$class[PV_ID,] <- 1:ncol(pv$class)
  res    <- pv.DEedgeR(pv,g1,!g1,"1","2",bSubControl=bMinus,bFullLibrarySize=bFullLib,bNormOnly=TRUE)
  return(res$samples$norm.factors)
}

pv.stripDBA <- function(conrec) {
  if(!is.null(conrec$edgeR)) {
    conrec$edgeR  <- pv.stripEdgeR(conrec$edgeR)
  }
  
  if(!is.null(conrec$DESeq2)) {
    conrec$DESeq2 <- pv.stripDESeq2(conrec$DESeq2)
  }
  
  return(conrec)
}


pv.normLibsize <- function(pv, score) {
  
  if (score == PV_SCORE_READS_MINUS_FULL ||
      score == PV_SCORE_READS_MINUS_EFFECTIVE) {
    minus <- TRUE
  } else {
    minus <- FALSE
  }
  
  counts <- matrix(0, nrow(pv$peaks[[1]]), length(pv$peaks))
  for (i in 1:length(pv$peaks)) {
    if (minus) {
      counts[, i] <- pv$peaks[[i]]$Reads - pv$peaks[[i]]$cReads
    } else {
      counts[, i] <- pv$peaks[[i]]$Reads
    }
  }
  
  if (score == PV_SCORE_READS_MINUS_FULL ||
      score == PV_SCORE_READS_FULL) {
    libsize <- as.numeric(pv$class[PV_READS, ])
  } else {
    libsize <- colSums(counts)
  }
  
  normfacs <- libsize / min(libsize)
  counts <- t(t(counts) / normfacs)
  
  return(counts)
}

pv.designAnalysis <- function(pv, method, tagwise=TRUE, bParallel=TRUE) {
  
  reslist <- NULL
  if(bParallel && (pv$config$parallelPackage > 0)) {
    pv <- dba.parallel(pv)
    reslist <- dba.parallel.lapply(pv$config,NULL,method,
                                   pv.designAnalysis_parallel,pv,
                                   tagwise)
  } else {
    reslist <- NULL
    for(meth in method) {
      reslist <- pv.listadd(reslist,pv.designAnalysis_parallel(meth,pv,tagwise))
    }
  }
  
  for(i in 1:length(method)) {
    if(method[i] == DBA_EDGER) {
      pv$edgeR <- reslist[[i]]
    } else if(method[i] == DBA_DESEQ2) {
      pv$DESeq2 <- reslist[[i]]
    }
  }
  
  return(pv)
}

pv.designAnalysis_parallel <- function(method,pv,tagwise) {
  
  if(method == DBA_EDGER) {
    res <- pv.edgeRdesign(pv,bSubControl=pv$norm$edgeR$bSubControl,
                          bTagwise=tagwise,
                          existing=pv$edgeR,
                          filter=pv$norm$edgeR$filter.val,
                          filterFun=pv$norm$edgeR$filter.fun)
  } else if(method == DBA_DESEQ2) {
    res <- pv.DESeq2design(pv,bSubControl=pv$norm$DESeq2$bSubControl,
                           existing=pv$DESeq2,
                           filter=pv$norm$DESeq2$filter.val,
                           filterFun=pv$norm$DESeq2$filter.fun)
  }
  return(res)
}

pv.getMethod <- function(str) {   
  if (str == "DBA_EDGER" || str == DBA_EDGER) {
    ret=DBA_EDGER
  } else if (str == "DBA_EDGER_GLM" || str == DBA_EDGER_GLM) {
    ret=DBA_EDGER_GLM  
  } else  if (str == "DBA_DESEQ2" || str == DBA_DESEQ2) {
    ret=DBA_DESEQ2 
  } else ret <- NULL
  return(ret)
}

