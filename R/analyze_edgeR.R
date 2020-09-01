pv.DEedgeR <- function(pv,group1,group2,label1="Group 1",label2="Group 2",blockList=NULL,
                       bSubControl=FALSE,bFullLibrarySize=FALSE,bTagwise=TRUE,bGLM=TRUE,bNormOnly=FALSE,
                       bWeighting=FALSE, filter=0,filterFun=max) {
  
  fdebug('Enter pv.DEedgeR')
  
  #require(edgeR)
  
  res <- pv.DEinit(pv,group1,group2,label1,label2,method='edgeR',
                   bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,
                   filter=filter, filterFun=filterFun)
  
  # if(is.null(pv$norm$edgeR))  {
  # message("edgeR: re-calc norm.factors")
  res <- edgeR::calcNormFactors(res,method="TMM",doWeighting=bWeighting)
  fdebug(sprintf('calcNormFactors: %f',res$counts[7,1]))
  # }  else {
  #   message("edgeR: DONT re-calc norm.factors (use  dba.normalize())")
  # }
  
  if(bNormOnly) {
    return(res)	
  }
  
  if(is.null(blockList)) {
    fdebug('blockList is NULL')
  } else {
    fdebug('blockList is not NULL')
  }
  fdebug('pv.DEedgeR: check for blocking factor')
  if(is.null(blockList)) {
    fdebug('pv.DEedgeR: NO blocking factor')
    res <- edgeR::estimateCommonDisp(res)
    fdebug(sprintf('estimateCommonDisp: %f',res$counts[7,1]))
    if(bGLM){
      res$design <- model.matrix(~res$samples$group)
      if(bTagwise) {
        res <- edgeR::estimateGLMCommonDisp(res,res$design)
        res <- edgeR::estimateGLMTagwiseDisp(res,res$design)
      } else {
        res <- edgeR::estimateGLMCommonDisp(res,res$design)
      }
      res$GLM <- edgeR::glmFit(res,res$design)
      res$LRT <- edgeR::glmLRT(res$GLM,2)
    } else {
      if(bTagwise){
        res <- edgeR::estimateTagwiseDisp(res,prior.df=50,trend="none")
        #res <- edgeR::estimateTagwiseDisp(res,prior.n=getPriorN(res),trend="movingave")
        res$db     <- edgeR::exactTest(res,dispersion='tagwise')
      } else {
        res$db     <- edgeR::exactTest(res,dispersion='common')	
      }
    }
    fdebug(sprintf('Fit and test: %f',res$counts[7,1]))
    fdebug('pv.DEedgeR: estimateTagwiseDisp complete')
    
    fdebug(sprintf('pv.DEedgeR: exactTest complete:%s-%s',res$db$comparison[1],res$db$comparison[2]))
    #res$db$fdr <- edgeR::topTags(res$db,nrow(res$db$counts))
  } else {
    fdebug('pv.DEedgeR: BLOCKING FACTOR')
    
    targets <- pv.blockFactors(pv,group1,group2,label1,label2,blockList)
    if(is.null(targets)){
      return(res)	
    }
    res$samples <- data.frame(cbind(targets,res$samples[,2:3]))
    
    attr <-  blockList[[1]]$attribute
    if(attr=='Replicate') {
      res$designmatrix <- model.matrix(~ Replicate + group,data = targets)
    } else if(attr=='Tissue') {
      res$designmatrix <- model.matrix(~ Tissue + group,data = targets)
    } else if(attr=='Factor') {
      res$designmatrix <- model.matrix(~ Factor + group,data = targets)
    } else if(attr=='Condition') {
      res$designmatrix <- model.matrix(~ Condition + group,data = targets)
    } else if(attr=='Caller') {
      res$designmatrix <- model.matrix(~ Caller + group,data = targets)
    } else if(attr=='Treatment') {
      res$designmatrix <- model.matrix(~ Treatment + group,data = targets)
    } else if(attr=='Block') {
      res$designmatrix <- model.matrix(~ Block + group,data = targets)
    } else {
      warning('Unsupported blocking attribute: ',attr,call.=FALSE)
      return(NULL)	
    }
    message('edgeR multi-factor analysis.')
    
    # if(is.null(pv$norm$edgeR))  {
    # message("edgeR: re-calc norm.factors")
    # res <- edgeR::calcNormFactors(res,method="TMM",doWeighting=bWeighting)
    # fdebug(sprintf('calcNormFactors: %f',res$counts[7,1]))
    # }  else {
    #   message("edgeR: DONT re-calc norm.factors (use  dba.normalize())")
    # }
    
    res <- edgeR::estimateGLMCommonDisp(res,res$designmatrix)
    if(bTagwise) {
      res <- edgeR::estimateGLMTagwiseDisp(res,res$designmatrix)
    }
    res$GLM <- edgeR::glmFit(res,res$designmatrix)
    res$LRT <- edgeR::glmLRT(res$GLM,ncol(res$designmatrix))
    res$counts <- NULL	 
    #res$fdr <- edgeR::topTags(res$LRT,nrow(res$counts))
  }
  
  res$bSubControl      <- bSubControl
  res$bFullLibrarySize <- bFullLibrarySize
  
  fdebug(sprintf('Exit pv.DEedgeR: %f',res$counts[7,1]))
  return(res)	
  
}

pv.DEedgeR_parallel <- function(contrast,pv,blockList,bSubControl,
                                bFullLibrarySize,bTagwise,bGLM,
                                filter=0,filterFun=max) {
  crec <- pv$contrasts[[contrast]]
  if(!is.null(blockList)) {
    blockList <- crec$blocklist
  }
  if(!is.null(crec$contrast)){
    res <- pv.edgeRContrast(pv,crec)
  } else {
    res <- pv.DEedgeR(pv,crec$group1,crec$group2,crec$name1,crec$name2,
                      blockList=blockList,
                      bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,
                      bTagwise=bTagwise,bGLM=bGLM,
                      filter=filter,filterFun=filterFun)
    
    fdebug(sprintf('Exit pv.DEedgeR_parallel: %f',res$counts[7,1]))
    pv.gc()
  }
  return(res)
}


pv.allDEedgeR <- function(pv,block,bFullLibrarySize=FALSE,bParallel=FALSE,bSubControl=FALSE,
                          bTagwise=TRUE,bGLM=FALSE,filter=filter,filterFun=filterFun) {
  
  fdebug('ENTER pv.allDEedgeR')
  #require(edgeR)
  
  if(is.null(pv$contrasts)) {
    if(missing(block)) {
      pv$contrasts <- pv.contrast(pv)
    } else {
      pv$contrasts <- pv.contrast(pv,block=block)
    }
  }
  
  if(bParallel) {
    pv <- dba.parallel(pv)
    jobs <- NULL
    blocks <- NULL
  }
  
  fdebug('pv.allDEedgeR: for each contrast')
  reslist <- NULL
  
  if(bParallel && (pv$config$parallelPackage > 0)) {
    params <- dba.parallel.params(pv$config,c('pv.DEedgeR_parallel','pv.DEedgeR',
                                              'pv.DEinit','calcNormFactors',
                                              'edgeR::estimateCommonDisp',
                                              'edgeR::estimateTagwiseDisp',
                                              'edgeR::estimateGLMCommonDisp',
                                              'edgeR::estimateGLMTagwiseDisp',
                                              'edgeR::exactTest','edgeR::topTags',
                                              'edgeR::glmFit','edgeR::glmLRT'))
    reslist <- dba.parallel.lapply(pv$config,params,1:length(pv$contrasts),pv.DEedgeR_parallel,pv,
                                   NULL,bSubControl,bFullLibrarySize,bTagwise,bGLM=bGLM,
                                   filter=filter,filterFun=filterFun)
    
    fdebug(sprintf('Return from parallel call to pv.DEedgeR_parallel: %f',reslist[[1]]$counts[7,1]))
    
    blist <- NULL
    for(i in 1:length(pv$contrasts)) {
      if(!is.null(pv$contrasts[[i]]$blocklist)) {
        blist <- c(blist,i)
      }
    }
    if(length(blist > 0)) {
      bres <-  dba.parallel.lapply(pv$config,params,blist,pv.DEedgeR_parallel,pv,
                                   TRUE,bSubControl,bFullLibrarySize,bTagwise,
                                   filter=filter,filterFun=filterFun)
      for(i in 1:length(blist)) {
        reslist[[blist[i]]]$block <- bres[[i]]
      }    
    }     
  } else { #SERIAL
    for(i in 1:length(pv$contrast)) { 	
      if(!is.null(pv$contrasts[[i]]$contrast)){
        res <- pv.edgeRContrast(pv,pv$contrasts[[i]])
      } else {
        res <- pv.DEedgeR(pv,pv$contrasts[[i]]$group1,pv$contrasts[[i]]$group2,
                          pv$contrasts[[i]]$name1,pv$contrasts[[i]]$name2,
                          bSubControl=bSubControl,
                          bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise,
                          bGLM=bGLM,filter=filter,filterFun=filterFun)
        
        if(!is.null(pv$contrasts[[i]]$blocklist)) {
          res$block <- pv.DEedgeR(pv,pv$contrasts[[i]]$group1,pv$contrasts[[i]]$group2,
                                  pv$contrasts[[i]]$name1,pv$contrasts[[i]]$name2,
                                  pv$contrasts[[i]]$blocklist,
                                  bSubControl=bSubControl,
                                  bFullLibrarySize=bFullLibrarySize,
                                  bTagwise=bTagwise,filter=filter,filterFun=filterFun)   
        }
      }
      reslist <- pv.listadd(reslist,res)   
    }
  }
  
  fdebug(sprintf('Exit pv.allDEedgeR: %f',reslist[[1]]$counts[7,1]))
  return(reslist)
}

pv.edgeRdesign <- function(pv,
                           bSubControl=TRUE,
                           bTagwise=TRUE, bWeighting=FALSE,
                           existing=NULL,
                           filter=0, filterFun=max){
  
  if (!requireNamespace("edgeR",quietly=TRUE)) {
    stop("Package edgeR not installed",call.=FALSE)
  }
  
  #Can we avoid re fitting model?
  if(!is.null(existing$DEdata)) {
    if(bTagwise != existing$bTagwise) {
      existing <- NULL
    } else if(bSubControl != existing$bSubControl) {
      existing <- NULL
    } 
    # else if(bFullLibrarySize != existing$bFullLibrarySize) {
    #   existing <- NULL
    # }
  }
  
  res <- existing
  
  if(is.null(existing$DEdata)) {
    
    res$DEdata <- pv.DEinitedgeR(pv,
                                 bSubControl=bSubControl,
                                 # bFullLibrarySize=bFullLibrarySize,
                                 filter=filter, filterFun=filterFun)
    
    if(is.null(pv$norm$edgeR))  {
      # message("edgeR: re-calc norm.factors")
      res$DEdata <- edgeR::calcNormFactors(res$DEdata,method="TMM",
                                           doWeighting=bWeighting)
      fdebug(sprintf('calcNormFactors: %f',res$counts[7,1]))
    }  else {
      # message("edgeR: DONT re-calc norm.factors (use  dba.normalize())")
    }
    
    design <- pv$design
    design.matrix <- pv.edgeRdesign.matrix(pv,design)
    
    res$DEdata <- edgeR::estimateGLMTrendedDisp(res$DEdata,design=design.matrix)
    if(bTagwise) {
      res$DEdata <- edgeR::estimateGLMTagwiseDisp(res$DEdata,design=design.matrix)
    }
    res$DEdata <- edgeR::glmQLFit(res$DEdata,design=design.matrix)
    
    res$bSubControl      <- bSubControl
    # res$bFullLibrarySize <- bFullLibrarySize
    res$bTagwise         <- bTagwise
    res$design           <- design.matrix
    
  } else {
    res <- existing
  }
  
  return(res)
}

pv.edgeRdesign.matrix <- function(pv,design){
  meta <- pv.getMeta(pv)
  design.matrix <- model.matrix(formula(design),data=meta)
  return(design.matrix)
}

pv.DEinitedgeR <- function(pv,
                           bSubControl=FALSE,
                           bFullLibrarySize,
                           bRawCounts=FALSE,
                           filter=0,filterFun=max,
                           numReads) {
  
  if (!requireNamespace("edgeR",quietly=TRUE)) {
    stop("Package edgeR not installed",call.=FALSE)
  }    
  
  srcmask <- pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts")
  if(sum(srcmask)==0) {
    stop('Error: no samples present with read counts',call.=FALSE)
  }
  
  counts  <- pv.get_reads(pv, srcmask, bSubControl=bSubControl,numReads=numReads)
  
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
  
  if(!is.null(pv$norm$edgeR)) {
    # message("edgeR: Using lib.sizes from dba.normalize()")
    libsize <- pv$norm$edgeR$lib.size
  } else {
    if(missing(bFullLibrarySize)) {
      stop('Internal error: bFullLibrarySize missing.')
    }
    # message("edgeR: NO lib.sizes from dba.normalize()")
    if(!is(bFullLibrarySize,"logical")) {
      libsize <- bFullLibrarySize
    } else {
      if(bFullLibrarySize) {
        libsize <- as.numeric(pv$class[PV_READS,srcmask])
      } else {
        libsize <- colSums(counts)
      }
    }
  }
  
  offsets <- NULL
  if(!is.null(pv$norm$edgeR)) {
    # message("edgeR: Using norm.factors from dba.normalize()")
    normfacs <- pv$norm$edgeR$norm.facs
    if(pv$norm$edgeR$norm.method == PV_NORM_OFFSETS) {
      if(!is.null(pv$norm$offsets$offsets)) {
        # message("edgeR: use offsets")
        offsets <- assay(pv$norm$offsets$offsets, "offsets")
      } else {
        stop('Internal error: no offsets available.',call.=FALSE)
      }
    }
  } else {
    # message("edgeR: NO norm.factors from dba.normalize()")
    normfacs <- NULL
  }
  
  meta <- pv.getMeta(pv)
  res <- #suppressMessages(
    edgeR::DGEList(counts,lib.size=libsize, norm.factors=normfacs,
                   samples=meta,genes=as.character(1:nrow(counts)))
  #)
  
  if(!is.null(offsets)) {
    # res <- edgeR::scaleOffset(res,offsets)
    res$offset <- offsets
  }
  
  return(res)
}


pv.edgeRContrast <- function(pv, contrast, shrink=TRUE) {
  
  if(is.null(pv$edgeR)) {
    stop("Can not test contrast: model has not been run",call.=FALSE)
  }
  
  res <-  de <- coeffs <- NULL
  
  coeffs <- pv.getCoeffs(contrast,pv)
  
  if(is.null(coeffs)) {
    stop("edgeR: unsupported contrast.",call.=FALSE)
  }
  
  de <- edgeR::glmQLFTest(pv$edgeR$DEdata,contrast=coeffs)
  
  if(is.null(de)) {
    stop("edgeR: unsupported contrast.",call.=FALSE)
  }
  
  res$de <- edgeR::topTags(de,nrow(de))$table
  
  res$de$PValue[is.na(res$de$PValue)]=1
  res$de$FDR[is.na(res$de$FDR)]=1
  res$de <- res$de[order(res$de$FDR),]
  res$de <- cbind(as.numeric(rownames(res$de)),res$de$PValue,
                  res$de$FDR,res$de$logFC, res$de$logCPM)
  colnames(res$de) <- c("id","pval","padj","fold","logCPM")
  
  rownames(res$de) <- res$de[,1]
  res$de <- data.frame(res$de)
  res$coeffs <- coeffs
  
  return(res)
}

pv.getCoeffs <- function(contrast,pv) {
  
  coef1 <- coef2 <- coeffs <- NULL
  if(contrast$contrastType=="bycolumn") {
    coeffs <- contrast$contrast
  } else if(contrast$contrastType=="simple") {
    coef1 <- pv.getCoef(contrast$contrast, pv)
  } else if(contrast$contrastType=="byname") {
    coef1 <- pv.getCoef(contrast$contrast, pv)
  } else if(is(contrast$contrast,"list")) {
    coef1 <- pv.getCoef(contrast$contrast[[1]], pv)
    if(contrast$contrastType=="byresults2") {
      coef2 <- pv.getCoef(contrast$contrast[[2]], pv)
    }
  } 
  
  if(!is.null(coef1)) {
    if(!is.null(coef2)) {
      coeffs <- coef1 - coef2
    } else {
      coeffs <- coef1
    }
  } 
  
  return(coeffs)
}

pv.getCoef <- function(contrast, pv) {
  
  coefs <- pv$DESeq2$names
  
  if(is(contrast,"vector")) {
    vals <- c(contrast[1],contrast[2],"vs",contrast[3])
    contrast <- paste(contrast[1],contrast[2],"vs",
                      contrast[3],sep="_")
  } else {
    vals  <- strsplit(contrast,"_")[[1]]
  }
  
  res   <- rep(0,length(coefs))
  
  coef  <- pv.matchCoefs(contrast,coefs)
  
  if(is.na(coef)) {
    contrast <- paste(vals[1],vals[4],vals[3],vals[2],sep="_")
    coef <- pv.matchCoefs(contrast,coefs)
    if(!is.na(coef)) {
      res[coef] <- -1
      return(res)
    }
  } else {
    res[coef] <- 1
    return(res)
  }
  
  coefs <- colnames(pv$edgeR$design)
  num <- paste(vals[1],vals[2],sep="")
  den <- paste(vals[1],vals[4],sep="")
  num <- match(num,coefs)
  den <- match(den,coefs)
  
  
  if(is.na(num) && is.na(den)) {
    return(NULL)
  }
  
  if(!is.na(num)) {
    res[num] <- 1
  }
  
  if(!is.na(den)) {
    res[den] <- -1
  }
  
  return(res)
}

pv.matchCoefs <- function(contrast,coefs) {
  
  coef <- match(contrast, coefs)
  if(!is.na(coef)) {
    return(coef)
  }
  
  coef <- match(make.names(contrast), coefs)
  
  return(coef)
  
}

pv.checkIntercept <- function(pv, facname, valname) {
  fac <- match(facname, names(pv$meta))
  if(pv$meta[[fac]][1] ==  valname) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

pv.normTMM <- function(pv,bMinus=TRUE,bFullLib=FALSE,bCPM=FALSE,bReturnFactors=FALSE){
  
  if(length(pv$peaks)<2) {
    warning('Unable to TMM normalize -- not enough peaksets',call.=FALSE)
    return(pv)	
  }
  
  vColors <- pv.colsv
  
  g1     <- rep(FALSE,length(pv$peaks))
  g1[1]  <- TRUE
  
  savenames <- pv$class[PV_ID,]
  pv$class[PV_ID,] <- 1:ncol(pv$class)
  res    <- pv.DEedgeR(pv,g1,!g1,"1","2",bSubControl=bMinus,bFullLibrarySize=bFullLib,bNormOnly=TRUE)
  #res    <- estimateCommonDisp(res)
  counts <- res$counts
  sizes  <- res$samples$lib.size * res$samples$norm.factors
  
  if(bReturnFactors) {
    return(list(libsize=res$samples$lib.size,normfactors=res$samples$norm.factors))
  }
  
  counts <- t(t(counts)/sizes)
  
  if(bCPM) {
    counts <- counts * 1E06     
  } else {
    counts <- counts * mean(res$samples$lib.size)
  }
  
  colnames(counts) <- savenames
  return(counts)
  
}
pv.stripEdgeR <- function(erec) {
  if(!is.null(erec)) {
    
    if(is(erec,"DGEList")) {
      erec$counts     <- NULL
      erec$pseudo.alt <- NULL
      erec$conc       <- NULL
      #erec$genes      <- NULL
      erec$all.zeros  <- NULL
      erec$tagwise.dispersion <- NULL
    }
    
    if(!is.null(erec$GLM)) {
      erec$GLM <- pv.stripEdgeRGLM(erec$GLM)
    }
    if(!is.null(erec$LRT)) {
      erec$LRT <- pv.stripEdgeRLRT(erec$LRT)
    }
    if(!is.null(erec$block)) {
      erec$block <- pv.stripEdgeR(erec$block)
    }
    if(!is.null(erec$db)) {
      erec$db <- pv.stripEdgeR(erec$db)
    }      
  }
  return(erec)
}

pv.stripEdgeRGLM <- function(grec) {
  if(!is.null(grec)) {
    grec$counts       <- NULL
    grec$fitted.values= NULL
    grec$offset       <- NULL
    grec$coefficients <- NULL
    grec$deviance     <- NULL
    grec$df.residual  <- NULL
    grec$abundance    <- NULL
    #grec$genes        <- NULL
  }
  return(grec)
}

pv.stripEdgeRLRT <- function(lrec) {
  if(!is.null(lrec)) {
    lrec <- pv.stripEdgeR(lrec)
    lrec <- pv.stripEdgeRGLM(lrec)
    lrec$GLM <- pv.stripEdgeRGLM(lrec$GLM)
    
    lrec$dispersion.used   <- NULL
  }
  return(lrec)
}      