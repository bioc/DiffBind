pv.DBAreport <- function(pv,contrast=1,method='edgeR',th=0.05,bUsePval=FALSE,
                         bCalled=FALSE,bCounts=FALSE,bCalledDetail=FALSE,
                         file,initString='reports/DBA',bNormalized=TRUE,
                         ext="csv",minFold=0,bSupressWarning=FALSE,
                         bFlip=FALSE, precision=2:3) {
  
  
  if(length(contrast)>1){
    stop('Can only specify one contrast unless requesting a report-based DBA.',
         call.=FALSE)
  }
  
  if(contrast > length(pv$contrasts)) {
    stop('Specified contrast number is greater than number of contrasts',
         call.=FALSE)
    return(NULL)
  }
  con <- pv$contrasts[[contrast]]
  
  if(is.null(con$group1)) {
    groups <- FALSE
  } else {
    groups <- TRUE
  } 
  
  group1 <- con$group1
  name1  <- con$name1
  group2 <- con$group2
  name2  <- con$name2
  facs   <- 1:(sum(group1)+sum(group2))
  if(bFlip) {
    group1 <- con$group2
    name1  <- con$name2
    group2 <- con$group1
    name2  <- con$name1
    facs   <- (sum(group2)+1):length(facs)
    facs   <- c(facs,1:sum(group2))
  }
  
  if(method=='edgeR' || method=='edgeRGLM'){
    if(is.null(con$edgeR) || is(con$edgeR,"try-error")) {
      stop('edgeR analysis has not been run for this contrast',call.=FALSE)
      return(NULL)
    }
    
    filter <- pv$norm$edgeR$filter.val
    filterFun <- pv$norm$edgeR$filter.fun
    
    normalized <- FALSE
    if(!is.null(con$contrast)) {
      
      if(is.null(pv$edgeR$DEdata)) {
        
        message('edgeR analysis missing -- re-running...')
        
        if(is.null(pv$config$edgeR$bTagwise)) {
          bTagwise <- TRUE
        } else {
          bTagwise <- pv$config$edgeR$bTagwise 
        }
        
        edger <- pv.edgeRdesign(pv,bSubControl=pv$norm$edgeR$bSubControl,
                                bTagwise=bTagwise,
                                existing=pv$edgeR,
                                filter=pv$norm$edgeR$filter.val,
                                filterFun=pv$norm$edgeR$filter.fun)
        pv$edgeR$DEdata <- edger$DEdata
      }
      
      counts <- pv.edgeRCounts(pv,method,bNormalized)
      normalized <- TRUE
      
      if(groups) {
        facs <-  c(which(group1),which(group2))
        counts <- counts[,facs]
      } else {
        facs <- 1:ncol(counts)
      }
      
    } else {
      if(is.null(con$edgeR$counts)) {
        counts <- pv.DEinit(pv,group1,group2,name1,name2,
                            method='edgeR',
                            bSubControl=con$edgeR$bSubControl,
                            bFullLibrarySize=con$edgeR$bFullLibrarySize,
                            bRawCounts=TRUE,
                            filter=filter, filterFun=filterFun)
      } else {
        counts <- con$edgeR$counts	
      }
    }
    if("id" %in% names(con$edgeR$de)) {
      siteCol <- 1
      pvCol   <- 2
      fdrCol  <- 3
      data <- con$edgeR$de   
      if(is.null(pv$edgeR$DEdata$pseudo.lib.size)) {
        sizes <-  pv$edgeR$DEdata$samples$norm.factors[facs]
        pls <- 1
      } else {
        sizes <- pv$edgeR$DEdata$samples$lib.size[facs] * 
          pv$edgeR$DEdata$samples$norm.factors[facs]
        pls <- pv$edgeR$DEdata$pseudo.lib.size[facs]
      }
    } else {
      if(!is.null(con$edgeR$LRT)) {
        siteCol <- 1
        pvCol   <- 5
        fdrCol  <- 6
        data <- edgeR::topTags(con$edgeR$LRT,nrow(counts))$table   	  	
      } else {
        siteCol <- 1
        pvCol   <- 4
        fdrCol  <- 5
        data <- edgeR::topTags(con$edgeR$db,nrow(counts))$table
      }
      sizes <- con$edgeR$samples$lib.size[facs] * 
        con$edgeR$samples$norm.factors[facs]
      pls <- con$edgeR$pseudo.lib.size
    }
    if(bNormalized && !normalized){
      counts <- t(t(counts)/sizes)
      counts <- counts * pls
    } 
  } else if(method=='edgeRlm'){
    
    filter <- pv$norm$edgeR$filter.val
    filterFun <- pv$norm$edgeR$filter.fun
    
    if(is.null(con$edgeR$counts)) {
      counts <- pv.DEinit(pv,group1,group2,name1,name2,
                          method='edgeR',
                          bSubControl=con$edgeR$block$bSubControl,
                          bFullLibrarySize=con$edgeR$block$bFullLibrarySize,
                          bRawCounts=TRUE,
                          filter=filter, filterFun=filterFun)
    } else {
      counts <- con$edgeR$counts	
    }
    
    siteCol <- 1
    pvCol   <- 5  
    fdrCol  <- 6
    data <- edgeR::topTags(con$edgeR$block$LRT,nrow(counts))$table
    
    if(bNormalized){
      sizes <- con$edgeR$samples$lib.size[facs] * 
        con$edgeR$samples$norm.factors[facs]
      counts <- t(t(counts)/sizes)
      counts <- counts * con$edgeR$pseudo.lib.size
    } 
  } else if (method=='DESeq2' || method=='DESeq2Block') {
    if (!requireNamespace("DESeq2",quietly=TRUE)) {
      stop("Package DESeq2 not installed",call.=FALSE)
    }
    if(is.null(con$DESeq2) || is(con$DESeq2,"try-error")) {
      stop('DESeq2 analysis has not been run for this contrast',call.=FALSE) 
      return(NULL) 
    }
    
    filter <- pv$norm$DESeq2$filter.val
    filterFun <- pv$norm$DESeq2$filter.fun
    
    siteCol <- 1
    pvCol   <- 2
    fdrCol <-  3
    
    if(!is.null(con$contrast)) {
      
      if(is.null(pv$DESeq2$DEdata)) {
        message('DESeq2 analysis missing -- re-running...')
        deseq2 <- pv.DESeq2design(pv,
                                  bSubControl=pv$norm$DESeq2$bSubControl,
                                  existing=pv$DESeq2,
                                  filter=pv$norm$DESeq2$filter.val,
                                  filterFun=pv$norm$DESeq2$filter.fun)
        pv$DESeq2$DEdata <- deseq2$DEdata
      }
      
      counts <- DESeq2::counts(pv$DESeq2$DEdata, normalized = bNormalized)
      
      data     <- con$DESeq2$de
      normfacs <- pv$DESeq2$facs
      
      if(groups) {
        counts   <- cbind(counts[,group1],counts[,group2])
        normfacs <- c(normfacs[group1],normfacs[group2])
      }
      
    } else {
      # counts <- pv.DEinit(pv,group1,group2,name1,name2,
      #                     method='DESeq2',
      #                     bSubControl=con$DESeq2$bSubControl,
      #                     bFullLibrarySize=con$DESeq2$bFullLibrarySize,
      #                     bRawCounts=TRUE,
      #                     filter=filter, filterFun=filterFun)
      
      counts <- pv.countsMA(pv, method=DBA_DESEQ2, contrast=NULL, 
                            bNormalized=bNormalized, bCountsOnly=TRUE, filter=0,
                            noSub=!pv$norm$DESeq2$bSubControl)
      counts <- cbind(counts[,group1],counts[,group2])
      
      
      if(method=='DESeq2Block') {
        data <- con$DESeq2$block$de
        normfacs <- con$DESeq2$block$facs[facs]
      } else {
        data <- con$DESeq2$de
        normfacs <- con$DESeq2$facs[facs]
      }
      # if(bNormalized){
      #   counts <- t(t(counts)/normfacs)
      # }
    }
  } else {
    stop('Unknown DE method: ',method,call.=FALSE)
    return(NULL)
  }
  
  if(bUsePval) {
    thCol <- pvCol	
  } else {
    thCol <- fdrCol	
  }
  
  x <- which(is.na(data[,pvCol]))
  if(length(x)>0){
    data[x,pvCol] <- 1
  }
  x <- which(is.na(data[,fdrCol]))
  if(length(x)>0){
    data[x,fdrCol] <- 1
  }
  
  keep <-  data[,thCol] <= th
  sites <- as.numeric(data[keep,siteCol])
  sites <- match(sites,as.integer(rownames(counts)))
  if(sum(keep)==0) {
    if(!bSupressWarning) {
      warning('No sites above threshold',call.=FALSE)
    }
    return(NULL)
  } else if(sum(keep)==1) {
    cnames <- colnames(counts)
    counts <- matrix(counts[sites,],nrow=1,ncol=ncol(counts))
    colnames(counts) <- cnames
    rownames(counts) <- sites
  } else {
    counts <- counts[sites,]
  }
  
  siteids <- as.numeric(rownames(counts))[1:length(sites)]
  
  if(groups) {
    if(length(sites)==1) {
      conc <- log2(mean(counts))
      if(sum(group1)>1) {
        con1 <- log2(mean(counts[,1:sum(group1)]))
      } else {
        con1 <- log2(counts[,1])
      }
      if(sum(group2)>1) {
        con2 <- log2(mean(counts[,(sum(group1)+1):ncol(counts)]))
      } else {
        con2 <- log2(counts[,sum(group1)+1])
      }        
    } else {
      conc <- log2(apply(counts,1,mean))
      if(sum(group1)>1) {
        con1 <- log2(apply(counts[,1:sum(group1)],1,mean))
      } else {
        con1 <- log2(counts[,1])
      }
      if(sum(group2)>1) {
        con2 <- log2(apply(counts[,(sum(group1)+1):ncol(counts)],1,mean))
      } else {
        con2 <- log2(counts[,sum(group1)+1])
      }
    }
    con1[con1<0] <- 0
    con2[con2<0] <- 0
    conc[conc<0] <- 0
    
    if(!bNormalized || is.null(data$fold)) {
      fold <- con1 - con2
    } else {
      fold <- data$fold[keep]
      if(bFlip == TRUE) {
        fold <- fold * -1
      }
    }

    data <- cbind(pv.getsites(pv,siteids),conc,con1,con2,fold,data[keep,c(pvCol,fdrCol)])
    conc1 <- sprintf('Conc_%s',name1)
    conc2 <- sprintf('Conc_%s',name2)
    colnames(data) <- c('Chr','Start','End','Conc',conc1,conc2,'Fold','p-value','FDR')
    
  } else {
    conc <- log2(apply(counts,1,mean))
    conc[conc<0] <- 0
    fold <- data$fold[keep]
    data <- cbind(pv.getsites(pv,siteids),conc,fold,data[keep,c(pvCol,fdrCol)])
    colnames(data) <- c('Chr','Start','End','Conc','Fold','p-value','FDR')
  }
  
  if(bCounts) {
    counts <- round(counts,2)
    if(groups) {
      colnames(counts) <- c(pv$class[PV_ID,group1],pv$class[PV_ID,group2])
    } else {
      colnames(counts) <- pv$class[PV_ID,]
    }
    if(length(sites)>1){
      data <- cbind(data,counts)
    } else {
      dnames <- colnames(data)
      cnames <- colnames(counts)
      data <- cbind(data,matrix(counts,1,ncol(counts)))
      colnames(data) <- c(dnames,cnames)
    }
  }
  
  if(groups) {
    
    if(bCalled & !is.null(pv$called)) {
      Called1 <- apply(pv$called[siteids,group1],1,sum)
      Called2 <- apply(pv$called[siteids,group2],1,sum)
      data <- cbind(data,Called1,Called2)
    }
    
    if(bCalledDetail & !is.null(pv$called)) {
      newd <- pv$called[siteids,c(which(group1),which(group2))]
      newd[newd==1] <- '+'
      newd[newd==0] <- '-'
      colnames(newd) <- c(pv$class[PV_ID,group1],pv$class[PV_ID,group2])
      data <- cbind(data,newd)
    }
  }
  
  if(minFold>0) {
    data <- data[abs(data$Fold)>=minFold,]
  }
  
  if(bUsePval) {
    data <- data[order(data$'p-value'),]
  } else {
    data <- data[order(data$FDR,data$'p-value'),]
  }
  
  if(length(precision==1) && precision[1]>0) {
    precision <- c(precision,precision)
  }
  
  if(length(precision)==2) {
    if(groups) {
      data[,4:7] <- round(data[,4:7],precision[1])
      data[,8:9] <- signif(data[,8:9],precision[2])
    } else {
      data[,4:5] <- round(data[,4:5],precision[1])
      data[,6:7] <- signif(data[,6:7],precision[2])
    }
  }
  
  if(!missing(file)) {
    if(!is.null(initString)) {
      initString <- sprintf("%s_",initString)
    } else {
      initString <- ""
    }
    if(is.null(file)) {
      file=sprintf("%s%s_vs_%s_%s.%s",initString,name1,name2,method,ext)
    } else {
      file=sprintf("%s%s.%s",initString,file,ext)
    }
    write.csv(data,row.names=FALSE,file=file)
  }
  
  return(data)
  
}

pv.getsites <- function(pv,sites){
  if(is.logical(sites)) {
    sites <- which(sites)
  }
  siteNum <- as.integer(sites)
  idx   <- match(siteNum,rownames(pv$binding))
  sites <- data.frame(pv$binding[idx,1:3])
  if(length(idx)==1) {
    sites <- data.frame(t(sites))
    rownames(sites) <- siteNum
  }
  sites[,1] <- pv$chrmap[sites[,1]]
  return(sites)
}

pv.getPlotData <- 
  function(pv,attributes=PV_GROUP,contrast=1,method=DBA_DESEQ2,th=0.05,
           bUsePval=FALSE,bNormalized=TRUE,report,
           bPCA=FALSE,bLog=TRUE,minval,maxval,mask,fold=0,
           bFlip=FALSE, precision=2:3) {
    
    if(contrast > length(pv$contrasts)) {
      stop('Specified contrast number is greater than number of contrasts',call.=FALSE)
      return(NULL)
    }
    
    con <- pv$contrasts[[contrast]]
    
    if(is.null(con$group1)) {
      pv <- pv.getPlotDataComplexContrast(pv=pv,attributes-attributes, contrast=contrast,
                                          method=method, bUsePval=bUsePval,mask=mask,
                                          bNormalized=bNormalized,report=report,
                                          bLog=bLog,minval=minval,maxval=maxval,
                                          fold=fold, precision=precision)
      return(pv)
    }
    
    if(missing(report)) {
      report <- pv.DBAreport(pv,contrast=contrast,method=method,th=th,bUsePval=bUsePval,
                             bNormalized=bNormalized,bCounts=TRUE,bSupressWarning=TRUE,
                             minFold=fold,bFlip=bFlip,precision=precision)
      if(is.null(report)) {
        stop('Unable to plot -- no sites within threshold',call.=FALSE)	
      }
    } else {
      report <- report[abs(report$Fold)>=fold,]
    }
    
    if(nrow(report)==1) {
      stop("Only one site to plot -- need at least 2!",call.=FALSE)   
    }
    
    if(!missing(mask)){
      if (!is.logical(mask)) {
        if (max(mask) > length(pv$peaks)) {
          stop("Invalid sample number in mask.",call.=FALSE)
        }
        temp <- rep(FALSE, length(pv$peaks))
        temp[mask] <- TRUE
        mask <- temp
      }
      if(!bFlip) {
        group1 <- con$group1 & mask
        group2 <- con$group2 & mask
      } else {
        group1 <- con$group2 & mask
        group2 <- con$group1 & mask
      }
      sites <- as.numeric(rownames(report))
      extra <- mask & !(group1 | group2)
      allsamps <- c(which(group1), which(group2), which(extra))
      numsamps <- length(allsamps)
      domap <- matrix(0,length(sites),0)
      if(sum(group1)) {     
        domap <- cbind(domap,pv$binding[sites,3+which(group1)])
      }
      if(sum(group2)) {
        domap <- cbind(domap,pv$binding[sites,3+which(group2)])
      }
      if(sum(extra)) {
        domap <- cbind(domap,pv$binding[sites,3+which(extra)])
      }
      rownames(domap) <- rownames(report)
      colnames(domap) <- pv$class[PV_ID,allsamps] 
      con$group1 <- group1
      con$group2 <- group2
      peaks <- pv$peaks[allsamps]
      for(i in 1:length(peaks)){
        peaks[[i]] <- peaks[[i]][sites,]
      }  
    } else {
      if(!bFlip) {
        group1 <- con$group1 
        group2 <- con$group2 
      } else {
        group1 <- con$group2 
        group2 <- con$group1 
      }
      allsamps <- c(which(group1),which(group2))
      extra <- rep(FALSE,ncol(pv$class)) 
      repcols <- colnames(report)
      numsamps <- sum(group1)+sum(group2)
      if(length(repcols) < (numsamps+9)) {
        stop('Report does not have count data, re-run dba.report with bCounts=TRUE',call.=FALSE)
      }
      first <- 10
      if(repcols[10]=="Called1") {
        if(length(repcols) < (numsamps+11)) {
          stop('Report does not have count data, re-run dba.report with bCounts=TRUE',call.=FALSE)
        }
        first <- 12	 
      }
      domap <- report[,first:(first+numsamps-1)]
      group1 <- rep(FALSE,numsamps)
      group2 <- rep(TRUE,numsamps)
      group1[1:sum(con$group1)] <- TRUE
      group2[1:sum(con$group1)] <- FALSE
      con$group1 <- group1
      con$group2 <- group2
      
      sites <- as.numeric(rownames(report))
      peaks <- pv$peaks[allsamps]
      for(i in 1:length(peaks)){
        peaks[[i]] <- peaks[[i]][sites,]
      }  
    }
    
    if(bLog) {
      domap[domap<=0]=1
      domap <- log2(domap)
      if(missing(minval)) {
        minval <- 0
      } else {
        minval <- max(0,minval)
      }
    }
    
    if(!missing(minval)) {
      domap[domap< minval]= minval
    }
    if(!missing(maxval)) {
      domap[domap>maxval] <- maxval
    }
    
    pv$binding <- cbind(report[,1:3],domap)
    pv$class <- pv$class[,allsamps]
    pv$peaks <- peaks
    pv$called <- pv$called[sites,allsamps]
    pv$merged <- pv$binding[,1:3]
    pv$totalMerged <- nrow(pv$merged)
    pv$contrasts <- list(pv$contrasts[[contrast]])
    pv$contrasts[[1]]$group1 <- rep(FALSE,ncol(pv$class))
    pv$contrasts[[1]]$group2 <- rep(FALSE,ncol(pv$class))
    if(sum(con$group1)) {
      pv$contrasts[[1]]$group1[1:sum(con$group1)] <- TRUE
    }
    if(sum(con$group2)) {
      pv$contrasts[[1]]$group2[(sum(con$group1)+1):(sum(con$group1)+sum(con$group2))] <- TRUE 
    }
    
    if(bPCA)  {
      
      #pv$pc <- princomp(domap,cor=bPCAcor)
      
      if(attributes[1] == PV_GROUP) {
        pv$class[PV_ID,] <- c( rep(con$name1,sum(con$group1)), rep(con$name2,sum(con$group2)), rep("other",sum(extra)) )
      }
    } 
    
    return(pv)     
  }

pv.getPlotDataComplexContrast <- 
  function(pv,attributes=PV_GROUP,contrast=1,method=DBA_DESEQ2,th=0.05,
           bUsePval=FALSE,mask,bNormalized=TRUE,report,
           bLog=TRUE,minval,maxval,fold=0,
           precision=2:3) {
    
    if(contrast > length(pv$contrasts)) {
      stop('Specified contrast number is greater than number of contrasts',call.=FALSE)
      return(NULL)
    }
    
    con <- pv$contrasts[[contrast]]
    
    if(missing(report)) {
      report <- pv.DBAreport(pv,contrast=contrast,method=method,th=th,bUsePval=bUsePval,
                             bNormalized=bNormalized,bCounts=TRUE,bSupressWarning=TRUE,
                             minFold=fold,precision=precision)
      if(is.null(report)) {
        stop('Unable to plot -- no sites within threshold',call.=FALSE)	
      }
    } else {
      report <- report[abs(report$Fold)>=fold,]
    }
    
    if(nrow(report)==1) {
      stop("Only one site to plot -- need at least 2!",call.=FALSE)   
    }
    
    sites <- as.numeric(rownames(report))
    
    if(missing(mask)){
      numsamps <- ncol(pv$class)
      allsamps <- 1:numsamps
    } else {
      if (!is.logical(mask)) {
        if (max(mask) > length(pv$peaks)) {
          stop("Invalid sample number in mask.",call.=FALSE)
        }
        temp <- rep(FALSE, length(pv$peaks))
        temp[mask] <- TRUE
        mask <- temp
      }
      
      allsamps <- which(mask)
      numsamps <- length(allsamps)
      
      repmask <- rep(TRUE,ncol(report))
      repmask[(length(repmask)-length(mask)+1):length(repmask)] <- mask
      report <- report[,repmask]
    } 
    
    repcols <- colnames(report)
    if(length(repcols) < (numsamps+7)) {
      stop('Report does not have count data, re-run dba.report with bCounts=TRUE',call.=FALSE)
    }
    
    first <- match("FDR",repcols) + 1
    if(repcols[first]=="Called1") {
      if(length(repcols) < (numsamps+first-1)) {
        stop('Report does not have count data, re-run dba.report with bCounts=TRUE',call.=FALSE)
      }
      first <- first+2	 
    }
    
    domap <- report[,first:(first+numsamps-1)]
    rownames(domap) <- rownames(report)
    colnames(domap) <- pv$class[PV_ID,allsamps]
    
    peaks <- pv$peaks[allsamps]
    for(i in 1:length(peaks)){
      peaks[[i]] <- peaks[[i]][sites,]
    }  
    
    if(bLog) {
      domap[domap<=0]=1
      domap <- log2(domap)
      if(missing(minval)) {
        minval <- 0
      } else {
        minval <- max(0,minval)
      }
    }
    
    if(!missing(minval)) {
      domap[domap< minval]= minval
    }
    if(!missing(maxval)) {
      domap[domap>maxval] <- maxval
    }
    
    pv$binding <- cbind(report[,1:3],domap)
    pv$class <- pv$class[,allsamps]
    pv$peaks <- peaks
    pv$called <- pv$called[sites,allsamps]
    pv$merged <- pv$binding[,1:3]
    pv$totalMerged <- nrow(pv$merged)
    pv$contrasts <- list(pv$contrasts[[contrast]])
    
    return(pv)     
  }

pv.resultsDBA <- function(DBA,contrasts,methods=DBA$config$AnalysisMethod,
                          th=0.05,bUsePval=FALSE,fold=0,
                          bDB=TRUE,bNotDB=FALSE,bUp=FALSE,bDown=FALSE,bAll=TRUE,
                          bFlip=FALSE) {
  
  if(missing(contrasts)) {
    contrasts <- 1:length(DBA$contrasts)
  }
  
  res <- NULL
  for(contrast in contrasts) {
    for(method in methods) {
      res <- pv.doResults(res,DBA,contrast,method,th,bUsePval,fold,
                          bDB=bDB,bNotDB=bNotDB,bUp=bUp,bDown=bDown,bAll=bAll,
                          bFlip=bFlip)		
    }
  }
  
  if(is.null(res)) {
    warning('No valid contrasts/methods specified.',call.=FALSE)	
  }
  
  res$config <- DBA$config
  res$config$id        <- "Contrast"
  res$config$factor    <- "DB"
  res$config$tissue    <- "Direction"
  res$config$condition <- "Method"
  res$config$treatment <- "Block"
  
  res <- pv.model(res, minOverlap=1)
  res <- pv.ResetScores(res, ones=FALSE)
  res$resultsObject <- TRUE
  res$score <- DBA_SCORE_FOLD
  
  class(res) <- "DBA"
  
  return(res)	
  
}

pv.doResults <- function(res,DBA,contrast,method,th,bUsePval,
                         fold=0,bDB=TRUE,bNotDB=FALSE,bAll=FALSE,bUp=FALSE,bDown=FALSE,
                         bFlip=FALSE) {
  
  if(method=='edgeR' || method=="edgeRGLM") {
    if(is.null(DBA$contrasts[[contrast]]$edgeR) && is.null(DBA$edgeR)) {
      return(res)	
    }
    methname <- "edgeR"
    block <- ""
  }
  if (method=='edgeRlm') {
    if(is.null(DBA$contrasts[[contrast]]$edgeR$block)) {
      return(res)	
    }
    methname="edgeR"
    block="block"	
  }
  if(method=='DESeq2' || method=='DESeq2GLM') {
    if(is.null(DBA$contrasts[[contrast]]$DESeq2) && is.null(DBA$DEseq2)) {
      return(res)  
    }      
    methname <- "DESeq2"
    block <- ""
  }
  if (method=='DESeq2Block') {
    if(is.null(DBA$contrasts[[contrast]]$DESeq2$block)) {
      return(res)	
    }           
    methname="DESeq2"
    block="block"	
  }
  
  if(bNotDB) {
    useth <- 1
  } else {
    useth <- th
  }
  
  rep <- suppressWarnings(dba.report(DBA,contrast=contrast,method=method,
                                     th=useth,bUsePval=bUsePval,fold=0,
                                     bFlip=bFlip,
                                     DataType=DBA_DATA_FRAME))
  
  if(is.null(rep)) {
    return(res)	
  }
  
  rep <- rep[order(as.integer(rownames(rep))),]
  
  id <- pv.getContrastString(DBA$contrast[[contrast]],bFlip) 
  
  if(bUsePval) {
    scores <- rep$"p-value"
  } else {
    scores <- rep$FDR
  }    
  
  db <- (scores <= th) & (abs(rep$Fold) >= fold)
  up <- rep$Fold >= 0
  peaks <- cbind(rep[,1:3],rep$Fold,rep[,4:ncol(rep)])
  
  if(bDB) {
    if(bAll) {
      if(sum(db)) {
        res <- pv.peakset(res,peaks=peaks[db,],sampID=id,factor="DB",tissue="All",
                          condition=methname,treatment=block)	
      }
    }
    if(bUp) {
      if(sum(db&up)) {	
        res <- pv.peakset(res,peaks=peaks[db&up,],sampID=id, factor ="DB", tissue ="Gain",
                          condition=methname,treatment=block)	
      }
    }
    if(bDown) {
      if(sum(db&!up)) {
        res <- pv.peakset(res,peaks=peaks[db&!up,],sampID=id, factor ="DB", tissue ="Loss",
                          condition=methname,treatment=block)	
      }
    }		
  }
  if(bNotDB) {
    if(bAll) {
      if(sum(!db)) {
        res <- pv.peakset(res,peaks=peaks[!db,],sampID=id, factor ="!DB", tissue ="All",
                          condition=methname,treatment=block)	
      }
    }
    if(bUp) {
      if(sum(!db&up)) {
        res <- pv.peakset(res,peaks=peaks[!db&up,],sampID=id, factor ="!DB", tissue ="Gain",
                          condition=methname,treatment=block)	
      }
    }
    if(bDown) {
      if(sum(!db&!up)) {
        res <- pv.peakset(res,peaks=peaks[!db&!up,],sampID=id, factor ="!DB", tissue ="Loss",
                          condition=methname,treatment=block)	
      }
    }		
  }    
  
  return(res)	
}

pv.ResetScores <- function(pv, ones=FALSE){
  
  return(pv)
  
  binding <- pv$binding[,1:3]
  binding[,1] <- pv$chrmap[binding[,1]]
  binding <- GRanges(data.frame(binding))
  GenomicRanges::mcols(binding) <- pv$binding[,4:ncol(pv$binding)]
  
  for(i in 1:length(pv$peaks)) {
    peaks <- GRanges(data.frame(pv$peaks[[i]]))
    matches <- which(binding %over% peaks)
    keep <- which(peaks %over% binding)
    if(length(matches) != length(keep)) {
      return(pv)
    }
    scores <- pv$peaks[[i]]$Score
    if(ones==TRUE) {
      scores[scores<1] <- -1
      scores[scores>1] <-  1
    }
    GenomicRanges::mcols(binding)[matches,i] <- scores[keep]
  }
  
  pv$binding[,4:ncol(pv$binding)] <- as.matrix(GenomicRanges::mcols(binding))
  return(pv)
}

PV_SCORE_FOLD             <- "score_fold"
PV_SCORE_CONCENTRATION    <- "score_conc"
PV_SCORE_CONC_NUMERATOR   <- "score_num"
PV_SCORE_CONC_DENOMINATOR <- "score_denom"
PV_SCORE_PVAL             <- "score_pval"
PV_SCORE_FDR              <- "score_fdr"

pv.ResetResultScores <- function(pv, score) {
  
  if(score == pv$score) {
    return(pv)
  }
  
  if(score == PV_SCORE_FOLD) {
    pv <- pv.changeScore(pv, 4)
  } else if(score == PV_SCORE_CONCENTRATION) {
    pv <- pv.changeScore(pv, 1)      
  } else if(score == PV_SCORE_CONC_NUMERATOR) {
    pv <- pv.changeScore(pv, 2)      
  } else if(score == PV_SCORE_CONC_DENOMINATOR) {
    pv <- pv.changeScore(pv, 3)      
  } else if(score == PV_SCORE_PVAL) {
    pv <- pv.changeScore(pv, 5)      
  } else if(score == PV_SCORE_FDR) {
    pv <- pv.changeScore(pv, 6)      
  } else {
    message("Invalid score")
    return(pv)
  }
  
  config <- pv$config
  pv <- pv.model(pv, minOverlap=pv$minOverlap) 
  pv$config <- config
  pv$resultsObject <- TRUE
  pv$score <- score
  
  class(pv) <- "DBA"
  
  return(pv)
}

pv.changeScore <- function(pv, colnum=1) {
  for(i in 1:length(pv$peaks)) {
    colnum <- 4+colnum
    pv$peaks[[i]]$Score <- pv$peaks[[i]][,colnum]
  }
  return(pv)
}

