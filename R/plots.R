pv.DBAplotMA <- function(pv,contrast,method='edgeR',bMA=TRUE,bXY=FALSE,th=0.05,
                         bUsePval=FALSE,fold=0,facname="",bNormalized=TRUE,
                         cex=.15,bSignificant=TRUE, bSmooth=TRUE, bFlip=FALSE,
                         xrange, yrange, bLoess=TRUE, ...) {
  
  noreport <- FALSE
  if(missing(contrast)){
    contrast <- 1:length(pv$contrasts)
  } else if(is(contrast,"list")){
    noreport <- TRUE
    group1 <- contrast[[1]]
    name1 <- names(contrast)[1]
    if(!is(group1,"logical")) {
      g1 <- rep(FALSE,length(pv$peaks))
      g1[group1] <- TRUE
      group1 <- g1
    }
    if(length(contrast) == 2) {
      group2 <- contrast[[2]]
      name2 <- names(contrast)[2]
      if(!is(group2,"logical")) {
        g1 <- rep(FALSE,length(pv$peaks))
        g1[group2] <- TRUE
        group2 <- g1
      }
    } else if (length(contrast) == 1) {
      group2 <- !group1
      name2  <- paste("!",name1,sep="")
    } else {
      stop("If contrast is a list, if muct be of length 1 or 2",call.=FALSE)
    }
    contrast <- contrast[1]
  } else if(max(contrast) > length(pv$contrasts)) {
    stop('Specified contrast number is greater than number of contrasts',call.=FALSE)
    return(NULL)
  }
  
  
  plotfun <- plot
  if (bSmooth) {
    plotfun <- smoothScatter
  }
  
  numSites <- nrow(pv$binding)
  
  for(con in 1:length(contrast)) {
    if(noreport) {
      conrec <- NULL
      conrec$group1 <- group1
      conrec$group2 <- group2
      conrec$name1  <- name1
      conrec$name2  <- name2
    } else {
      conrec <- pv$contrasts[[contrast[con]]]
      name1 <- conrec$name1
      name2 <- conrec$name2
    }
    if(bFlip) {
      name1 <- conrec$name2
      name2 <- conrec$name1   
    }
    for(meth in method) {
      if(noreport) {
        res <- pv.countsMA(pv, meth, conrec, bNormalized)
      } else {
        res <- pv.DBAreport(pv,contrast=contrast[con],method=meth,bUsePval=TRUE,
                            th=100,bNormalized=bNormalized,bFlip=bFlip,precision=0)
      }
      if(!is.null(res)) {
        if(noreport) {
          bSignificant <- FALSE
          idx <- rep(FALSE, length(res$Fold))
        } else {
          if(bUsePval) {
            idx <- res$"p-value" <= th
            tstr <- "p"
          } else {
            idx <- res$FDR <= th
            tstr <- "FDR"
          }
          idx <- idx & (abs(res$Fold) >= fold)
        }
        if(bMA){
          if(missing(xrange)) {
            xmin  <- floor(min(res$Conc))
            xmax  <- ceiling(max(res$Conc))
          } else {
            if (length(xrange) != 2) {
              stop("xrange must be vector of two numbers",call.=FALSE)
            }
            xmin <- xrange[1]
            xmax <- xrange[2]
          }
          if(missing(yrange)) {
            ymin  <- floor(min(res$Fold))
            ymax  <- ceiling(max(res$Fold))
          } else {
            if (length(yrange) != 2) {
              stop("yrange must be vector of two numbers",call.=FALSE)
            }
            ymin <- yrange[1]
            ymax <- yrange[2]
          }
          
          constr <- pv.getContrastString(conrec)
          
          if(is(res,"list")) {
            mainstr <- constr
          } else {
            if(facname=="") {
              mainstr <- sprintf('%s (%s %s < %1.3f)', 
                                 constr,sum(idx),tstr,th) 
            } else {
              mainstr <- sprintf('%s: %s (%s %s < %1.3f)', 
                                 facname, constr,sum(idx),tstr,th)               
            }
          }
          
          if(is.null(conrec$group1)) {
            ylabstr <- sprintf('log Fold Change')
          } else {
            ylabstr <- sprintf('log Fold Change: %s - %s',name1,name2)
          }
          
          if(!is.na(match("main",names(list(...))))){
            mainstr <- NULL
          }
          
          if(bSmooth | !bSignificant) {
            plotfun(res$Conc,res$Fold,pch=20,cex=cex,col=crukBlue,
                    xaxp=c(xmin,xmax,xmax-xmin),xlim=c(xmin,xmax),
                    xlab='log concentration',
                    yaxp=c(ymin,ymax,(ymax-ymin)),ylim=c(ymin,ymax),
                    ylab=ylabstr,main=mainstr,...)              	
          } else {
            plotfun(res$Conc[!idx],res$Fold[!idx],pch=20,cex=cex, col=crukBlue,
                    xaxp=c(xmin,xmax,xmax-xmin),xlim=c(xmin,xmax),
                    xlab='log concentration',
                    yaxp=c(ymin,ymax,(ymax-ymin)),ylim=c(ymin,ymax),
                    ylab=ylabstr, main=mainstr,...)
          }
          if(bSignificant) {
            points(res$Conc[idx],res$Fold[idx],pch=20,cex=cex,col=crukMagenta)
          }
          abline(h=0,col='dodgerblue')
          if(bLoess) {
            lfit <- limma::loessFit(y=res$Fold,x=res$Conc)
            o <- order(res$Conc)
            lines(res$Conc[o],lfit$fitted[o],col="red")
          } 
        }
      }
      if(bXY){
        if(is.null(conrec$group1)) {
          stop('Can not plot XY for complex contrast.',call.=FALSE)
        }
        if(missing(xrange)) {
          if(is(res,"list")) {
            xmin <- floor(min(res$Con1))
            xmax <- ceiling(max(res$Con1))
          } else {
            xmin  <- floor(min(res[,5]))
            xmax  <- ceiling(max(res[,5]))
          }
        } else {
          if (length(xrange) != 2) {
            stop("xrange must be vector of two numbers",call.=FALSE)
          }
          xmin <- xrange[1]
          xmax <- xrange[2]
        }
        if(missing(yrange)) {
          if(is(res,"list")) {
            ymin <- floor(min(res$Con2))
            ymax <- ceiling(max(res$Con2))
          } else {
            ymin  <- floor(min(res[,6]))
            ymax  <- ceiling(max(res[,6]))
          }
        } else {
          if (length(yrange) != 2) {
            stop("yrange must be vector of two numbers",call.=FALSE)
          }
          ymin <- yrange[1]
          ymax <- yrange[2]
        }
        xymin <- min(xmin,ymin)
        xymin <- max(xymin,0)
        xymax <- max(xmax,ymax)
        
        constr <- pv.getContrastString(conrec)
        
        if(is(res,"list")) {
          mainstr <- constr
          xvals <- res$Con1
          yvals <- res$Con2
        } else {
          xvals <- res[,5]
          yvals <- res[,6]
          if(facname=="") {
            mainstr <- sprintf('%s (%s %s < %1.3f)', 
                               constr,sum(idx),tstr,th) 
          } else {
            mainstr <- sprintf('%s: %s (%s %s < %1.3f)', 
                               facname, constr,sum(idx),tstr,th)               
          }
        }
        
        plotfun(yvals[!idx], xvals[!idx], pch=20,cex=cex,col=crukBlue,
                xaxp=c(xymin,xymax,xymax-xymin),xlim=c(xymin,xymax),
                xlab=sprintf('log concentration :%s',name2),
                yaxp=c(xymin,xymax,(xymax-xymin)),ylim=c(xymin,xymax),
                ylab=sprintf('log concentration :%s',name1),
                main=mainstr,...)
        points(yvals[idx],xvals[idx],pch=20,cex=cex,col=crukMagenta)
        abline(0,1,col='red')
        
      }
    }
  }      	
}	

pv.countsMA <- function(pv, method, contrast, bNormalized=FALSE) {
  
  if(bNormalized) {
    
    design <- pv$design
    if(is.null(design)) {
      pv$design <- "~1"
    }
    
    if(method == DBA_DESEQ2) {
      
      if(is.null(pv$norm$DESeq2)) {
        pv <- dba.normalize(pv, method=DBA_DESEQ2, background=FALSE)
      }
      
      if(is.null(pv$DESeq2$DEdata)) {
        dedata <- pv.DEinitDESeq2(pv,
                                  bSubControl=pv$norm$DESeq2$bSubControl,
                                  bFullLibrarySize=pv$norm$DESeq2$lib.sizes,
                                  filter=pv$norm$DESeq2$filter.val,
                                  filterFun=pv$norm$DESeq2$filter.fun) 
      } else {
        dedata <- pv$DESeq2$DEdata
      }
      counts <- DESeq2::counts(dedata, normalized = bNormalized)
      
    } else if(method == DBA_EDGER) {
      
      if(is.null(pv$norm$edgeR)) {
        pv <- dba.normalize(pv, method=DBA_EDGER, background=FALSE)
      }
      
      if(is.null(pv$edgeR$DEdata)) {
        dedata <- pv.DEinitedgeR(pv,
                                 bSubControl=pv$norm$edgeR$bSubControl,
                                 filter=pv$norm$edgeR$filter.val,
                                 filterFun=pv$norm$edgeR$filter.fun)
      } else {
        dedata <-pv$norm$edgeR 
      }
      counts <- pv.edgeRCounts(pv,method,bNormalized)
    }
    counts <- cbind(counts[,contrast$group1], counts[,contrast$group2])
    pv$design <- design
  } else {
    counts <- pv.DEinit(pv,contrast$group1, contrast$group2, bRawCounts=TRUE)
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

Legend <- Fold <- NULL
pv.DBAplotVolcano <- function(pv,contrast,method='edgeR', th=0.05,
                              bUsePval=FALSE,fold=0,facname="",
                              bLabels=FALSE,maxLabels=50,
                              dotSize=1,bSignificant=TRUE, bFlip=FALSE,
                              xrange,yrange) {
  
  if(missing(contrast)){
    contrast <- 1:length(pv$contrasts)
  } else {
    if(contrast > length(pv$contrasts)) {
      stop('Specified contrast number is greater than number of contrasts',call.=FALSE)
      return(NULL)
    }
  }
  
  for(con in 1:length(contrast)) {
    conrec <- pv$contrasts[[contrast[con]]]
    name1 <- conrec$name1
    name2 <- conrec$name2
    if(bFlip) {
      name1 <- conrec$name2
      name2 <- conrec$name1   
    }
    for(meth in method) {
      res <- pv.DBAreport(pv,contrast=contrast[con],method=meth,bUsePval=TRUE,
                          th=100,bNormalized=TRUE,bFlip=bFlip,precision=0)
      
      if(!is.null(res)) {
        if(bUsePval) {
          vals <- res$"p-value" 
          idx  <- vals <= th
          tstr <- "p"
          res = mutate(res,
                       Legend=ifelse(res$"p-value"<=th,
                                     sprintf(" p-val<=%1.2f",th),
                                     sprintf(" p-val >%1.2f",th)))
        } else {
          vals <- res$FDR 
          idx  <- vals <= th
          tstr <- "FDR"
          res = mutate(res,
                       Legend=ifelse(res$FDR<th,
                                     sprintf(" FDR<=%1.2f",th),
                                     sprintf(" FDR >%1.2f",th)))
        }
        
        res$Legend[idx & abs(res$Fold) < fold] <- 
          sprintf("abs(Fold)<%1.2f",2^fold)
        idx <- idx & abs(res$Fold) >= fold 
        
        sigSites <- res[idx,]
        rownames(sigSites) <- 1:sum(idx)
        
        res <- cbind(0,res)
        colnames(res)[1] <- "SiteNum"
        res[idx,1] <- 1:sum(idx)
        
        constr <- pv.getContrastString(conrec)
        plotTitle <- sprintf('%s Contrast: %s  [%s %s<=%1.3f',
                             facname, constr,sum(idx),tstr,th)
        if(fold>0) {
          plotTitle <- sprintf("%s & abs(Fold)>=%1.2f]",
                               plotTitle, 2^fold)
        } else {
          plotTitle <- sprintf("%s]",plotTitle)
        }
        if(is.null(conrec$name2)) {
          xLabel <- "log Fold Change"
        } else {
          xLabel <- sprintf('log Fold Change [log2(%s) - log2(%s)]',name1,name2)
        }
        
        yLabel <- sprintf("-log10(%s)",tstr)
        
        p <- ggplot(res,aes(Fold,-log10(vals))) +
          geom_point(aes(col=Legend),size=dotSize) +
          scale_color_manual(values=c(crukBlue,crukMagenta,crukGrey)) + 
          labs(title=plotTitle,x=xLabel,y=yLabel)
        
        if(bLabels) {
          maxLabels <- min(sum(idx),maxLabels)
          if(maxLabels > 0) {
            xx <-  which(idx)[1:maxLabels]
            p <- p + geom_text_repel(data=sigSites[1:maxLabels,],
                                     aes(x=Fold, 
                                         y = -log10(vals[xx]),
                                         label=rownames(sigSites)[1:maxLabels]))
          }
        }
        plot(p)
      }
    }
  }
  return(sigSites[,-10])
}

## pv.plotHeatmap -- draw a heatmap using a subset of binding sites
PV_ONLYA <- 3
PV_ONLYB <- 4
PV_INALL <- 5
PV_COR   <- 6
PV_OLAP  <- 7
PV_TOTAL <- 0
pv.plotHeatmap <-
  function(pv,numSites = 1000,attributes = pv$attributes,mask,sites,contrast,
           overlaps, olmask, olPlot = PV_COR,divVal,RowAttributes,ColAttributes,rowSideCols,colSideCols,
           bTop = TRUE,minval,maxval,bReorder = FALSE,ColScheme =
             "Greens",distMeth = "pearson",...) {
    pv <- pv.check(pv)
    
    if (missing(mask)) {
      mask <- rep(TRUE,ncol(pv$class))
    } else if (is.null(mask)) {
      mask <- rep(TRUE,ncol(pv$class))
    }
    
    ocm <- NULL
    if (!missing(overlaps)) {
      cres  <- overlaps
      if (!missing(olmask)) {
        cres <- cres[olmask,]
      }
      if (is.null(nrow(cres))) {
        cres <- matrix(cres,1,length(cres))
      }
      #cres <- pv.overlapToLabels(pv,cres,attributes)
      if (missing(divVal)) {
        ocm <- pv.occ2matrix(cres,olPlot)
      } else {
        ocm <- pv.occ2matrix(cres,olPlot,divVal)
      }
      labels <- pv.nums2labels(pv,rownames(ocm),attributes)
      rowlab <- collab <- labels
      rownames(ocm) <- labels
      colnames(ocm) <- labels
      domap <- ocm
      
    } else {
      if (missing(sites)) {
        sites <- 1:nrow(pv$binding)
        numSites <- min(length(sites),numSites)
      } else {
        if (sum(sites) <= length(sites)) {
          numSites <- min(sum(sites),numSites)
          sites <- which(sites)
        } else {
          numSites <- length(sites)
        }
      }
      
      if (bTop == TRUE) {
        sites <- sites[1:numSites]
      } else {
        tsites <- length(sites)
        sites <- sites[(tsites - numSites + 1):tsites]
      }
      colnames <- NULL
      for (i in 1:ncol(pv$class)) {
        colnames <-
          c(colnames,pv.namestrings(pv$class[attributes,i])$tstring)
      }
      domap <- matrix(0.1,length(sites),sum(mask))
      for (i in 1:ncol(domap)) {
        domap[,i] <- as.numeric(pv$binding[sites,c(FALSE,FALSE,FALSE,mask)][,i])
      }
      rowlab <- ""
      collab <- colnames[mask]
    }
    
    if (!missing(minval)) {
      domap[domap < minval] <- minval
    }
    if (!missing(maxval)) {
      domap[domap > maxval] <- maxval
    }
    
    if (missing(overlaps)) {
      for (i in 1:ncol(domap)) {
        if (sum(domap[,i] == 0) == nrow(domap)) {
          domap[1,i] <- 0.000001
        }
      }
    }
    
    domap <- domap[rowSums(domap)!=0,]
    
    if (length(ColScheme) == 1) {
      cols <- colorRampPalette(brewer.pal(9,ColScheme))(256)
    } else
      cols <- ColScheme
    
    if (missing(rowSideCols)) {
      rowSideCols <- pv.colsv
    }
    rowatts <- NULL
    rowcols <- 0
    if (missing(RowAttributes)) {
      if (!missing(contrast)) {
        if (sum(pv$contrasts[[contrast]]$group1) ||
            sum(pv$contrasts[[contrast]]$group2)) {
          rowatts <-
            pv.attributematrix(pv,mask,contrast = contrast,PV_GROUP,rowSideCols)
          rowcols <- length(unique(as.vector(rowatts)))
        }
      }
    } else if (!is.null(RowAttributes)) {
      rowatts <-
        pv.attributematrix(
          pv,mask,contrast = contrast,RowAttributes,rowSideCols,bReverse = TRUE
        )
      rowcols <- length(unique(as.vector(rowatts)))
    }
    
    if (missing(colSideCols)) {
      colSideCols <- pv.colsv
    }
    if (missing(ColAttributes)) {
      colatts <-
        pv.attributematrix(pv,mask,contrast = contrast,NULL,colSideCols,bAddGroup =
                             is.null(ocm))
      colcols <- length(unique(as.vector(colatts)))
    } else if (!is.null(ColAttributes)) {
      colatts <-
        pv.attributematrix(pv,mask,contrast = contrast,ColAttributes,colSideCols)
      colcols <- length(unique(as.vector(colatts)))
    } else {
      colatts <- NULL
      colcols <- 0
    }
    
    if (is.null(ocm)) {
      if (!missing(RowAttributes)) {
        warning("Row color bars not permitted for peak score heatmaps.",
                call.=FALSE)
      }
      heatmap.3(
        domap,labCol = collab,col = cols,trace = "none",labRow = rowlab,
        KeyValueName = "Score",
        distfun = function(x) amap::Dist(x,method = distMeth),
        ColSideColors = colatts,NumColSideColors = colcols,
        ...
      )
    } else {
      res <-
        heatmap.3(
          domap,labCol = collab,col = cols,trace = "none",labRow = rowlab, 
          KeyValueName ="Correlation",
          distfun = function(x) amap::Dist(x,method = distMeth),
          symm = TRUE,revC = TRUE,Colv = TRUE,
          RowSideColors = rowatts,ColSideColors = colatts,
          NumRowSideColors = rowcols,NumColSideColors = colcols,
          ...
        )
      if (bReorder) {
        if (length(unique(rownames(ocm))) == nrow(ocm)) {
          ocm <- pv.reorderM(ocm,res$rowDendrogram)
        } else {
          warning("Unable to re-order returned correlation matrix as labels are non-unique",call.=FALSE)
        }
      }
      ocm <- signif(ocm,2) 
      return(ocm)
    }
  }



## pv.plotPCA -- 3D plot of PCA
pv.plotPCA <-
  function(pv,attributes = PV_ID,second,third,fourth,size,mask,
           numSites,sites,cor = FALSE,comps=1:3, b3D = TRUE,vColors,
           label = NULL,bLog = TRUE,labelSize = .8,labelCols = "black",...) {
    
    pv <- pv.check(pv)
    
    if(length(comps)==1) {
      c1 <- comps
      c2 <- comps+1
      c3 <- comps+2
    } else{
      c1 <- comps[1]
      c2 <- comps[2]         
      if(length(comps)==2) {
        c3 <- comps[2]+1
      } else {
        c3 <- comps[3]
      }
    }
    compnums <- c(c1,c2,c3)
    
    class  <- attributes[1]
    if (length(attributes) > 1) {
      second <- attributes[2]
      if (length(attributes) > 2) {
        third <- attributes[3]
      }
      if (length(attributes) > 3) {
        fourth <- attributes[4]
      }
    }
    
    if (missing(sites))
      sites <- NULL
    
    if (missing(numSites)) {
      numSites <- nrow(pv$binding)
    }
    
    config <- pv$config
    
    if (!missing(mask) || !missing(numSites)) {
      if (missing(mask)) {
        mask <- rep(TRUE,ncol(pv$class))
      }
      pv <-
        pv.pcmask(pv,numSites,mask,sites,cor = cor,bLog = bLog)
    }
    pv$config <- config
    pc <- pv$pc
    
    if (is.null(pc)) {
      warning("Unable to perform PCA. Make sure there aren't fewer sites than there are samples.",
              call. = FALSE)
      return(NULL)
    }
    
    #if(!is.null(pv$mask)) {
    #   classes <- pv$class[,which(pv$mask)]
    #} else {
    classes <- pv$class[,mask]
    #}
    
    if (max(class) > nrow(classes)) {
      return(FALSE)
    }
    
    vr <- rep(0,length(pc$sdev))
    for (i in 1:length(vr)) {
      vr[i] <- pc$sdev[i] ^ 2
    }
    
    if (b3D) {
      #startComp <- 1
      pvar <- sum(vr[compnums[1:3]]) / sum(vr) * 100
    } else {
      pvar <- sum(vr[compnums[1:2]]) / sum(vr) * 100
    }
    c1p <- vr[c1] / sum(vr) * 100
    c2p <- vr[c2] / sum(vr) * 100
    c3p <- vr[c3] / sum(vr) * 100
    
    if (!missing(second)) {
      if (!missing(third)) {
        if (!missing(fourth)) {
          classvec <-
            sprintf("%s:%s:%s:%s",classes[class,],classes[second,],classes[third,],classes[fourth,])
          thetitle <-
            sprintf(
              "PCA: %s+%s+%s+%s",pv.attname(class,pv),
              pv.attname(second,pv),
              pv.attname(third,pv),
              pv.attname(fourth,pv),pvar
            )
        } else {
          classvec <-
            sprintf("%s:%s:%s",classes[class,],classes[second,],classes[third,])
          thetitle <-
            sprintf(
              "PCA: %s+%s+%s",pv.attname(class,pv),
              pv.attname(second,pv),
              pv.attname(third,pv),pvar
            )
        }
      } else {
        classvec <- sprintf("%s:%s",classes[class,],classes[second,])
        thetitle <- sprintf("PCA: %s+%s",pv.attname(class,pv),
                            pv.attname(second,pv),pvar)
      }
    } else {
      classvec <- classes[class,]
      thetitle <- sprintf("PCA: %s",pv.attname(class,pv),pvar)
    }
    
    if (!is.null(label)) {
      addlabels <- classes[label,]
    } else
      addlabels <- NULL
    
    numsamps <- ncol(classes)
    if (numsamps <= 10) {
      sval <- 2.25
    } else if (numsamps <= 25) {
      sval <- 1.75
    } else {
      sval <- 1.5
    }
    if (!missing(size)) {
      if (!is.null(size)) {
        sval <- size
      }
    }
    
    if (missing(vColors)) {
      vColors <- pv.colsv
    }
    
    if (b3D) {
      if (requireNamespace("rgl",quietly = TRUE)) {
        rgl::plot3d(
          pc$x[,compnums[c(1,3,2)]],
          col = pv.colorv(classvec,vColors),type = 's',size = sval,
          xlab = sprintf('PC #%d [%2.0f%%]',c1,c1p),
          ylab = sprintf('PC #%d [%2.0f%%]',c3,c3p),
          zlab = sprintf('PC #%d [%2.0f%%]',c2,c2p),
          aspect = c(1,1,1),main = thetitle,...
        )
        uclass <- unique(classvec)
        p <- matrix(vColors[1:length(uclass)],length(uclass),1)
        rownames(p) <- uclass
        colnames(p) <- "Legend"
        for (i in 1:nrow(p)) {
          if (p[i] == crukBlue)
            p[i] <- "crukBlue"
          if (p[i] == crukMagenta)
            p[i] <- "crukMagenta"
          if (p[i] == crukCyan)
            p[i] <- "crukCyan"
        }
        print(p)
      } else {
        warning("Package rgl not installed",call.=FALSE)
        p <-
          pv.doPCAplot(pc,classvec,c1,c2,sval,vColors,thetitle,c1p,c2p,addlabels,...)
      }
    } else {
      if (!missing(size)) {
        #sval <- sval + .5
      }
      p <-
        pv.doPCAplot(
          pc,classvec,c1,c2,sval,vColors,thetitle,c1p,c2p,
          addlabels,labelSize,labelCols,...
        )
    }
    return(p)
  }

pv.doPCAplot <-
  function(pc,classvec,c1,c2,sval,vColors,thetitle,c1p,c2p,
           addlabels = NULL,labelSize = .8,labelCols = "black",...) {
    
    plotData <- as.data.frame(pc$x[,c(c1,c2)])
    colnames(plotData) <- c("PC1","PC2")
    p <- xyplot(
      PC2 ~ PC1,
      #groups=classvec,
      data = plotData,
      pch = 16, cex = sval,aspect = 1,
      col = pv.colorv(classvec,vColors),
      xlab = sprintf('Principal Component #%d [%2.0f%%]',c1,c1p),
      ylab = sprintf('Principal Component #%d [%2.0f%%]',c2,c2p),
      main = thetitle,
      key = list(
        space = "right",
        rect = list(col = unique(pv.colorv(classvec,vColors))[1:length(unique(classvec))]),
        text = list(unique(classvec)),
        rep = FALSE
      ),
      ...
    )
    
    if (!is.null(addlabels)) {
      if (length(labelCols) > 1) {
        newcols <- rep("",length(addlabels))
        colvals <- unique(addlabels)
        for (i in 1:length(colvals)) {
          newcols[addlabels %in% colvals[i]] <- labelCols[i]
        }
        labelCols <- newcols
      }
      p <- update(
        p, panel = function(x, y, ...) {
          lattice::panel.xyplot(x, y, ...);
          lattice::ltext(
            x = x, y = y, labels = as.character(addlabels), pos = 1,
            offset = 1, cex = labelSize,col = labelCols
          )
        }
      )
    }
    print(p)
  }

## pv.plotBoxplot -- Boxplots
pv.plotBoxplot <-
  function(DBA, contrast, method=DBA_DESEQ2, th=0.05, bUsePval=FALSE, bNormalized=TRUE, 
           attribute=DBA_GROUP, mask,
           bAll, bAllIncreased, bAllDecreased, bDB, bDBIncreased, bDBDecreased,
           pvalMethod=wilcox.test,  bReversePos=FALSE, attribOrder, vColors, 
           varwidth=TRUE, notch=TRUE, ...) {
    
    if (missing(bAll) &&
        missing(bAllIncreased) && missing(bAllDecreased)) {
      bMissingAll <- TRUE
    } else
      bMissingAll <- FALSE
    
    if (missing(bDB) &&
        missing(bDBIncreased) && missing(bDBDecreased)) {
      bMissingDB <- TRUE
    } else
      bMissingDB <- FALSE
    
    if (missing(contrast)) {
      if (bMissingAll) {
        bAll          <- TRUE
        bAllIncreased <- FALSE
        bAllDecreased <- FALSE
        bDB           <- FALSE
        bDBIncreased  <- FALSE
        bDBDecreased  <- FALSE
      }
    } else {
      if (bMissingAll && bMissingDB) {
        bAll          <- FALSE
        bAllIncreased <- FALSE
        bAllDecreased <- FALSE
        bDB           <- TRUE
        bDBIncreased  <- FALSE
        bDBDecreased  <- FALSE
      }
    }
    
    if (missing(vColors)) {
      vColors <- pv.colsv
      #vColors <- vColors[2:length(vColors)]
    }
    
    if (attribute == DBA_GROUP &&
        !is.null(DBA$contrasts[[contrast]]$group1)) {
      numPlots <- 2
      cols <- vColors[1:2]
      groups <-
        list(DBA$class[PV_ID,DBA$contrasts[[contrast]]$group1],DBA$class[PV_ID,DBA$contrasts[[contrast]]$group2])
      names <-
        c(DBA$contrasts[[contrast]]$name1,DBA$contrasts[[contrast]]$name2)
      if (!missing(attribOrder)) {
        if (attribOrder[1] == 2 & attribOrder[2] == 1) {
          groups <-
            list(DBA$class[PV_ID,DBA$contrasts[[contrast]]$group2],DBA$class[PV_ID,DBA$contrasts[[contrast]]$group1])
          names <-
            c(DBA$contrasts[[contrast]]$name2,DBA$contrasts[[contrast]]$name1)
        }
      }
      if(bReversePos) {
        groups <- list(groups[[2]],groups[[1]]) 
        names  <- list(names[[2]], names[[1]])
      }
    } else {
      if(attribute == DBA_GROUP) {
        attribute <- pv.attributePCA(DBA)
      }
      if(!is.null(DBA$contrasts[[contrast]]$group1)) {
        samples <-
          which(DBA$contrasts[[contrast]]$group1 |
                  DBA$contrasts[[contrast]]$group2)
      } else {
        if(missing(mask)) {
          samples <- 1:ncol(DBA$class)
        } else {
          samples <- which(mask)
        }
      }
      classes <- DBA$class[, samples]
      names <- unique(classes[attribute,])
      numPlots <- length(names)
      if (!missing(attribOrder)) {
        if (length(attribOrder < numPlots)) {
          neworder <- 1:numPlots
          neworder[1:length(attribOrder)] <- attribOrder
          attribOrder <- neworder
        }
        names  <- names[attribOrder[1:numPlots]]
      }
      cols <- vColors[1:numPlots]
      groups <- NULL
      for (name in names) {
        groups <-
          pv.listadd(groups,classes[PV_ID, classes[attribute,] == name])
      }
    }
    
    subtitle <- FALSE
    
    if (bAll | bAllIncreased | bAllDecreased) {
      report <-
        pv.DBAreport(
          DBA, contrast=contrast, method=method,th=1,
          bNormalized=bNormalized,bCounts = TRUE, precision=0
        )
      if (bReversePos) {
        increase <- report$Fold < 0
      } else {
        increase <- report$Fold > 0
      }
      if (bUsePval) {
        DB <- report$p <= th
      } else {
        DB <- report$FDR <= th
      }
    } else {
      report <- NULL
    }
    
    toplot <- NULL
    vcols  <- NULL
    vnames <- NULL
    
    if (bAll) {
      res       <- lapply(groups,pv.box,report)
      for (i in 1:length(res)) {
        names(res)[i] <- sprintf("%s",names[i])
      }
      toplot   <- c(toplot,res)
      vnames   <- c(vnames,names)
      vcols    <- c(vcols,cols)
    }
    
    if (bAllIncreased) {
      res       <- lapply(groups,pv.box,report[increase,])
      for (i in 1:length(res)) {
        names(res)[i] <- sprintf("%s+",names[i])
      }
      toplot   <- c(toplot,res)
      vnames   <- c(vnames,rep("+",numPlots))
      vcols    <- c(vcols,cols)
      subtitle <- TRUE
    }
    
    if (bAllDecreased) {
      res       <- lapply(groups,pv.box,report[!increase,])
      for (i in 1:length(res)) {
        names(res)[i] <- sprintf("%s-",names[i])
      }
      toplot   <- c(toplot,res)
      vnames   <- c(vnames,rep("-",numPlots))
      vcols    <- c(vcols,cols)
      subtitle <- TRUE
    }
    
    
    if (is.null(report)) {
      report <-
        pv.DBAreport(
          DBA, contrast=contrast, method=method,th=th,bUsePval=bUsePval,
          bNormalized=bNormalized,bCounts=TRUE, precision=0
        )
    } else {
      report <- report[DB,]
    }
    
    if (nrow(report) == 1) {
      stop('Need more than one DB site for boxplot',call.=FALSE)
    }
    
    if(bReversePos) {
      increase <- report$Fold < 0
    } else {
      increase <- report$Fold > 0
    }
    posgroup <- names[[1]]
    neggroup <- names[[2]]
    
    if (bDB) {
      res       <- lapply(groups,pv.box,report)
      for (i in 1:length(res)) {
        names(res)[i] <- sprintf("%s.DB",names[i])
      }
      toplot   <- c(toplot,res)
      vnames   <- c(vnames,names)
      vcols    <- c(vcols,cols)
    }
    
    if (bDBIncreased) {
      res <- lapply(groups,pv.box,report[increase,])
      for (i in 1:length(res)) {
        names(res)[i] <- sprintf("%s.DB+",names[i])
      }
      toplot <- c(toplot,res)
      vnames <- c(vnames,rep("+",numPlots))
      vcols  <- c(vcols,cols)
      subtitle <- TRUE
    }
    
    if (bDBDecreased) {
      res <- lapply(groups,pv.box,report[!increase,])
      for (i in 1:length(res)) {
        names(res)[i] <- sprintf("%s.DB-",names[i])
      }
      toplot <- c(toplot,res)
      vnames <- c(vnames,rep("-",numPlots))
      vcols <- c(vcols,cols)
      subtitle <- TRUE
    }
    
    if (bNormalized == TRUE) {
      ystr <- "log2 normalized reads in binding sites"
    } else {
      ystr <- "log2 reads in binding sites"
    }
    
    if (subtitle == TRUE && !is.null(DBA$contrasts[[contrast]]$group1)) {
      subt <-
        sprintf(
          "+ indicates sites with increased affinity in %s\n- indicates sites with increased affinity in %s",
          posgroup,neggroup
        )
    } else {
      subt <- ""
    }
    
    toplot <- lapply(toplot,pv.infToNA)
    mainstr <- pv.getContrastString(DBA$contrasts[[contrast]])
    boxplot(
      toplot,notch=notch, varwidth=varwidth,
      col=vcols,names=vnames,main=mainstr,
      sub=subt,ylab=ystr
    )
    
    if (!is.null(pvalMethod)) {
      pvals <- matrix(1,length(toplot),length(toplot))
      rownames(pvals) <- names(toplot)
      colnames(pvals) <- names(toplot)
      for (i in 1:(length(toplot) - 1)) {
        for (j in (i + 1):length(toplot)) {
          if (length(toplot[[i]]) == length(toplot[[j]])) {
            pvals[i,j] <-
              pvalMethod(toplot[[i]],toplot[[j]],paired=TRUE)$p.value
          } else {
            pvals[i,j] <-
              pvalMethod(toplot[[i]],toplot[[j]],paired=FALSE)$p.value
          }
          pvals[j,i] <- pvals[i,j]
        }
      }
    } else {
      pvals <- NULL
    }
    pvals <- signif(pvals,3)
    return(pvals)
  }

pv.box <- function(ids,report) {
  idx <- match(ids,colnames(report))
  res <- log2(apply(report[,idx],1,mean))
  return(res)
}

pv.infToNA <- function(vals) {
  if(min(vals)==-Inf) {
    vals[vals==-Inf] <- NA
  }
  return(vals)
}

## pv.plotVenn -- draw venn diagrams
pv.plotVenn <-
  function(ovrec,label1 = "A",label2 = "B",label3 = "C",label4 = "D",
           main = "",sub = "") {
    if (length(ovrec) == 3) {
      pv.venn2(ovrec,label1,label2,main,sub)
    }
    
    if (length(ovrec) == 7) {
      pv.venn3(ovrec,label1,label2,label3,main,sub)
    }
    
    if (length(ovrec) == 15) {
      pv.venn4(ovrec,label1,label2,label3,label4,main,sub)
    }
  }


pv.venn2 <- function(olaps,l1,l2,main="",sub="") {
  
  counts <- c(nrow(olaps$onlyA),
              nrow(olaps$onlyB),
              nrow(olaps$inAll))
  names(counts) <- c("A","B","A_B")		
  counts <- list(counts)
  vennPlot(counts,setlabels=c(l1,l2),mysub=sub,mymain=main)
  
}

pv.venn3 <- function(olaps,l1,l2,l3,main="",sub="") {
  
  counts <- c(nrow(olaps$onlyA),
              nrow(olaps$onlyB),
              nrow(olaps$onlyC),	
              nrow(olaps$notC),
              nrow(olaps$notB),
              nrow(olaps$notA),
              nrow(olaps$inAll))
  names(counts) <- c("A","B","C","A_B","A_C","B_C","A_B_C")		
  counts <- list(counts)
  vennPlot(counts,setlabels=c(l1,l2,l3),mysub=sub,mymain=main)
  
}


pv.venn4 <- function(olaps,l1,l2,l3,l4,main="",sub="") {
  
  counts <- c(nrow(olaps$onlyA),
              nrow(olaps$onlyB),
              nrow(olaps$onlyC),	
              nrow(olaps$onlyD),
              nrow(olaps$AandB),
              nrow(olaps$AandC),
              nrow(olaps$AandD),
              nrow(olaps$BandC),
              nrow(olaps$BandD),
              nrow(olaps$CandD),
              nrow(olaps$notD),
              nrow(olaps$notC),
              nrow(olaps$notB),
              nrow(olaps$notA),
              nrow(olaps$inAll))
  names(counts) <- c("A","B","C","D","A_B","A_C","A_D","B_C","B_D","C_D",
                     "A_B_C","A_B_D","A_C_D","B_C_D","A_B_C_D")		
  counts <- list(counts)
  
  vennPlot(counts,setlabels=c(l1,l2,l3,l4),mysub=sub,mymain=main)
  
}


crukMagenta=rgb(236,0,140,maxColorValue=255) #CRUK Magenta
crukBlue=rgb(46,0,139,maxColorValue=255) #CRUK Blue
crukGrey='grey'#rgb(200,201,199,maxColorValue=255) #CRUK grey
crukCyan=rgb(0,182,237,maxColorValue=255) #CRUK light blue

pv.colsv <- c(crukBlue, crukMagenta, crukCyan,crukGrey,
              "lightgreen", "orange","rosybrown",
              "black","red","dodgerblue","darkgreen",
              "yellow","grey50","purple3",
              "sienna","limegreen","lightblue",
              "violet","seagreen1",
              "lavender","olivedrab")

pv.colorv <- function(classes,cols=pv.colsv){
  
  colv <- rep(0,length(classes))
  uv <- unique(classes)
  for(i in 1:length(uv)){
    colv[classes==uv[i]] <- cols[i]
  }
  return(colv)	
} 



