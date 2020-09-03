PV_NORM_LIB            <- "lib"
PV_NORM_TMM            <- "TMM" 
PV_NORM_RLE            <- "RLE"
PV_NORM_DEFAULT        <- "default"
PV_NORM_NATIVE         <- "native"
PV_NORM_SPIKEIN        <- "spike-in"
PV_NORM_USER           <- "user"
PV_NORM_OFFSETS        <- "offsets"
PV_NORM_OFFSETS_ADJUST <- "adjust offsets"

PV_LIBSIZE_DEFAULT     <- "default"
PV_LIBSIZE_FULL        <- "full"
PV_LIBSIZE_PEAKREADS   <- "RiP"
PV_LIBSIZE_CHRREADS    <- "chrs"
PV_LIBSIZE_USER        <- "user"

PV_OFFSETS_LOESS       <- "loess"
PV_OFFSETS_USER        <- "user"

PV_BACKGROUND_BINSIZE  <- 15000
PV_NORM_LIBFUN         <- mean

pv.normalize <- function(pv, 
                         method    = DBA_ALL_METHODS,
                         libSizes  = PV_LIBSIZE_FULL,
                         normalize = PV_NORM_DEFAULT,
                         bSubControl = is.null(pv$greylist),
                         filter=0, filterFun=max, 
                         background=FALSE, 
                         offsets=FALSE,
                         libFun=PV_NORM_LIBFUN,
                         bRetrieve=FALSE, ...) {
  
  if(bRetrieve==TRUE) {
    res <- pv.normalizeRetrieve(pv, method)
    return(res)
  }
  
  if(all(method == DBA_ALL_METHODS)) {
    
    if(background[1] != FALSE || libSizes==DBA_LIBSIZE_CHRREADS) {
      pv$norm$background <- pv.getBackground(pv,background)
    }
    
    for(method in c(DBA_DESEQ2, DBA_EDGER)) {
      pv <- pv.normalize(pv, method=method,
                         libSizes=libSizes,normalize=normalize,
                         bSubControl=bSubControl,
                         filter=filter, filterFun=filterFun, 
                         background=background, 
                         offsets=offsets,
                         libFun=libFun, bRetrieve=FALSE)
    }
    return(pv)
  } else if(method == DBA_DESEQ2) {
    pv$norm$DESeq2 <- NULL
  } else if(method == DBA_EDGER) {
    pv$norm$edgeR  <- NULL
  } else {
    stop('Invalid method.',call.=FALSE)
  }
  
  norm <- NULL
  norm$bSubControl <- bSubControl
  norm$filter.val <- filter
  if(filter > 0) {
    norm$filter.fun <- filterFun
  } 
  norm$lib.method <- libSizes
  norm$background <- FALSE
  
  if(background[1] != FALSE || libSizes==DBA_LIBSIZE_CHRREADS) {
    norm.background <- pv.getBackground(pv,background) 
    binned   <- norm.background$binned
    bin.size <- norm.background$bin.size
  } else {
    binned <- pv$norm$background$binned
    bin.size <- pv$norm$background$bin.size
  }
  
  if(length(libSizes) == ncol(pv$samples)) {
    norm$lib.calc  <- "User supplied"
    norm$lib.sizes <- libSizes
    norm$lib.method <- PV_LIBSIZE_USER
  } else if (length(libSizes) == 1) {
    if(libSizes==PV_LIBSIZE_FULL) {
      norm$lib.calc  <- "Full library"
      norm$lib.sizes <- as.numeric(pv$class[PV_READS,]) 
    } else if(libSizes==PV_LIBSIZE_PEAKREADS) {
      norm$lib.calc  <- "Reads in peaks"
      norm$lib.sizes <- pv.readsInPeaks(pv, bSubControl=bSubControl,
                                        filter=filter, filterFun=filterFun)
    }  else if(libSizes==PV_LIBSIZE_CHRREADS) {
      norm$lib.calc  <- "Reads in Background"
      norm$lib.sizes <- binned$totals
    } else {
      stop('Invalid libSizes',call.=FALSE)      
    }
  } else {
    stop('libSizes invalid length',call.=FALSE)
  }
  
  doOffsets <- FALSE
  if(is(offsets,"logical")) {
    if(offsets) {
      doOffsets <- TRUE
    }
  } else if(is(offsets,"matrix")) {
    doOffsets <- TRUE
  } else if(is(offsets,"SummarizedExperiment")) {
    if("offsets" %in% names(assays(offsets))) {
      offsets <- assay(offsets,"offsets")
      doOffsets <- TRUE      
    } else {
      stop("No assay named offsets in passes SummarizedExperiment",call.=FALSE)
    }
  }  else {
    stop("offsets must be a logical, matrix, or SummarizedExperiment",call.=FALSE)
  }
  if(doOffsets) {
    if(normalize!=DBA_NORM_OFFSETS_ADJUST) {
      norm$norm.method <- DBA_NORM_OFFSETS
    } else {
      norm$norm.method <- DBA_NORM_OFFSETS_ADJUST
    }
    pv <- pv.normalizeOffsets(pv, offsets=offsets, norm,
                              method=method, bSubControl=bSubControl,
                              filter=filter, filterFun=filterFun,
                              libFun=libFun)  
    return(pv)
  }
  
  if(length(normalize) == ncol(pv$samples)) {
    norm$norm.calc <- "User supplied"
    norm$norm.facs <- normalize
    norm$norm.method <- PV_NORM_USER
  } else if (length(normalize) == 1) {
    if(normalize == PV_NORM_DEFAULT) {
      if(method == DBA_EDGER) {
        normalize <- PV_NORM_TMM
      } else if (method == DBA_DESEQ2) {
        if(norm$lib.calc != "Reads in peaks") {
          normalize <- PV_NORM_LIB        
        } else {
          normalize <- PV_NORM_RLE                  
        }
      }
    } else if (normalize == PV_NORM_NATIVE) {
      if(method == DBA_EDGER) {
        normalize <- PV_NORM_TMM
      } else if (method == DBA_DESEQ2) {
        normalize <- PV_NORM_RLE                  
      }
    }
    
    norm$norm.method <- normalize
    if(normalize == PV_NORM_LIB) {
      norm$norm.calc  <- "Library size"
      norm <- pv.normfacsLIB(pv, norm=norm, method=method,
                             libFun=libFun)
    } else if(normalize == PV_NORM_TMM) {
      norm$norm.calc  <- "edgeR/TMM"
      norm <- pv.normfacsTMM(pv,norm=norm,method=method,
                             bSubControl=bSubControl,
                             filter=filter, filterFun=filterFun,
                             libFun=libFun, binned=binned,
                             background=background)
      
    } else if(normalize == PV_NORM_RLE) {
      norm$norm.calc  <- "DESeq2/RLE"
      norm <- pv.normfacsRLE(pv,norm=norm,method=method,
                             bSubControl=bSubControl,
                             filter=filter, filterFun=filterFun,
                             libFun=libFun, binned=binned,
                             background=background)
    } else {
      stop('Invalid normalization',call.=FALSE)      
    }
  } else {
    stop('normalize invalid length',call.=FALSE)
  }
  
  if(method==DBA_DESEQ2) {
    pv$norm$DESeq2 <- norm
    pv$DESeq2$DEdata <- NULL
    
  } else if(method==DBA_EDGER) {
    pv$norm$edgeR <- norm
    pv$edgeR$DEdata <- NULL
  } 
  
  if(!is.null(binned)){
    pv$norm$background$binned   <- binned
    pv$norm$background$bin.size <- bin.size
  }
  
  return(pv)
}

pv.readsInPeaks <- function(pv, bSubControl=bSubControl,
                            filter=filter, filterFun=filterFun) {
  
  counts <- pv.DEinitedgeR(pv,bSubControl=bSubControl,
                           bRawCounts=TRUE,
                           filter=filter,filterFun=filterFun)
  return(colSums(counts))
}

pv.normfacsTMM <- function(pv,norm,method,bSubControl=FALSE,
                           filter=0, filterFun=max,libFun=PV_NORM_LIBFUN,
                           binned=NULL, background=FALSE) {
  
  if(background[1] != FALSE) {
    if(!is.null(binned)) {# TMM on Background bins
      binned$totals <- norm$lib.sizes
      norm$norm.facs <- csaw::normFactors(binned,method="TMM", se.out=FALSE)
    }  else {
      stop("No binned counts for background TMM",call.=FALSE)
    }
    norm$background <- TRUE
  } else { # edgeR TMM
    
    edger <- pv.DEinitedgeR(pv,bSubControl=bSubControl,
                            bFullLibrarySize=norm$lib.sizes,
                            filter=filter,filterFun=filterFun)
    edger <- edgeR::calcNormFactors(edger,method="TMM", 
                                    lib.size=norm$lib.sizes, doWeighting=FALSE)
    norm$norm.facs <- edger$samples$norm.factors
    
  }
  
  if(method==DBA_DESEQ2) { # Convert from edgeR to DESeq2 factors
    norm$norm.facs <- pv.edgeRtoDESeq2norm(norm, libFun=libFun)
  }
  
  return(norm)
}

pv.normfacsRLE <- function(pv,norm,method,bSubControl=FALSE,
                           filter=0, filterFun=max,libFun=PV_NORM_LIBFUN,
                           binned=NULL, background=FALSE){
  
  if(background[1] != FALSE) {
    if(!is.null(binned)) {# RLE on Background bins
      binned$totals <- norm$lib.sizes
      norm$norm.facs <- csaw::normFactors(binned,method="RLE", se.out=FALSE)
      if(method==DBA_DESEQ2) {
        norm$norm.facs <- pv.edgeRtoDESeq2norm(norm, libFun=libFun)
      }
    } else {
      stop("No binned counts for background RLE",call.=FALSE)
    }
    norm$background <- TRUE
  } else if(method == DBA_DESEQ2 && norm$lib.calc == "Reads in peaks") {
    
    deseq <- pv.DEinitDESeq2(pv,bSubControl=bSubControl,
                             bFullLibrarySize=FALSE,
                             filter=filter,filterFun=filterFun)
    norm$norm.facs <- DESeq2::sizeFactors(deseq)
  } else {
    edger <- pv.DEinitedgeR(pv,bSubControl=bSubControl,
                            bFullLibrarySize=norm$lib.sizes,
                            filter=filter,filterFun=filterFun)
    edger <- edgeR::calcNormFactors(edger,method="RLE",
                                    lib.size=norm$lib.sizes, doWeighting=FALSE)
    norm$norm.facs <- edger$samples$norm.factors
    
    if(method==DBA_DESEQ2) {
      norm$norm.facs <- pv.edgeRtoDESeq2norm(norm, libFun=libFun)
    }
  }
  
  return(norm)
}

pv.normfacsLIB <- function(pv, norm=norm, method=method,
                           libFun=libFun){
  
  norm$lib.fun <- libFun
  
  if(method==DBA_DESEQ2) {
    norm$norm.facs  <- norm$lib.sizes/libFun(norm$lib.sizes)
  } else if(method==DBA_EDGER) {
    norm$norm.facs  <- 1/(norm$lib.sizes/libFun(norm$lib.sizes))
    norm$norm.facs <- pv.makeProd1(norm$norm.facs)
  }
  
  return(norm)
}

pv.makeProd1 <- function(f){
  return(f/exp(mean(log(f))))
}

pv.edgeRtoDESeq2norm <- function(norm, libFun) {
  efflib <- norm$norm.facs * norm$lib.sizes
  norm.facs <- efflib / libFun(norm$lib.sizes)
  return(norm.facs)
}

pv.getBackground <- function(pv,background=PV_BACKGROUND_BINSIZE) {
  
  if (!requireNamespace("csaw",quietly=TRUE)) {
    stop("Package csaw not installed",call.=FALSE)
  }
  
  if(is(background,"logical")) {
    if(!is.null(pv$norm$background)) {
      background <- pv$norm$background$bin.size
    } else {
      background <- PV_BACKGROUND_BINSIZE
    }
  } 
  
  binned <- NULL
  
  if(!is.null(pv$norm$background)) {
    if(pv$norm$background$bin.size==background) {
      binned <- pv$norm$background$binned
    }
  }
  
  if(is.null(binned)) {
    message("Generating background bins with csaw...")
    rParams <- pv.readParams(pv)
    if(pv$config$RunParallel) {
      if(pv$config$parallelPackage == DBA_PARALLEL_MULTICORE) {
        if(is.null(pv$config$cores)) {
          cores <- BiocParallel::multicoreWorkers()
        } else {
          cores <- pv$config$cores
        }
        mcparam <- BiocParallel::MulticoreParam(workers=cores)
      }
    } else {
      mcparam <- BiocParallel::SerialParam()
    }
    binned <- suppressWarnings(suppressMessages(
      csaw::windowCounts(pv$class[PV_BAMREADS,], bin=TRUE,
                         width=background, param=rParams,
                         BPPARAM=mcparam)))
  }
  
  return(list(binned=binned,bin.size=background))
}

pv.readParams <- function(pv) {
  pe <- "none"
  if(!is.null(pv$config$singleEnd)) {
    if(!pv$config$singleEnd) {
      pe <- "both"
    }
  }
  minq <- 15
  if(!is.null(pv$config$minQCth)) {
    minq <- pv$config$minQCth
  }
  rp <- csaw::readParam(pe=pe,minq=minq,restrict=pv$chrmap)
  
  return(rp)
}

pv.normalizeOffsets <- function(pv, offsets=offsets, norm,
                                method=method, bSubControl=bSubControl,
                                filter=filter, filterFun=filterFun, libFun,
                                ...) {
  
  if (!requireNamespace("csaw",quietly=TRUE)) {
    stop("Package csaw not installed",call.=FALSE)
  }
  
  counts <- pv.DEinitedgeR(pv,bSubControl=bSubControl,
                           bRawCounts=TRUE,
                           filter=filter,filterFun=filterFun)
  
  offset.meth <- PV_OFFSETS_LOESS
  
  if(is(offsets,"matrix")) {
    if(sum(dim(counts) != dim(offsets))) {
      stop("offsets matrix must have same dimensions as binding matrix (after filtering).",
           call.=FALSE)
    }
    offset.meth <- PV_OFFSETS_USER
    offsets <- SummarizedExperiment(list(offsets=offsets))
  } else if(!is.null(pv$norm$offsets)) {
    if (pv$norm$offsets$offset.method == PV_OFFSETS_LOESS) {
      offsets <- pv$norm$offsets$offset
    }
  }
  if(!is(offsets,"SummarizedExperiment")) {
    counts <- SummarizedExperiment(list(counts=counts))
    counts$totals <- norm$lib.sizes
    offsets <- csaw::normOffsets(counts, se.out=FALSE, ...)
    offsets <- SummarizedExperiment(list(offsets=offsets))
    counts$totals <- norm$lib.sizes
  }
  
  if(offset.meth == PV_OFFSETS_LOESS) {
    norm$norm.method <- PV_NORM_OFFSETS_ADJUST
  }
  
  norm$norm.calc   <- "Use offsets"
  norm$lib.fun     <- libFun
  
  pv <- pv.setNorm(pv, norm, method)
  
  pv$norm$offsets$offsets <- offsets
  pv$norm$offsets$offset.method <- offset.meth
  
  return(pv)
}

pv.offsetsAdjust <- function(pv, offsets, deobj) {
  offsets <- edgeR::scaleOffset(assay(deobj),offsets)
  offsets <- offsets / exp(rowMeans(log(offsets)))
  libs <- pv$norm$DESeq2$lib.sizes 
  nf <- libs/pv$norm$DESeq2$lib.fun(libs)
  nfs <- matrix(nf, nrow(offsets), ncol(offsets), byrow=TRUE)
  offsets <- offsets * nfs 
  return(offsets)
}

pv.getNorm <- function(pv,method) {
  if(method==DBA_EDGER) {
    return(pv$norm$edgeR)
  } else if (method==DBA_DESEQ2) {
    return(pv$norm$DESeq2)
  } else {
    stop("Internal error: Invalid method.")
  }
}

pv.setNorm <- function(pv,norm,method) {
  if(method==DBA_EDGER) {
    pv$norm$edgeR <- norm
    pv$edgeR$DEdata <- NULL
  } else if (method==DBA_DESEQ2) {
    pv$norm$DESeq2 <- norm
    pv$DESeq2$DEdata <- NULL
  } else {
    stop("Internal error: Invalid method.")
  }
  return(pv)
}

pv.reNormalize <- function(pv) {
  
  if(is.null(pv$norm)) {
    return(pv)
  }
  
  if(!is.null(pv$norm$offsets)){
    if(pv$norm$offsets$offset.method == PV_OFFSETS_USER) {
      warning("Re-run dba.normalize() to add user-supplied offsets.",
              call.=FALSE)
      pv$norm$offsets$offsets <- NULL
    }
  }
  
  if(!is.null(pv$norm$DESeq2) || !is.null(pv$norm$edgeR)) {
    message("Re-normalizing...")
    pv$norm$DESeq2 <- pv.doRenormalize(pv,DBA_DESEQ2)
    pv$norm$edgeR  <- pv.doRenormalize(pv,DBA_EDGER)
  }
  
  return(pv)
}

pv.doRenormalize <- function(pv, method) {
  
  if(method == DBA_EDGER) {
    norm <- pv$norm$edgeR
  }
  
  if(method == DBA_DESEQ2) {
    norm <- pv$norm$DESeq2
  }
  
  if(is.null(norm)) {
    return(NULL)
  }
  
  if(norm$lib.method==PV_LIBSIZE_USER) {
    libsizes <- norm$lib.sizes
  } else {
    libsizes <- norm$lib.method
  }
  
  if(norm$norm.method==PV_LIBSIZE_USER) {
    normfacs <- norm$norm.facs
  } else {
    normfacs <- norm$norm.method
  }
  
  if(is.null(norm$filter.fun)) {
    filterfun <- max
  } else {
    filterfun <- norm$filter.fun
  }
  
  if(is.null(norm$lib.fun)) {
    libfun <- PV_NORM_LIBFUN
  } else {
    libfun <- norm$lib.fun
  }
  
  offsets <- FALSE
  if(norm$method == PV_NORM_OFFSETS) {
    if(is.null(pv$norm$offsets)) {
      offsets <- TRUE
    } else if(pv$norm$offsets$offset.method == PV_OFFSETS_USER) {
      pv <- pv.setNorm(pv, NULL, method)
      return(pv)
    } else {
      offsets <- TRUE
    }
  }
  
  pv <- pv.normalize(pv,
                     method    = method,
                     libSizes  = libsizes,
                     normalize = normfacs,
                     bSubControl = norm$bSubControl,
                     filter=norm$filter.val, filterFun=filterfun, libFun=libfun, 
                     background=norm$background, offsets=offsets)
  return(pv)
}

pv.edgeRCounts <- function(pv,method,bNormalized=TRUE) {
  
  if(is.null(pv$edgeR$DEdata)) {
    stop('No edgeR data object -- re-run dba.analyze with method=DBA_EDGER',
         call.=FALSE)
  } 
  
  counts <- pv$edgeR$DEdata$counts
  
  if(bNormalized) {
    if(is.null(pv$edgeR$DEdata$offset)) {
      samples <- pv$edgeR$DEdata$samples
      eff.lib <- samples$lib.size * samples$norm.factors
      counts <- edgeR::cpm(counts, normalized.lib.sizes=TRUE, 
                           lib.size=pv$edgeR$DEdata$samples)
    } else {
      counts <- edgeR::cpm(counts, normalized.lib.sizes=TRUE, 
                           lib.size=pv$edgeR$DEdata$samples$lib.size,
                           offset=pv$edgeR$DEdata$offset)
    }
  }
  
  return(counts)
}

pv.normalizeRetrieve <- function(pv, method) {
  res <- NULL
  
  if(is.null(pv$norm)) {
    return(NULL)
  }
  
  if(all(method==DBA_DESEQ2)) {
    res <- pv.formatNorm(pv$norm$DESeq2)
  }
  if(all(method==DBA_EDGER)) {
    res <- pv.formatNorm(pv$norm$edgeR)
  }
  if(all(method==DBA_ALL_METHODS)) {
    res$edgeR      <- pv.formatNorm(pv$norm$edgeR)
    res$DESeq2     <- pv.formatNorm(pv$norm$DESeq2)
    res$background <- pv$norm$background
    res$offsets    <- pv$norm$offsets
  }    
  
  return(res)
}

pv.formatNorm <- function(norm) {
  res <- NULL
  
  if(norm$background) {
    res$background <- norm$background
  }
  
  res$norm.method    <- norm$norm.method
  res$norm.factors   <- norm$norm.facs
  res$lib.method     <- norm$lib.method
  res$lib.sizes      <- norm$lib.sizes
  
  if(norm$bSubControl) {
    res$control.subtract <- norm$bSubControl
  }
  if(norm$filter.val > 0) {
    res$filter.value <- norm$filter.val
  }
  
  return(res)
}


