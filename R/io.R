pv.save <- function(DBAobject,file='model',dir='Robjects',pre='pv_',ext='RData',
                    compress=TRUE,ascii=FALSE) {
   fn <- sprintf('%s/%s%s.%s',dir,pre,file,ext)
   
   if(is(compress,"logical")) {
      if(compress==TRUE) {
         compress_level <- 9
      } 
   } else {
      compress_level <- compress
      compress <- TRUE
   }
   
   save(DBAobject,file=fn,compress=compress,
        compression_level=compress_level,ascii=ascii)
   return(fn)
}

pv.load <- function(file='model',dir='Robjects',pre='pv_',ext='RData') {
   DBAobject <- NULL
   pv <- NULL
   load(sprintf('%s/%s%s.%s',dir,pre,file,ext))
   if(is.null(DBAobject)) {
      DBAobject <- pv
   } else {
      if(!is.null(DBAobject$config$Version1)) {
         if(as.numeric(DBAobject$config$Version1) < 3) {
            if(as.numeric(DBAobject$config$Version2) < 99) {
               DBAobject <- pv.loadPre3(DBAobject)
            }
         }
      }
   }
   return(DBAobject)
}


pv.loadPre3 <- function(pv) {
   
   # Make sure it is a DBA object
   if(!is(pv,"DBA")) {
      class(pv) <- "DBA"
   }
   
   # If it already has normalization don't do anything
   if(!is.null(pv$norm)) {
      return(pv)
   }
   
   # If no counts available don't do anything
   srcmask <- pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts")
   if(sum(srcmask)==0) {
      return(pv)
   }
   
   # Set minimum count to 1 (default is 0 in 3.0)
   if (is.null(pv$minCount)) {
      pv$minCount <- 1
   }
   
   # If an analysis has been run,  check for normalization parameters
   bSubControlD <- bSubControlE <- TRUE
   bFullLibrarySizeD <- bFullLibrarySizeE <- TRUE
   con <- pv$contrasts[[1]]
   if(!is.null(con)) {
      if(!is.null(con$DESeq2)) {
         bSubControlD <- con$DESeq2$bSubControl
         bFullLibrarySizeD <- con$DESeq2$bFullLibrarySize 
      }
      if(!is.null(con$edgeR)) {
         bSubControlE <- con$edgeR$bSubControl
         bFullLibrarySizeE <- con$edgeR$bFullLibrarySize 
      }
   }
   # Set mode to pre-2.0
   if(is.null(pv$design)) {
      pv$design <- FALSE
   }
   
   # Normalize
   filtval <- 0
   filtFun <- max
   if(!is.null(pv$filter)) {
      filtval <- pv$filter
   }
   if(!is.null(pv$filterFun)) {
      filtFun <- pv$filterFun
   }
   if(bFullLibrarySizeD) {
      pv <- dba.normalize(pv, method = DBA_DESEQ2, 
                          normalize = DBA_NORM_LIB,
                          library = DBA_LIBSIZE_FULL, 
                          bSubControl = bSubControlD, 
                          filter = filtval,filterFun = filtFun)
   } else {
      pv <- dba.normalize(pv, method = DBA_DESEQ2, 
                          normalize = DBA_NORM_RLE,
                          library = DBA_LIBSIZE_PEAKREADS, 
                          bSubControl = bSubControlD,
                          filter = filtval,filterFun = filtFun)      
   }
   
   if(bFullLibrarySizeE) {
      pv <- dba.normalize(pv, method = DBA_EDGER, 
                          normalize = DBA_NORM_TMM,
                          library = DBA_LIBSIZE_FULL, 
                          bSubControl = bSubControlE, 
                          filter = filtval,filterFun = filtFun)
   } else {
      pv <- dba.normalize(pv, method = DBA_EDGER,
                          normalize = DBA_NORM_RLE,
                          library = DBA_LIBSIZE_PEAKREADS, 
                          bSubControl = bSubControlE,
                          filter = filtval,filterFun = filtFun)      
   }
   return(pv)
}


## pv.writePeakset --- write out vectorized peaks as a bed file for external 
pv.writePeakset <- function(pv,fname,peaks,numCols=4){
   
   if(missing(peaks)) {
      peaks <- rep(T,nrow(pv$binding))
   } else {
      if(sum(class(peaks)=='logical')) {
         peaks <- which(peaks)[1]	
      }	
   }
   
   if(missing(fname)) {
      fname <- NULL
   }
   
   if(sum((class(peaks)=='numeric')) || sum((class(peaks)=='integer'))) {
      peaks=pv$peaks[[peaks]]
   }           
   
   if(!is.null(dim(peaks))) {
      if(is(peaks[1,1],"character")) {
         bed <- pv.do_peaks2bed(peaks,NULL,fname,numCols=numCols)
      } else {
         bed <- pv.do_peaks2bed(peaks,pv$chrmap,fname,numCols=numCols)
      }
   } else {
      bed <- pv.do_peaks2bed(pv$binding,pv$chrmap,fname,numCols=ncol(pv$binding))
   }
   
   return(bed)
}
