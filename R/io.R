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
   }
   return(DBAobject)
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
