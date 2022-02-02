#####################################
## pv_helper.R -- support for pv   ##
## 20 October 2009                 ##
## 3 February 2011 -- packaged     ##
## Rory Stark                      ##
## Cancer Research UK              ##
#####################################

#########################
## pv HELPER FUNCTIONS ##
#########################

pv.check <- function(pv,bCheckEmpty=FALSE,bCheckSort=TRUE,bDoVectors=TRUE) {
   
   if(missing(pv)) {
      stop('DBA object missing!',call.=FALSE)
   }
   
   if(is.null(pv)) {
      return(NULL)	
   }
   
   if(is(pv,"ChIPQCexperiment")){
      saveqc <- pv
      pv <- pv@DBA
      x <- "DBA object placeholder"
      class(x) <- "DBA"
      saveqc@DBA <- x
      pv$ChIPQCobj <- saveqc
   }
   
   if(is.null(pv$peaks)) {
      return(pv)	
   }
   
   if(bCheckEmpty && sum(sapply(pv$peaks,nrow))==0) {
      stop('DBA object has no peaks.',call.=FALSE)
   }
   if(!is.null(pv$vectors)) {
      if(!is.null(pv$allvectors)) {
         pv$binding     <- pv$vectors
         pv$totalMerged <- nrow(pv$allvectors)
         pv$merged      <- pv$allvectors[,1:3]
      } else {
         pv$binding     <- pv$vectors
         pv$merged      <- pv$vectors[,1:3]
         pv$totalMerged <- nrow(pv$merged)
      }
      pv$vectors    <- NULL
      pv$allvectors <- NULL
   }
   
   if(is.null(pv$binding)) {
      if(is.null(pv$minOverlap)) {
         minOverlap=2
      } else {
         minOverlap <- pv$minOverlap
      }
      if(bDoVectors) {
         contrasts  <- pv$contrasts
         called     <- pv$called
         allcalled  <- pv$allcalled
         attributes <- pv$attributes
         for(i in 1:length(pv$peaks)) {
            if(is.factor(pv$peaks[[i]][,1])) {
               pv$peaks[[i]][,1] <- as.character(pv$peaks[[i]][,1])
            }
         }
         pv <- pv.vectors(pv,minOverlap=minOverlap,bAllSame=pv.allSame(pv),
                          merge=is.null(pv$merged))
         pv$contrasts  <- contrasts
         pv$called     <- called
         pv$allcalled  <- allcalled
         pv$attributes <- attributes
      }
   }
   
   pv <- pv.checkCalled(pv)
   
   if(is.null(pv$config$DataType)) {
      pv$config$DataType=DBA_DATA_DEFAULT
      if(!is.null(pv$config$RangedData)) {
         if(pv$config$RangedData==FALSE) {
            pv$config$DataType <- DBA_DATA_FRAME   	
         } 
      }
   }
   
   if(is.null(pv$config$bCorPlot)) {
      pv$config$bCorPlot <- FALSE
   }
   
   if(is(pv$attributes,"function")) {
      pv$attributes <- NULL
   }
   
   if(is.null(pv$config$th)){
      pv$config$th <- 0.05
   }
   if(is.null(pv$config$bUsePval)){
      pv$config$bUsePval <- FALSE
   }   
   
   if(nrow(pv$class)<PV_TREATMENT) {
      pv$class <- rbind(pv$class,'')
      rownames(pv$class)[PV_TREATMENT]='Treatment'	
   }
   
   pv$config <- as.list(pv$config)
   
   if (is.null(pv$config$mapQCth)) {
      pv$config$mapQCth <- 15   
   }
   
   if (is.null(pv$config$fragmentSize)) {
      pv$config$fragmentSize <- 125
   }   
   
   if(bCheckSort && is.unsorted(pv$chrmap)) {
      message("Older version of DBA object. Recommend running DBA <- dba(DBA).")
      message('Converting...')
      pv <- dba(pv)
   }
   
   if(nrow(pv$class) < PV_SPIKEIN) {
      pv$class <- rbind(pv$class, Spikein=NA)
   }
   
   return(pv)
}

pv.version <- function(pv,v1,v2,v3){
   
   warn <- FALSE
   if(is.null(pv$config$Version1)) {
      warn <- TRUE   	
   } else {
      if(pv$config$Version1 != v1) {
         warn <- TRUE   	
      }	
   }
   if(is.null(pv$config$Version2)) {
      warn <- TRUE   	
   } else {
      if(pv$config$Version2 != v2) {
         warn <- TRUE   	
      }	
   }
   
   if(warn) {
      warning('Loading DBA object from a previous version -- updating...',call.=FALSE)
   }
   
   if(!is.null(pv$contrasts)) {
      for(i in 1:length(pv$contrasts)) {
         if(!is.null(pv$contrasts[[i]]$DESeq) & is.null(pv$contrasts[[i]]$DESeq1)
            & is.null(pv$contrasts[[i]]$DESeq2)) {
            pv$contrasts[[i]]$DESeq1 <- pv$contrasts[[i]]$DESeq
            pv$contrasts[[i]]$DESeq  <- NULL
         }
      }   
   }
   
   if(nrow(pv$class)<PV_TREATMENT) {
      pv$class <- rbind(pv$class,'')
      rownames(pv$class)[PV_TREATMENT]='Treatment'	
   }
   
   if(is.unsorted(pv$chrmap)) {
      pv <- dba(pv)
   }
   
   pv$config$Version1 <- v1
   pv$config$Version2 <- v2
   pv$config$Version3 <- v3   	
   
   return(pv)
}

checkQCobj <- function(resQC,res) {
   ## add better checks that objects are in sync
   sampnames <- unique(c(unique(res$class[PV_ID,]),unique(res$class[PV_CONTROL,])))
   sampnames <- sampnames[!is.na(sampnames)]
   sampnames <- sampnames[sampnames!=""]
   if (sum(sampnames %in% names(resQC@Samples))==length(sampnames)) {
      if(length(sampnames) == length(names(resQC@Samples))) {
         resQC@DBA <- res
      } else resQC=NULL
   } else resQC= NULL
   if(is.null(resQC)) {
      warning("ChIPQexperiment out of sync -- returning new DBA object",call.=FALSE)
      return(res)
   } else {
      return(resQC)
   }
}

pv.whichPeaksets <- function(pv,mask) {
   
   if(missing(mask)) {
      warning('mask required',call.=FALSE)
      return(NULL)
   }
   if(is.null(mask)) {
      warning('mask required',call.=FALSE)
      return(NULL)
   }
   if(is(mask,'logical')) {
      mask <- which(mask)
   }
   
   A <- mask[1]
   B <- mask[2]
   if(length(mask) >= 3) {
      C <- mask[3]
   } else {
      C <- NULL
   }
   if(length(mask) >= 4) {
      D <- mask[4]
   } else {
      D <- NULL
   }   
   
   return(list(A=A,B=B,C=C,D=D))
}

pv.listadd <- function(a,b){
   b <- list(b)
   if (is.null(a)) return(b)
   return(c(a,b))
}

pv.listaddto <- function(a,b){
   if (is.null(a)) return(b)
   return(c(a,b))
}

fdebug <- function(str,file='debug.txt'){
   
   PV_DEBUG=FALSE
   
   if(PV_DEBUG == FALSE){
      return()
   }
   
   #write(sprintf('%s\n',str),file=file,append=TRUE)
   
}

pv.peaksort <- function(peaks,chrmap) {
   if(is.character(peaks[1,1])) {
      chrs <- TRUE
      if(missing(chrmap)) {
         chrmap <- sort(unique(peaks[,1]))
      }
      peaks[,1] <- match(peaks[,1],chrmap)
   } else chrs=FALSE
   
   if(nrow(peaks) < 2) {
      return(peaks)
   }
   
   o <- peakOrder(peaks[,1],as.integer(peaks[,2]),as.integer(peaks[,3]))
   peaks <- peaks[o,]
   if(chrs) {
      peaks[,1] <- chrmap[peaks[,1]]
   }
   pv.gc()
   return(peaks)
}

pv.contrast2 <- function(peaks,A,B,v1,v2){
   
   allpeaks <- v1 | v2
   inAll    <- v1 & v2
   onlyA    <- v1 & !v2
   onlyB    <- !v1 & v2
   
   res <- list(onlyA=data.frame(peaks[onlyA,c(1:3,(3+A))]),
               onlyB=data.frame(peaks[onlyB,c(1:3,(3+B))]),
               inAll=data.frame(peaks[inAll,c(1:3,(3+A),(3+B))]))
   
   for(i in 1:3){
      if(is.null(nrow(res[[i]]))){
         res[[i]] <- data.frame(t(res[[i]]))
      }
   }
   
   colnames(res[[1]]) <- c("chr","start","end","score")	
   colnames(res[[2]]) <- c("chr","start","end","score")	
   colnames(res[[3]]) <- c("chr","start","end","scoreA","scoreB")	
   
   return(res)   	
}

pv.contrast3 <- function(peaks,A,B,C,v1,v2,v3){
   
   allpeaks <- v1 | v2 | v3
   
   inAll  <-  v1 &  v2 &  v3
   onlyA  <-  v1 & !v2 & !v3
   onlyB  <- !v1 &  v2 & !v3 
   onlyC  <- !v1 & !v2 &  v3
   
   notA   <- !v1 &  (v2 & v3)
   notB   <- !v2 &  (v1 & v3)
   notC   <- !v3 &  (v1 & v2)    
   
   res <- list(onlyA=data.frame(peaks[onlyA,c(1:3,(3+A))]),
               onlyB=data.frame(peaks[onlyB,c(1:3,(3+B))]),
               onlyC=data.frame(peaks[onlyC,c(1:3,(3+C))]),
               notA =data.frame(peaks[notA,c(1:3,(3+B),(3+C))]),              
               notB =data.frame(peaks[notB,c(1:3,(3+A),(3+C))]), 
               notC =data.frame(peaks[notC,c(1:3,(3+A),(3+B))]), 
               inAll=data.frame(peaks[inAll,c(1:3,(3+A),(3+B),(3+C))]))
   
   for(i in 1:7){
      if(is.null(nrow(res[[i]]))){
         res[[i]] <- data.frame(t(res[[i]]))
      }
   }	
   
   colnames(res[[1]]) <- c("chr","start","end","score")	
   colnames(res[[2]]) <- c("chr","start","end","score")	
   colnames(res[[3]]) <- c("chr","start","end","score")	
   colnames(res[[4]]) <- c("chr","start","end","scoreB","scoreC")
   colnames(res[[5]]) <- c("chr","start","end","scoreA","scoreC")
   colnames(res[[6]]) <- c("chr","start","end","scoreA","scoreB")
   colnames(res[[7]]) <- c("chr","start","end","scoreA","scoreB","scoreC")
   
   return(res)
}

pv.contrast4 <- function(peaks,A,B,C,D,v1,v2,v3,v4){
   
   allpeaks <- v1 | v2 | v3 | v4
   
   inAll  <-  v1 &  v2 &  v3 & v4
   onlyA  <-  v1 & !v2 & !v3 & !v4
   onlyB  <- !v1 &  v2 & !v3 & !v4
   onlyC  <- !v1 & !v2 &  v3 & !v4
   onlyD  <- !v1 & !v2 & !v3 & v4
   
   notA   <- !v1 &  (v2 & v3 & v4)
   notB   <- !v2 &  (v1 & v3 & v4)
   notC   <- !v3 &  (v1 & v2 & v4)
   notD   <- !v4 &  (v1 & v2 & v3)      
   
   AandB  <- (v1 & v2) & !(v3 | v4)
   AandC  <- (v1 & v3) & !(v2 | v4)
   AandD  <- (v1 & v4) & !(v2 | v3)
   BandC  <- (v2 & v3) & !(v1 | v4)
   BandD  <- (v2 & v4) & !(v1 | v3)
   CandD  <- (v3 & v4) & !(v1 | v2)     
   
   res <- list(onlyA=data.frame(peaks[onlyA,c(1:3,(3+A))]),
               onlyB=data.frame(peaks[onlyB,c(1:3,(3+B))]),
               onlyC=data.frame(peaks[onlyC,c(1:3,(3+C))]),
               onlyD=data.frame(peaks[onlyD,c(1:3,(3+D))]),
               AandB=data.frame(peaks[AandB,c(1:3,(3+A),(3+B))]),
               AandC=data.frame(peaks[AandC,c(1:3,(3+A),(3+C))]),
               AandD=data.frame(peaks[AandD,c(1:3,(3+A),(3+D))]),
               BandC=data.frame(peaks[BandC,c(1:3,(3+B),(3+C))]),
               BandD=data.frame(peaks[BandD,c(1:3,(3+B),(3+D))]),
               CandD=data.frame(peaks[CandD,c(1:3,(3+C),(3+D))]),
               notA =data.frame(peaks[notA, c(1:3,(3+B),(3+C),(3+D))]),              
               notB =data.frame(peaks[notB, c(1:3,(3+A),(3+C),(3+D))]), 
               notC =data.frame(peaks[notC, c(1:3,(3+A),(3+B),(3+D))]),
               notD =data.frame(peaks[notD, c(1:3,(3+A),(3+B),(3+C))]), 
               inAll=data.frame(peaks[inAll,c(1:3,(3+A),(3+B),(3+C),(3+D))]))
   
   for(i in 1:15){
      if(is.null(nrow(res[[i]]))){
         res[[i]] <- data.frame(t(res[[i]]))
      }
   }	
   
   colnames(res[[1]])  <- c("chr","start","end","score")	
   colnames(res[[2]])  <- c("chr","start","end","score")	
   colnames(res[[3]])  <- c("chr","start","end","score")	
   colnames(res[[4]])  <- c("chr","start","end","score")
   colnames(res[[5]])  <- c("chr","start","end","scoreA","scoreB")
   colnames(res[[6]])  <- c("chr","start","end","scoreA","scoreC")
   colnames(res[[7]])  <- c("chr","start","end","scoreA","scoreD")
   colnames(res[[8]])  <- c("chr","start","end","scoreB","scoreC")   
   colnames(res[[9]])  <- c("chr","start","end","scoreB","scoreD") 
   colnames(res[[10]]) <- c("chr","start","end","scoreC","scoreD")
   colnames(res[[11]]) <- c("chr","start","end","scoreB","scoreC","scoreD")    
   colnames(res[[12]]) <- c("chr","start","end","scoreA","scoreC","scoreD")       
   colnames(res[[13]]) <- c("chr","start","end","scoreA","scoreB","scoreD")    
   colnames(res[[14]]) <- c("chr","start","end","scoreA","scoreB","scoreC")      
   colnames(res[[15]]) <- c("chr","start","end","scoreA","scoreB","scoreC","scoreD")
   
   return(res)
}


pv.analysis <- function(pv,attributes=pv$attributes,bPCA=TRUE,distMeth="pearson") {
   
   #require(amap)
   
   peaks  <- pv$binding 
   values <- matrix(as.numeric(as.matrix(peaks[,4:ncol(peaks)])),nrow(peaks),ncol(peaks)-3)
   
   if(sum(pv$binding[,4:ncol(pv$binding)] != 1) == 0) {
      return(pv)
   }
   
   if(sum(pv$binding[,4:ncol(pv$binding)] != -1) == 0) {
      return(pv)
   }
   
   if(nrow(peaks) <= length(pv$peaks) ) {
      bPCA <- FALSE
   }
   #values <- unique(values)
   cnames=NULL
   for(i in 1:ncol(pv$class)) {
      cnames <- c(cnames,pv.namestrings(pv$class[attributes,i])$tstring)
   }
   colnames(values) <- cnames
   x <- apply(values,1,pv.howmany)
   values <- values[x>1,]
   #pv$hc <- hclust(Dist(t(values),method=distMeth))
   if(bPCA) pv$pc <- princomp(values)
   #pv$values <- values
   
   return(pv)  
}

pv.howmany <- function(vals){
   return(sum(vals>0))	
}


pv.readPeaks <- function(peaks,peak.format,skipLines=0){
   if(peak.format == "macs") {
      peaks <- pv.macs(peaks)
   } 
   else if(peak.format == "bayes") {
      peaks <- pv.bayes(peaks)
   }
   else if(peak.format == "swembl") {
      peaks <- pv.swembl(peaks)
   } 
   else if(peak.format == "raw") {
      peaks <- pv.readbed(peaks,checkHeader=TRUE)
   } 
   else if(peak.format == "fp4") {
      peaks <-  pv.readbed(peaks,1)
   } 
   else if(peak.format == "bed") {
      peaks <-  pv.readbed(peaks,checkHeader=TRUE)
   } 
   else if(peak.format == "tpic") {
      peaks <-  pv.tpic(peaks)
   } 
   else if(peak.format == "sicer") {
      peaks <-  pv.sicer(peaks)
   }    
   else if(peak.format == "narrow") {
      peaks <-  pv.readbed(peaks,skipLines)
   } 
   #else if(peak.format == "raw") {
   #   peaks <-  pv.readbed(peaks,skipLines)
   #}  
   else if(peak.format == "csv") {
      peaks <-  pv.csv(peaks)
   }  else if(peak.format == "report") {
      peaks <-  pv.csv(peaks)
   } else {
      peaks <-  pv.readbed(peaks,skipLines)      
   }     
}

pv.defaultScoreCol <- function(peak.format){
   if(is.null(peak.format)) {
      return(0)
   }
   if(peak.format == "macs") {
      val <- 7
   } 
   else if(peak.format == "bayes") {
      val <- 0
   }
   else if(peak.format == "swembl") {
      val <- 4
   } 
   else if(peak.format == "raw") {
      val <- 4
   } 
   else if(peak.format == "fp4") {
      val <- 5
   } 
   else if(peak.format == "bed") {
      val <- 5
   } 
   else if(peak.format == "tpic") {
      val <- 0
   } 
   else if(peak.format == "sicer") {
      val <- 7
   }    
   else if(peak.format == "narrow") {
      val <- 8
   } 
   else if(peak.format == "raw") {
      val <- 4
   } else if(peak.format == "csv") {
      val <- 4
   } else if(peak.format == "report") {
      val <- 9
   } else {
      val <- 4 	
   } 
   return(val)    
}


FDRth=100
pv.macs <- function(fn){
   data <- read.table(fn,blank.lines.skip=TRUE,header=TRUE)
   res  <- pv.peaksort(data)
   return(res)
}

pv.swembl <- function(fn){
   data <- read.table(fn,skip=14)
   res  <- pv.peaksort(data)
   return(res) 
}

hasHeader <- function(fn) {
   line <- readLines(fn,n=1)
   flds <- strsplit(line,"\\s")
   f2 <- suppressWarnings(as.numeric(flds[[1]][2]))
   f3 <- suppressWarnings(as.numeric(flds[[1]][3]))
   return(is.na(f2) || is.na(f3))
}

pv.readbed <- function(fn,skipnum=0,checkHeader=FALSE){
   if (checkHeader) {
      skipnum <- if (hasHeader(fn)) 1 else 0
   }
   data <- read.table(fn,skip=skipnum)
   res  <- pv.peaksort(data)
   if(ncol(res)==3) {
      res=cbind(res,1)   
   }
   return(res)
}

pv.bayes <- function(fn){
   data <- read.table(fn)
   #idx <- data[,4]>0.5
   #data <- data[idx,]
   res  <- pv.peaksort(data)
   return(res)
}

pv.tpic <- function(fn){
   data <- read.table(fn)
   res  <- pv.peaksort(data)
   return(cbind(res,1))
}

pv.sicer <- function(fn){
   data <- read.table(fn)
   res  <- pv.peaksort(data)
   return(res[,c(1:3,7)]) 
}

# pv.sourcedata <- function(fn,maxval){
#    data <- read.table(fn)
#    vals <- data[,6]
#    if(!missing(maxval)) {
#       vals[vals>maxval]=maxval
#    }
#    #vals <- vals/100
#    vals <- log2(vals)
#    vals[vals<0]=0
#    data <- cbind(data[,1:3],vals,data[,4:5])
#    data <- pv.peaksort(data)
#    return(data)
# }

pv.csv <- function(fn){
   data <- read.csv(fn)
   res  <- pv.peaksort(data)
   if(ncol(res)==3) {
      res=cbind(res,1)	
   }
   return(res)
}

pv.checkCounts <- function(pv,pv2=NULL) {
   callers <- unique(pv$class[PV_CALLER,])
   if(length(callers)>1) {
      return(FALSE)
   }
   if (callers != "counts") {
      return(FALSE)
   }
   if(!is.null(pv2)) {
      if(!pv.checkCounts(pv2)) {
         return(FALSE)
      }
      if(nrow(pv$binding)!=nrow(pv2$binding)) {
         return(FALSE)
      }
   }
   return(TRUE)
}

pv.peakset_all <- function(pv, addpv, minOverlap) {
   
   
   allcounts <- pv.checkCounts(pv,addpv)
   
   for(i in 1:length(addpv$peaks)) {
      
      message(addpv$class[PV_ID,i],' ',
              addpv$class[PV_TISSUE,i],' ',
              addpv$class[PV_FACTOR,i],' ',
              addpv$class[PV_CONDITION,i],' ',
              addpv$class[PV_TREATMENT,i],' ',
              addpv$class[PV_REPLICATE,i],' ',
              addpv$class[PV_CALLER,i])
      
      pv <- pv.peakset(pv,peaks=addpv$peaks[[i]],
                       sampID      = addpv$class[PV_ID,i],
                       tissue      = addpv$class[PV_TISSUE,i],
                       factor      = addpv$class[PV_FACTOR,i],
                       condition   = addpv$class[PV_CONDITION,i],
                       treatment   = addpv$class[PV_TREATMENT,i],
                       replicate   = addpv$class[PV_REPLICATE,i],
                       control     = addpv$class[PV_CONTROL,i],
                       peak.caller = addpv$class[PV_CALLER,i],
                       reads       = addpv$class[PV_READS,i],
                       consensus   = addpv$class[PV_CONSENSUS,i],
                       readBam     = addpv$class[PV_BAMREADS,i],
                       controlBam  = addpv$class[PV_BAMCONTROL,i],
                       spikein     = addpv$class[PV_SPIKEIN,i]
      )
   }
   
   if(minOverlap>0 && minOverlap<1) {
      minOverlap <- ceiling(length(pv$peaks) * minOverlap)	
   }
   if(allcounts) {
      pv <- pv.vectors(pv,minOverlap=minOverlap,
                       attributes=attributes,
                       bAllSame=TRUE)
   } else {
      pv <- dba(pv, minOverlap=minOverlap)
   }
   
   return(pv)
   
}

pv.minOverlap <- function(vec,minval){
   res <- sum(vec != -1)
   if(minval==0) {
      if (length(vec) == res) {
         return(TRUE)
      } else {
         return(FALSE)
      }
   } else {
      if(res >= minval) {
         return(TRUE)
      } else {
         return(FALSE)
      }
   }	
}

pv.countOverlap <- function(vec,minval= -1){
   res <- sum(vec > minval)
   return(res)
}

pv.domean <- function(vals){
   return(mean(vals[vals>0]))	
}

pv.do_peaks2bed <- function(peaks,chrmap=NULL,fn,numCols=4) {
   peaks <- data.frame(peaks)
   if(!is.null(chrmap)) {
      peaks[,1] <- chrmap[peaks[,1]]
   }
   if(numCols>0) {
      numCols <- max(3,numCols)
      mcols <- min(numCols, ncol(peaks))
   } else {
      mcols <- ncol(peaks)
   }
   
   if(!is.null(fn)) {
      ds <- options("scipen")
      options(scipen=8)
      write.table(peaks[,1:mcols],file=fn,#sprintf("%s.bed",fn), 
                  quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
      options(scipen=ds$scipen)
   }
   return(peaks[,1:mcols])
}

pv.catstr <- function(strvec){
   unq <- unique(strvec)
   if(length(unq) == 1){
      return(unq)
   }
   str <- unq[1]
   for(i in 2:length(unq)){
      str <- sprintf("%s-%s",str,unq[i])	
   }
   return(str)	
}

pv.namestrings <- function(crec1,crec2,crec3,crec4) {
   s1 <- NULL
   s2 <- NULL
   s3 <- NULL
   s4 <- NULL
   t1 <- NULL
   if(missing(crec2)) {
      crec2=crec1
   }
   if(missing(crec3)) {
      crec3=crec1   
   }
   if(missing(crec4)) {
      crec4=crec1
   }
   for(i in 1:length(crec1)) {
      if( (crec1[i]==crec2[i]) && (crec1[i]==crec3[i]) && (crec1[i]==crec4[i])) {
         t1 <- pv.addstr(t1,crec1[i])
      } else {
         s1 <- pv.addstr(s1,crec1[i])
         s2 <- pv.addstr(s2,crec2[i])
         s3 <- pv.addstr(s3,crec3[i])
         s4 <- pv.addstr(s4,crec4[i])
      }
   }
   if(is.null(s1)) { s1 <- ""}
   if(is.null(s2)) { s2 <- ""}
   if(is.null(s3)) { s3 <- ""}
   if(is.null(s4)) { s4 <- ""}
   if(is.null(t1)) { t1 <- ""}
   
   return(list(n1=s1,n2=s2,n3=s3,n4=s4,tstring=t1))	
}

pv.addstr <- function(s1,a1) {
   if(is.null(a1)) {
      return(s1)
   } else if (a1 == "") {
      return(s1)
   }
   
   if(a1=="T") {
      a1 <- "TRUE"
   }
   if(a1=="F") {
      a1 <- "FALSE"
   }
   
   if(is.null(s1)){
      s1 <- a1
   } else {
      s1 <- sprintf("%s:%s",s1,a1)
   }
   return(s1)	
}

pv.merge <- function(allpeaks,peaks=NULL,classes,maxgap=-1,
                     useExternal=TRUE,useC=TRUE, defVal=0){
   
   if(!useC) {
      stop('pv.dovectors called with useC set to FALSE.')	
   }
   
   allpeaks <- as.data.frame(allpeaks)
   
   chrmap <- sort(unique(allpeaks[,1]))
   allpeaks[,1] <- match(allpeaks[,1],chrmap)
   colnames(allpeaks) <- c("CHR","START","END")
   allpeaks <- pv.peaksort(allpeaks)
   allpeaks$CHR   <- as.integer(allpeaks$CHR)
   allpeaks$START <- as.integer(allpeaks$START)
   allpeaks$END   <- as.integer(allpeaks$END)
   maxGap <- as.integer(maxgap)
   merged <- mergePeaks(allpeaks,maxGap)
   
   included     <- matrix(defVal,nrow(merged),length(peaks))
   result       <- matrix(defVal,nrow(merged),length(peaks)+3)
   result[,1:3] <- as.matrix(merged[,1:3])
   def <- rep(defVal,nrow(merged))
   if(length(peaks)>0) {
      for(i in 1:length(peaks)) {
         peakset <- peaks[[i]][,1:4]
         peakset[,1] <- match(peakset[,1],chrmap)
         if(is.unsorted(unique(peakset[,1]))) {
            peakset <- pv.peaksort(peakset)
         }
         res <- mergeScores(merged,def,peakset,TRUE)
         result[,i+3] <- res$score
         included[,i] <- res$included
      }
      #included <- pv.called(merged, chrmap, peaks)
   }
   
   colnames(result) <- rep("",ncol(result))
   colnames(result)[1:3] <- c("CHR","START","END")
   if(ncol(result)>3) {
      colnames(result)[4:ncol(result)] <- colnames(classes)
   }
   
   pv.gc()
   return(list(merged=result,included=included, chrmap=chrmap))
}

pv.called <- function(merged, chrmap, peaks) {
   if(!is(merged,"GRanges")) {
      colnames(merged) <- c("chr","start","end")
      merged[,1] <- chrmap[merged[,1]]
      merged <- GRanges(merged)
   }
   called <- NULL
   for(i in 1:length(peaks)) {
      called <- cbind(called, merged %over% GRanges(peaks[[i]]))
   }
   called[called==TRUE] <- 1
   return(called)
}

pv.CalledMasks <- function(pv,newpv,master) {
   master <- cbind(master[,1:3],1)
   spare <- pv.peakset(pv,master,peak.caller='raw',scoreCol=4,bLowerScoreBetter=FALSE)
   spare <- pv.model(spare)
   res   <- pv.list(spare,spare$masks$counts)
   masternum <- length(spare$peaks)
   resl <- NULL
   for(i in 1:nrow(res)) {
      sampl <- pv.matchingSamples(res[i,1],spare$class)
      sampvec <- rep(FALSE,nrow(master))
      for(samp in sampl) {
         sampvec <- sampvec | pv.whichCalled(spare,samp,masternum)
      }
      resl <- pv.listadd(resl,sampvec)      
   }
   names(resl) <- res$ID
   return(resl)
}

pv.matchingSamples <- function(id,classes) {
   res <- which(!classes[PV_CALLER,] %in% "counts" & classes[PV_ID,] %in% id)
   return(res)  
}

pv.whichCalled <- function(pv,called,master,minVal=-1) {
   called <- pv$binding[,called+3] > minVal
   master <- pv$binding[,master+3] > minVal
   res <- called[master]     	
   return(res)
}


## pv.pairs -- compare all pairs of peaksets, rank by % overlap
pv.pairs <- function(pv,mask,bPlot=FALSE,attributes=pv$attributes,bAllVecs=TRUE,
                     CorMethod="pearson",bCorOnly=FALSE,bNonZeroCors=FALSE,minVal=0,bFixConstantVecs=TRUE) {
   
   if(missing(mask)) {
      mask=rep(TRUE,ncol(pv$class))
      peaks <- pv$peaks
   } else {
      peaks <- NULL
      for(i in 1:length(mask)){
         if(mask[i]) {
            peaks <- pv.listadd(peaks,pv$peaks[[i]])
         }
      }
   }
   
   tmp <- NULL
   if(bAllVecs==TRUE) {
      if(pv$totalMerged != nrow(pv$binding)) {
         pv <- pv.vectors(pv,minOverlap=1)  
      }
   }
   
   tmp$binding    <- pv$binding[,c(TRUE,TRUE,TRUE,mask)]
   tmp$class      <- pv$class[,mask]
   tmp$peaks      <- peaks
   tmp$chrmap     <- pv$chrmap
   tmp$totalMerged <- pv$totalMerged
   tmp$called      <- pv$called[,mask]
   if(!is.null(tmp$allcalled)) {
      tmp$allcalled      <- pv$allcalled[,mask]
   } else {
      tmp$allcalled <- NULL
   }
   
   numSets <- sum(mask)
   resm <- NULL
   resl <- NULL
   cvecs <- NULL
   for(first in 1:(numSets-1)) {
      if(bCorOnly==FALSE){
         #cat(".")
      }
      for(second in (first+1):numSets) {
         if(!bCorOnly) {
            res  <- pv.overlap(tmp,mask=c(first,second),minVal=minVal)
            resl <- pv.listadd(resl,res)
            inall <- nrow(res$inAll)
            onlya <- nrow(res$onlyA)
            onlyb <- nrow(res$onlyB)
            prop <- inall / (inall + onlya + onlyb)
         } else {
            inall <- 0
            onlya <- 0
            onlyb <- 0
            prop  <- 0
         }
         v1 <- tmp$binding[,first+3]
         v2 <- tmp$binding[,second+3]
         if(bNonZeroCors) {
            zeros <- v1==0 & v2==0
            v1 <- v1[!zeros]
            v2 <- v2[!zeros]
         }
         if(bFixConstantVecs) {
            if(sd(v1)==0) {
               #fval <- v1[1]
               #v1[1] <- fval+.000001
               #v1[2] <- fval-.000001
               cvecs <- c(cvecs,colnames(tmp$binding)[first+3])	
            }
            if(sd(v2)==0) {
               #fval <- v2[1]
               #v2[1] <- fval+.000001
               #v2[2] <- fval-.000001
               cvecs <- c(cvecs,colnames(tmp$binding)[second+3])	
            }		
         }
         if(sd(v1) && sd(v2)) {
            corval <- cor(v1,v2,method=CorMethod)
         } else corval <- 0
         resm <- rbind(resm,c(which(mask)[first],which(mask)[second],
                              onlya,onlyb,inall,corval,prop))
      }
   }
   if(bCorOnly==FALSE) {
      #cat("\n")
   }
   cvecs <- unique(cvecs)
   if(!is.null(cvecs)) {
      cvecs <- sort(cvecs)
      for(cv in cvecs) {
         #          warning(sprintf('Scores for peakset %s are all the same -- correlations set to zero.',cv),call.=FALSE)	
      }	
   }
   
   o <- order(resm[,7],decreasing=TRUE)
   colnames(resm) <- c("A","B","onlyA","onlyB","inAll","Cor","Overlap")
   
   if(bPlot & !bCorOnly){
      warning('Plotting in pv.occupancy unsupported',call.=FALSE)
      #for(i in 1:length(o)) {
      #	 recnum <- o[i]
      #   pv.PlotContrast(pv,resl[[recnum]],resm[recnum,1],resm[recnum,2],attributes=attributes)
      #}
   }
   
   return(resm[o,])	
}

pv.overlapToLabels <- function(pv,overlap,labelatts=PV_ID) {
   for(i in 1:nrow(overlap)) {
      overlap[i,1] <- pv.namestrings(pv$class[labelatts,as.numeric(overlap[i,1])])$tstring
      overlap[i,2] <- pv.namestrings(pv$class[labelatts,as.numeric(overlap[i,2])])$tstring
   }
   #overlap <- overlap[order(overlap[,1],overlap[,2]),]
   return(data.frame(overlap))
}

pv.orderfacs <- function(facvec,decreasing=FALSE) {
   res <- order(as.numeric(as.character(facvec)),decreasing=decreasing)
   return(res)	
}

pv.normalizeScores <- function(peaks,pCol,zeroVal=-1,bLog=FALSE,bDensity=FALSE){
   if(bDensity) {
      width   <- peaks[,3] - peaks[,2]
      width[width==0]=1
      density <- peaks[,pCol]/width 
      res <- density/max(density)
   } else {
      res <- peaks[,pCol]/max(peaks[,pCol])
   }
   if(bLog) {
      res <- log2(res)
      x <- res == -Inf
      res[x] <- 1
      res <- res - min(res)
      res[x] <- zeroVal
   } else {
      res[res == 0] <- zeroVal
   }
   return(res)
}

pv.activefun <- function(x){
   if(sum(x>0)>0){
      return(TRUE)
   } else {
      return(FALSE)
   }	
}

pv.addrow <- function(x,a,y){
   if(is.null(y)) return(x)
   nm <- nrow(y)
   if(is.null(nm)) return(rbind(x,a))
   if(nm == 0)     return(x)
   ncl <- length(a)
   return(rbind(x,matrix(a,nm,ncl,byrow=TRUE)))
}

pv.reorderM <- function(ocm,dgram) {
   lab <- rownames(ocm)	
   newlab <- rev(rapply(dgram,pv.dval))
   neword <- match(newlab,lab)
   newocm <- ocm[neword,]
   newocm <- newocm[,neword]
   return(newocm)
}
pv.dval <- function(dgram) {
   att <- attributes(dgram)
   if(!is.null(att$leaf)) {
      if(!is.null(att$label)) {
         if(att$leaf==TRUE) {
            return(att$label)
         }
      }
   }
}

pv.pcmask <- function(pv,numSites, mask, sites,removeComps,cor=FALSE,bLog=TRUE){
   
   if(missing(numSites)) numSites <- nrow(pv$binding)
   if(is.null(numSites)) numSites <- nrow(pv$binding)  
   numSites <- min(numSites,nrow(pv$binding))
   
   if(missing(sites)) sites <- 1:numSites
   if(is.null(sites)) sites <- 1:numSites
   
   if(missing(mask)) mask <- rep(TRUE,ncol(pv$class))
   for(i in which(mask)) {
      if(nrow(pv$peaks[[i]])==0) {
         mask[i] <- FALSE
      }
   }
   if(sum(mask)<2) {
      stop('Need at least two samples for PCA.',call.=FALSE)
   }
   
   res <- NULL   
   res$class <- pv$class
   pv$values <- pv$binding[sites,c(FALSE,FALSE,FALSE,mask)]
   active   <- apply(pv$values,1,pv.activefun)
   numSites <- min(numSites,sum(active))
   
   pv$values <- pv$values[active,][1:numSites,]
   
   if(!missing(removeComps)) {
      pv$values <- pv.removeComp(pv$values,numRemove=removeComps)
   }
   
   if(bLog) {
      if(max(pv$values)>1) {
         pv$values[pv$values<=1] <- 1
         pv$values <- log2(pv$values)
      }
   }
   
   if(nrow(pv$values) >= sum(mask)) {
      #res$pc <- prcomp(pv$values) #,cor=cor)
      res$pc <- prcomp(t(pv$values))
   }
   res$mask <- mask
   
   return(res)
}


pv.Signal2Noise <- function(pv) {
   sns <- rep("",length(pv$peaks))
   for(i in 1:length(sns)) {
      rip <- sum(pv$peaks[[i]]$Reads)
      treads <- as.integer(pv$class[PV_READS,i])
      if(treads) {
         sn  <- (rip-sum(pv$peaks[[i]]$Reads==1))/treads
         sns[i] <- sprintf("%1.2f",sn)
      }	
   }
   if(sum(as.numeric(sns) > 0,na.rm=TRUE)) {
      return(sns)  	
   } else {
      return(NULL)
   }
}


pv.peaksetCounts <- function(pv=NULL,peaks=NULL,counts,
                             sampID="",tissue="",factor="",condition="",treatment="",replicate) {
   
   
   if(!is.null(pv$peaks)) {
      if(sum(!pv$class[PV_CALLER,] %in% "counts")) {
         stop("DBA object can only have count peaksets",call.=FALSE)	
      }	
   }
   
   IDs <- "counts"
   froms <- NULL
   tos   <- NULL
   if(is.null(dim(counts)) && length(counts)==1) { # filename
      counts <- read.table(counts,as.is=TRUE,fill=TRUE)
      counts <- counts[!is.na(counts[,2]),]
      if(ncol(counts) == 3) {
         if(sum(is.na(counts[,3]))==nrow(counts)) {
            counts <- counts[,1:2]
         }
      }
      #annotation <- as.character(counts[,1])
      #counts[,1] <- 1:nrow(counts)
   }
   if(!is.vector(counts)) {   
      numcols <- ncol(counts)
      if(numcols>1) {
         if(numcols == 2) {
            numcounts <- nrow(counts)
            footer <- match("Assigned",counts[,1])
            if(!is.na(footer)) {
               numcounts <- footer-1
            }
            IDs <- counts[1:numcounts,1]
            counts <- counts[1:numcounts,2]
            #annotation <- annotation[1:numcounts]
         } else if(numcols == 4) {
            IDs    <- counts[,1]
            froms  <- counts[,2]
            tos    <- counts[,3]
            counts <- counts[,4]	
         } else {
            stop("Counts must have 1, 2, or 4 columns.",call.=FALSE)	
         }
      } else counts <- counts[,1]
   }
   
   if(!is.null(dim(counts))) {
      stop('counts must be vector of counts',call.=FALSE)	
   }
   numcounts <- length(counts)
   #if(length(peaks)>0) {
   #   warning('Specified peaks ignored',call.=FALSE)
   #}
   if(is.null(froms)) {
      froms <- tos <- 1:numcounts	
   }
   if(is.null(peaks)) {
      peaks <- data.frame(cbind(IDs,froms,tos))
   } else {
      peaks <- data.frame(peaks[,1:3])
   }
   
   peaks <- cbind(peaks,counts,counts,rep(0,numcounts),counts,rep(0,numcounts),rep(0,numcounts))
   colnames(peaks) <- c("Chr","Start","End", "Score", "Score","RPKM", "Reads","cRPKM","cReads")
   
   res <- dba.peakset(pv,
                      peaks       = peaks,
                      sampID      = sampID,
                      tissue      = tissue,
                      factor      = factor,
                      condition   = condition,
                      treatment   = treatment,
                      consensus   = TRUE,
                      peak.caller = "counts",
                      reads       = NA,
                      replicate   = replicate,
                      bMerge      = FALSE)
   
   res$class[PV_READS,length(res$peaks)] <- sum(counts)
   
   if(nrow(peaks) != nrow(res$peaks[[1]])) {
      stop('Mismatch in number of intervals',call.=FALSE)
   }
   
   if(!is.null(pv$rownames)) {
      if(sum(pv$rownames %in% rownames)!=length(rownames)) {
         stop('Mismatch in interval annotations',call.=FALSE)           
      }
   }
   
   #     if(sum(as.integer(res$peaks[[1]][,1]) != as.integer(res$peaks[[length(res$peaks)]][,1]))) {
   #         stop("Mismatch in ID",call.=FALSE)	
   #     }
   if(sum(res$peaks[[1]][,1] != res$peaks[[length(res$peaks)]][,1])) {
      stop("Mismatch in ID",call.=FALSE)    
   }
   
   if(sum(res$peaks[[1]][,2] != res$peaks[[length(res$peaks)]][,2])) {
      stop("Mismatch in interval start",call.=FALSE)	
   }
   if(sum(res$peaks[[1]][,3] != res$peaks[[length(res$peaks)]][,3])) {
      stop("Mismatch in interval end",call.=FALSE)	
   }
   
   #res$annotation <- annotation
   
   return(res)
}

pv.DBA2SummarizedExperiment <- function(DBA, bAssays=TRUE, report) {
   if (!missing(report)) {
      peaks <- report[,1:9]
   } else {
      peaks <- pv.writePeakset(DBA, peaks=DBA$binding, numCols=3)
   }
   rnames <- rownames(peaks)
   peaks <- pv.peaks2DataType(peaks,DBA_DATA_GRANGES)
   meta  <- t(DBA$class)
   meta[meta==""]=NA
   meta <- meta[,apply(meta,2,function(x) {sum(is.na(x))!=length(x)})]
   meta <- DataFrame(meta[,-1])
   counts <- as.matrix(DBA$binding[,4:ncol(DBA$binding)])
   assays <- SimpleList(scores=counts)
   if(bAssays) {
      extra <- pv.assaysFromPeaks(DBA) 
      if(!is.null(extra)){
         for(newassay in extra) {
            rownames(newassay) <- rnames
            assays <- pv.listadd(assays,newassay)
         }   
         names(assays) <- c("scores",names(extra))
      }  
   }
   res <- SummarizedExperiment(assays    = assays,
                               rowRanges = peaks,
                               colData   = meta)
   colnames(res) <- colnames(DBA$class)                           
   return(res)                           
}

pv.assaysFromPeaks <- function(DBA) {
   peaks <- DBA$peaks
   numpeaks  <- nrow(peaks[[1]])
   numfields <- ncol(peaks[[1]])
   if(numfields <= 4) {
      return(NULL)
   }
   fields  <- colnames(peaks[[1]])
   samples <- length(peaks)
   assays <- NULL
   for(i in 5:numfields) {
      toadd <- matrix(0,numpeaks,samples)
      rownames(toadd) <- 1:numpeaks
      colnames(toadd) <- colnames(DBA$class)
      assays <- pv.listadd(assays,toadd)
   }   
   names(assays) <- fields[5:numfields]
   for(sampnum in 1:samples) {
      peak <- peaks[[sampnum]]
      if(numpeaks!=nrow(peak) || numfields!=ncol(peak)) {
         return(NULL)
      }
      for(i in 5:numfields) {
         assays[[i-4]][,sampnum] <- peak[,i]
      }
   }
   
   if(!is.null(DBA$called)) {
      if(nrow(DBA$called)==nrow(assays[[1]])){
         assays <- c(assays,Called=list(DBA$called))
      }
   }
   
   return(assays)
}

stripSpaces <- function(data) {
   data[is.na(data)]=""
   cnames <- names(data)
   for (cname in cnames) {
      stripped <- trimws(data[,cname])
      changed <- !(stripped == data[,cname])
      if (sum(changed) > 0) {
         rows <- seq(1,length(stripped))
         rownums <- rows[changed]
         orig <- paste(stripped[changed],sep=",",collapse=",")
         warning(sprintf("Removed white space from %s in column %s (%s %s)",
                         orig,cname,
                         if (length(rownums)==1) "row" else "rows",
                         paste(rownums,sep=",",collapse=",")),call.=FALSE)
         data[,cname] <- stripped
      }
   }
   return(data)
}

pv.gc <- function(force=FALSE){
   if(force) {
      gc(verbose=FALSE)
   }
}

pv.overlaps <- function(pv,minOverlap) {
   if(is.null(pv$called)) {
      return(NULL)
   }
   if(!is.null(pv$allcalled)) {
      pv$called <- pv$allcalled
   }
   overlaps <- apply(pv$called,1,sum)>=minOverlap
   return(overlaps)
}

pv.allSame <- function(pv) {
   callers <- unique(pv$class[DBA_CALLER,])
   if((length(callers)==1) & (callers[1]=='counts')) {
      if(length(unique(sapply(pv$peaks,nrow)))==1) {
         return(TRUE)
      }
   }
   return(FALSE)
}

pv.makeGRanges <- function(data, chrmap) {
   data <- data.frame(data)
   colnames(data)[1:3] <- c("chr","start","end")
   data[,1] <- chrmap[data[,1]]
   data <- GRanges(data)
   return(data)
}

pv.checkCalled <- function(pv){
   
   if(!is.null(pv$called)) {
      if(nrow(pv$called) != nrow(pv$binding)) {
         if(is.null(pv$allcalled)) {
            if(nrow(pv$called) != nrow(pv$merged)) {
               pv$called <- NULL
               return(pv)
            }
         }
      }
   }
   
   if(nrow(pv$merged) != nrow(pv$binding)) {
      if (!is.null(pv$called) && is.null(pv$allcalled)) {
         if(nrow(pv$called) == nrow(pv$merged)) {
            pv$allcalled <- pv$called
            pv$called <- pv$allcalled[pv.makeGRanges(pv$merged,pv$chrmap) %over% 
                                         pv.makeGRanges(pv$binding, pv$chrmap),]
         }
      }
   }
   
   if(!is.null(pv$allcalled)) {
      if(nrow(pv$allcalled) != nrow(pv$merged)) {
         pv$allcalled <- NULL
      }
   }
   
   if(!is.null(pv$called)) {
      if(nrow(pv$called) != nrow(pv$binding)) {
         pv$called <- NULL
      }
   }
   
   return(pv)
}

