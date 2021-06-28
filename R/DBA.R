#############################################
## DBA.R -- Differential Binding Analysis  ##
## 20 April 2011                           ##
## Rory Stark                              ##
## Cancer Research UK                      ##
#############################################

## dba	            Construct a dba object

## dba.peakset	    Add a peakset to a dba object
## dba.overlap	    Compute binding site overlaps
## dba.blacklist    apply/generate blacklists and greylists

## dba.count  	    Count reads in binding sites
## dba.normalize    Normalize dataset

## dba.contrast	    Establish contrast(s) for analysis
## dna.analyze  	  Execute affinity analysis
## dba.report	      Generate report for a contrast analysis

## dba.plotClust	  Cluster dendrogram plo
## dba.plotHeatmap	Heatmap plot
## dba.plotPCA	    Principal Components plot
## dba.plotBox	    Boxplots
## dba.plotMA	      MA/scatter plot
## dba.plotVenn	    2, 3, or 4-way Venn diagram plot
## dba.plotVolcano	Volcano plots
## dba.plotProfile  Profile heatmaps

## dba.show	        List dba metadata
## dba.mask	        Mask samples or sites 

## dba.save	        Save dba object
## dba.load	        Load dba object

### NOTE: DBA is a front-end to a package formerly known as pv, with most DBA
### functions being simple pass-throughs to pv functions

##########################
### CONSTANT VARIABLES ###
##########################
DBA_GROUP     <- PV_GROUP
DBA_ID        <- PV_ID 
DBA_TISSUE    <- PV_TISSUE 
DBA_FACTOR    <- PV_FACTOR
DBA_CONDITION <- PV_CONDITION
DBA_TREATMENT <- PV_TREATMENT
DBA_CONSENSUS <- PV_CONSENSUS
DBA_CALLER    <- PV_CALLER
DBA_CONTROL   <- PV_CONTROL
DBA_READS     <- PV_READS
DBA_REPLICATE <- PV_REPLICATE
DBA_INTERVALS <- PV_INTERVALS
DBA_FRIP      <- PV_SN_RATIO
DBA_ALL_ATTRIBUTES <- c(DBA_ID,DBA_TISSUE,DBA_FACTOR,
                        DBA_CONDITION,DBA_TREATMENT,
                        DBA_REPLICATE,DBA_CALLER)

DBA_EDGER_BLOCK   <- 'edgeRlm'
DBA_EDGER_GLM     <- 'edgeRGLM'
DBA_EDGER         <- DBA_EDGER_GLM

DBA_DESEQ2        <- 'DESeq2'
DBA_DESEQ2_BLOCK  <- 'DESeq2Block'

DBA_ALL_METHODS <- c(DBA_EDGER,DBA_DESEQ2)
DBA_ALL_BLOCK   <- c(DBA_EDGER_BLOCK,DBA_DESEQ2_BLOCK)
DBA_ALL_METHODS_BLOCK <- c(DBA_ALL_METHODS, DBA_ALL_BLOCK)

DBA_DATA_FRAME                    <- 0
DBA_DATA_RANGEDDATA               <- 1
DBA_DATA_GRANGES                  <- 2
DBA_DATA_SUMMARIZED_EXPERIMENT    <- 3
DBA_DATA_DBAOBJECT                <- 4
DBA_DATA_DEFAULT                  <- DBA_DATA_GRANGES

#########################################################
## dba -- construct DBA object, e.g. from sample sheet ##
#########################################################
dba <- function(DBA,mask, minOverlap=2,
                sampleSheet="dba_samples.csv", 
                config=data.frame(AnalysisMethod=DBA_DESEQ2,th=0.05,
                                  DataType=DBA_DATA_GRANGES, RunParallel=TRUE, 
                                  minQCth=15, fragmentSize=125, 
                                  bCorPlot=FALSE, reportInit="DBA", 
                                  bUsePval=FALSE, design=TRUE,
                                  doBlacklist=TRUE, doGreylist=TRUE),
                peakCaller="raw", peakFormat, scoreCol, bLowerScoreBetter, 
                filter, skipLines=0, bAddCallerConsensus=FALSE, 
                bRemoveM=TRUE, bRemoveRandom=TRUE, bSummarizedExperiment=FALSE,
                attributes, dir) 
{
  if(!missing(DBA)){
    if(is(DBA,"character")) {
      stop("DBA object is a character string; perhaps meant to be argument \'sampleSheet\'?")	
    }
    DBA <- pv.check(DBA,bCheckSort=FALSE)	
  } 
  
  res <- pv.model(DBA, mask=mask, minOverlap=minOverlap, samplesheet=sampleSheet, 
                  config=config, caller=peakCaller, format=peakFormat, 
                  scorecol=scoreCol,bLowerBetter=bLowerScoreBetter,
                  skipLines=skipLines,bAddCallerConsensus=bAddCallerConsensus, 
                  bRemoveM=bRemoveM, bRemoveRandom=bRemoveRandom,
                  filter=filter, attributes=attributes, dir)
  
  res$contrasts=NULL
  
  if(nrow(res$peaks[[1]])>0) {
    if(sum(res$peaks[[1]]$Score<0)>0) {
      res <- pv.ResetScores(res,ones=FALSE)
    }
  }
  
  if(is.null(res$config$DataType)) {
    res$config$DataType <- DBA_DATA_DEFAULT
  }
  
  if(is.null(res$config$parallelPackage)){
    res$config$parallelPackage=DBA_PARALLEL_MULTICORE
  }
  if(is.null(res$config$RunParallel)){
    res$config$RunParallel <- TRUE
  }
  if(is.null(res$config$AnalysisMethod)){
    res$config$AnalysisMethod=DBA_DESEQ2
  }
  if(is.null(res$config$bCorPlot)){
    res$config$bCorPlot <- FALSE
  } 
  if(is.null(res$config$th)){
    res$config$th=0.05
  }
  if(is.null(res$config$bUsePval)){
    res$config$bUsePval <- FALSE
  }
  
  if(is.null(res$config$cores)){
    res$config$cores <- parallel::detectCores(logical=FALSE)
  }
  
  if(is.null(res$config$doBlacklist)){
    res$config$doBlacklist <- TRUE
  }
  if(is.null(res$config$doGreylist)){
    res$config$doGreylist <- TRUE
  }
  
  if(missing(DBA)){
    DBA <- NULL
  } 
  if(is.null(DBA$config$reportInit)){
    res$config$reportInit <- "DBA"
  } else {
    res$config$reportInit <- DBA$config$reportInit
  }
  if(!is(res,"DBA")){
    class(res) <- "DBA"
  }
  
  if(res$config$bCorPlot) {
    try(dba.plotHeatmap(res),silent=TRUE)      
  }
  
  if(bSummarizedExperiment) {
    res <- pv.DBA2SummarizedExperiment(res)
  } else {
    if(!is.null(DBA$ChIPQCobj)) {
      #          resQC <- DBA$ChIPQCobj
      #          resQC@DBA <- res
      #          res <- resQC
      warning('Returning new DBA object (not ChIPQCexperiment object)')
    }
    res$ChIPQCobj <- NULL   
  }
  
  gc(verbose=FALSE)
  
  return(res)                 
}

###############################################
## dba.peakset -- add a peakset to the model ##
###############################################

dba.peakset <- function(DBA=NULL, peaks, sampID, tissue, factor, 
                        condition, treatment, replicate,control, peak.caller, 
                        peak.format, reads=0, consensus=FALSE, 
                        bamReads, bamControl, spikein,
                        scoreCol, bLowerScoreBetter, filter, counts, 
                        bRemoveM=TRUE, bRemoveRandom=TRUE,
                        minOverlap=2, bMerge=TRUE,
                        bRetrieve=FALSE, writeFile, numCols=4,
                        DataType=DBA$config$DataType)
{
  res <- NULL
  
  if(!missing(peaks)){
    if(!is(peaks,"DBA")) {
      peaks <- pv.DataType2Peaks(peaks)
    } else {
      peaks <- pv.check(peaks,bDoVectors=FALSE)	
    }
  }
  
  if(missing(DataType)) {
    DataType <- DBA_DATA_DEFAULT
  }
  
  if(bRetrieve==TRUE || !missing(writeFile)) { ## RETRIEVE/WRITE PEAKSETS
    
    if(missing(writeFile)) {
      writeFile <- NULL
    }
    
    if(missing(peaks)) {
      if(!is.null(DBA)) {
        DBA <- pv.check(DBA,bCheckEmpty=TRUE,bDoVectors=FALSE)
      } else {
        stop("DBA object is NULL; can't retrieve peakset.")	
      }	
    } else {
      if(is.vector(peaks)) {
        if(is.null(DBA)) {
          stop("DBA object is NULL; can't retrieve peakset.")	
        }	
        if(length(peaks) > 1) {
          if(minOverlap > length(peaks)) {
            stop('minOverlap is greater than number of specified peaks.')	
          } else {
            saveCorPlot <- DBA$config$bCorPlot
            DBA$config$bCorPlot <- FALSE
            DBA <- dba(DBA,mask=peaks,minOverlap=minOverlap)
            peaks <- NULL
            DBA$config$bCorPlot<- DBA$config$bCorPlot
          }
        }      
      }	
    }
    
    res <- pv.writePeakset(DBA, fname=writeFile, peaks=peaks, numCols=numCols)     
    
    if(DataType!=DBA_DATA_FRAME) {
      res <- pv.peaks2DataType(res,DataType)
    }    
    
  } else { ## ADD PEAKSET(S)
    
    bCheckS <- TRUE
    bDoV <- TRUE
    
    if(!is.null(DBA)) {
      DBA <- pv.check(DBA,bDoVectors=FALSE)
    }
    if(!missing(peaks)) {
      if(is(peaks,"DBA")) {
        res <- pv.peakset_all(DBA, addpv=peaks, minOverlap=minOverlap)
        bCheckS <- FALSE
        bDoV <- FALSE
      }
    }
    if(is.null(res)) {
      if(!missing(consensus) && pv.isConsensus(DBA)) {
        if(consensus != TRUE) {
          stop("DBA object is already formed from a consensus peakset!")
        }
      }
      res <- pv.peakset(DBA, peaks=peaks, 
                        sampID=sampID, tissue=tissue, factor=factor,
                        condition=condition,treatment=treatment,
                        replicate=replicate,control=control,
                        peak.caller=peak.caller, peak.format=peak.format, 
                        reads=reads, consensus=consensus, 
                        readBam=bamReads, controlBam=bamControl,
                        scoreCol=scoreCol, bLowerScoreBetter=bLowerScoreBetter, 
                        bRemoveM=bRemoveM, bRemoveRandom=bRemoveRandom,
                        minOverlap=minOverlap, filter=filter, counts=counts, 
                        spikein=spikein)
    }
    
    if(!is(res,"DBA")) {
      class(res) <- "DBA"
    }
    
    if(is.null(res$config$DataType)) {
      res$config$DataType <- DBA_DATA_DEFAULT
    }
    if(is.null(res$config$parallelPackage)){
      res$config$parallelPackage <- DBA_PARALLEL_MULTICORE
    }
    if(is.null(res$config$RunParallel)){
      res$config$RunParallel <- TRUE
    }
    if(is.null(res$config$th)){
      res$config$th <- 0.05
    }
    if(is.null(res$config$bUsePval)){
      res$config$bUsePval <- FALSE
    }
    
    if(is.null(DBA$config$reportInit)){
      res$config$reportInit <- "DBA"
    } else {
      res$config$reportInit <- DBA$config$reportInit
    }
    if(is.null(res$config$AnalysisMethod)){
      res$config$AnalysisMethod <- DBA_DESEQ2
    }
    if(is.null(res$config$bCorPlot)){
      res$config$bCorPlot <- FALSE
    } 
    
    if(is.null(res$config$cores)){
      res$config$cores <- parallel::detectCores(logical=FALSE)
    }
    
    if(is.null(res$config$doBlacklist)){
      res$config$doBlacklist <- TRUE
    }
    if(is.null(res$config$doGreylist)){
      res$config$doGreylist <- TRUE
    }
    
    if(bMerge) {
      res <- pv.check(res,bCheckSort=bCheckS,bDoVectors=bDoV)
      gc(verbose=FALSE)
    }
    
    if(!is.null(DBA$ChIPQCobj)) {
      #          resQC <- DBA$ChIPQCobj
      #          resQC@DBA <- res
      #          res <- resQC
      warning('Returning new DBA object (not ChIPQCexperiment object)')
    }   
  }   
  
  numpeaks <- length(res$peaks)
  if(numpeaks > 1) {
    if(is.null(res$called)) {
      res <- dba(res)
    } else 
      if(numpeaks != ncol(res$called)) {
        res <- dba(res)
      }
  }
  
  return(res)                       
  
}                      

##################################################
## dba.overlap -- compute binding site overlaps ##
##################################################

DBA_OLAP_PEAKS <- 1 # Return list of peaksets (common/unique peaks) 
DBA_OLAP_ALL   <- 2 # Return overlap report with statstics for peakset pairs
DBA_OLAP_RATE  <- 3 # Return vector of number of peaks in overlap for all values of minOverlap

DBA_OLAP  <- PV_OLAP
DBA_COR   <- PV_COR
DBA_INALL <- PV_INALL

dba.overlap <- function(DBA, mask, mode=DBA_OLAP_PEAKS, 
                        contrast, method=DBA$config$AnalysisMethod, 
                        th=DBA$config$th, bUsePval=DBA$config$bUsePval, report,
                        byAttribute, bCorOnly=TRUE, CorMethod="pearson", 
                        DataType=DBA$config$DataType)
{                      
  DBA <- pv.check(DBA,bCheckEmpty=TRUE)
  
  if( (mode == DBA_OLAP_ALL) | (!missing(contrast)) | (!missing(report)) ) {
    
    if( (!missing(contrast)) | (!missing(report)) ) {
      
      if(missing(report)) {
        report   <- dba.report(DBA,method=method, contrast=contrast,th=th,
                               bUsePval=bUsePval,DataType=DBA_DATA_FRAME)
      } else {
        if(!is(report,"data.frame")) {
          stop('Class not supported for report parameter. Call dba.report with DataType=DBA_DATA_FRAME.')	
        }
        if(!missing(contrast)) {
          DBA <- pv.getOverlapData(DBA,contrast,report)
        }
      }
      
      sites <- as.numeric(rownames(report))
      
      if(missing(mask)) {
        if(missing(contrast)) {
          mask <- 1:length(DBA$peaks)
        } else {
          mask  <- DBA$contrasts[[contrast]]$group1 | DBA$contrasts[[contrast]]$group2
        }   
      }  else if (length(mask)==1) {
        mask <- 1:length(DBA$peaks)         
      }
      
      res <- pv.occupancy(DBA, mask=mask, sites=sites, byAttribute=byAttribute, 
                          Sort='cor', CorMethod=CorMethod,
                          bCorOnly=bCorOnly)         
      
    } else {
      
      res <- pv.occupancy(DBA, mask=mask, byAttribute=byAttribute, 
                          Sort='cor', CorMethod=CorMethod,
                          bCorOnly=bCorOnly)
    } 
    
  }  else if(mode == DBA_OLAP_RATE) {
    
    res <- pv.overlapRate(DBA,mask=mask)
    
  }  else {
    
    res <- pv.overlap(DBA,mask=mask)
    
    if(DataType!=DBA_DATA_FRAME) {
      res <- lapply(res,pv.peaks2DataType,DataType)
    }   
    if(DataType == DBA_DATA_GRANGES) {
      for(i in 1:length(res)) {
        if(!is.null(res[[i]])) {
          if(is.null(res[[i]]$score)) {
            res[[i]]$score <- rowMeans(data.frame(mcols(res[[i]])), 
                                            na.rm=TRUE)
          }
        } else {
          res[[i]] <- GRanges(NULL)
        }
      }
      res <- GRangesList(res)
    }
  }                       
  
  return(res)
}
   
##############################################################
## dba.blacklist -- apply/generate blacklists and greylists ##
##############################################################

DBA_BLACKLIST_HG19    <- PV_BLACKLIST_HG19
DBA_BLACKLIST_HG38    <- PV_BLACKLIST_HG38 
DBA_BLACKLIST_GRCH37  <- PV_BLACKLIST_GRCH37 
DBA_BLACKLIST_GRCH38  <- PV_BLACKLIST_GRCH38 
DBA_BLACKLIST_MM9     <- PV_BLACKLIST_MM9
DBA_BLACKLIST_MM10    <- PV_BLACKLIST_MM10   
DBA_BLACKLIST_CE10    <- PV_BLACKLIST_CE10  
DBA_BLACKLIST_CE11    <- PV_BLACKLIST_CE11  
DBA_BLACKLIST_DM3     <- PV_BLACKLIST_DM3   
DBA_BLACKLIST_DM6     <- PV_BLACKLIST_DM6   
DBA_BLACKLIST         <- "blacklist"
DBA_GREYLIST          <- "greylist"
DBA_BLACKLISTED_PEAKS <- "blacklisted peaks"

dba.blacklist <- function(DBA, blacklist=DBA$config$doBlacklist, 
                          greylist=DBA$config$doGreylist, 
                          Retrieve, cores=DBA$config$cores) {
  
  DBA <- pv.check(DBA)
  
  if(is.null(blacklist)) {
    blacklist <- TRUE
  }
  
  if(is.null(greylist)) {
    if(pv.noControls(DBA$class["bamControl",])) {
      greylist <- FALSE
    } else {
      greylist <- TRUE
    }
  }
  
  res <- pv.BlackGreyList(DBA=DBA, blacklist=blacklist, greylist=greylist,
                          Retrieve=Retrieve, cores=cores)
  
  return(res)
  
}

###############################################  
## dba.count -- count reads in binding sites ##
###############################################   

DBA_SCORE_NORMALIZED          <- PV_SCORE_NORMALIZED
DBA_SCORE_RPKM                <- PV_SCORE_RPKM
DBA_SCORE_RPKM_FOLD           <- PV_SCORE_RPKM_FOLD
DBA_SCORE_RPKM_MINUS          <- PV_SCORE_RPKM_MINUS
DBA_SCORE_READS               <- PV_SCORE_READS
DBA_SCORE_CONTROL_READS       <- PV_SCORE_CONTROL_READS
DBA_SCORE_READS_FOLD          <- PV_SCORE_READS_FOLD
DBA_SCORE_READS_MINUS         <- PV_SCORE_READS_MINUS
DBA_SCORE_READS_FULL          <- PV_SCORE_READS_FULL
DBA_SCORE_READS_EFFECTIVE     <- PV_SCORE_READS_EFFECTIVE
DBA_SCORE_READS_MINUS_FULL    <- PV_SCORE_READS_MINUS_FULL
DBA_SCORE_READS_MINUS_EFFECTIVE <- PV_SCORE_READS_MINUS_EFFECTIVE
DBA_SCORE_TMM_MINUS_FULL      <- PV_SCORE_TMM_MINUS_FULL
DBA_SCORE_TMM_MINUS_EFFECTIVE <- PV_SCORE_TMM_MINUS_EFFECTIVE
DBA_SCORE_TMM_READS_FULL      <- PV_SCORE_TMM_READS_FULL
DBA_SCORE_TMM_READS_EFFECTIVE <- PV_SCORE_TMM_READS_EFFECTIVE
DBA_SCORE_TMM_MINUS_FULL_CPM      <- PV_SCORE_TMM_MINUS_FULL_CPM
DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM <- PV_SCORE_TMM_MINUS_EFFECTIVE_CPM
DBA_SCORE_TMM_READS_FULL_CPM      <- PV_SCORE_TMM_READS_FULL_CPM
DBA_SCORE_TMM_READS_EFFECTIVE_CPM <- PV_SCORE_TMM_READS_EFFECTIVE_CPM
DBA_SCORE_SUMMIT              <- PV_SCORE_SUMMIT
DBA_SCORE_SUMMIT_ADJ          <- PV_SCORE_SUMMIT_ADJ
DBA_SCORE_SUMMIT_POS          <- PV_SCORE_SUMMIT_POS

DBA_READS_DEFAULT <- PV_READS_DEFAULT
DBA_READS_BAM     <- PV_READS_BAM
DBA_READS_BED     <- PV_READS_BED

DBA_SCORE_FOLD             <- PV_SCORE_FOLD 
DBA_SCORE_CONCENTRATION    <- PV_SCORE_CONCENTRATION  
DBA_SCORE_CONC_NUMERATOR   <- PV_SCORE_CONC_NUMERATOR 
DBA_SCORE_CONC_DENOMINATOR <- PV_SCORE_CONC_DENOMINATOR
DBA_SCORE_PVAL             <- PV_SCORE_PVAL    
DBA_SCORE_FDR              <- PV_SCORE_FDR  

dba.count <- function(DBA, peaks, minOverlap=2, score=DBA_SCORE_NORMALIZED,
                      fragmentSize=DBA$config$fragmentSize, summits=200, 
                      filter=1, bRemoveDuplicates=FALSE,
                      bScaleControl=TRUE, bSubControl = is.null(DBA$greylist), 
                      mapQCth=DBA$config$mapQCth, filterFun=max, minCount=0, 
                      bLog=FALSE, bUseSummarizeOverlaps=TRUE,  
                      readFormat=DBA_READS_DEFAULT,bParallel=DBA$config$RunParallel) 
{
  DBA <- pv.check(DBA, missing(peaks))
  
  if(!is.null(DBA$resultObject)) {
    if(DBA$resultObject == TRUE) {
      return(pv.ResetResultScores(DBA,score))
    }
  }
  
  if(minOverlap > length(DBA$peaks) && missing(peaks)) {
    minOverlap <- length(DBA$peaks)
  }
  
  if(is.null(DBA$config$mergeOverlap)) {
    maxGap <- as.integer(-1)
  } else {
    maxGap <- as.integer(-DBA$config$mergeOverlap)
  }
  
  bUseLast <- FALSE
  
  res <- NULL
  resetContrasts <- TRUE
  
  if(is.logical(summits)) {
    if(summits==FALSE) {
      DBA$summits=NULL
      summits <- -1
    } else {
      summits <- 0
    }
  }
  
  if(!missing(peaks)) {
    if(is.null(peaks) && missing(summits)) { # don't recenter if summits not specified
      summits <- -1
    } else if(is.null(peaks) && summits == 0) { # don't recenter if summit amount not specified
      summits <- -1
    }
  }
  
  if( (summits > -1) && missing(peaks) && !is.null(DBA$summits)) {
    if(DBA$summits == 0 ) {
      peaks=NULL
    } else {
      if(summits != DBA$summits) {
        stop('Can not change value of summits. Re-run from peaks.')
      } else {
        warning('No action taken, returning passed object...')
        res <- DBA
      }
    }
  }
  
  if(!missing(peaks)) {
    if(is.null(peaks)) {
      if(!is.null(DBA$minCount)) {
        if(DBA$minCount != minCount) {
          warning("Unable to reset minCount without re-counting (peaks can not be NULL).",
                  call.=FALSE)
        } 
      } else {
        DBA$minCount <- minCount 
      }
    }
  }
  
  if(!missing(peaks) || length(filter)>1) {
    if(is.null(peaks) || length(filter)>1) {
      callers <- unique(DBA$class[DBA_CALLER,])
      if((length(callers)==1) & (callers[1]=='counts')) {
        DBA <- pv.check(DBA)
        if(summits != -1) {
          if(summits>0) {
            res <- pv.Recenter(DBA,summits,1:length(DBA$peaks),DBA$called)
            res <- pv.counts(res,peaks=res$merged,
                             defaultScore=score, bLog=bLog, 
                             insertLength=fragmentSize, bOnlyCounts=TRUE,
                             bCalledMasks=TRUE, minMaxval=filter,
                             bParallel=bParallel, bUseLast=bUseLast,
                             bWithoutDupes=bRemoveDuplicates,
                             bScaleControl=bScaleControl,filterFun=filterFun,
                             bLowMem=bUseSummarizeOverlaps,
                             readFormat=readFormat,summits=0,
                             minMappingQuality=mapQCth,bRecentered=TRUE,
                             minCount=minCount, bSubControl=bSubControl, 
                             maxGap=maxGap)
            res$summits <- summits
            res$norm <- DBA$norm
            res <- pv.reNormalize(res)
          } 
        } else {
          if(length(filter)>1) {
            DBA <- pv.setScore(DBA,score=score,bLog=bLog,minMaxval=0,
                               filterFun=filterFun)
            res <- pv.filterRate(DBA,filter,filterFun=filterFun)	
          } else {
            res <- pv.setScore(DBA,score=score,bLog=bLog,minMaxval=filter,
                               filterFun=filterFun)
          }
          resetContrasts=FALSE	
        }
      } else {
        stop('DBA object must contain only counts')	
      }	
    } else {
      peaks <- pv.DataType2Peaks(peaks)
    }
  }
  
  if(is.null(res)) {
    res <- pv.counts(DBA, peaks=peaks, minOverlap=minOverlap, 
                     defaultScore=score, bLog=bLog, insertLength=fragmentSize, 
                     bOnlyCounts=TRUE, bCalledMasks=TRUE, minMaxval=filter,
                     bParallel=bParallel, bUseLast=bUseLast,
                     bWithoutDupes=bRemoveDuplicates,bScaleControl=bScaleControl,
                     filterFun=filterFun,
                     bLowMem=bUseSummarizeOverlaps,readFormat=readFormat,
                     summits=summits,
                     minMappingQuality=mapQCth,
                     minCount=minCount, bSubControl=bSubControl,
                     maxGap=maxGap)
    if(summits != -1) {
      res$summits <- summits
    }
    res$norm <- DBA$norm
    res <- pv.reNormalize(res)
  }
  if(resetContrasts && length(res$contrasts)>0) {
    for(i in 1:length(res$contrasts)) {
      res$contrasts[[i]]$edgeR <- NULL
      res$contrasts[[i]]$DESeq <- NULL         	
    }
  }
  
  if(DBA$config$bCorPlot){
    try(dba.plotHeatmap(res,correlations=TRUE),silent=TRUE)
  }
  
  if(!is(res,"DBA")) {
    class(res) <- "DBA"
  }
  
  if(!is.null(DBA$ChIPQCobj)) {
    res <- checkQCobj(DBA$ChIPQCobj,res)
  }
  
  gc(verbose=FALSE)
  
  return(res)
}

########################################  
## dba.normalize -- normalize dataset ##
######################################## 

DBA_LIBSIZE_DEFAULT    <- PV_LIBSIZE_DEFAULT
DBA_LIBSIZE_FULL       <- PV_LIBSIZE_FULL
DBA_LIBSIZE_PEAKREADS  <- PV_LIBSIZE_PEAKREADS
DBA_LIBSIZE_BACKGROUND <- PV_LIBSIZE_CHRREADS
DBA_LIBSIZE_USER       <- PV_LIBSIZE_USER

DBA_NORM_DEFAULT        <- PV_NORM_DEFAULT
DBA_NORM_LIB            <- PV_NORM_LIB
DBA_NORM_TMM            <- PV_NORM_TMM 
DBA_NORM_RLE            <- PV_NORM_RLE
DBA_NORM_NATIVE         <- PV_NORM_NATIVE
DBA_NORM_SPIKEIN        <- PV_NORM_SPIKEIN
DBA_NORM_USER           <- PV_NORM_USER
DBA_NORM_OFFSETS        <- PV_NORM_OFFSETS
DBA_NORM_OFFSETS_ADJUST <- PV_NORM_OFFSETS_ADJUST

DBA_OFFSETS_LOESS     <- PV_OFFSETS_LOESS 
DBA_OFFSETS_USER      <- PV_OFFSETS_USER      

dba.normalize <- function(DBA, method = DBA$config$AnalysisMethod,
                          normalize   = DBA_NORM_DEFAULT,
                          library     = DBA_LIBSIZE_DEFAULT, 
                          background=FALSE, spikein=FALSE, offsets=FALSE, 
                          libFun=mean, bRetrieve=FALSE, ...) {
  
  DBA <- pv.check(DBA,TRUE)
  
  ### SET DEFAULTS ###
  pre3 <- FALSE
  if(is.null(DBA$design) && !is.null(DBA$contrasts)) {
    pre3 <- TRUE
  }
  if(!is.null(DBA$config$design)) {
    if(DBA$config$design == FALSE) {
      pre3 <- TRUE
    }
  }
  
  if(pre3 == TRUE) {
    # defaults are pre-3.0
    deflib  <- DBA_LIBSIZE_FULL
    defnorm <- DBA_NORM_DEFAULT
    defback <- FALSE   
  } else {
    # defaults are 3.0+    
    deflib  <-  DBA_LIBSIZE_BACKGROUND
    if(is(background,"logical")) {
      if(background == FALSE) {
        deflib  <- DBA_LIBSIZE_FULL
      }
    } 
    defnorm <- DBA_NORM_LIB
    defback <- FALSE    
  }
  
  if(length(library)==1 && library[1]==DBA_LIBSIZE_DEFAULT) {
    library <- deflib
  }
  
  if(length(normalize)==1 && normalize[1]==DBA_NORM_DEFAULT) {
    normalize <- defnorm
  }
  
  if(missing(background)) {
    background <- defback
  }
  
  res <- pv.normalize(DBA, method=method, 
                      libSizes=library, normalize=normalize,
                      libFun=libFun,
                      background=background, spikein=spikein, offsets=offsets, 
                      bRetrieve=bRetrieve, filter=0, ...)
  
  if(bRetrieve == FALSE) {
    if(res$score == DBA_SCORE_NORMALIZED) {
      res <- pv.doResetScore(res)
    }
  }
  
  return(res)
}

########################################################
## dba.contrast -- establish contrast(s) for analysis ##
########################################################

dba.contrast <- function(DBA, design=missing(group1), contrast,
                         group1, group2=!group1, name1, name2,
                         minMembers=3, block, bNot = FALSE, bComplex = FALSE,
                         categories=c(DBA_TISSUE,DBA_FACTOR,DBA_CONDITION,DBA_TREATMENT), 
                         bGetCoefficients=FALSE, reorderMeta)
{
  if(minMembers < 2) {
    stop('minMembers must be at least 2. Use of replicates strongly advised.')	
  }
  
  DBA <- pv.check(DBA,TRUE)
  
  if(missing(design) && design) {
    if(!is.null(DBA$config$design)) {
      design <- DBA$config$design 
    }
  }
  
  res <- pv.contrast(DBA, group1=group1, group2=group2, name1=name1, name2=name2,
                     minMembers=minMembers, categories=categories,block=block, 
                     bMulti = bComplex, bNot=bNot,
                     design=design, contrast=contrast, 
                     bGetNames=bGetCoefficients, reorderMeta=reorderMeta)
  
  if(!bGetCoefficients) {
    if(!is(res,"DBA")) {
      class(res) <- "DBA"
    }
    
    if(!is.null(DBA$ChIPQCobj)) {
      res <- checkQCobj(DBA$ChIPQCobj,res)
    }
  }
  
  return(res)                       	
}

###################################################################
##  -- perform differential binding affinity analysis ##
###################################################################

dba.analyze <- function(DBA, method=DBA$config$AnalysisMethod, design,
                        bBlacklist=DBA$config$doBlacklist,
                        bGreylist=DBA$config$doGreylist,
                        bRetrieveAnalysis=FALSE, bReduceObjects=TRUE, 
                        bParallel=DBA$config$RunParallel)
{
  
  if(!is(DBA,"DBA") && !is(DBA,"list")) {
    if(is(DBA,"character") || is(DBA,"data.frame")) {
      message("Loading sample sheet...")
      DBA <- dba(sampleSheet = DBA) 
    }
  }
  
  DBA <- pv.check(DBA,bCheckEmpty=TRUE)
  
  if(bRetrieveAnalysis!=FALSE) {
    if(bRetrieveAnalysis==TRUE || bRetrieveAnalysis==DBA_DESEQ2) {
      if(!is.null(DBA$DESeq2$DEdata)){
        return(DBA$DESeq2$DEdata)
      }
    }
    if (bRetrieveAnalysis==TRUE || bRetrieveAnalysis==DBA_EDGER){
      if(!is.null(DBA$edgeR$DEdata)) {
        return(DBA$edgeR$DEdata)
      } 
    }
    if(is.null(DBA$design)) {
      stop("Unable to return DEObject: design must be present, and analysis run.")
    } else {
      stop("Unable to return DEObject: run dba.analyze() with bRetrieveAnalysis=FALSE first.")
    }
  }
  
  if(is.null(bBlacklist)) {
    bBlacklist <- TRUE
  }
  if(is.null(bGreylist)) {
    bGreylist <- TRUE
  }
  
  doblacklist <- dogreylist <-  FALSE
  if(is.null(DBA$blacklist)) {
    if(bBlacklist == TRUE) {
      doblacklist <- TRUE
    }
  }
  
  if(is.null(DBA$greylist)) {
    if(bGreylist == TRUE) {
      if(!pv.noControls(DBA$class["bamControl",])) {
        dogreylist <- TRUE
      }
    } 
  }
  
  res <- NULL
  
  if(doblacklist || dogreylist) {
    message("Applying Blacklist/Greylists...")
    res <- tryCatch(
      dba.blacklist(DBA, blacklist=doblacklist, greylist=dogreylist),
      error=function(x){message("Blacklist error: ",x);NULL}
    )
    if(is.null(res)){
      message("Unable to apply Blacklist/Greylist.")
      return(DBA)
    } else {
      DBA <- res
    }
  }
  
  if(sum(DBA$class[DBA_CALLER,] %in% c("source","counts"))==0) {
    message("Forming consensus peakset and counting reads...")
    res <- tryCatch(
      dba.count(DBA, bParallel=bParallel),
      error=function(x){message("Count error: ",x);NULL}
    ) 
    if(is.null(res)){
      message("Unable to count overlapping reads.")
      return(DBA)
    } else {
      DBA <- res
    }
  }
  
  
  if (DBA_EDGER %in% method) {
    if(is.null(DBA$norm$edgeR)) {
      message("Normalize edgeR with defaults...")
      res <- tryCatch(
        dba.normalize(DBA, method=DBA_EDGER),
        error=function(x){message("Normalize error: ",x);NULL}
      )
      if(is.null(res)){
        message("Unable to normalize datset with edgeR.")
        return(DBA)
      } else {
        DBA <- res
      }
    }
  }
  
  if (DBA_DESEQ2 %in% method) {
    if(is.null(DBA$norm$DESeq2)) {
      message("Normalize DESeq2 with defaults...")
      res <- tryCatch(
        dba.normalize(DBA, method=DBA_DESEQ2),
        error=function(x){message("Normalize error: ",x);NULL}
      )
      if(is.null(res)){
        message("Unable to normalize datset with DESeq2.")
        return(DBA)
      } else {
        DBA <- res
      }
    }
  }
  
  if(!missing(design)) {
    message("Setting design...")
    res <- tryCatch(
      dba.contrast(DBA, design=design),
      error=function(x){message("Contrast error: ",x);NULL}
    )
    if(is.null(res)){
      message("Unable to set design.")
      return(DBA)
    } else {
      DBA <- res
    }
  }
  
  if(is.null(DBA$contrasts)) {
    if(missing(design) && is.null(DBA$design)) {
      message("Forming default model design and contrast(s)...")
    } else {
      message("Setting default contrast(s)...")      
    }
    res <- tryCatch(
      dba.contrast(DBA),
      error=function(x){message("Contrast error: ",x);NULL}
    )
    if(is.null(res)){
      message("Unable to establish model design and contrasts(s). 
Check that here are at least two groups each with at least three samples, 
and that the default model is of full rank, and call dba.contrast() explicitly.")
      return(DBA)
    } else {
      DBA <- res
    }
  }
  
  if(is.null(DBA$config$edgeR$bTagwise)) {
    bTagwise <- TRUE
  } else {
    bTagwise <- DBA$config$edgeR$bTagwise 
  }
  
  message("Analyzing...")
  res <- tryCatch(
    pv.DBA(DBA, method,bTagwise=bTagwise,
           minMembers=3, bParallel=bParallel),
    error=function(x){message("Analyze error: ",x);NULL}
  )
  if(is.null(res)){
    message("Unable to complete analysis.")
    return(DBA)
  } 
  
  if(bReduceObjects) {
    if(!is.null(res$contrasts)) {
      res$contrasts <- lapply(res$contrasts,pv.stripDBA)
    }      	
  }
  
  if(DBA$config$bCorPlot){
    warn <- TRUE
    rep <- pv.DBAreport(res,contrast=1,method=method[1],th=0.05,bSupressWarning=TRUE)
    if(!is.null(rep)) {
      if(!is.null(dim(rep))) {
        if(nrow(rep)>1) {
          warn=FALSE
          x <- try(dba.plotHeatmap(res,contrast=1,method=method[1],
                                   correlations=TRUE),silent=TRUE)
        }	
      }
    }
    if(warn) {
      warning('No correlation heatmap plotted -- contrast 1 does not have enough differentially bound sites.')	
    }
  }
  
  if(!is(res,"DBA")) {
    class(res) <- "DBA"
  }
  
  if(!is.null(DBA$ChIPQCobj)) {
    res <- checkQCobj(DBA$ChIPQCobj,res)
  }
  
  gc(verbose=FALSE)
  
  return(res)
}

###########################################################
## dba.report -- generate report for a contrast analysis ##
###########################################################

dba.report <- function(DBA, contrast, method=DBA$config$AnalysisMethod, 
                       th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                       fold=0, bNormalized=TRUE, bFlip=FALSE, precision,
                       bCalled=FALSE, bCounts=FALSE, bCalledDetail=FALSE,
                       bDB, bNotDB, bAll=TRUE, bGain=FALSE, bLoss=FALSE,
                       file,initString=DBA$config$reportInit,ext='csv',
                       DataType=DBA$config$DataType) 
  
{
  
  DBA <- pv.check(DBA,bCheckEmpty=TRUE) 
  
  if(DataType==DBA_DATA_SUMMARIZED_EXPERIMENT) {
    bCounts <- TRUE
  }
  
  if(!missing(bDB) | !missing(bNotDB) | bGain | bLoss) {
    if(missing(bDB)) {
      if(bGain | bLoss) {
        bDB <- TRUE
      } else {
        bDB <- FALSE
      }
    }
    if(missing(bNotDB)) {
      bNotDB <- FALSE
    }
    message("Generating report-based DBA object...")
    res <- pv.resultsDBA(DBA,contrasts=contrast,methods=method,
                         th=th,bUsePval=bUsePval,fold=fold,
                         bDB=bDB,bNotDB=bNotDB,bUp=bGain,bDown=bLoss,bAll=bAll,
                         bFlip=bFlip)
    res$resultObject <- TRUE
    return(res)                    
  }
  
  if(missing(contrast)) {
    contrast=1
  }
  
  if(missing(precision)) {
    if(DataType==DBA_DATA_SUMMARIZED_EXPERIMENT) {
      precision <- 0
    } else {
      precision <- 2:3
    }
  }
  
  res <- pv.DBAreport(pv=DBA,contrast=contrast,method=method,th=th,bUsePval=bUsePval,
                      bCalled=bCalled,bCounts=bCounts,bCalledDetail=bCalledDetail,
                      file=file,initString=initString,bNormalized=bNormalized,
                      ext=ext,minFold=fold,
                      bFlip=bFlip, precision=precision) 
  
  if(DataType==DBA_DATA_SUMMARIZED_EXPERIMENT) {
    DBA <- pv.getPlotData(DBA,contrast=contrast,report=res,
                          method=method,th=th,bUsePval=bUsePval,bNormalized=TRUE,
                          bFlip=bFlip, precision=0)
    res <- pv.DBA2SummarizedExperiment(DBA,report=res)
    return(res)
  }
  
  if(DataType!=DBA_DATA_FRAME) {
    res <- pv.peaks2DataType(res,DataType)
  }
  
  return(res)	
  
}                      

################################################
## dba.plotHeatmap -- Heatmap with clustering ##
################################################

dba.plotHeatmap <- function(DBA, attributes=DBA$attributes, maxSites=1000,
                            minval, maxval,
                            contrast, method=DBA$config$AnalysisMethod, 
                            th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                            report, score, bLog=TRUE, mask, sites, sortFun=sd, 
                            correlations=TRUE, olPlot=DBA_COR, 
                            ColAttributes, RowAttributes, colSideCols, 
                            rowSideCols=colSideCols,
                            margin=10, colScheme="Greens", distMethod="pearson",
                            ...)
{
  DBA <- pv.check(DBA,bCheckEmpty=TRUE)
  
  if( (missing(contrast) || !missing(mask)) && !missing(score) ) {
    saveCorPlot <- DBA$config$bCorPlot
    DBA$config$bCorPlot <- FALSE
    DBA <- dba.count(DBA,peaks=NULL,score=score)	
    DBA$config$bCorPlot <- saveCorPlot
  }
  
  #mask <- pv.setMask(DBA,mask,contrast)
  
  if(!missing(contrast)) {
    if(!missing(report)) {
      report <- pv.DataType2Peaks(report)
    }   	
    DBA <- pv.getPlotData(DBA,attributes=attributes,contrast=contrast,report=report,
                          method=method,th=th,bUsePval=bUsePval,bNormalized=TRUE,
                          bPCA=FALSE,bLog=FALSE,minval=minval,maxval=maxval,mask=mask)
    contrast <- 1
    mask <- NULL                     
  }
  
  if(bLog) {
    vectors <- DBA$binding[,4:ncol(DBA$binding)]
    vectors[is.na(vectors)] <- 0
    if(max(vectors) > 1) { # must have positive counts to do log
      vectors[vectors<1]=1
      if(missing(minval)) {
        minval <- 0
      } else {
        minval <- max(0,minval)
      }
      vectors <- log2(vectors)
      DBA$binding[,4:ncol(DBA$binding)] <- vectors
    }
  }
  
  if(length(correlations)==1 & ((correlations[1] == DBA_OLAP_ALL) | 
                                (correlations[1] == TRUE)))  {
    if(nrow(DBA$binding)>1) {
      if(!missing(sites)) {
        if(is.logical(sites)) {
          sites <- which(sites)	
        }	   
      } 	
      correlations <- pv.occupancy(DBA, mask=mask, sites=sites, Sort='cor', 
                                   bCorOnly=TRUE,CorMethod=distMethod)
    } else {
      warning('No correlation heatmap plotted -- contrast does not have enough differentially bound sites.')	
      return(NULL)   	     	
    }
  }
  
  if(correlations[1]!=FALSE) {
    res <- pv.plotHeatmap(DBA, attributes=attributes, overlaps=correlations, 
                          olPlot=olPlot, mask=mask,
                          ColScheme=colScheme, distMeth=distMethod, 
                          bReorder=TRUE, contrast=contrast,
                          RowAttributes=RowAttributes,ColAttributes=ColAttributes,
                          rowSideCols=rowSideCols,colSideCols=colSideCols,
                          minval=minval, maxval=maxval, margins=c(margin,margin),
                          ...)
  } else {
    
    if(!missing(contrast)) {
      if(nrow(DBA$binding)<2) { 	
        warning('No heatmap plotted -- contrast does not have enough differentially bound sites.')	
        return(NULL)   	     	
      }
      res <- pv.plotHeatmap(DBA, numSites=maxSites, attributes=attributes,
                            contrast=contrast,
                            RowAttributes=RowAttributes,ColAttributes=ColAttributes,
                            rowSideCols=rowSideCols,colSideCols=colSideCols,
                            ColScheme=colScheme, distMeth=distMethod, 
                            margins=c(margin,margin), ...)
      maxSites <- min(maxSites, nrow(DBA$binding))
      res <- data.frame(DBA$binding[1:maxSites,][res$rowInd,c(1:3,3+res$colInd)])
      if(!is.character(res[1,1])) {
        res[,1] <- DBA$chrmap[res[,1]]
      }
      res <- pv.peaks2DataType(res,datatype=DBA_DATA_GRANGES) 
    } else {
      
      if(is(sortFun,"function")) {
        savevecs <- DBA$binding
        DBA <- pv.sort(DBA, sortFun, mask=mask)
      } else savevecs <- NULL
      
      res <- pv.plotHeatmap(DBA, numSites=maxSites, attributes=attributes, 
                            mask=mask, sites=sites,
                            RowAttributes=RowAttributes,ColAttributes=ColAttributes,
                            rowSideCols=rowSideCols,colSideCols=colSideCols,
                            minval=minval, maxval=maxval, ColScheme=colScheme, 
                            distMeth=distMethod, 
                            margins=c(margin,margin),...)
      maxSites <- min(maxSites, nrow(DBA$binding))
      res <- data.frame(DBA$binding[1:maxSites,][res$rowInd,c(1:3,3+res$colInd)])
      if(!is.character(res[1,1])) {
        res[,1] <- DBA$chrmap[res[,1]]
      }
      res <- pv.peaks2DataType(res,datatype=DBA_DATA_GRANGES) 
      if(!is.null(savevecs)) {
        DBA$binding <- savevecs
      }
    }
  }
  
  invisible(res)	
}

#######################################################
## dba.plotPCA -- Principal Components Analysis plot ##
#######################################################

dba.plotPCA <- function(DBA, attributes, minval, maxval,
                        contrast, method=DBA$config$AnalysisMethod, 
                        th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                        report, score, bLog=TRUE, mask, sites, label, cor=FALSE,
                        b3D=FALSE, vColors, dotSize, labelSize, labelCols, 
                        components=1:3, ...)
  
{
  DBA <- pv.check(DBA,bCheckEmpty=TRUE)
  
  mask <- pv.setMask(DBA,mask,contrast)
  
  if(missing(contrast) && !missing(score)) {
    saveCorPlot <- DBA$config$bCorPlot
    DBA$config$bCorPlot <- FALSE
    DBA <- dba.count(DBA,peaks=NULL,score=score)	
    DBA$config$bCorPlot <- saveCorPlot
  } else if (!missing(score)) {
    warning('score parameter ignored when contrast is specified')	
  }  
  
  if(missing(label)) {
    label <- NULL
  } else {
    if(missing(labelSize)) {
      labelSize=.8
    }
    if(missing(labelCols)) {
      labelCols="black"
    }
  }
  
  if(!missing(contrast)){
    if(missing(attributes)) {
      attributes <- DBA_GROUP
    }
    if(missing(dotSize)) {
      dotSize <- NULL
    }
    if(!missing(report)) {
      report <- pv.DataType2Peaks(report)
    }  	  
    DBA <- pv.getPlotData(DBA,attributes=attributes,contrast=contrast,report=report,
                          method=method,th=th,bUsePval=bUsePval,bNormalized=TRUE,
                          bPCA=TRUE,minval=minval,maxval=maxval,mask=mask,bLog=FALSE)                     
    if(attributes[1] == PV_GROUP) {
      if( length(unique(DBA$class[PV_ID,])) != ncol(DBA$class) ){
        attributes <- PV_ID
      } else {
        attributes <- pv.attributePCA(DBA)
      }
    }
    res <- pv.plotPCA(DBA,attributes=attributes,size=dotSize,cor=cor,
                      b3D=b3D,vColors=vColors,label=label,bLog=bLog,
                      labelSize=labelSize,labelCols=labelCols,
                      comps=components,...)
  } else {
    if(missing(attributes)) {
      attributes <- pv.attributePCA(DBA)
    }
    
    res <- pv.plotPCA(DBA, attributes=attributes, size=dotSize, mask=mask, 
                      sites=sites, b3D=b3D, vColors=vColors, label=label,
                      bLog=bLog,labelSize=labelSize,labelCols=labelCols,
                      comps=components,...)  
  }
  
  invisible(res)	
}

#############################
## dba.plotBox --Boxplots  ##
#############################

dba.plotBox <- function(DBA, contrast=1, method=DBA$config$AnalysisMethod, 
                        th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                        bNormalized=TRUE, attribute=DBA_GROUP, mask,
                        bAll=FALSE, bAllIncreased=FALSE, bAllDecreased=FALSE, 
                        bDB=TRUE, bDBIncreased=TRUE, bDBDecreased=TRUE,
                        pvalMethod=wilcox.test,  bReversePos=FALSE, attribOrder, 
                        vColors, varwidth=TRUE, notch=TRUE, ...) 
  
{
  DBA <- pv.check(DBA,bCheckEmpty=TRUE)
  
  if(contrast > length(DBA$contrasts)) {
    stop('Supplied contrast greater than number of contrasts')	
  }
  
  res <- pv.plotBoxplot(DBA, contrast=contrast, method=method, th=th, 
                        bUsePval=bUsePval, bNormalized=bNormalized,
                        attribute=attribute,mask=mask,bAll=bAll, 
                        bAllIncreased=bAllIncreased, bAllDecreased=bAllDecreased, 
                        bDB=bDB, bDBIncreased=bDBIncreased, 
                        bDBDecreased=bDBDecreased,
                        pvalMethod=pvalMethod,  bReversePos=bReversePos, 
                        attribOrder=attribOrder, vColors=vColors, 
                        varwidth=varwidth, notch=notch, ...)
  
  invisible(res)	
}

#########################################
## dba.plotMA -- MA or XY scatter plot ##
#########################################

dba.plotMA <- function(DBA, contrast=1, method=DBA$config$AnalysisMethod, 
                       th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                       fold=0, bNormalized=TRUE,
                       factor="", bFlip=FALSE, bXY=FALSE, dotSize=.45, 
                       bSignificant=TRUE, highlight=NULL,
                       bSmooth=TRUE, bLoess=TRUE, xrange, yrange, ...)
  
{
  DBA <- pv.check(DBA,bCheckEmpty=TRUE)
  
  res <- pv.DBAplotMA(DBA, contrast=contrast, method=method, bMA=!bXY, bXY=bXY,
                      th=th, bUsePval=bUsePval, fold=fold,
                      facname=factor, bNormalized=bNormalized, cex=dotSize, 
                      bSignificant=bSignificant, bSmooth=bSmooth,bFlip=bFlip,
                      highlight=highlight, xrange=xrange, yrange=yrange,
                      bLoess=bLoess,...)
  
  invisible(res)
}

#####################################
## dba.plotVolcano -- Volcano plot ##
#####################################

dba.plotVolcano <- function(DBA, contrast=1, method=DBA$config$AnalysisMethod, 
                            th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                            fold=0, factor="", bFlip=FALSE, 
                            bLabels=FALSE, maxLabels=50, dotSize=1)
  
{
  DBA <- pv.check(DBA,bCheckEmpty=TRUE)
  
  res <- pv.DBAplotVolcano(DBA, contrast=contrast, method=method, 
                           th=th, bUsePval=bUsePval, fold=fold,
                           facname=factor, dotSize=dotSize,
                           bFlip=bFlip, bLabels=bLabels, maxLabels=maxLabels)
  
  invisible(pv.peaks2DataType(res,datatype=DBA_DATA_GRANGES))
}


#############################################
## dba.plotProfile -- Profile heatmap plot ##
#############################################

dba.plotProfile <- function(Object, samples, sites,
                            scores="Score", labels, # annotate=TRUE, 
                            normalize=TRUE, merge=DBA_REPLICATE,
                            maxSites=1000, absScores=TRUE, 
                            doPlot=is(Object,"profileplyr"),
                            ...)
{
  Object <- pv.check(Object,bCheckEmpty=TRUE)
  
  res <- pv.plotProfile(pv=Object, mask=samples, sites=sites, maxSites=maxSites, 
                        scores=scores, annotate=FALSE, labels=labels,
                        normalize=normalize,merge=merge, absScores=absScores,
                        doPlot=doPlot, returnVal="profileplyr",
                        ...) 
  
  invisible(res)
}


####################################################
## dba.plotVenn -- Venn diagram plots of overlaps ##
####################################################

dba.plotVenn <- function(DBA, mask, overlaps, label1, label2, label3, label4,
                         main, sub, 
                         contrast, method=DBA$config$AnalysisMethod, 
                         th=DBA$config$th, bUsePval=DBA$config$bUsePval,
                         bDB=TRUE, bNotDB, bAll=TRUE, bGain=FALSE, bLoss=FALSE,
                         labelAttributes, DataType=DBA$config$DataType)
{
  DBA <- pv.check(DBA,bCheckEmpty=TRUE)
  
  newDBA <- NULL
  
  if (!missing(overlaps)){
    
    if(missing(label1)) {
      label1 <- "A"
    }
    if(missing(label2)) {
      label2 <- "B"
    }
    if(missing(label3)) {
      label3 <- "C"
    }
    if(missing(label4)) {
      label4 <- "D"
    }
    
    overlaps <- lapply(overlaps,pv.DataType2Peaks)
    
  } else if(!missing(contrast)){
    if(max(contrast)>length(DBA$contrasts)) {
      stop('Contrast greater than number of contrasts.')   
    }
    newDBA <- dba.report(DBA,contrast=contrast,method=method,th=th,bUsePval=bUsePval,
                         bDB=bDB, bNotDB=bNotDB,bAll=bAll, bGain=bGain, bLoss=bLoss)
    
    if(is.null(newDBA$peaks)){
      stop('No peaksets meet specified criteria.')   
    }
    if(missing(mask)) {
      mask <- 1:length(newDBA$peaks)
      if(length(mask)>4){
        stop('Too many peaksets meet specified criteria.')
      }
      if(length(mask)==1) {
        stop('Only one peakset meets specified criteria.')
      }
    } else {
      if(is.logical(mask)) {
        if(length(mask)!=length(newDBA$peaks)) {
          stop('Logical mask doe not have same number of elements as there are peaksets.')
        }
        if(sum(mask)>4) {
          stop('Too many peaksets in mask.')
        } else if(length(mask)>4){
          stop('Too many peaksets in mask.')
        }
        mask <- which(mask)
      }
      if(length(mask)>length(newDBA$peaks) | max(mask)>length(newDBA$peaks) ) {
        stop('Peakset specified in mask is out of range.')
      }
    }
    
    overlaps <- dba.overlap(newDBA,mask,mode=DBA_OLAP_PEAKS,DataType=DBA_DATA_FRAME)
    
    res <- pv.whichPeaksets(newDBA,mask)
    
    if(missing(labelAttributes)) {
      labelAttributes=c(DBA_ID,DBA_FACTOR,DBA_TISSUE,DBA_CONDITION,DBA_TREATMENT)  
    }
    crec <- matrix(newDBA$class[labelAttributes,mask],length(labelAttributes),length(mask))
    if(length(mask)==2) {
      labels <- pv.namestrings(crec[,1],crec[,2])
    } else if (length(mask)==3) {
      labels <- pv.namestrings(crec[,1],crec[,2],crec[,3])         
    } else {
      labels <- pv.namestrings(crec[,1],crec[,2],crec[,3],crec[,4])         
    }
    
    if(missing(label1)) {
      label1 <- labels$n1
    }
    if(missing(label2)) {
      label2 <- labels$n2
    }
    if(missing(label3)) {
      label3 <- labels$n3
    }
    if(missing(label4)) {
      label4 <- labels$n4
    }
    
    if (missing(sub)) {
      sub <- labels$tstring
    }
    
  } else if(!missing(mask)) {
    
    if(is.logical(mask)) {
      if(length(mask)!=length(DBA$peaks)) {
        stop('Logical mask does not have same number of elements as there are peaksets.')
      }
      if(sum(mask)>4) {
        stop('Too many peaksets in mask.')
      }
      mask <- which(mask)
    } else if(length(mask)>4){
      stop('Too many peaksets in mask.')
    }      
    
    overlaps <- dba.overlap(DBA,mask,mode=DBA_OLAP_PEAKS,DataType=DBA_DATA_FRAME)
    
    res <- pv.whichPeaksets(DBA,mask)
    if(missing(labelAttributes)) {
      if(pv.checkValue(DBA$resultObject,TRUE)) {
        labelAttributes <- c(DBA_ID,DBA_FACTOR,DBA_TISSUE,DBA_CONDITION,DBA_TREATMENT)  
      } else {
        labelAttributes <- DBA_ID
      }
    }
    crec <- matrix(DBA$class[labelAttributes,mask],length(labelAttributes),length(mask))
    if(length(mask)==2) {
      labels <- pv.namestrings(crec[,1],crec[,2])
    } else if (length(mask)==3) {
      labels <- pv.namestrings(crec[,1],crec[,2],crec[,3])         
    } else {
      labels <- pv.namestrings(crec[,1],crec[,2],crec[,3],crec[,4])         
    }
    if(missing(label1)) {
      label1 <- labels$n1
    }
    if(missing(label2)) {
      label2 <- labels$n2
    }
    if(missing(label3)) {
      label3 <- labels$n3
    }
    if(missing(label4)) {
      label4 <- labels$n4
    }
    if (missing(sub)) {
      sub <- labels$tstring
    }
  } else {
    stop("Must specify one of mask, overlaps, or contrast.")
  }
  
  if(missing(main)) {
    main <- "Binding Site Overlaps"
  }
  if (missing(sub)) {
    sub <- ""
  }
  
  pv.plotVenn(overlaps,label1=label1,label2=label2,label3=label3,label4=label4,
              main,sub)
  
  if(DataType == DBA_DATA_DBAOBJECT) {
    if(!is.null(newDBA)) {
      return(newDBA)
    } else {
      warning('No DBA object to return.')
    }
  } else {
    overlaps <- lapply(overlaps, pv.peaks2DataType, DataType)
    if(DataType == DBA_DATA_GRANGES) {
      for(i in 1:length(overlaps)) {
        if(!is.null(overlaps[[i]])) {
          if(is.null(overlaps[[i]]$score)) {
            overlaps[[i]]$score <- rowMeans(data.frame(mcols(overlaps[[i]])), 
                                            na.rm=TRUE)
          }
        }
        else {
          overlaps[[i]] <- GRanges(NULL)
        }
      }
      overlaps <- GRangesList(overlaps)
    }
    invisible(overlaps)
  }   
}

###################################
## dba.show -- List DBA metadata ##
###################################

dba.show <- function(DBA, mask, attributes, bContrasts=FALSE, bDesign=FALSE,
                     th=DBA$config$th) 
{
  DBA <- pv.check(DBA)
  
  res <- pv.list(DBA, mask=mask, bContrasts=bContrasts, bDesign=bDesign,
                 attributes=attributes, th=th)
  
  return(res)
}

#######################################
## dba.mask -- mask samples or sites ##
#######################################

dba.mask <- function(DBA, attribute, value, combine='or', mask, merge='or', 
                     bApply=FALSE,peakset, minValue=-1)
  
{
  DBA <- pv.check(DBA)
  
  if(missing(peakset)) {
    
    res <- pv.mask(DBA, attribute=attribute, value=value, combine=combine,
                   mask=mask, merge=merge, bApply=bApply)
    
  } else {
    
    res <- pv.whichSites(DBA, pnum=peakset, combine=combine, minVal=minValue)
    
  }
  
  return(res)
  
}

###################################
## dba.save -- save a DBA object ##
###################################

dba.save <- function(DBA, file='DBA', dir='.', pre='dba_', ext='RData',
                     bRemoveAnalysis=FALSE, bRemoveBackground=FALSE, 
                     bCompress=FALSE)
{
  
  if(is(DBA,"ChIPQCexperiment")) {
    saveChIPQC <- DBA
    DBA <- DBA@DBA
  } else saveChIPQC <- NULL
  
  if(nrow(DBA$class)<DBA_TREATMENT) {
    DBA$class <- rbind(DBA$class,'')
    rownames(DBA$class)[DBA_TREATMENT]='Treatment'	
  }
  
  
  DBA$binding <- NULL	
  DBA$values  <- NULL
  DBA$hc      <- NULL
  DBA$pc      <- NULL
  
  DBA$config$lsfInit      <- NULL
  DBA$config$parallelInit <- NULL
  DBA$config$initFun      <- NULL
  DBA$config$paramFun     <- NULL
  DBA$config$addjobFun    <- NULL
  DBA$config$lapplyFun    <- NULL
  DBA$config$wait4jobsFun <- NULL
  DBA$config$parallelInit <- NULL
  
  version <- strsplit(as.character(packageVersion("DiffBind")),".",fixed=TRUE)[[1]]
  
  DBA$config$Version1 <- version[1]
  DBA$config$Version2 <- version[2]
  DBA$config$Version3 <- version[3]
  
  if(bRemoveAnalysis) {
    DBA$DESeq2$DEdata <- NULL
    DBA$edgeR$DEdata <- NULL
  }
  
  if(bRemoveBackground) {
    DBA$norm$background$binned <- NULL
  }
  if(!is.null(saveChIPQC)) {
    saveChIPQC@DBA <- DBA
    DBA <- saveChIPQC
  }
  res <- pv.save(DBA,file=file ,
                 dir=dir, pre=pre, ext=ext,
                 compress=bCompress)
  
  return(res)
} 

###################################
## dba.load -- load a DBA object ##
###################################

dba.load <- function(file='DBA', dir='.', pre='dba_', ext='RData')
{
  
  res <- pv.load(file=file, dir=dir, pre=pre, ext=ext)
  
  if(is(res,"ChIPQCexperiment")) {
    saveChIPQC <- res
    res <- res@DBA
  } else saveChIPQC <- NULL
  
  if(!is.null(res$vectors)) {
    if(!is.null(res$allvectors)) {
      res$binding     <- res$vectors
      res$totalMerged <- nrow(res$allvectors)
      res$merged      <- res$allvectors[,1:3]
    } else {
      res$binding     <- res$vectors
      res$merged      <- res$vectors[,1:3]
      res$totalMerged <- nrow(res$merged)
    }
    res$vectors    <- NULL
    res$allvectors <- NULL
  }
  
  if(is.null(res$binding)) {
    if(is.null(res$minOverlap)) {
      minOverlap=2
    } else {
      minOverlap <- res$minOverlap
    }
    contrasts  <- res$contrasts
    called     <- res$called
    attributes <- res$attributes
    for(i in 1:length(res$peaks)) {
      if(is.factor(res$peaks[[i]][,1])) {
        res$peaks[[i]][,1] <- as.character(res$peaks[[i]][,1])
      }
    }
    res <- pv.vectors(res,minOverlap=minOverlap,bAllSame=pv.allSame(res),
                      merge=is.null(res$merged))
    res$contrasts  <- contrasts
    res$called     <- called
    res$attributes <- attributes
  }
  
  if(is.null(res$config$DataType)) {
    res$config$DataType=DBA_DATA_DEFAULT
    if(!is.null(res$config$RangedData)) {
      if(res$config$RangedData==FALSE) {
        res$config$DataType=DBA_DATA_FRAME   	
      } 
    }
  }
  
  if(is.null(res$config$savePrefix)) {
    res$config$savePrefix <- "dba_"
  }
  
  if(is.null(res$config$saveExt)) {
    res$config$saveExt <- "RData"
  }
  
  if(is.null(res$config$reportInit)) {
    res$config$reportInit <- "reports/DBA"
  }
  if(is.null(res$config$AnalysisMethod)){
    res$config$AnalysisMethod=DBA_DESEQ2
  }
  
  if(is.null(res$config$bCorPlot)){
    res$config$bCorPlot=FALSE
  }
  
  if(is(res$attributes,"function")) {
    res$attributes <- NULL
  }
  
  if(is.null(res$config$th)){
    res$config$th=0.05
  }
  if(is.null(res$config$bUsePval)){
    res$config$bUsePval=FALSE
  }   
  
  
  
  res$config$lsfInit      <- NULL
  res$config$parallelInit <- NULL
  res$config$initFun      <- NULL
  res$config$paramFun     <- NULL
  res$config$addjobFun    <- NULL
  res$config$lapplyFun    <- NULL
  res$config$wait4jobsFun <- NULL
  res$config$parallelInit <- NULL
  res$config$multicoreInit <- NULL
  
  res$config <- as.list(res$config)
  
  if (is.null(res$config$mapQCth)) {
    res$config$mapQCth=15   
  }
  
  if (is.null(res$config$fragmentSize)) {
    res$config$fragmentSize=125
  }   
  version <- strsplit(as.character(packageVersion("DiffBind")),".",fixed=TRUE)[[1]]
  res <- pv.version(res,version[1],version[2], version[3])
  
  if(!is(res,"DBA")) {
    class(res) <- "DBA"
  }
  
  if(!is.null(saveChIPQC)) {
    saveChIPQC@DBA <- res
    res <- saveChIPQC
  }
  gc(verbose=FALSE)
  return(res)
} 

##########################
## DBA object functions ##
##########################

print.DBA <- function(x,...){
  x <- pv.check(x)
  cat(sprintf("%s:\n",summary(x)))
  toshow <- dba.show(x)
  
  if(nrow(toshow) > 1) {
    if(length(unique(toshow$Reads))==1) {
      delcol <- which(colnames(toshow) %in% 'Reads')
      toshow <- toshow[,-delcol]
    }
    if(length(unique(toshow$Intervals))==1) {
      delcol <- which(colnames(toshow) %in% 'Intervals')
      toshow <- toshow[,-delcol]
    }
    if(length(unique(toshow$Caller))==1) {
      delcol <- which(colnames(toshow) %in% 'Caller')
      toshow <- toshow[,-delcol]
    }
  } else {
    toshow <-  toshow[,which(toshow[1,] != "")]
    toshow <-  toshow[,which(toshow[1,] != 0)]
  }
  print(toshow)
  
  if(!is.null(x$contrasts)){
    cat("\n")
    if(!is.null(x$design)) {
      cat(sprintf("Design: [%s] | ",dba.show(x,bDesign=TRUE)))
    }
    if(length(x$contrasts) == 1) {
      cat(sprintf("1 Contrast:") )
    } else {
      cat(sprintf("%d Contrasts:",length(x$contrasts)))
    }
    cat("\n")
    print(dba.show(x,bContrasts=TRUE,th=x$config$th))
  } else {
    if(!is.null(x$design)) {
      cat(sprintf("\nDesign: [%s]\n",dba.show(x,bDesign=TRUE)))
    }
  }
}

summary.DBA <- function(object,...) {
  if(is.null(object$binding)) {
    cat('Run dba first\n')
    return()
  }
  res <- sprintf("%d Samples, %d sites in matrix",
                 length(object$peaks),nrow(object$binding))
  if(nrow(object$binding) != object$totalMerged) {
    res <- sprintf("%s (%d total)",res,object$totalMerged)
  }
  return(res)
}

plot.DBA <- function(x,...){
  DBA <- pv.check(x)
  res <- dba.plotHeatmap(x,...)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(" >>> DiffBind 3.2.")
}

