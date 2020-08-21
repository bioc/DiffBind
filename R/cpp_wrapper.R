################################################
## cpp_wrapper.R -- wrap C/C++ functions in R ##
## first created: 2014 03 12                  ##
## Gord Brown                                 ##
## CR-UK Cambridge Institute                  ##
################################################

################################################
## Entry Points to C/C++ Routines             ##
################################################

## cpp_count_reads -- count reads on intervals

cpp_count_reads <- function(bamfile,insertLength,fileType,bufferSize,
                            intervals,bWithoutDupes,summits,minMappingQual=0,
                            minVal=1) {
  icount <- length(intervals[[1]])
  counts <- vector(mode="integer",length=icount)
  if (summits >= 0) { # RJS 2/7/2020: change from !missing(summits)) {
    summits.vec <- vector(mode="integer",length=icount)
    heights.vec <- vector(mode="integer",length=icount)
    bSummits = TRUE
  } else {
    summits.vec <- vector()
    heights.vec <- vector()
    bSummits = FALSE
  }
  bamfile <- path.expand(bamfile)
  libsize <- .Call("croi_count_reads",PACKAGE='DiffBind',
                   bamfile,
                   as.integer(insertLength),
                   as.integer(fileType),
                   as.integer(bufferSize),
                   as.integer(minMappingQual),
                   as.character(intervals[[1]]),
                   as.integer(intervals[[2]]),
                   as.integer(intervals[[3]]),
                   as.integer(icount),
                   as.logical(bWithoutDupes),
                   as.logical(bSummits),
                   counts,
                   summits.vec,
                   heights.vec)
  
  # RJS 10 July 2020 -- Added minVal parameter  
  mincounts <- counts
  mincounts[mincounts<minVal] <- minVal
  
  widths = intervals[,3] - intervals[,2]
  rpkm = (mincounts/(widths/1000))/(libsize/1E6)
  
  result <- list(counts=counts,rpkm=rpkm,libsize=libsize)
  if (bSummits==TRUE) {
    result$summits <- summits.vec;
    result$heights <- heights.vec;
  }
  return(result)
}
