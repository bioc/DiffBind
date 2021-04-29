\name{DBA object methods}
\alias{DBA-object}
\alias{DBA.object}
\alias{print.DBA}
\alias{summary.DBA}
\alias{plot.DBA}
\alias{DiffBind-config}
\alias{DBA-config}
\alias{DiffBind.config}
\alias{DBA.config}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standard S3 methods for DBA object
}
\description{
Standard S3 methods for DBA object.
}
\usage{
\method{print}{DBA}(x, \dots)
\method{summary}{DBA}(object, \dots)
\method{plot}{DBA}(x, \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
DBA object
}
\item{object}{
DBA object
}

\item{\dots}{
Arguments passed on to parent methods
}
}
\details{
S3 methods for DBA object from the \code{\link{DiffBind}} package.

DBA objects are usually constructed using the \code{\link{dba}} function.

There are a number of internal parameters that can be set, 
and defaults overridden, by setting \code{DBA$config} options:

\itemize{
\item {\code{DBA$config$AnalysisMethod:}}	
  {either \code{DBA_DESEQ2} or \code{DBA_EDGER}.}
  
\item {\code{DBA$config$th:}}   
  {default threshold for reporting and plotting analysis results.}
  
\item {\code{DBA$config$DataType:}}	
  {default class for peaks and reports
(\code{DBA_DATA_GRANGES, DBA_DATA_RANGEDDATA, or DBA_DATA_FRAME}).}

\item {\code{DBA$config$RunParallel:}}	
  {logical indicating if counting and analysis 
  operations should be run in parallel using multicore by default.}

\item {\code{DBA$config$cores:}}	
  {number of cores to use when performing multi-core parallel processing.}
  
\item {\code{DBA$config$minQCth:}}   
  {numeric, for filtering reads based on mapping 
quality score; only reads with a mapping quality score
greater than or equal to this will be counted.}

\item {\code{DBA$config$fragmentSize:}} 
  {numeric indicating mean fragment size for single-end counting.
Reads will be extended to this length before counting overlaps. 
May be a vector of lengths, one for each sample.}

\item {\code{DBA$config$bCorPlot:}}
  {logical indicating that a correlation heatmap 
should be plotted automatically}

\item {\code{DBA$config$ReportInit:}}  
  {string to append to the beginning of saved report file names.}
  
\item {\code{DBA$config$bUsePval:}}   
  {logical, default indicating whether to use FDR 
(\code{FALSE}) or p-values (\code{TRUE}).}

\item {\code{DBA$config$doBlacklist:}}  
  {logical, whether to attempt to find and apply
a blacklist if none is present when running \code{\link{dba.analyze}}}.

\item {\code{DBA$config$doGreylist}}   
  {logical, whether to attempt to generate and apply 
a greylist if none is present when running \code{\link{dba.analyze}}}.

\item {\code{DBA$config$DataType}}
The class of object for returned reports and peaksets:
\itemize{
\item \code{\link{DBA_DATA_GRANGES}}
\item \code{\link{DBA_DATA_RANGEDDATA}}
\item \code{\link{DBA_DATA_FRAME}}
\item \code{\link{DBA_DATA_SUMMARIZED_EXPERIMENT}}
}

\item {\code{DBA$config$mergeOverlap:}}	
  {The overlap (in basepairs) between peaks to merge when generating 
  a consensus peakset.
  A positive valuecontrols how many basepairs peaks 
  must overlap to be merged,
  while a negative value will result in non-overlapping peaks to be merged, 
  
  If absent, the default value of \code{1} will result in any peaks
  overlapping by at least one basepair to be merged into a single interval.}

\item {\code{DBA$config$design:}}  
  {When calling \code{\link{dba.contrast}}, if \code{design} parameter is
  missing, this will be used as the value for that parameter.}
  
\item {\code{DBA$config$edgeR$bTagwise:}}  
  {\code{logical} indicating if \code{edgeR::estimateGLMTagwiseDisp} should
  be called when performing an \code{edgeR} analysis. 
  If absent the default is \code{TRUE}, so setting this to \code{FALSE}
  prevents the tagwise dispersion estimate form being calculated.
}  

\item {\code{DBA$config$DESeq2$fitType:}}  
  {\code{logical} indicating the \code{fitType}
  to be used in \code{DESeq2::estimateDispersions} 
  when performing a \code{DESeq2} analysis. 
  If absent the default is \code{local}.
}  

\item {\code{DBA$config$savePrefix:}}  
  {When calling \code{\link{dba.save}} or \code{\link{dba.load}}, 
  this value (if present)
  will override the default value for the \code{pre} parameter.
  }
  
\item {\code{DBA$config$saveExt:}}  
  {When calling \code{\link{dba.save}} or \code{\link{dba.load}},
  this value (if present)
  will override the default value for the \code{ext} parameter.
  }  
  
\item {\code{DBA$config$greylist.pval:}}  
  {pvalue cutoff to use when generating a greylist
  using \code{GreyListChIP::calcThreshold}.
  If missing, the default is \code{0.999}
  } 
  
\item {\code{DBA$config$saveExt:}}  
  {When calling \code{\link{dba.save}}, this value (if present)
  will override the default value for the \code{ext} parameter.
  } 
  
\item {\code{DBA$config$yieldSize:}} 
  {yieldSize indicating how many reads to process at one time; default is 5000000. 
The lower this value, the less memory will be used, but the more time it will take to 
complete the count operation.}

\item {\code{DBA$config$intersectMode:}}
  {mode indicating which overlap algorithm to use; 
  default is \code{"IntersectionNotEmpty"}}
  
\item {\code{DBA$config$singleEnd:}} 
  {\code{logical} indicating if reads are single end; 
  if \code{NULL}, status will be automatically detected.}
  
\item {\code{DBA$config$fragments:}}
  {\code{logical} indicating how unmatched reads are counted; 
  default is \code{FALSE}.}

\item {\code{DBA$config$scanbamparam:}} 
  {\code{ScanBamParam} object to pass to \code{\link{summarizeOverlaps}}. 
  If present, \code{bRemoveDuplicates} is ignored.}
  
\item {\code{DBA$config$pp.style:}} 
  {Sets \code{style} parameter for \code{profileplyr::BamBigwig_to_chipProfile}
  when calling \code{\link{dba.plotProfile}.}}

\item {\code{DBA$config$pp.nOfWindows:}} 
  {Sets \code{nOfWindow} parameter for 
  \code{profileplyr::BamBigwig_to_chipProfile}
  when calling \code{\link{dba.plotProfile}.}}
  
\item {\code{DBA$config$bin_size:}} 
  {Sets \code{bin_size} parameter for 
  \code{profileplyr::BamBigwig_to_chipProfile}
  when calling \code{\link{dba.plotProfile}.}}

\item {\code{DBA$config$distanceAround:}} 
  {Sets \code{distanceAround} parameter for 
  \code{profileplyr::BamBigwig_to_chipProfile}
  when calling \code{\link{dba.plotProfile}.}}
  
\item {\code{DBA$config$distanceUp:}} 
  {Sets \code{distanceUp} parameter for 
  \code{profileplyr::BamBigwig_to_chipProfile}
  when calling \code{\link{dba.plotProfile}.}}
  
\item {\code{DBA$config$distanceDown:}} 
  {Sets \code{distanceDown} parameter for 
  \code{profileplyr::BamBigwig_to_chipProfile}
  when calling \code{\link{dba.plotProfile}.}}
  
\item {\code{DBA$config$id:}}  
  {\code{character} string to use to replace \code{"ID"}
  when displaying a \code{DBA} object (\code{\link{dba.show}})}
  
\item {\code{DBA$config$group:}}  
  {\code{character} string to use to replace \code{"Group"}
  when displaying a \code{DBA} object (\code{\link{dba.show}})}

\item {\code{DBA$config$tissue:}}  
  {\code{character} string to use to replace \code{"Tissue"}
  when displaying a \code{DBA} object (\code{\link{dba.show}})}

\item {\code{DBA$config$factor:}}  
  {\code{character} string to use to replace \code{"Factor"}
  when displaying a \code{DBA} object (\code{\link{dba.show}})}
  
\item {\code{DBA$config$condition:}}  
  {\code{character} string to use to replace \code{"Condition"}
  when displaying a \code{DBA} object (\code{\link{dba.show}})}

\item {\code{DBA$config$treatment:}}  
  {\code{character} string to use to replace \code{"Treatment"}
  when displaying a \code{DBA} object (\code{\link{dba.show}})}

\item {\code{DBA$config$replicate:}}  
  {\code{character} string to use to replace \code{"Replicate"}
  when displaying a \code{DBA} object (\code{\link{dba.show}})}
  
\item {\code{DBA$config$caller:}}  
  {\code{character} string to use to replace \code{"Caller"}
  when displaying a \code{DBA} object (\code{\link{dba.show}})}

\item {\code{DBA$config$reads:}}  
  {\code{character} string to use to replace \code{"Reads"}
  when displaying a \code{DBA} object (\code{\link{dba.show}})}
  
}
}

%\value{
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(tamoxifen_peaks)
tamoxifen
data(tamoxifen_counts)
tamoxifen
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
