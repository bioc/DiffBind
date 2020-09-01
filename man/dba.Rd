\name{dba}
\alias{dba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct a DBA object
}
\description{
Constructs a new DBA object from a sample sheet,  or based on an existing DBA object
}
\usage{
dba(DBA,mask, minOverlap=2,
    sampleSheet="dba_samples.csv", 
    config=data.frame(AnalysisMethod=DBA_DESEQ2,th=0.05,
                      DataType=DBA_DATA_GRANGES, RunParallel=TRUE, 
                      minQCth=15, fragmentSize=125, 
                      bCorPlot=FALSE, reportInit="DBA",
                      bUsePval=FALSE, design=TRUE),
    peakCaller="raw", peakFormat, scoreCol, bLowerScoreBetter, 
    filter, skipLines=0, 
    bAddCallerConsensus=FALSE, 
    bRemoveM=TRUE, bRemoveRandom=TRUE, 
    bSummarizedExperiment=FALSE,
    attributes, dir)  
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Optional or defaulted arguments}{}

  \item{DBA}{
existing DBA object -- if present, will return a fully-constructed DBA object based on the passed one,
using criteria specified in the \code{mask} and/or \code{minOverlap} parameters. 
If missing, will create a new DBA object based on the \code{sampleSheet}.
}
\item{mask}{
logical or numerical vector indicating which peaksets to include in the resulting model 
if basing DBA object on an existing one. 
See  \code{\link{dba.mask}}.
}
\item{minOverlap}{
only include peaks in at least this many peaksets in the main binding matrix
if basing DBA object on an existing one. 
If \code{minOverlap} is between zero and one, peak will be included from at least this proportion of peaksets.
}
\item{sampleSheet}{
data frame containing sample sheet, or file name of sample sheet to load (ignored if DBA is specified). 
Columns names in sample sheet may include:
\itemize{
\item {\code{SampleID:}}	{Identifier string for sample. 
Must be unique for each sample.}
\item {\code{Tissue:}}	{Identifier string for tissue type}
\item {\code{Factor:}}	{Identifier string for factor}
\item {\code{Condition:}}	{Identifier string for condition}
\item {\code{Treatment:}}	{Identifier string for treatment}
\item {\code{Replicate:}}	{Replicate number of sample}
\item {\code{bamReads:}}	{file path for bam file containing aligned reads for ChIP sample}
\item {\code{bamControl:}}	{file path for bam file containing aligned reads for control sample}
\item {\code{ControlID:}}	{Identifier string for control sample}
\item {\code{Peaks:}}	{path for file containing peaks for sample. format determined by PeakCaller field or caller parameter}
\item {\code{PeakCaller:}}	{Identifier string for peak caller used. 
If Peaks is not a bed file, this will determine how the Peaks file is parsed. 
If missing, will use default peak caller specified in caller parameter. Possible values:}
\itemize{
\item {\dQuote{raw}:}	{text file file; peak score is in fourth column}
\item {\dQuote{bed}:}	{.bed file; peak score is in fifth column}
\item {\dQuote{narrow}:}  {default peak.format: narrowPeaks file}
\item {\dQuote{macs}:}	{MACS .xls file}
\item {\dQuote{swembl}:}	{SWEMBL .peaks file}
\item {\dQuote{bayes}:}	{bayesPeak file}
\item {\dQuote{peakset}:}	{peakset written out using pv.writepeakset}
\item {\dQuote{fp4}:}	{FindPeaks v4}
}
\item {\code{PeakFormat:}}	{string indicating format for peak files; see PeakCaller and  \code{\link{dba.peakset}}}
\item {\code{ScoreCol:}}	{column in peak files that contains peak scores}
\item {\code{LowerBetter:}}	{logical indicating that lower scores signify better peaks}
\item{\code{Counts:}} {file path for externally computed read counts; see \code{\link{dba.peakset}}
(\code{counts} parameter)}
}
For sample sheets loaded from a file, the accepted formats are comma-separated values 
(column headers, followed by one line per sample), or Excel-formatted spreadsheets (\code{.xls} 
or \code{.xlsx} extension). 
Leading and trailing white space will be removed from all values, with a warning.
}
  \item{config}{
data frame containing configuration options, 
or file name of config file to load when constructing a new DBA object from a sample sheet. 
\code{NULL} indicates no config file. 
Relevant fields include:
\itemize{
\item {\code{AnalysisMethod:}}	{either \code{DBA_DESEQ2} or \code{DBA_EDGER}.}
\item {\code{th:}}   {default threshold for reporting and plotting analysis results.}
\item {\code{DataType:}}	{default class for peaks and reports (\code{DBA_DATA_GRANGES, DBA_DATA_RANGEDDATA, or DBA_DATA_FRAME}).}
\item {\code{RunParallel:}}	{logical indicating if counting and analysis operations should be run in parallel using multicore by default.}
\item {\code{minQCth:}}   {numeric, for filtering reads based on mapping quality score; only reads with a mapping quality score
greater than or equal to this will be counted.}
\item {\code{fragmentSize:}} {numeric with mean fragment size.
Reads will be extended to this length before counting overlaps. 
May be a vector of lengths, one for each sample.}
\item {\code{bCorPlot:}} {logical indicating that a correlation heatmap should be plotted automatically}
\item {\code{ReportInit:}}   {string to append to the beginning of saved report file names.}
\item {\code{bUsePval:}}   {logical, default indicating whether to use FDR (\code{FALSE}) or p-values (\code{TRUE}).}
}
}
  \item{peakCaller}{
if a \code{sampleSheet} is specified, the default peak caller that will be used 
if the \code{PeakCaller} column is absent.
}
 \item{peakFormat}{
if a \code{sampleSheet} is specified, the default peak file format that will be used 
if the \code{PeakFormat} column is absent.
}
 \item{scoreCol}{
if a \code{sampleSheet} is specified, the default column  in the peak files that will be used 
for scoring if the \code{ScoreCol} column is absent.
}
 \item{bLowerScoreBetter}{
if a \code{sampleSheet} is specified, the  sort order for peak scores 
if the \code{LowerBetter} column is absent.
}
 \item{filter}{
if a \code{sampleSheet} is specified, a filter value if the \code{Filter} column is absent.
Peaks with scores lower than this value 
(or higher if \code{ bLowerScoreBetter} or \code{LowerBetter} is \code{TRUE}) will be removed.
}
  \item{skipLines}{
if a \code{sampleSheet} is specified, the number of lines (ie header lines)
at the beginning of each peak file to skip.
}
  \item{bAddCallerConsensus}{
add a consensus peakset for each sample with more than one peakset 
(i.e. different peak callers) when constructing a new DBA object from a \code{sampleSheet}.
}
  \item{bRemoveM}{
logical indicating whether to remove peaks on chrM (mitochondria) 
when constructing a new DBA object from a sample sheet.
}
\item{bRemoveRandom}{
logical indicating whether to remove peaks on chrN_random when 
constructing a new DBA object from a sample sheet.
}

\item{bSummarizedExperiment}{
logical indicating whether to return resulting object as a \code{\link{SummarizedExperiment}}.
}
\item{bCorPlot}{
logical indicating that a correlation heatmap should be plotted before returning. If \code{DBA} is \code{NULL} (a new DBA object is being created), and \code{bCorPlot} is missing, then this will take the default value (\code{FALSE}). However if \code{DBA} is \code{NULL} (a new DBA object is being created), and \code{bCorPlot} is specified, then the specified value will become the default value of \code{bCorPlot} for the resultant \code{DBA} object.
}
\item{attributes}{
vector of attributes to use subsequently as defaults when generating labels in plotting functions:
\itemize{ 
\item \code{DBA_ID}
\item \code{DBA_TISSUE}
\item \code{DBA_FACTOR}
\item \code{DBA_CONDITION}
\item \code{DBA_TREATMENT}
\item \code{DBA_REPLICATE}
\item \code{DBA_CONSENSUS}
\item \code{DBA_CALLER}
\item \code{DBA_CONTROL}
}
}
\item{dir}{
Directory path. 
If supplied, files referenced in the \code{sampleSheet} will have
this path prepended. 
Applies to \code{PeakFiles}, \code{bamReads}, and \code{bamControl}, if present.
If \code{sampleSheet} is a filepath, this will prepended to that as well.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

MODE: Construct a new DBA object from a samplesheet:

\code{dba(sampleSheet, config,
    bAddCallerConsensus, bRemoveM, bRemoveRandom,
    attributes)}

MODE: Construct a DBA object based on an existing one:

\code{dba(DBA, mask, attributes)}

MODE: Convert a DBA object to a SummarizedExperiment object:

\code{dba(DBA, bSummarizedExperiment=TRUE)}
}
\value{
DBA object
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark and Gordon Brown
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{dba.peakset}},  \code{\link{dba.show}}
}

\examples{
# Create DBA object from a samplesheet
\dontrun{
basedir <- system.file("extra", package="DiffBind")
tamoxifen <- dba(sampleSheet="tamoxifen.csv", dir=basedir)
tamoxifen

tamoxifen <- dba(sampleSheet="tamoxifen_allfields.csv")
tamoxifen

tamoxifen <- dba(sampleSheet="tamoxifen_allfields.csv",config="config.csv")
tamoxifen
}

#Create a DBA object with a subset of samples
data(tamoxifen_peaks)
Responsive <- dba(tamoxifen,tamoxifen$masks$Responsive)
Responsive

# change peak caller but leave peak format the same
basedir <- system.file("extra", package="DiffBind")
tamoxifen <- dba(sampleSheet="tamoxifen.csv", dir=basedir,
                 peakCaller="macs", peakFormat="raw", scoreCol=5 )
dba.show(tamoxifen, attributes=c(DBA_TISSUE,DBA_CONDITION,DBA_REPLICATE,DBA_CALLER))

# Convert DBA object to SummarizedExperiment
data(tamoxifen_counts)
sset <- dba(tamoxifen,bSummarizedExperiment=TRUE)
sset
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
