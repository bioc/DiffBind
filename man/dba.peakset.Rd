\name{dba.peakset}
\alias{dba.peakset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add a peakset to, or retrieve a peakset from, a DBA object
}
\description{
Adds a peakset to, or retrieves a peakset from, a DBA object
}
\usage{
dba.peakset(DBA=NULL, peaks, sampID, tissue, factor, condition, treatment, replicate,
            control, peak.caller, peak.format, reads=0, consensus=FALSE, 
            bamReads, bamControl, spikein,
            scoreCol, bLowerScoreBetter, filter, counts,
            bRemoveM=TRUE, bRemoveRandom=TRUE,
            minOverlap=2, bMerge=TRUE,
            bRetrieve=FALSE, writeFile, numCols=4,
            DataType=DBA$config$DataType)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


\item{DBA}{
DBA object. Required unless creating a new DBA object by adding an initial peakset.
}

\item{peaks}{
When adding a specified peakset: set of peaks, either a \code{\link{GRanges}} object,
or a peak dataframe or matrix (chr,start,end,score), or a filename where the peaks are stored.

When adding a consensus peakset: a sample mask or vector of peakset numbers to include in the consensus. 
If missing or \code{NULL}, a consensus is derived from all peaksets present in the model. 
See  \code{\link{dba.mask}}, or \code{\link{dba.show}} to get peakset numbers.

When adding and empty peakset (zero peaks), set \code{peaks=NA}.

When adding a set of consensus peaksets: a sample mask or vector of peakset numbers. 
Sample sets will  be derived only from subsets of these peaksets.

When adding all the peaks from one DBA object to another: a DBA object. 
In this case, the only other parameter to have an effect is \code{minOverlap}.

When retrieving and/or writing a peakset: either a \code{\link{GRanges}}, or a peak dataframe or 
matrix (chr,start,end,score), or a peakset number; 
if \code{NULL}, retrieves/writes the full binding matrix.
}

%\subsection{Optional/defaulted arguments}{}

\item{sampID}{
ID string for the peakset being added; if missing, one is assigned 
(a serial number for a new peakset, or a concatenation of IDs for a consensus peakset).
Must be unique for each sample.
}
\item{tissue}{
tissue name for the peakset being added; if missing, one is assigned 
for a consensus peakset (a concatenation of tissues).
}
  \item{factor}{
factor name for the peakset being added; if missing, one is assigned 
for a consensus peakset (a concatenation of factors).
}
  \item{condition}{
condition name for the peakset being added; 
if missing, one is assigned for a consensus peakset (a concatenation of conditions).
}
  \item{treatment}{
treatment name for the peakset being added;
if missing, one is assigned for a consensus peakset (a concatenation of treatment).
}
  \item{replicate}{
replicate number for the peakset being added; 
if missing, one is assigned for a consensus peakset (a concatenation of replicate numbers).
}
  \item{control}{
control name for the peakset being added; 
if missing, one is assigned for a consensus peakset (a concatenation of control names).
}
\item{peak.caller}{
peak caller name string. If peaks is specified as a file, and peak.format is missing, a default fie format for the caller will be used (see peak.format).  Supported values:
\itemize{
\item {\dQuote{raw}:}	       {default peak.format: raw text file}
\item {\dQuote{bed}:}	       {default peak.format: bed file}
\item {\dQuote{narrow}:}  {default peak.format: narrowPeaks file}
\item {\dQuote{macs}:}     {default peak.format: MACS .xls file}
\item {\dQuote{bayes}:}    {default peak.format: bayesPeak file}
\item {\dQuote{tpic}:}	        {default peak.format: TPIC file}
\item {\dQuote{sicer}:}      {default peak.format: SICER file}
\item {\dQuote{fp4}:}	        {default peak.format: FindPeaks v4 file}
\item {\dQuote{swembl}:} {default peak.format: SWEMBL file}
\item {\dQuote{csv}:} {default peak.format: comma separated value file}
\item {\dQuote{report}:} {default peak.format: csv file saved via \code{\link{dba.report}}}
}
When adding a consensus peakset, a default value (a concatenation of peak caller names) is assigned if this is missing. 
}
\item{peak.format}{
peak format string. If specified, overrides the default file format for the specified peak caller. Supported formats (with default score column):
\itemize{
\item {\dQuote{raw}:}	  {raw text file file; scoreCol=4}
\item {\dQuote{bed}:}	  {bed file; scoreCol=5}
\item {\dQuote{narrow}:}  {narrowPeaks file; scoreCol=8}
\item {\dQuote{macs}:}  {MACS .xls file; scoreCol=7}
\item {\dQuote{bayes}:}  {bayesPeak file; scoreCol=4, filter=0.5}
\item {\dQuote{tpic}:}	  {TPIC file; scoreCol=0 (all scores=1)}
\item {\dQuote{sicer}:}  {SICER file; scoreCol=7}
\item {\dQuote{fp4}:}	  {FindPeaks v4 file; scoreCol=5}
\item {\dQuote{swembl}:}	{SWEMBL file; scoreCol=4}
\item {\dQuote{csv}:}	{csv file; scoreCol=4}
\item {\dQuote{report}:}	{report file; scoreCol=9, bLowerScoreBetter=T}
}
}
\item{reads}{
total number of ChIPed library reads for the peakset being added.
}
 \item{consensus}{
either the logical value of the consensus attribute when adding a specific peakset 
(set to \code{TRUE} for consensus peaksets generated by  \code{\link{dba.peakset}}), 
or a metadata attribute or vector of attributes when generating a set of consensus peaksets. 
In the latter case, a consensus peakset will be added for each set of samples that have 
the same values for the specified attributes. 
Alternatively, attributes may be specified proceeded by a negative sign, 
in which case a consensus peakset will be added for each set of samples 
that differ only in their values for those attributes. 
See examples.  Allowable attributes:
\itemize{
\item \code{DBA_TISSUE; -DBA_TISSUE}
\item \code{DBA_FACTOR; -DBA_FACTOR}
\item \code{DBA_CONDITION; -DBA_CONDITION}
\item \code{DBA_TREATMENT; -DBA_TREATMENT}
\item \code{DBA_REPLICATE; -DBA_REPLICATE}
\item \code{DBA_CALLER; -DBA_CALLER}
}
}
\item{bamReads}{
file path of the BAM/BED file containing the aligned reads for the 
peakset being added.
}
\item{bamControl}{
file path of the BAM/BED file containing the aligned reads for the control 
used for the peakset being added.
}
\item{spikein}{
file path of the BAM/BED file containing the aligned reads for the spike-ins
used for the peakset being added.
}
  \item{scoreCol}{
peak column to normalize to 0...1 scale when adding a peakset; 
0 indicates no normalization
}
  \item{bLowerScoreBetter}{
Logical indicating that lower scores indicate higher confidence peaks; 
default is that higher scores indicate better peaks.
}
  \item{filter}{
Numeric indicating a filter value for peaks. 
If present, any peaks with a score less than this value 
(or higher if \code{bLowerScoreBetter==TRUE}) will be removed from the peakset. 
}
  \item{counts}{
Used for adding externally computed peak counts. 
Can be a filename or a dataframe. 
Can consist of a single column (or vector) with the counts, or two columns, 
with an ID for each interval in the first column and the counts in the second column,
or four columns (chr, start, end, counts). 
When \code{counts} is specified, \code{peaks} and related parameters are ignored, 
and all peaksets in the DBA object must be specified in this way, 
all with exactly the same number of intervals.
}
  \item{bRemoveM}{
logical indicating whether to remove peaks on chrM when adding a peakset
}
  \item{bRemoveRandom}{
logical indicating whether to remove peaks on chrN_random when adding a peakset
}
  \item{minOverlap}{
the minimum number of peaksets a peak must be in to be included when adding a consensus peakset.
When retrieving, if the peaks parameter is 
a vector (logical mask or vector of peakset numbers), a binding matrix will be retrieved 
including all peaks in at least this many peaksets. 
If \code{minOverlap} is between zero and one, 
peak will be included from at least this proportion of peaksets.
}
  \item{bMerge}{
logical indicating whether global binding matrix should be compiled after adding the peakset. 
When adding several peaksets via successive calls to  \code{\link{dba.peakset}},
it may be more efficient to set this parameter to \code{FALSE} and call \code{dba(DBA)} 
after all of the peaksets have been added.
}
  \item{bRetrieve}{
logical indicating that a peakset is being retrieved and/or written, not added.
}
 \item{writeFile}{
file to write retrieved peakset.
}
 \item{numCols}{
number of columns to include when writing out peakset.
First four columns are chr, start, end, score; 
the remainder are maintained from the original peakset. 
Ignored when writing out complete binding matrix.
}
\item{DataType}{
The class of object for returned peaksets:
\itemize{
\item \code{DBA_DATA_GRANGES}
\item \code{DBA_DATA_FRAME}
}
Can be set as default behavior by setting \code{DBA$config$DataType}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

MODE: Add a specified peakset:

\code{dba.peakset(DBA=NULL, peaks, sampID, tissue, factor, condition, replicate,
                       control, peak.caller, reads, consensus, bamReads, bamControl,
                       normCol, bRemoveM, bRemoveRandom)}

MODE: Add a consensus peakset (derived from overlapping peaks in peaksets already present):

\code{dba.peakset(DBA, peaks, minOverlap)}

MODE: Add a sets of consensus peaksets bases on sample sets that share or differ in specified attributes

\code{dba.peakset(DBA, peaks, consensus, minOverlap)}

MODE: Retrieve a peakset:

\code{dba.peakset(DBA, peaks, bRetrieve=T)}

MODE: Write a peakset out to a file:

\code{dba.peakset(DBA, peaks, bRetrieve=T,  writeFile, numCols)}
}
\value{
DBA object when adding a peakset.
Peakset matrix or \code{\link{GRanges}} object when retrieving and/or writing a peakset.
}
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

\seealso{
to add peaksets using a sample sheet, see  \code{\link{dba}}.

\code{$config$} options are described in \code{\link{DBA.config}}.
}

\examples{

# create a new DBA object by adding three peaksets
mcf7 <- dba.peakset(NULL,
                   peaks=system.file("extra/peaks/MCF7_ER_1.bed.gz", package="DiffBind"),
                   peak.caller="bed", sampID="MCF7.1",tissue="MCF7",
                   factor="ER",condition="Responsive",replicate=1)
mcf7 <- dba.peakset(mcf7,
                   peaks=system.file("extra/peaks/MCF7_ER_2.bed.gz", package="DiffBind"),    
                   peak.caller="bed", sampID="MCF7.2",tissue="MCF7",
                   factor="ER",condition="Responsive",replicate=2)
mcf7 <- dba.peakset(mcf7,
                   peaks=system.file("extra/peaks/MCF7_ER_3.bed.gz", package="DiffBind"),      
                   peak.caller="bed", sampID="MCF7.3",tissue="MCF7",
                   factor="ER",condition="Responsive",replicate=3)
mcf7

#retrieve peaks that are in all three peaksets
mcf7.consensus <- dba.peakset(mcf7, 1:3, minOverlap=3, bRetrieve=TRUE)
mcf7.consensus

#add a consensus peakset -- peaks in all three replicates
mcf7 <- dba.peakset(mcf7, 1:3, minOverlap=3,sampID="MCF7_3of3")
mcf7

#add consensus peaksets for all sample types by combining replicates
data(tamoxifen_peaks)
tamoxifen <- dba.peakset(tamoxifen,consensus = -DBA_REPLICATE)
dba.show(tamoxifen,mask=tamoxifen$masks$Consensus)

#add consensus peaksets for all sample types by (same tissue and condition) 
data(tamoxifen_peaks)
tamoxifen <- dba.peakset(tamoxifen,consensus = c(DBA_TISSUE,DBA_CONDITION))
dba.show(tamoxifen,mask=tamoxifen$masks$Consensus)
dba.plotVenn(tamoxifen,tamoxifen$masks$Responsive & tamoxifen$masks$Consensus)

#create consensus peaksets from sample type consensuses for Responsive and Resistant sample groups
tamoxifen <- dba.peakset(tamoxifen,peaks=tamoxifen$masks$Consensus,consensus=DBA_CONDITION)
dba.show(tamoxifen,mask=tamoxifen$masks$Consensus)
dba.plotVenn(tamoxifen,17:18)
 
#retrieve the consensus peakset as GRanges object
mcf7.consensus <- dba.peakset(mcf7,mcf7$masks$Consensus,bRetrieve=TRUE)
mcf7.consensus

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
