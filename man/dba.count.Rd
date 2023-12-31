\name{dba.count}
\alias{dba.count}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count reads in binding site intervals
}
\description{
Counts reads in binding site intervals.  Files must be one of bam,
bed and gzip-compressed bed.  File suffixes must be ".bam", ".bed", or ".bed.gz"
respectively.
}
\usage{
dba.count(DBA, peaks, minOverlap=2, score=DBA_SCORE_NORMALIZED,
          fragmentSize=DBA$config$fragmentSize, 
          summits=200, filter=1, bRemoveDuplicates=FALSE, bScaleControl=TRUE,
          bSubControl = is.null(DBA$greylist),
          mapQCth=DBA$config$mapQCth, filterFun=max, minCount=0, 
          bLog=FALSE, bUseSummarizeOverlaps=TRUE, 
          readFormat=DBA_READS_DEFAULT, bParallel=DBA$config$RunParallel) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


 \item{DBA}{
DBA object
}
 \item{peaks}{
If GRanges, RangedData, dataframe, or matrix, this parameter contains the intervals to use for counting. If character string, it specifies a file containing the intervals to use (with the first three columns specifying chromosome, startpos, endpos).If missing or a mask, generates a consensus peakset using minOverlap parameter (after applying the mask if present). If NULL, the \code{score}, \code{filter}, and \code{summits} parameters are honored, updating the global binding matrix without re-counting in the cases of \code{score} and \code{filter}, and only counting after re-centering in the case of \code{summits}.
}
\item{minOverlap}{
only include peaks in at least this many peaksets when generating consensus peakset (i.e. when peaks parameter is missing).  If minOverlap is between zero and one, peak will be included from at least this proportion of peaksets.
}
\item{score}{
which score to use in the binding affinity matrix.
Note that all raw read counts are maintained for use by \code{\link{dba.analyze}}, regardless of how this is set.
One of:
\tabular{ll}{
DBA_SCORE_NORMALIZED \tab normalized reads, as set by \code{\link{dba.normalize}}\cr
DBA_SCORE_READS \tab raw read count for interval using only reads from ChIP\cr
DBA_SCORE_CONTROL_READS \tab raw read count for interval using only reads from Control\cr
DBA_SCORE_READS_FOLD \tab raw read count for interval from ChIP divided by read count for interval from control\cr
DBA_SCORE_READS_MINUS \tab raw read count for interval from ChIP minus read count for interval from control\cr
%DBA_SCORE_READS_FULL \tab read count for interval from ChIP normalized based on Full Library size\cr
%DBA_SCORE_READS_MINUS_FULL \tab read count for interval from ChIP minus read count for interval from control normalized based on Full Library size\cr
%DBA_SCORE_READS_EFFECTIVE \tab read count for interval from ChIP normalized based on Effective Library size\cr
%DBA_SCORE_READS_MINUS_EFFECTIVE \tab read count for interval from ChIP minus read count for interval from control normalized based on Effective Library size\cr
DBA_SCORE_RPKM \tab RPKM for interval using only reads from ChIP\cr
DBA_SCORE_RPKM_FOLD \tab RPKM for interval from ChIP divided by RPKM for interval from control\cr
DBA_SCORE_RPKM_MINUS \tab RPKM for interval from ChIP minus RPKM for interval from control\cr
%DBA_SCORE_TMM_READS_FULL \tab TMM normalized (using edgeR), using ChIP read counts and Full Library size\cr
%DBA_SCORE_TMM_READS_EFFECTIVE \tab TMM normalized (using edgeR), using ChIP read counts and Effective Library size\cr
%DBA_SCORE_TMM_MINUS_FULL \tab TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Full Library size\cr
%DBA_SCORE_TMM_MINUS_EFFECTIVE \tab TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Effective Library size\cr
%DBA_SCORE_TMM_READS_FULL_CPM \tab same as \code{DBA_SCORE_TMM_READS_FULL}, but reported in counts-per-million.\cr
%DBA_SCORE_TMM_READS_EFFECTIVE_CPM \tab same as \code{DBA_SCORE_TMM_READS_EFFECTIVE}, but reported in counts-per-million.\cr
%DBA_SCORE_TMM_MINUS_FULL_CPM \tab same as \code{DBA_SCORE_TMM_MINUS_FULL}, but reported in counts-per-million.\cr
%DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM \tab same as \code{DBA_SCORE_TMM_MINUS_EFFECTIVE}, but reported in counts-per-million.\cr
DBA_SCORE_SUMMIT \tab summit height (maximum read pileup value)\cr
DBA_SCORE_SUMMIT_ADJ \tab summit height (maximum read pileup value), normalized to relative library size\cr
DBA_SCORE_SUMMIT_POS \tab summit position (location of maximum read pileup)\cr
}

If \code{DBA} is a report-based object, the allowable scores are:
\tabular{ll}{
DBA_SCORE_FOLD  \tab log2 Fold Change\cr        
DBA_SCORE_CONCENTRATION \tab mean concentration (log2)\cr      
DBA_SCORE_CONC_NUMERATOR \tab mean concentration (log2) of first group in contrast\cr
DBA_SCORE_CONC_DENOMINATOR \tab mean concentration (log2) of second group in contrast\cr
DBA_SCORE_PVAL  \tab p-value\cr            
DBA_SCORE_FDR  \tab FDR\cr  
}
}

\item{fragmentSize}{
This value will be used as the length of the reads. 
Each read will be extended from its endpoint along the appropriate strand by this many bases. 
If set to zero, the read size indicated in the BAM/BED file will be used.
\code{fragmentSize} may also be a vector of values, one for each ChIP sample plus one for 
each unique Control library.
}
\item{summits}{
unless set to \code{FALSE}, summit heights (read pileup) and locations will be calculated for each peak. 
The values can retrieved using \code{\link{dba.peakset}}. 
The summits can also be used as a read score in the global binding matrix (see \code{score}). 

If the value of \code{summits} is \code{TRUE} (or \code{0}), 
the summits will be calculated but the peaksets will be unaffected.  
If the value is greater than zero, all consensus peaks will be re-centered around a consensus summit, 
with the value of \code{summits} indicating how many base pairs to include upstream and downstream 
of the summit (so all consensus peaks will be of the same width, namely \code{2 * summits + 1}). 

Note that if \code{summits} is greater than zero, the counting procedure will 
take twice as long.
}

\item{filter}{
value to use for filtering intervals with low read counts. 
The \code{filterFun} will be applied to the counts for each interval, 
and if it returns a value below the \code{filter} value, 
the interval will be removed from further analysis.
If \code{peaks} is \code{NULL}, will remove sites from existing DBA object without recounting. 
If filter is a vector of values, \code{dba.count} will return a vector of the same length, 
indicating how many intervals will be retained for each specified \code{filter} level.

NB: the filtering will be based on RPKM values. 
If \code{bSubControl} is \code{FALSE}, this is the RPKM value of the read counts
(equivalent to \code{score=DBA_SCORE_RPKM}.
If \code{bSubControl} is \code{TRUE}, this is the RPKM value of the 
control counts subtracted from the RPKM of the read counts
(equivalent to \code{score=DBA_SCORE_RPKM_MINUS}).

}
\item{bRemoveDuplicates}{
logical indicating if duplicate reads (ones that map to exactly the same genomic position) should be removed.
If \code{TRUE}, any location where multiple reads map will be counted as a single read. 
Note that if \code{bLowMem} is set, duplicates needs to have been already marked 
in all of the BAM files.
The built-in counting code may not correctly handle certain cases when the 
\code{bRemoveDuplicates} parameter is set to \code{TRUE}. 
These cases include paired-end data and datasets where the read length may differ within a single BAM file. 
In these cases, see the \code{bUseSummarizeOverlaps} parameter.
}
\item{bScaleControl}{
logical indicating if the Control reads should be scaled based on relative library sizes. 
If \code{TRUE}, and there are more reads in the Control library than in the ChIP library, 
the number of Control reads for each peak will be multiplied by a scaling factor 
determined by dividing the total number of reads in the ChIP library by 
the total number of reads in the Control library. 
If this value is not an integer, the number of Control reads for each peak  
will be the next highest integer.
}

\item{bSubControl}{
logical indicating whether Control read counts are subtracted for 
each site in each sample.
If there are more overlapping control reads than ChIP reads, 
the count will be set to the \code{minCount} value specified
when \code{\link{dba.count}} was called, or zero if no value is specified.

If \code{bSubControl} is not explicitly specified,
it will be set to \code{TRUE} unless a greylist has been applied
(see \code{\link{dba.blacklist}}).

}

\item{mapQCth}{
for filtering by mapping quality (mapqc). 
Only alignments with mapping scores of at least this value will be included. 
Only applicable for bam files when \code{bUseSummarizeOverlaps=FALSE} 
(setting \code{DBA$config$scanbamparam} appropriately to filter on quality scores when using \code{\link{summarizeOverlaps}}.)
}

\item{filterFun}{
function that will be invoked for each interval with a
vector of scores for each sample. 
Returns a score that will be evaluated against the \code{filter} value 
(only intervals with a score 
at least as high as \code{filter} will be kept).
Default is \code{max}, indicating that at least one sample 
should have a score of at least \code{filter}; 
other useful values include \code{sum} (indicating that all the scores added together should be at least \code{filter}) and \code{mean}
(setting a minimum mean normalized count level). 
Users can supply their own function as well.
}
\item{minCount}{
minimum read count value. 
Any interval with fewer than this many overlapping reads will
be set to have this count.
Also applies to scores.
}

\item{bLog}{
logical indicating whether log2 of score should be used (only applies to 
DBA_SCORE_RPKM_FOLD and DBA_SCORE_READS_FOLD).
}

\item{bUseSummarizeOverlaps}{
logical indicating that \code{\link{summarizeOverlaps}} should be used for
counting instead of the built-in counting code. 
This option is slower but uses the more standard counting function.
If \code{TRUE}, all read files must be BAM (.bam extension), with associated index files 
(.bam.bai extension). 
The \code{fragmentSize} parameter must absent. 

See notes for when the \code{bRemoveDuplicates} parameter is set to \code{TRUE}, 
where the built-in counting code may not correctly handle certain cases and 
\code{bUseSummarizeOverlaps} should be set to \code{TRUE}.

Five additional parameters for \code{\link{summarizeOverlaps}} may be specified in \code{DBA$config}:
\tabular{ll}{
\code{DBA$config$yieldSize} \tab yieldSize indicating how many reads to process at one time; default is 5000000. 
The lower this value, the less memory will be used, but the more time it will take to 
complete the count operation.\cr
\code{DBA$config$intersectMode} \tab mode indicating which overlap algorithm to use; default is \code{"IntersectionNotEmpty"}\cr 
\code{DBA$config$singleEnd} \tab logical indicating if reads are single end; 
if NULL, status will be automatically detected.\cr 
\code{DBA$config$fragments} \tab logical indicating how unmatched reads are counted; default is \code{FALSE}\cr 
\code{DBA$config$inter.feature} \tab logical indicating the setting
for the \code{inter.feature} parameter; default is \code{TRUE}\cr 
\code{DBA$config$scanbamparam} \tab \code{ScanBamParam} object to pass to \code{\link{summarizeOverlaps}}. 
If present, \code{bRemoveDuplicates} is ignored.\cr 
}
}
\item{readFormat}{
Specify the file type of the read files, over-riding the file extension. Possible values:
\tabular{ll}{
\code{DBA_READS_DEFAULT} \tab use file extension (.bam, .bed, .bed.gz) to determine file type\cr
\code{DBA_READS_BAM} \tab assume the file type is BAM, regardless of the file extension\cr
\code{DBA_READS_BED} \tab assume the file type is BED (or zipped BED), regardless of the file extension.
}
Note that if \code{readFormat} is anything other than \code{DBA_READS_DEFAULT}, 
all the read files must be of the same file type.
}
\item{bParallel}{
if \code{TRUE}, use multicore to get counts for each read file in parallel
}
}

\value{
DBA object with binding affinity matrix based on read count scores.
}

%\note{
%The consensus set will be determine
%}

\author{
Rory Stark and Gordon Brown
}
 
 \seealso{
 \code{\link{dba.analyze}}
 }

\examples{
# These won't run unless you have the reads available in a BAM or BED file
data(tamoxifen_peaks)
\dontrun{tamoxifen <- dba.count(tamoxifen)}


# Count using a peakset made up of only peaks in all responsive MCF7 replicates
data(tamoxifen_peaks)
mcf7Common <- dba.overlap(tamoxifen,tamoxifen$masks$MCF7&tamoxifen$masks$Responsive)
\dontrun{tamoxifen <- dba.count(tamoxifen,peaks=mcf7Common$inAll)}
tamoxifen

#First make consensus peaksets from each set of replicates, 
#then derive master consensus set for counting from those
data(tamoxifen_peaks)
tamoxifen <- dba.peakset(tamoxifen,consensus = -DBA_REPLICATE)
\dontrun{tamoxifen <- dba.count(tamoxifen, peaks=tamoxifen$masks$Consensus)}
tamoxifen

# Change binding affinity scores
data(tamoxifen_counts)
dba.peakset(tamoxifen, bRetrieve=TRUE) # default: DBA_SCORE_NORMALIZED
tamoxifen <- dba.count(tamoxifen,peaks=NULL,score=DBA_SCORE_READS)
dba.peakset(tamoxifen, bRetrieve=TRUE)
tamoxifen <- dba.count(tamoxifen,peaks=NULL,score=DBA_SCORE_RPKM_MINUS)
dba.peakset(tamoxifen, bRetrieve=TRUE)

# Plot effect of a range of filter values and then apply filter 
data(tamoxifen_counts)
rate.max <- dba.count(tamoxifen, peaks=NULL, filter=0:250)
rate.sum <- dba.count(tamoxifen, peaks=NULL, filter=0:250,filterFun=sum)
plot(0:250,rate.max/rate.max[1],type='l',xlab="Filter Value",ylab="Proportion Retained Sites")
lines(0:250,rate.sum/rate.sum[1],col=2)
tamoxifen <- dba.count(tamoxifen,peaks=NULL,filter=125,filterFun=sum)
tamoxifen

# Calculate summits
data(tamoxifen_counts) 
# pre-counted with summits=250 or 501bp intervals
as.numeric(dba.show(tamoxifen)$FRiP)
\dontrun{tamoxifen <- dba.count(tamoxifen,peaks=NULL,summits=50)}
# re-counted with summits=50 or 101bp intervals
as.numeric(dba.show(tamoxifen)$FRiP)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
