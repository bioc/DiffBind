\name{DiffBind -- DBA global constant variables}
\alias{DBA_ID}
\alias{DBA_FACTOR}
\alias{DBA_TISSUE}
\alias{DBA_CONDITION}
\alias{DBA_TREATMENT}
\alias{DBA_REPLICATE}
\alias{DBA_CALLER}
\alias{DBA_CONSENSUS}
\alias{DBA_CONTROL}
\alias{DBA_READS}
\alias{DBA_ALL_ATTRIBUTES}
\alias{DBA_INTERVALS}
\alias{DBA_FRIP}
\alias{DBA_GROUP}
\alias{DBA_OLAP_PEAKS}
\alias{DBA_OLAP_ALL}
\alias{DBA_OLAP_RATE}
\alias{DBA_COR}
\alias{DBA_OLAP}
\alias{DBA_INALL}
\alias{DBA_SCORE_READS}
\alias{DBA_SCORE_CONTROL_READS}
\alias{DBA_SCORE_READS_MINUS}
\alias{DBA_SCORE_READS_FULL}
\alias{DBA_SCORE_READS_MINUS_FULL}
\alias{DBA_SCORE_READS_EFFECTIVE}
\alias{DBA_SCORE_READS_MINUS_EFFECTIVE}
\alias{DBA_SCORE_READS_FOLD}
\alias{DBA_SCORE_RPKM}
\alias{DBA_SCORE_NORMALIZED}
\alias{DBA_SCORE_RPKM_FOLD}
\alias{DBA_SCORE_TMM_READS_FULL}
\alias{DBA_SCORE_TMM_READS_EFFECTIVE}
\alias{DBA_SCORE_TMM_MINUS_FULL}
\alias{DBA_SCORE_TMM_MINUS_EFFECTIVE}
\alias{DBA_SCORE_TMM_READS_FULL_CPM}
\alias{DBA_SCORE_TMM_READS_EFFECTIVE_CPM}
\alias{DBA_SCORE_TMM_MINUS_FULL_CPM}
\alias{DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM}
\alias{DBA_SCORE_SUMMIT}
\alias{DBA_SCORE_SUMMIT_ADJ}
\alias{DBA_SCORE_SUMMIT_POS}
\alias{DBA_READS_DEFAULT}
\alias{DBA_READS_BAM}
\alias{DBA_READS_BED}
\alias{DBA_EDGER}
\alias{DBA_DESEQ2}
\alias{DBA_EDGER_BLOCK}
\alias{DBA_DESEQ2_BLOCK}
\alias{DBA_EDGER_GLM}
\alias{DBA_ALL_METHODS}
\alias{DBA_ALL_BLOCK}
\alias{DBA_ALL_METHODS_BLOCK}
\alias{DBA_DATA_FRAME}
\alias{DBA_DATA_GRANGES}
\alias{DBA_DATA_RANGEDDATA}
\alias{DBA_DATA_SUMMARIZED_EXPERIMENT}
\alias{DBA_DATA_DBAOBJECT}
\alias{DBA_BLACKLIST_CE10}
\alias{DBA_BLACKLIST_CE11}
\alias{DBA_BLACKLIST_DM3}
\alias{DBA_BLACKLIST_DM6}
\alias{DBA_BLACKLIST_GRCH37} 
\alias{DBA_BLACKLIST_GRCH38}
\alias{DBA_BLACKLIST_HG19}
\alias{DBA_BLACKLIST_HG38}
\alias{DBA_BLACKLIST_MM9}
\alias{DBA_BLACKLIST_MM10} 
\alias{DBA_BLACKLIST}
\alias{DBA_GREYLIST}
\alias{DBA_BLACKLISTED_PEAKS}
\alias{DBA_LIBSIZE_DEFAULT}
\alias{DBA_LIBSIZE_FULL}
\alias{DBA_LIBSIZE_PEAKREADS}
\alias{DBA_LIBSIZE_BACKGROUND}
\alias{DBA_LIBSIZE_USER}
\alias{DBA_NORM_DEFAULT}
\alias{DBA_NORM_NATIVE}
\alias{DBA_NORM_LIB}
\alias{DBA_NORM_TMM}
\alias{DBA_NORM_RLE}
\alias{DBA_NORM_SPIKEIN}
\alias{DBA_NORM_USER}
\alias{DBA_NORM_OFFSETS}
\alias{DBA_NORM_OFFSETS_ADJUST}
\alias{DBA_OFFSETS_LOESS}
\alias{DBA_OFFSETS_USER}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constant variables used in DiffBind package
}
\description{
Constant variables used in DiffBind package
}
\usage{
DBA_ID
DBA_FACTOR
DBA_TISSUE
DBA_CONDITION
DBA_TREATMENT
DBA_REPLICATE
DBA_CALLER
DBA_CONSENSUS
DBA_CONTROL
DBA_READS
DBA_ALL_ATTRIBUTES

DBA_INTERVALS
DBA_FRIP

DBA_GROUP

DBA_OLAP_PEAKS
DBA_OLAP_ALL
DBA_OLAP_RATE

DBA_COR
DBA_OLAP
DBA_INALL

DBA_SCORE_READS
DBA_SCORE_NORMALIZED
DBA_SCORE_CONTROL_READS
DBA_SCORE_READS_MINUS
DBA_SCORE_READS_FULL
DBA_SCORE_READS_MINUS_FULL
DBA_SCORE_READS_EFFECTIVE
DBA_SCORE_READS_MINUS_EFFECTIVE
DBA_SCORE_READS_FOLD
DBA_SCORE_RPKM
DBA_SCORE_RPKM_FOLD
DBA_SCORE_TMM_READS_FULL
DBA_SCORE_TMM_READS_EFFECTIVE
DBA_SCORE_TMM_MINUS_FULL
DBA_SCORE_TMM_MINUS_EFFECTIVE
DBA_SCORE_TMM_READS_FULL_CPM
DBA_SCORE_TMM_READS_EFFECTIVE_CPM
DBA_SCORE_TMM_MINUS_FULL_CPM
DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM
DBA_SCORE_SUMMIT
DBA_SCORE_SUMMIT_ADJ
DBA_SCORE_SUMMIT_POS

DBA_READS_DEFAULT
DBA_READS_BAM
DBA_READS_BED

DBA_EDGER
DBA_DESEQ2
DBA_EDGER_BLOCK
DBA_DESEQ2_BLOCK
DBA_EDGER_GLM
DBA_ALL_METHODS
DBA_ALL_BLOCK
DBA_ALL_METHODS_BLOCK

DBA_DATA_FRAME
DBA_DATA_GRANGES
DBA_DATA_RANGEDDATA
DBA_DATA_SUMMARIZED_EXPERIMENT
DBA_DATA_DBAOBJECT

DBA_BLACKLIST_CE10
DBA_BLACKLIST_CE11
DBA_BLACKLIST_DM3
DBA_BLACKLIST_DM6
DBA_BLACKLIST_GRCH37
DBA_BLACKLIST_GRCH38
DBA_BLACKLIST_HG19
DBA_BLACKLIST_HG38
DBA_BLACKLIST_MM9
DBA_BLACKLIST_MM10
DBA_BLACKLIST
DBA_GREYLIST
DBA_BLACKLISTED_PEAKS

DBA_LIBSIZE_DEFAULT
DBA_LIBSIZE_FULL
DBA_LIBSIZE_PEAKREADS
DBA_LIBSIZE_BACKGROUND
DBA_LIBSIZE_USER
DBA_NORM_DEFAULT
DBA_NORM_NATIVE
DBA_NORM_LIB
DBA_NORM_TMM
DBA_NORM_RLE
DBA_NORM_SPIKEIN
DBA_NORM_USER
DBA_NORM_OFFSETS
DBA_NORM_OFFSETS_ADJUST
DBA_OFFSETS_LOESS
DBA_OFFSETS_USER

}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{DBA_ID}{
DBA peakset metadata: Peakset ID
}
\item{DBA_FACTOR}{
DBA peakset metadata: Factor
}
\item{DBA_TISSUE}{
DBA peakset metadata: Tissue
}
\item{DBA_CONDITION}{
DBA peakset metadata: Condition
}
\item{DBA_TREATMENT}{
DBA peakset metadata: Treatment
}
\item{DBA_REPLICATE}{
DBA peakset metadata: Replicate
}
\item{DBA_CALLER}{
DBA peakset metadata: Peak Caller
}
\item{DBA_CONSENSUS}{
DBA peakset metadata: Is this a consensus peakset?
}
\item{DBA_CONTROL}{
DBA peakset metadata: ID of Control sample
}
\item{DBA_READS}{
Number of reads counted in BAM file.
}
\item{DBA_ALL_ATTRIBUTES}{
DBA peakset metadata: all attributes that can be used in certain plot labels (cf \code{\link{dba.plotVenn}}), equivalent to \code{c(DBA_ID, DBA_TISSUE, DBA_FACTOR, DBA_CONDITION,DBA_TREATMENT, DBA_REPLICATE, DBA_CALLER)} 
}
\item{DBA_INTERVALS}{
DBA peakset metadata: Number of intervals in peakset
}
\item{DBA_FRIP}{
DBA peakset metadata: Fraction of Reads in Peaks (number of reads in intervals divided by total number of reads in library)
}

\item{DBA_GROUP}{
DBA peakset metadata: color PCA plot using contras groups
}

\item{DBA_OLAP_PEAKS}{
dba.overlap mode: return overlapping/unique peaksets
}
\item{DBA_OLAP_ALL}{
dba.overlap mode: return report of correlations/overlaps for each pair of samples
}
\item{DBA_OLAP_RATE}{
dba.overlap mode: return overlap rates
}
\item{DBA_COR}{
When plotting a heatmap from an overlap record, use the correlation value.
}
\item{DBA_OLAP}{
When plotting a heatmap from an overlap record, use the percentage overlap value.
}
\item{DBA_INALL}{
When plotting a heatmap from an overlap record, use the number of peaks in common to both samples.
}


\item{DBA_SCORE_READS}{
dba.count score is number of reads in ChIP
}
\item{DBA_SCORE_CONTROL_READS}{
dba.count score is number of reads in Control
}
\item{DBA_SCORE_READS_FOLD}{
dba.count score is number of reads in ChIP divided by number of reads in Control
}
\item{DBA_SCORE_READS_MINUS}{
dba.count score is number of reads in ChIP minus  number of reads in Control
}

\item{DBA_SCORE_READS_FULL}{
dba.count score is normalized ChIP read counts, using Full Library size
}
\item{DBA_SCORE_READS_MINUS_FULL}{
dba.count score is normalized ChIP read counts minus Control read counts, using Full Library size
}
\item{DBA_SCORE_READS_EFFECTIVE}{
dba.count score is normalized ChIP read counts, using Effective Library size
}
\item{DBA_SCORE_READS_MINUS_EFFECTIVE}{
dba.count score is normalized ChIP read counts minus Control read counts, using Effective Library size
}

\item{DBA_SCORE_NORMALIZED}{
dba.count score is normalized reads
}

\item{DBA_SCORE_RPKM}{
dba.count score is RPKM of ChIP 
}
\item{DBA_SCORE_RPKM_FOLD}{
dba.count score is RPKM of ChIP divided by RPKM of Control
}
\item{DBA_SCORE_TMM_READS_FULL}{
dba.count score is TMM normalized (using edgeR), using ChIP read counts and Full Library size
}
\item{DBA_SCORE_TMM_READS_EFFECTIVE}{
dba.count score is TMM normalized (using edgeR), using ChIP read counts and Effective Library size
}
\item{DBA_SCORE_TMM_MINUS_FULL}{
dba.count score is TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Full Library size
}
\item{DBA_SCORE_TMM_MINUS_EFFECTIVE}{
dba.count score is TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Effective Library size
}

\item{DBA_SCORE_TMM_READS_FULL_CPM}{
dba.count score is TMM normalized (using edgeR), using ChIP read counts and Full Library size, reported in counts-per-million.
}
\item{DBA_SCORE_TMM_READS_EFFECTIVE_CPM}{
dba.count score is TMM normalized (using edgeR), using ChIP read counts and Effective Library size, reported in counts-per-million.
}
\item{DBA_SCORE_TMM_MINUS_FULL_CPM}{
dba.count score is TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Full Library size, reported in counts-per-million.
}
\item{DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM}{
dba.count score is TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Effective Library size, reported in counts-per-million. 
}
\item{DBA_SCORE_SUMMIT}{
dba.count score is summit height (highest pile-up).
}
\item{DBA_SCORE_SUMMIT_ADJ}{
dba.count score is summit height (highest pile-up), adjusted for library size.
}
\item{DBA_SCORE_SUMMIT_POS}{
dba.count score is summit location (position of highest pile-up).
}
\item{DBA_READS_DEFAULT}{
When counting read files, use the file extension to determine the file type.
}
\item{DBA_READS_BAM}{
When counting read files, assume the file type is BAM, regardless of the file extension.
}
\item{DBA_READS_BED}{
When counting read files, assume the file type is BED (or zipped BED), regardless of the file extension.
}

\item{DBA_EDGER}{
differential analysis method: edgeR (default: DBA_EDGER_GLM)
}
\item{DBA_DESEQ2}{
differential analysis method: DESeq2 (using a single-factor GLM)
}
\item{DBA_EDGER_BLOCK}{
differential analysis method: edgeR with blocking factors (GLM)
}
\item{DBA_DESEQ2_BLOCK}{
differential analysis method: DESeq2 with blocking factors (GLM)
}

\item{DBA_EDGER_GLM}{
differential analysis method: use GLM in edgeR for two-group comparisons
}
\item{DBA_ALL_METHODS}{
use both analysis methods: \code{c(DBA_EDGER, DBA_DESEQ2)}
}
\item{DBA_ALL_BLOCK}{
report on block results for both analysis methods: \code{c(DBA_EDGER_BLOCK, DBA_DESEQ2_BLOCK)}
}
\item{DBA_ALL_METHODS_BLOCK}{
report on block results for all analysis methods, both blocked and unblocked: \code{c(DBA_ALL_METHODS, DBA_ALL_BLOCK)}
}

\item{DBA_DATA_GRANGES}{
Use GRanges class for peaksets and reports. This is the  default (DBA$config$DataType = DBA_DATA_GRANGES).
}
\item{DBA_DATA_RANGEDDATA}{
Use RangedData class for peaksets and reports. Can be set as default (DBA$config$DataType = DBA_DATA_RANGEDDATA).
}
\item{DBA_DATA_FRAME}{
Use data.frame class for peaksets and reports. Can be set as default (DBA$config$DataType = DBA_DATA_FRAME).
}
\item{DBA_DATA_SUMMARIZED_EXPERIMENT}{
Return report as a \code{\link{SummarizedExperiment}}.
}
\item{DBA_DATA_DBAOBJECT}{
Return a result-based DBA object from \code{\link{dba.plotVenn}}.
}

\item{DBA_BLACKLIST_HG19}{Homo sapiens 19 (chromosomes have "chr")}
\item{DBA_BLACKLIST_HG38}{Homo sapiens 38 (chromosomes have "chr")}  
\item{DBA_BLACKLIST_GRCH37}{Homo sapiens 37 (chromosomes are numbers)}
\item{DBA_BLACKLIST_GRCH38}{Homo sapiens 38 (chromosomes are numbers)}
\item{DBA_BLACKLIST_MM9}{Mus musculus 9 }
\item{DBA_BLACKLIST_MM10}{Mus musculus 10}
\item{DBA_BLACKLIST_CE10}{C. elegans 10}
\item{DBA_BLACKLIST_CE11}{C. elegans 11}
\item{DBA_BLACKLIST_DM3}{Drosophila melanogaster 3}
\item{DBA_BLACKLIST_DM6}{Drosophila melanogaster 6}
\item{DBA_BLACKLIST}{Retrieve blacklist}
\item{DBA_GREYLIST}{Retrieve greylist}
\item{DBA_BLACKLISTED_PEAKS}{Retrieve blacklisted peaks}

\item{DBA_LIBSIZE_DEFAULT}{Default library size
(\code{DBA_LIBSIZE_FULL} if no background, and \code{DBA_LIBSIZE_CHR} if
background present)}
\item{DBA_LIBSIZE_FULL}{Full library size (all reads in library)}
\item{DBA_LIBSIZE_PEAKREADS}{Library size is Reads in Peaks}
\item{DBA_LIBSIZE_BACKGROUND}{Library size is Reads in  Background}
\item{DBA_LIBSIZE_USER}{User supplied library sizes}
\item{DBA_NORM_DEFAULT}{Default normalization method}
\item{DBA_NORM_NATIVE}{"Native"" normalization method 
(TMM for \code{DBA_EDGER} and RLE for \code{DBA_DESEQ2})}
\item{DBA_NORM_LIB}{Normalize by library size only}
\item{DBA_NORM_TMM}{Normalize using TMM method (\code{edgeR})}
\item{DBA_NORM_RLE}{Normalize using RLE method (\code{DESeq2})}
\item{DBA_NORM_SPIKEIN}{Normalize based on spike-ins}
\item{DBA_NORM_USER}{User supplied normalization factors}
\item{DBA_NORM_OFFSETS}{Use offsets instead of normalization factors}
\item{DBA_NORM_OFFSETS_ADJUST}{Use offsets instead of normalization factors; 
adjust based on library size (DESeq) }
\item{DBA_OFFSETS_LOESS}{Compute offsets using loess fit}
\item{DBA_OFFSETS_USER}{Use offsetrs supplied by user}
}
%\details{
%%
%}
%\value{
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
\note{
Variables with ALL CAP names are used as constants within \code{DiffBind}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
