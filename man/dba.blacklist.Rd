\name{dba.blacklist}
\alias{dba.blacklist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply blacklists and/or greylists to peaks (and generate greylist)
}
\description{
Filters peak intervals that overlap a blacklist (from ENCODE or user supplied.) 
Filter peak intervals that overlap a greylist, either user supplied or generated based on 
aligned reads for control samples (e.g. Inputs).
}
\usage{
dba.blacklist(DBA, blacklist=TRUE, greylist=TRUE, 
              Retrieve, cores=DBA$config$cores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DBA}{
DBA object
}
  \item{blacklist}{
If not equal to \code{FALSE}, specifies that a blacklist should be 
applied to the peak intervals in the \code{DBA} object.

If equal to \code{TRUE}, the read bam files will be examined
to determine an appropriate reference genome.
If successful, and a blacklist is available for that genome, it
will be applied.

A user specified blacklist can be specified by setting this parameter to
a \code{\link{GRanges}} object containing the blacklisted regions.

Otherwise, this parameter may be set to one of the following constants, 
indicating which of the ENCODE blacklists should be applied:
 \itemize{
\item{\code{DBA_BLACKLIST_HG19}:} {Homo sapiens 19 (chromosomes have "chr")}
\item{\code{DBA_BLACKLIST_HG38}:} {Homo sapiens 38 (chromosomes have "chr")}  
\item{\code{DBA_BLACKLIST_GRCH37}:} {Homo sapiens 37 (chromosomes are numbers)}
\item{\code{DBA_BLACKLIST_GRCH38}:} {Homo sapiens 38 (chromosomes are numbers)}
\item{\code{DBA_BLACKLIST_MM9}:} {Mus musculus 9}
\item{\code{DBA_BLACKLIST_MM10}:} {Mus musculus 10}
\item{\code{DBA_BLACKLIST_CE10}:} {C. elegans 10}
\item{\code{DBA_BLACKLIST_CE11}:} {C. elegans 11}
\item{\code{DBA_BLACKLIST_DM3}:} {Drosophila melanogaster 3}
\item{\code{DBA_BLACKLIST_DM6}:} {Drosophila melanogaster 6}
 }
}
  \item{greylist}{
If not equal to \code{FALSE}, specifies that a greylist should be 
applied to the peak intervals in the \code{DBA} object.

If equal to \code{TRUE}, the control bam files (if present),
will be examined to determine an appropriate reference genome.
Genomes associated with a valid \code{BSgenome}
can be detected.
If successful, this genome will be used to generate greylists
for each available control 
(eg specified as \code{bamControls} in the sample sheet.)).

The \code{greylist} parameter can also be set explicitly to either
a valid \code{BSgenome} object, or a character string
with the name of a valid \code{BSgenome} object.
the following constants map to a subset of possible \code{BSgenome}
objects:

\itemize{
\item{\code{DBA_BLACKLIST_HG19}  :} {\code{ seqinfo from BSgenome.Hsapiens.UCSC.hg19}}
\item{\code{DBA_BLACKLIST_HG38}  :} {\code{ seqinfo from BSgenome.Hsapiens.UCSC.hg38}}  
\item{\code{DBA_BLACKLIST_GRCH38}:} {\code{ seqinfo from BSgenome.Hsapiens.NCBI.GRCh38}}
\item{\code{DBA_BLACKLIST_MM9}    :} {\code{ seqinfo from BSgenome.Mmusculus.UCSC.mm9}}
\item{\code{DBA_BLACKLIST_MM10}  :} {\code{ seqinfo from BSgenome.Mmusculus.UCSC.mm10}}
\item{\code{DBA_BLACKLIST_CE10}  :} {\code{ seqinfo from BSgenome.Celegans.UCSC.ce10}}
\item{\code{DBA_BLACKLIST_CE11}  :} {\code{ seqinfo from BSgenome.Celegans.UCSC.ce11}}
\item{\code{DBA_BLACKLIST_DM3}    :} {\code{ seqinfo from BSgenome.Dmelanogaster.UCSC.dm3}}
\item{\code{DBA_BLACKLIST_DM6}    :} {\code{ seqinfo from BSgenome.Dmelanogaster.UCSC.dm6}}
}     

A user specified greylist can also be specified by setting this parameter to
a \code{\link{GRanges}} object containing the greylisted regions.
It can also be a \code{list} with an element named \code{greylist$master},
which is a \code{\link{GRanges}} object containing the greylist to be applied.
}

\item{Retrieve}{

If present, some aspects of a previous run of the function is retrieved instead
of returning a \code{DBA} object.

If \code{Retrieve=DBA_BLACKLIST}, the blacklist, if present, is returned as a
\code{\link{GRanges}} object.

If \code{Retrieve=DBA_GREYLIST}, the greylist, if present, is returned.
If it was generated from more than one control, it will be returned as
a \code{list} object with the first element (named \code{$master}) a 
\code{\link{GRanges}} object containing the merged greylist, and the second
element (named \code{$controls}) being a \code{\link{GRangesList}}
with each element containing the greylist for one control

If \code{Retrieve=DBA_BLACKLISTED_PEAKS}, 
the excluded peaks for each sample will be returned in
a \code{\link{GRangesList}} object (with each element containing the
filtered peak intervals for each sample).
If counts are available for the peaks, this will include the following metadata columns:
\itemize{
\item{cReads:} {Number of control reads overlapping this interval}
\item{Reads:}  {Number of primary (ChIP) reads overlapping this interval}
\item{Score:}  {Read score calculated by \code{\link{dba.count}}}
}

Note that the if \code{Retrieve} is set,
\code{dba.blacklist} must have been previously run,
and all other parameters will be ignored.
}
  \item{cores}{
Parallel cores to use when running greylist generation.
}
}
\details{
This function is intended to filter peak intervals that fall in regions of the genome that
are known to be problematic for ChIP analysis.
Blacklists, which are derived for a reference genome and should be applied
for any experiments that use that reference,
are distinguished from greylists, which are derived on a per-experiment basis
using anomalous pileups in the control tracks (such as Inputs).

A core set of blacklists have been defined as part of the ENCODE project (see references).

Greylists can be generated using this function, which serves as a front-end to the \code{GreyListChIP} package. 
See the details of that package for more information on how it works.
Note that the \code{GreyListChIP} package can be utilized separately to generate greylists
with more fine-grained control, with the results passed back to \code{DiffBind} to filter peaks.

}
\value{
\code{DBA} object, with peaks filtered (unless \code{Retrieve} is specified.)
}
\references{
\itemize{
\item{Amemiya HM, Kundaje A, Boyle AP. The ENCODE blacklist: identification of problematic regions of the genome. Sci Rep. 2019 Dec; 9(1) 9354 DOI: 10.1038/s41598-019-45839-z
}
\item{
ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74. doi: 10.1038/nature11247.
}
\item{Brown, Gord. Generating Grey Lists from Input Libraries. Bioconductor. https://bioconductor.org/packages/release/bioc/html/GreyListChIP.html}
}
}
\author{
Rory Stark with thanks to Gord Brown
}
\note{

The \code{p} threshold can be altered by setting \code{DBA$config$greylist.pval}. 
The default is \code{0.999}. 
See \code{GreyListChIP::}\code{\link{calcThreshold}} for details.

Ideally, Blacklists and Greylists will be applied to the aligned reads prior
to calling peaks, as removing reads in anomalous regions
will yield better background noise models. 
Once greylists have been generated, peaks can be
re-called and read into \code{DiffBind}.
}

\seealso{
\code{GreyListChIP (\link{GreyList})}, \code{BSgenome}
}
\examples{

data(tamoxifen_peaks)
\dontrun{tamoxifen <- dba.blacklist(tamoxifen, blacklist=DBA_BLACKLIST_HG19,
                                    greylist="BSgenome.Hsapiens.UCSC.hg19")}
data(tamoxifen_greylist)
tamoxifen <- dba.blacklist(tamoxifen, blacklist=DBA_BLACKLIST_HG19,
                           greylist=tamoxifen.greylist$master)
dba.blacklist(tamoxifen,Retrieve=DBA_GREYLIST)

data(tamoxifen_counts)
tamoxifen <- dba.count(tamoxifen, peaks=NULL, score=DBA_SCORE_CONTROL_READS)
tamoxifen <- dba.blacklist(tamoxifen, blacklist=DBA_BLACKLIST_HG19,
                           greylist=tamoxifen.greylist$master)
blacklisted <- dba.blacklist(tamoxifen, Retrieve=DBA_BLACKLISTED_PEAKS)
mean(blacklisted[[1]]$cReads)
mean(dba.peakset(tamoxifen,peaks=1,bRetrieve=TRUE)$Score)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
