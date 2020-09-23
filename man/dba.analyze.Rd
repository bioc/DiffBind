\name{dba.analyze}
\alias{dba.analyze}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform differential binding affinity analysis
}
\description{
Performs differential binding affinity analysis. 
Performs default generation of a consensus peakset, read counting,
normalization, and setting up of contrasts if they have
not been specified.
}
\usage{
dba.analyze(DBA, method=DBA$config$AnalysisMethod, design,
            bBlacklist=DBA$config$doBlacklist,
            bGreylist=DBA$config$doGreylist,
            bRetrieveAnalysis=FALSE, bReduceObjects=TRUE, 
            bParallel=DBA$config$RunParallel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}

\item{DBA}{
Either a \code{DBA} object, or a sample sheet 
(either a \code{character} vector with the name of the
sample sheet, or a \code{data.frame} containing the experimental
metadata.

If no blacklist or greylists are included,
a call will be made to
\code{\link{dba.blacklist}} using defaults.
This can be skipped by setting the \code{bBlacklist} and/or
\code{bGreylist} parameters.

If no counts are included, a default consensus will be formed
and read counts computed via a call to
\code{\link{dba.count}} using defaults.

If no normalization has been specified,
the reads will be normalized via a call to
\code{\link{dba.normalize}} using defaults.

If no contrasts are specified (\code{DBA$contrast} is \code{NULL}),
default contrasts will be added via a call to \code{\link{dba.contrast}}
using defaults.
}

\item{method}{
Underlying method, or vector of methods, by which to analyze 
differential binding affinity. 

Supported methods:
\itemize{ 
  \item{\code{\link{DBA_EDGER}}} {use \code{edgeR} package for analysis}
  \item{\code{\link{DBA_DESEQ2}}} {use \code{DESeq2} package for analysis}
  \item{\code{\link{DBA_ALL_METHODS}}} {perform two analyses, 
using both \code{edgeR} and \code{DESeq2}}
}
}

\item{design}{
If present and a character string, will be used as the design formula
for the analysis, replacing any previously established design if present.

If \code{FALSE}, will complete analysis in pre-version 3 mode.

See \code{link{dba.contrast}}.
}

\item{bBlacklist}{
If \code{TRUE}, and no blacklist has been applied to the \code{DBA} object
(or when starting from a samplesheet), the read bam files will
be examined to determine the reference genome, and an appropriate blacklist
applied, if available.
See \code{link{dba.blacklist}}.
}

\item{bGreylist}{
If \code{TRUE}, and no greylist has been applied to the \code{DBA} object
(or when starting from a samplesheet), the control bam files, if present, will
be examined to determine the reference genome, greylists will be
computed for each, merged into a master greylist, and applied to the
peaksets.
See \code{link{dba.blacklist}}.
}

\item{bRetrieveAnalysis}{
If changed from \code{FALSE}, the underlying DE analysis object is returned
instead of running a new analysis.
Possible values for \code{bRetrieveAnalysis}:
\itemize{
  \item{\code{DBA_DESEQ2}} {Returns the \code{DESeq2} \code{\link{DESeqDataSet}}}.
  \item{\code{DBA_EDGER}} {Returns the \code{edgeR} \code{DGEList}}.  
  \item{\code{TRUE}} {Returns the \code{DESeq2} \code{\link{DESeqDataSet}}, 
  if present. If not, returns the \code{edgeR} \code{DGEList}, if present.}.
}
An analysis object will only be successfully returned if there is at
least one contrast utilizing an explicit design 
(see \code{\link{dba.contrast}}), and an analysis has been carried out.
}

\item{bReduceObjects}{
logical indicating whether strip the analysis objects of unnecessary
fields to save memory. 
If it is desired to use the \code{DBA$contrasts[[n]]$edgeR}
and/or \code{DBA$contrasts[[n]]$DESeq2} objects directly in the
\code{edgeR} and/or \code{DESeq2} packages, this should be set to \code{FALSE}.
}
\item{bParallel}{
logical indicating that the analyses is to be done in parallel using multicore 
(one process for each contrast for each method,
plus an additional process per method).
}
}
\details{
In general, prior to calling \code{dba.analyze}, \code{\link{dba.count}}
should have been run.
If no contrasts have been established prior to invoking \code{dba.analyze},
then the default set of contrasts will be added using (\code{\link{dba.contrast}}).

If no normalization parameters have been supplied by 
calling \code{\link{dba.normalize}},
default normalization parameters will be used.

See the DBA User Guide for more details on how 
the \code{edgeR} and \code{DESeq2} analyses are carried out.  
}
\value{
DBA object with results of analysis added to \code{DBA$contrasts}.

Alternatively, an analysis object (either a \code{\link{DESeqDataSet}}
or a \code{DGEList}) if \code{bRetrieveAnalysis} if not \code{FALSE}.
}

\author{
Rory Stark
}
\note{
If there is a blocking factor for the contrast(s) specified using a 
previous call to \code{\link{dba.contrast}} with \code{design=FALSE}, 
a multi-factor analysis will automatically be carried out in addition
to a single factor analysis.
}

\seealso{

\code{\link{dba.blacklist}}, \code{\link{dba.count}},
\code{\link{dba.contrast}}, \code{\link{dba.normalize}},
\code{\link{dba.report}}

}

\examples{
data(tamoxifen_counts)
dba.analyze(tamoxifen)

tamoxifen <- dba.analyze(tamoxifen, method=DBA_ALL_METHODS,
                         design="~Tissue + Condition")
dba.show(tamoxifen, bContrasts=TRUE)

dba.analyze(tamoxifen, bRetrieveAnalysis=TRUE)
edger.object <- dba.analyze(tamoxifen, bRetrieveAnalysis=DBA_EDGER)
class(edger.object)
}

