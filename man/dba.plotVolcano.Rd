\name{dba.plotVolcano}
\alias{dba.plotVolcano}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate volcano plots of differential binding analysis results
}
\description{
Generates volcano plots of differential binding analysis results.
}
\usage{
dba.plotVolcano(DBA, contrast=1, method=DBA$config$AnalysisMethod, 
                            th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                            fold=0, factor="", bFlip=FALSE, 
                            bLabels=FALSE, maxLabels=50, dotSize=1,
                            bReturnSites=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


\item{DBA}{
DBA object, on which  \code{\link{dba.analyze}} should have been successfully run.
}
%\subsection{Optional/defaulted arguments}{}
 \item{contrast}{
number of contrast to report on. 
See  \code{\link{dba.show}(DBA, bContrast=TRUE)} to get contrast numbers.
}
\item{method}{
 method or vector of methods to plot results for:
\itemize{ 
\item \code{\link{DBA_DESEQ2}}
\item \code{\link{DBA_DESEQ2_BLOCK}}
\item \code{\link{DBA_EDGER}}
\item \code{\link{DBA_EDGER_BLOCK}}
}
}


\item{th}{
significance threshold; sites with FDR (or p-values, see \code{bUsePval})
less than or equal to this value will be colored red in the plot
}
  \item{bUsePval}{
logical indicating whether to use FDR (\code{FALSE}) or p-value (\code{TRUE}) for thresholding.
}
  \item{fold}{
will only include sites with fold change greater than this 
as significant (colored red).

If \code{fold} is greater than zero, 
and an explicit design was used for the contrast,
the \code{p-value} and \code{FDR} will be re-calculated based on testing
for changes greater than the specified fold change.
For a \code{DESeq2} analysis, this involves including the \code{fold}
when calling \code{DESeq2::results}.
For a \code{edgeR} analysis, \code{edgeR::glmTreat} is used.
}

 \item{factor}{
string to be prepended to plot main title; e.g. factor name.
}

 \item{bFlip}{
logical indicating that order of groups in contrast should be "flipped", 
allowing control of which sample group
will have positive and which will have negative fold changes.
}

 \item{bLabels}{
logical indicating that labels should be drawn on the plot.
The labels are the site numbers, the row index in the (silently) returned set of
significant sites. The maximum number of sites can be specified 
using \code{maxLabels}.
}

 \item{maxLabels}{
The maximum number of labels to use in the plot.
Ignored if \code{bLabels=FALSE}.
}

  \item{dotSize}{
size of points on plot.
}
  \item{bReturnSites}{
If \code{TRUE}, silently returns the differential sites. 
If \code{FALSE}, the \code{ggplot} object is silently returned.
}
}

\details{
Makes a volcano plot.
}
\value{
Silently returns wither a \code{GRanges} object of the sites highlighted in red or a \code{ggplot} object.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
%\note{
%plots log2 values of read counts.
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{dba.analyze}},  \code{\link{dba.plotMA}}
}

\examples{
data(tamoxifen_analysis)

# default volcano plot
dba.plotVolcano(tamoxifen)

# only highlight significant sites with at least 3x Fold Change
sigSites <- dba.plotVolcano(tamoxifen, fold=log2(3))

# use labels to find outlier sites
sigSites <- dba.plotVolcano(tamoxifen, fold=log2(5), th=0.01, bLabels=TRUE)
sigSites
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
