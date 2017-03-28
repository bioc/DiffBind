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
                            fold=0, factor="", bFlip=FALSE, dotSize=1)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


\item{DBA}{
DBA object, on which  \code{\link{dba.analyze}} should have been successfully run.
}
%\subsection{Optional/defaulted arguments}{}
 \item{contrast}{
number of contrast to report on. See  \code{\link{dba.show}(DBA, bContrast=TRUE)} to get contrast numbers.
}
\item{method}{
 method or vector of methods to plot results for:
\itemize{ 
\item \link{DBA_DESEQ2}
\item \link{DBA_DESEQ2_BLOCK}
\item \link{DBA_EDGER}
\item \link{DBA_EDGER_BLOCK}
}
}


\item{th}{
significance threshold; sites with FDR (or p-values, see \code{bUsePval}) less than or equal to this value will be colored red in the plot
}
  \item{bUsePval}{
logical indicating whether to use FDR (FALSE) or p-value (TRUE) for thresholding.
}
  \item{fold}{
will only include sites with fold change greater than this as significant (colored red).
}

 \item{factor}{
string to be prepended to plot main title; e.g. factor name.
}

 \item{bFlip}{
logical indicating that order of groups in contrast should be "flipped", allowing control of which sample group
will have positive and which will have negative fold changes.
}


  \item{dotSize}{
size of points on plot.
}
}

%\details{
%}
%\value{
%}
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

# only highlight significant sites with at least 10x Fold Change
dba.plotVolcano(tamoxifen, fold=log2(10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
