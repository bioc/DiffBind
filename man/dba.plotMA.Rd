\name{dba.plotMA}
\alias{dba.plotMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate MA and scatter plots of differential binding analysis results
}
\description{
Generates MA and scatter plots of differential binding analysis results.
}
\usage{
dba.plotMA(DBA, contrast=1, method=DBA$config$AnalysisMethod, 
           th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
           fold=0, bNormalized=TRUE,
           factor="", bFlip=FALSE, bXY=FALSE, dotSize=.45, 
           bSignificant=TRUE, bSmooth=TRUE, 
           xrange, yrange, ...)
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
\item \code{\link{DBA_DESEQ2}}
\item \code{\link{DBA_DESEQ2_BLOCK}}
\item \code{\link{DBA_EDGER}}
\item \code{\link{DBA_EDGER_BLOCK}}
}
}


\item{th}{
significance threshold; all sites with FDR (or p-values, see \code{bUsePval}) 
less than or equal to this value will be colored red in the plot
}
  \item{bUsePval}{
logical indicating whether to use FDR (\code{FALSE}) or p-value (\code{TRUE}) for thresholding.
}
  \item{fold}{
will only include sites with fold change greater than this as significant (colored red).
}
  \item{bNormalized}{
logical indicating whether to plot normalized data using normalization factors computed by 
differential analysis method (\code{TRUE}) or raw read counts (\code{FALSE}).
}
 \item{factor}{
string to be prepended to plot main title; e.g. factor name.
}
 \item{bFlip}{
logical indicating that order of groups in contrast should be "flipped", 
allowing control of which sample group
will have positive and which will have negative fold changes.
}
\item{bXY}{
logical indicating whether to draw MA plot (\code{FALSE}) or XY scatter plot (\code{TRUE}).
}

  \item{dotSize}{
size of points on plot (\code{cex}).
}
  \item{bSignificant}{
Logical indicating if points corresponding to significantly differentially bound sites 
(based on \code{contrast}, \code{th}, \code{bUsePval}, and \code{fold} parameters) 
should be overlaid in red.
}
  \item{bSmooth}{
logical indicating that basic plot should be plotted using \code{\link{smoothScatter}}.
Note that overlaid significant sites will be not plotted using a smoothing function.
}

  \item{xrange}{
vector of length 2 containing the desired minimum and maximum concentrations to plot.
}

  \item{yrange}{
vector of length 2 containing the desired minimum and maximum fold changes to plot.
}

  \item{\dots}{
passed to underlying plotting functions.
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
 \code{\link{dba.analyze}}
}

\examples{
data(tamoxifen_analysis)

# default MA plot
dba.plotMA(tamoxifen)

#XY plots (with raw and normalized data)
par(mfrow=c(1,2))
dba.plotMA(tamoxifen,bXY=TRUE,bSmooth=FALSE,bNormalized=FALSE)
dba.plotMA(tamoxifen,bXY=TRUE,bSmooth=FALSE,bNormalized=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
