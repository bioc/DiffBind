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
           bSignificant=TRUE, bSmooth=TRUE, bLoess=TRUE,
           xrange, yrange, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


  \item{DBA}{
DBA object, on which  \code{\link{dba.analyze}} 
should have been successfully run.
}
%\subsection{Optional/defaulted arguments}{}
  \item{contrast}{
number of contrast to report on. 
See  \code{\link{dba.show}(DBA, bContrast=TRUE)} to get contrast numbers.

Alternatively, an MA plot can be generated without a specific contrast,
plotting one set of samples against another.
In this case, \code{contrast} should be a list on length one or two.
Each element of the list should be either a logical
sample mask, or a vector of sample numbers.
If the second set of samples is jot specified (list is length one),
all the samples other than those specified will be used for the second group.
The list elements should be named; these names will be used as labels
for the sample groups in the plot.
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
logical indicating whether to use FDR (\code{FALSE}) or p-value (\code{TRUE})
for thresholding.
}
  \item{fold}{
will only include sites with fold change greater than this as significant 
(colored red).
}
  \item{bNormalized}{
logical indicating whether to plot normalized data using normalization factors 
computed by differential analysis method (\code{TRUE}) 
or raw read counts (\code{FALSE}).
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
logical indicating whether to draw MA plot (\code{FALSE}) or 
XY scatter plot (\code{TRUE}).
}

  \item{dotSize}{
size of points on plot (\code{cex}).
}
  \item{bSignificant}{
Logical indicating if points corresponding to significantly 
differentially bound sites (based on \code{contrast},
\code{th}, \code{bUsePval}, and \code{fold} parameters) 
should be overlaid in red.
}
  \item{bSmooth}{
logical indicating that basic plot should be plotted using 
\code{\link{smoothScatter}}.
Note that overlaid significant sites will be not plotted 
using a smoothing function.
}
  \item{bLoess}{
logical indicating that a MA plot should include a fitted loess curve.
}

  \item{xrange}{
vector of length 2 containing the desired minimum and 
maximum concentrations to plot.
}

  \item{yrange}{
vector of length 2 containing the desired minimum and 
maximum fold changes to plot.
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

# Show different normalizations
tamoxifen <- dba.normalize(tamoxifen,method=DBA_ALL_METHODS, 
                           library=DBA_LIBSIZE_PEAKREADS, background=FALSE)
tamoxifen <- dba.analyze(tamoxifen, method=DBA_ALL_METHODS)

par(mfrow=c(3,2))
dba.plotMA(tamoxifen,th=0,bNormalized=FALSE,sub="NON-NORMALIZED")
dba.plotMA(tamoxifen,th=0,bNormalized=FALSE,sub="NON-NORMALIZED")

dba.plotMA(tamoxifen,method=DBA_DESEQ2,bNormalized=TRUE,
           sub="DESeq2_MRN-RiP")
dba.plotMA(tamoxifen,method=DBA_EDGER,bNormalized=TRUE,
           sub="edgeR_TMM-RiP")

tamoxifen <- dba.normalize(tamoxifen, method=DBA_ALL_METHODS,
                           normalize=DBA_NORM_LIB, background=FALSE) 
tamoxifen <- dba.analyze(tamoxifen,method=DBA_ALL_METHODS)

dba.plotMA(tamoxifen,method=DBA_DESEQ2,bNormalized=TRUE, 
           sub="DESeq2_LIB-FULL")
dba.plotMA(tamoxifen,method=DBA_EDGER,bNormalized=TRUE,
           sub="edgeR_LIB-FULL")
           
# MA plots of samples without a contrast
data(tamoxifen_counts)
par(mfrow=c(2,2))
dba.plotMA(tamoxifen,list(Resistant=tamoxifen$masks$Resistant,
                          Responsive=tamoxifen$masks$Responsive),
                          bNormalized=FALSE)
dba.plotMA(tamoxifen,list(MCF7=tamoxifen$masks$MCF7),
                          bNormalized=FALSE)
dba.plotMA(tamoxifen, list(Sample1=1), bNormalized=FALSE)
dba.plotMA(tamoxifen, list(Random=sample(1:11,5)), bNormalized=FALSE)

#XY plots (with raw and normalized data)
data(tamoxifen_analysis)
par(mfrow=c(1,2))
dba.plotMA(tamoxifen,bXY=TRUE,bSmooth=FALSE,bNormalized=FALSE, 
           sub="NON_NORMALIZED")
dba.plotMA(tamoxifen,bXY=TRUE,bSmooth=FALSE,bNormalized=TRUE,
           sub="DESeq2-MRN-Background")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
