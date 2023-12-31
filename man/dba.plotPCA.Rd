\name{dba.plotPCA}
\alias{dba.plotPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PCA plot
}
\description{
Principal Component Analysis plot
}
\usage{
dba.plotPCA(DBA, attributes, minval, maxval,
           contrast, method=DBA$config$AnalysisMethod, 
           th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
           report, score, bLog=TRUE, mask, sites, label, cor=FALSE,
           b3D=FALSE, vColors, dotSize, labelSize, labelCols, 
           components=1:3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}



\item{DBA}{
DBA object.
}

%\subsection{Optional/defaulted arguments}{}

\item{attributes}{
attribute or vector of attributes to use to color plotted points. 
Each unique combination of attribute values will be assigned a color. 
Chosen from:
\itemize{ 
\item \code{\link{DBA_GROUP}}
\item \code{\link{DBA_ID}}
\item \code{\link{DBA_TISSUE}}
\item \code{\link{DBA_FACTOR}}
\item \code{\link{DBA_CONDITION}}
\item \code{\link{DBA_TREATMENT}}
\item \code{\link{DBA_REPLICATE}}
\item \code{\link{DBA_CONSENSUS}}
\item \code{\link{DBA_CALLER}}
}
Note that \code{DBA_GROUP} is a special attribute which will result in samples from each group in a contrast
(if present) being colored separately.
}

\item{minval}{
Set all scores less than this to minval
}
\item{maxval}{
Set all scores greater than this to maxval
}

\item{contrast}{
number of contrast to use for PCA; if present, plots a PCA based on a differential binding affinity analysis
(see  \code{\link{dba.analyze}}). 
If \code{mask} is unspecified, only the samples in the contrast will be included. 
See  \code{\link{dba.show}(DBA, bContrast=T)} to get contrast numbers. 
If missing, uses scores in the main binding matrix.
}
\item{method}{
method used for analysis (used in conjunction with contrast):
\itemize{ 
\item \code{\link{DBA_DESEQ2}}
\item \code{\link{DBA_DESEQ2_BLOCK}}
\item \code{\link{DBA_EDGER}}
\item \code{\link{DBA_EDGER_BLOCK}}
}
}

\item{th}{
significance threshold; all sites with FDR (or p-values, see \code{bUsePval}) 
less than or equal to this value will be included in the PCA, subject to \code{maxVal}. 
Used in conjunction with \code{contrast}.
}
\item{bUsePval}{
if \code{TRUE}, uses p-value instead of FDR for thresholding. Used in conjunction with \code{contrast}.
}
\item{report}{
report (obtained from  \code{\link{dba.report}}) specifying the data to be used .
If this is present, the \code{method}, \code{th}, and \code{bUsePval} parameters are ignored.
}

\item{score}{
Score to use for count data.
Only used when plotting the global binding matrix (no \code{contrast} specified). 
One of:
\itemize{
\item \code{\link{DBA_SCORE_READS}}
\item \code{\link{DBA_SCORE_NORMALIZED}}
\item \code{\link{DBA_SCORE_CONTROL_READS}}
\item \code{\link{DBA_SCORE_READS_MINUS}}
%\item \code{\link{DBA_SCORE_READS_FULL}}
%\item \code{\link{DBA_SCORE_READS_EFFECTIVE}}
%\item \code{\link{DBA_SCORE_READS_MINUS_FULL}}
%\item \code{\link{DBA_SCORE_READS_MINUS_EFFECTIVE}}
\item \code{\link{DBA_SCORE_READS_FOLD}}
\item \code{\link{DBA_SCORE_RPKM}}
\item \code{\link{DBA_SCORE_RPKM_FOLD}}
\item \code{\link{DBA_SCORE_RPKM_MINUS}}
%\item \code{\link{DBA_SCORE_TMM_READS_FULL}}
%\item \code{\link{DBA_SCORE_TMM_READS_EFFECTIVE}}
%\item \code{\link{DBA_SCORE_TMM_MINUS_FULL}}
%\item \code{\link{DBA_SCORE_TMM_MINUS_EFFECTIVE}}
}
}
\item{bLog}{
Logical indicating that log2 values should be used. 
Only applicable to read count scores (not peak scores).
}
\item{mask}{
mask indicating a subset of peaksets to use when using global binding matrix scores.
If a \code{contrast} is specified, these peaksets will be included, 
but only the significantly differentially bound sites 
(using \code{th}, \code{bUsePval}, or \code{report}) will be included. 
See  \code{\link{dba.mask}}.
}
 \item{sites}{
logical vector indicating which sites to include in PCA. 
Only relevant when using global binding matrix (\code{contrast} is missing).
}
 \item{label}{
A metadata field to use as a label in 2D plots. 
The value for this field will be written directly on the plot near the dot for each sample. 
Values can be any of those valid for the \code{attributes} parameter.
}

 \item{cor}{
a logical value indicating whether the calculation should use the 
correlation matrix or the covariance matrix.
Passed into princomp.
}
\item{b3D}{
logical indicating that three principal components should be plotted (requires package \code{rgl}). 
If \code{FALSE}, the first two principal components are plotted.
}

\item{vColors}{
vector of custom colors; is absent, default colors will be used.
}
\item{dotSize}{
size of dots to plot; is absent, a default  will be calculated.
}

\item{labelSize}{
Scaling factor for labels if present. Default is 0.8.
}

\item{labelCols}{
Vector of colors to use for labels. Default is "black".
}

\item{components}{
Number(s) of the components to plot.
Can be a vector of two or three component numbers, or a single integer.
If an integer, that component, 
in addition to the succeeding one (\code{b3D=FALSE}) or two (\code{b3D=TRUE})
will be plotted.
}
\item{\dots}{
arguments passed to \code{\link{plot}} or \code{plot3d} (rgl).
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
MODE: PCA plot using significantly differentially bound sites:

\code{dba.plotPCA(DBA, attributes, minval, maxval,
                       contrast, method, th, bUsePval,
                       b3D=F, vColors, dotSize, ...)}

MODE: PCA plot using global binding matrix:

\code{dba.plotPCA(DBA, attributes, minval, maxval,
                       mask, sites,
                       b3D=F, vColors, dotSize, ...)}
}
\value{
\code{trellis} plot from \code{\link{lattice}} package; see \code{\link{xyplot}}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
\note{
uses rgl package for 3D plots (if available)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{dba.analyze}},  \code{\link{dba.plotHeatmap}}
}

\examples{
data(tamoxifen_peaks)

# peakcaller scores PCA
dba.plotPCA(tamoxifen)

# raw count correlation PCA
data(tamoxifen_analysis)
dba.plotPCA(tamoxifen)

#PCA based on normalized data for all sites
dba.plotPCA(tamoxifen,contrast=1,th=1)

#PCA based on DB sites only
p <- dba.plotPCA(tamoxifen,contrast=1)
p <- dba.plotPCA(tamoxifen,contrast=1,attributes=DBA_TISSUE)
p <- dba.plotPCA(tamoxifen,contrast=1,attributes=DBA_TISSUE,label=DBA_CONDITION)
p <- dba.plotPCA(tamoxifen,contrast=1,attributes=DBA_CONDITION,label=DBA_TISSUE)
p <- dba.plotPCA(tamoxifen,contrast=1,attributes=c(DBA_TISSUE,DBA_CONDITION),
                label=DBA_REPLICATE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
