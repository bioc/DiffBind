\name{dba.plotHeatmap}
\alias{dba.plotHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw a binding site heatmap
}
\description{
Draws a binding site heatmap
}
\usage{
dba.plotHeatmap(DBA, attributes=DBA$attributes, maxSites=1000, minval, maxval,
                contrast, method=DBA$config$AnalysisMethod, 
                th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                report, score, bLog=TRUE, mask, sites, sortFun=sd,
                correlations=TRUE, olPlot=DBA_COR, 
                ColAttributes,RowAttributes, colSideCols, rowSideCols=colSideCols,
                margin=10, colScheme="Greens", distMethod="pearson",
                ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{

%\subsection{Required arguments}{}

\item{DBA}{
DBA object. 
}

%\subsection{Optional/defaulted arguments}{}

\item{attributes}{
attribute or vector of attributes to use for column labels:
\itemize{ 
\item \code{\link{DBA_ID}}
\item \code{\link{DBA_TISSUE}}
\item \code{\link{DBA_FACTOR}}
\item \code{\link{DBA_CONDITION}}
\item \code{\link{DBA_TREATMENT}}
\item \code{\link{DBA_REPLICATE}}
\item \code{\link{DBA_CONSENSUS}}
\item \code{\link{DBA_CALLER}}
}
}

\item{maxSites}{
maximum number of binding sites to use in heatmap. 
Only used when not drawing a correlation heatmap (\code{correlations=FALSE})
}

\item{minval}{
Set all scores less than this to minval
}
\item{maxval}{
Set all scores greater than this to maxval
}

\item{contrast}{
number of contrast to report on; if present, draws a heatmap based on a differential binding affinity analysis
(see \code{\link{dba.analyze}}). 
Only significantly differentially bound sites will be used 
(subject to the \code{th} and \code{bUsePval} parameters). 
If mask is unspecified, only the samples in the contrast will be included. 
See  \code{\link{dba.show}(DBA, bContrast=TRUE)} to get contrast numbers. 
If missing, uses scores in the main binding matrix.
}
\item{method}{
analysis method (used in conjunction with contrast):
\itemize{ 
\item \code{\link{DBA_DESEQ2}}
\item \code{\link{DBA_DESEQ2_BLOCK}}
\item \code{\link{DBA_EDGER}}
\item \code{\link{DBA_EDGER_BLOCK}}
}
}

\item{th}{
significance threshold; all sites with FDR (or p-values, see b\code{UsePval})
less than or equal to this value will be included in the report (subject to \code{maxSites}). 
Used in conjunction with \code{contrast}.
}
\item{bUsePval}{
logical indicating whether to use FDR (\code{FALSE}) or p-value (\code{TRUE})
for thresholding. Used in conjunction with \code{contrast}.
}
\item{report}{
report (obtained from \code{\link{dba.report}} specifying the data to be used. 
If this is present, the \code{method}, \code{th}, and \code{bUsePval} parameters are ignored.
Used in conjunction with \code{contrast}.
}
\item{score}{
Score to use for count data. Only used when plotting the global binding matrix (no contrast specified). One of:
\itemize{
\item \code{\link{DBA_SCORE_NORMALIZED}}
\item \code{\link{DBA_SCORE_READS}}
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
Only applicable with read count scores (not peak scores).
}
\item{mask}{
mask indicating a subset of peaksets to use when using global binding matrix scores. 
If a \code{contrast} is specified, these peaksets will be included, 
but only the significantly differentially bound sites 
(using \code{th}, \code{bUsePva}l, and/or \code{report}) will be included.
}
  \item{sites}{
logical vector indicating which sites to include; first \code{maxSites} of these. 
Only relevant when using global binding matrix (\code{contrast} is missing).
}

 \item{sortFun}{
function taking a vector of scores and returning a single value.
Only relevant when using global binding matrix (\code{contrast} is missing).
If not equal to \code{FALSE}, the global binding matrix will be sorted (descending) on the results, 
and the first \code{maxSites} used in the heatmap. 
Recommended sort function options include \code{\link{sd}}, \code{\link{mean}}, \code{\link{median}}, \code{\link{min}}.
}

\item{correlations}{
logical indicating that a correlation heatmap should be plotted (\code{TRUE}). 
If \code{FALSE}, a binding heatmap of scores/reads is plotted. 
This parameter can also be set to a correlation record; 
see  \code{\link{dba.overlap}(mode=DBA_OLAP_ALL)}, in which case a correlation heatmap is plotted
based on the specified correlation record, using the statistic specified in \code{olPlot}.
}

 \item{olPlot}{
if correlations is specified as a dataframe returned by  \code{\link{dba.overlap}}, 
indicates which statistic to plot. 
One of:
\itemize{
\item{\code{\link{DBA_COR}}}  {Correlation}
\item{\code{\link{DBA_OLAP}}}  {Percentage overlap}
\item{\code{\link{DBA_INALL}}}  {number of peaks common to both samples}
}
}
\item{ColAttributes}{
Attribute or vector of attributes to plot for column color bars. 
If missing, all attributes with two or more unique non-NA values will be plotted. 
(For correlation heatmaps, \code{DBA_GROUP} will be plotted in the column color bar by default 
when a \code{contrast} is specified). 
A value of \code{NULL} indicates that no column color bar should be drawn. 
Allowable attribute values include:
\itemize{ 
\item \code{\link{DBA_GROUP}}
\item \code{\link{DBA_TISSUE}}
\item \code{\link{DBA_FACTOR}}
\item \code{\link{DBA_CONDITION}}
\item \code{\link{DBA_TREATMENT}}
\item \code{\link{DBA_REPLICATE}}
\item \code{\link{DBA_CALLER}}
}
}
\item{RowAttributes}{
Attribute or vector of attributes for row color bars. 
Row color bars are only allowed for correlation heatmaps.  
Same values as for \code{ColAttributes} parameter. 
Default is  to draw a row color bar only if a \code{contrast} is specified, 
in which case the plotted attribute is \code{DBA_GROUP} (if present).
}
\item{rowSideCols}{
Vector of colors to use in row color bars.
Uses default colors if missing. 
Can also be a list of color vectors.
}
\item{colSideCols}{
Vector of colors to use in column color bars.
Uses default colors if missing.
Can also be a list of color vectors.
}
  \item{margin}{
margin size of plot
}
  \item{colScheme}{
Color scheme; see \code{\link{colorRampPalette}} 
}
  \item{distMethod}{
distance method for clustering; see \code{\link{Dist}}
}

\item{\dots}{
passed on to \code{\link{heatmap.2}}, e.g. \code{scale} etc.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
MODE: Correlation Heatmap plot using statistics for global binding matrix:

\code{dba.plotHeatmap(DBA, attributes=DBA$attributes, minval, maxval,
                      correlations, olPlot,
                      colScheme="Greens", distMethod="pearson", ...)}
                      
MODE: Correlation Heatmap plot using statistics for significantly differentially bound sites:

\code{dba.plotHeatmap(DBA, attributes=DBA$attributes, minval, maxval,
                      contrast, method=DBA_DESEQ2, th=0.05, bUsePval=F, mask,
                      overlaps, olPlot=DBA_COR,
                      colScheme="Greens", distMethod="pearson", ...)}
                      
MODE: Binding heatmap plot using significantly differentially bound sites:

\code{dba.plotHeatmap(DBA, attributes, maxSites, minval, maxval,
                      contrast, method, th, bUsePval, correlations=FALSE,
                      colScheme, distMethod, ...)}


MODE: Binding heatmap plot using the global binding matrix:

\code{dba.plotHeatmap(DBA, attributes, maxSites, minval, maxval,
                      mask, sites, correlations=FALSE, sortFun,
                      colScheme, distMethod, ...)}
}
\value{
%%  ~Describe the value returned
if  correlations is not \code{FALSE}, the overlap/correlation matrix is returned.

if correlations is \code{FALSE}, the sites used in the heatmap are returned in a
\code{\link{GRanges}} object,
in the row order they appear (top to bottom).
The metadata contains a column for 
each sample (also in the order they are appear in the clustering plot), 
with the values being the actual plotted values.
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
 \code{\link{dba.overlap}}
}

\examples{
data(tamoxifen_peaks)
# peak overlap correlation heatmap
dba.plotHeatmap(tamoxifen)

data(tamoxifen_counts)
# counts correlation heatmap
dba.plotHeatmap(tamoxifen)

data(tamoxifen_analysis)
#correlation heatmap based on all normalized data
dba.plotHeatmap(tamoxifen,contrast=1,th=1)

#correlation heatmap based on DB sites only
dba.plotHeatmap(tamoxifen,contrast=1)

#binding heatmap based on DB sites
dba.plotHeatmap(tamoxifen,contrast=1,correlations=FALSE)

#binding heatmap based on 1,000 sites with highest variance
sites <- dba.plotHeatmap(tamoxifen,contrast=1,th=1,
                         correlations=FALSE,sortFun=var)
sites

data(tamoxifen_counts)
#Examples of  heatmaps using DB sites with different subsets of samples
#exclude T47D
tamoxifen <- dba.contrast(tamoxifen,design=FALSE,
                          group1=tamoxifen$masks$Resistant,
                          group2=c(3:5,10:11)) 
tamoxifen <- dba.analyze(tamoxifen)

# regular heatmaps with samples from two contrast groups only
dba.plotHeatmap(tamoxifen, contrast=1) 
#also include the T47D samples
dba.plotHeatmap(tamoxifen,contrast=1,mask=tamoxifen$masks$All) 

#correlation heatmap without MCF7 
plot(tamoxifen,contrast=1,mask=!tamoxifen$masks$MCF7) 

# binding heatmap using only the MCF7 samples
dba.plotHeatmap(tamoxifen,contrast=1,mask=tamoxifen$masks$MCF7,correlations=FALSE)  

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
