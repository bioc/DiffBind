\name{dba.plotBox}
\alias{dba.plotBox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Boxplots
}
\description{
Boxplots for read count distributions within differentially bound sites
}
\usage{
dba.plotBox(DBA, contrast=1, method=DBA$config$AnalysisMethod, 
            th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
            bNormalized=TRUE, attribute=DBA_GROUP, mask,
            bAll=FALSE, bAllIncreased=FALSE, bAllDecreased=FALSE, 
            bDB=TRUE, bDBIncreased=TRUE, bDBDecreased=TRUE,
            pvalMethod=wilcox.test,  bReversePos=FALSE, attribOrder, 
            vColors, varwidth=TRUE, notch=TRUE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}

\item{DBA}{
DBA object.
}

%\subsection{Optional/defaulted arguments}{}

\item{contrast}{
number of contrast to use for boxplot.
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
significance threshold; all sites with FDR (or p-values, see bUsePval) 
less than or equal to this value will be included in the boxplot.
}
\item{bUsePval}{
logical indicating whether to use FDR (FALSE) or p-value (TRUE) for thresholding.
}

\item{bNormalized}{
logical indicating that normalized data (using normalization factors computed by differential analysis method) 
should be plotted. \code{FALSE} uses raw count data.
}


\item{attribute}{
attribute to use for determining groups of samples. 
Default (\code{DBA_GROUP}) plots the two groups used in the contrast, if available. 
Possible values:
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
}

\item{mask}{
logical mask of samples to include when no groups are present.
}

\item{bAll}{
logical indicating if plot should include a set of boxplots using all counts, 
regardless of whether or not they pass the significance threshold.
}
\item{bAllIncreased}{
logical indicating if plot should include a set of boxplots using all counts that increase in affinity, regardless of whether or not they pass the significance threshold.
}
\item{bAllDecreased}{
logical indicating if plot should include a set of boxplots using all counts that decrease in affinity, regardless of whether or not they pass the significance threshold.
}
\item{bDB}{
logical indicating if plot should include a set of boxplots using all counts in significantly differentially bound sites (i.e. those that pass the significance threshold), 
regardless of whether they increase or decrease in affinity.
}
\item{bDBIncreased}{
logical indicating if plot should include a set of boxplots using all counts  in significantly differentially bound sites that increase in affinity.
}
\item{bDBDecreased}{
logical indicating if plot should include a set of boxplots using all counts  in significantly differentially bound sites that decrease in affinity.
}

\item{pvalMethod}{
method to use when computing matrix of p-values. 
If \code{NULL}, no matrix is computed, and NULL is returned; this may speed up processing if there are many boxplots.
}

\item{bReversePos}{
logical indicating if the default definition of positive affinity (higher affinity in the second group of the contrast) should be reversed 
(i.e. positive affinity is defined as being higher in the first group of the contrast).
}

\item{attribOrder}{
vector of group numbers used to change the order that groups are plotted. 
If \code{NULL}, default order is used (group order for \code{DBA_GROUP}, and the order the attribute values appear for other values of attribute).
}

\item{vColors}{
vector of custom colors; if absent, default colors will be used.
}
\item{varwidth}{
passed to boxplot
}
\item{notch}{
passed to boxplot
}
\item{\dots}{
other arguments passed to boxplot
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Draws a boxplot showing distributions of read counts for various groups of samples under various conditions.
In default mode, draws six boxes: one pair of boxes showing the distribution of read counts within all
significantly differentially bound sites (one box for each sample group), one pair of boxes showing the
distribution of read counts for significantly differentially bound sites that increase affinity in 
the second group, and a second pair of boxes showing the distribution of read counts 
for significantly differentially bound sites that have higher mean affinity in the first group. 
}
\value{
if \code{pvalMethod} is not \code{NULL}, returns a matrix of p-values indicating the significance of the difference between each pair of distributions.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
%\note{

%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(tamoxifen_analysis)

#default boxplot includes all DB sites, then divided into those increasing 
# affinity in each group
dba.plotBox(tamoxifen)

# plot non-normalized data for DB sites by tissue
# (changing order to place Resistant samples last)
dba.plotBox(tamoxifen, attribute=DBA_CONDITION, bDBIncreased=FALSE,
            bDBDecreased=FALSE, attribOrder=c(2,1), bNormalized=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
