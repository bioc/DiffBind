\name{dba.analyze}
\alias{dba.analyze}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform differential binding affinity analysis
}
\description{
Performs differential binding affinity analysis
}
\usage{
dba.analyze(DBA, method=DBA$config$AnalysisMethod, 
            bSubControl=TRUE, bFullLibrarySize=TRUE, bTagwise=TRUE,
            filter=0, filterFun=max,
            bCorPlot=DBA$config$bCorPlot,  
            bReduceObjects=TRUE, 
            bParallel=DBA$config$RunParallel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}

  \item{DBA}{
DBA object. If no contrasts are specified (DBA$contrast is NULL), default contrasts will be added via a call to \code{\link{dba.contrast}}.
}
%\subsection{Optional/defaulted arguments}{}
\item{method}{
method, or vector of methods, by which to analyze differential binding affinity. Supported methods:
\itemize{ 
\item \link{DBA_EDGER}
\item \link{DBA_DESEQ2}
}
also, for backward compatibility:
\itemize{
\item \link{DBA_DESEQ}
}
Additionally, if this value is set to \code{DBA_ALL_METHODS}, this is equivalent to specifying \code{c(DBA_EDGER,DBA_DESEQ2)}.
}

\item{bSubControl}{
logical indicating whether Control read counts are subtracted for each site in each sample before performing analysis.
}
\item{bFullLibrarySize}{
logical indicating if the full library size (total number of reads in BAM/SAM/BED file) for each sample is used for scaling normalization. If FALSE, the total number of reads present in the peaks for each sample is used (generally preferable if overall biding levels are expected to be similar between samples).
}
\item{bTagwise}{
logical indicating if dispersion should be calculated on a tagwise (or per-condition) basis. If there are only a very few members of each group in a contrast (e.g. no replicates), this should be set to FALSE.
}
\item{filter}{
value to use for filtering intervals with low read counts. 
Each contrast will be filtered separately.
The \code{filterFun} will be applied t each interval, and any scores below the \code{filter}
value will be removed prior to analysis.
}
\item{filterFun}{
function that will be invoked for each interval with a vector of scores for each sample. Returns a score that will be evaluated against the \code{filter} value (only intervals with a score at least as high as \code{filter} will be kept). Default is \code{max}, indicating that at least one sample should have a score of at least \code{filter}; other useful values include \code{sum} (indicating that all the scores added together should be at least \code{filter}) and \code{mean} (setting a minimum mean normalized count level). Users can supply their own function as well.
}
\item{bCorPlot}{
logical indicating whether to plot a correlation heatmap for the analyzed data (first contrast only). 
If no sites are significantly differentially bound using the default threholds, no heatmap will be plotted.
}
\item{bReduceObjects}{
logical indicating whether strip the analysis objects of unnecessary fields to save memory. If it is desired to used the DBA$contrasts[[n]]$edgeR and/or DBA$contrasts[[n]]$DESeq objects directly in the edgeR and/or DESeq packages, this should be set to FALSE.
}
\item{bParallel}{
logical indicating that the analyses is to be done in parallel using multicore (one process for each contrast for each method, plus an additional process per method).
}
}
\details{
See the DBA User Guide for more details on how the \code{edgeR} and \code{DESeq2} analyses are carried out.  
}
\value{
DBA object with results of analysis added to DBA$contrasts. 
}

\author{
Rory Stark
}
\note{
If there is a blocking factor for the contrast(s) specified using a previous call to \code{\link{dba.contrast}}, a multi-factor analysis will automatically be carried out in addition to a single factor analysis.
}


\seealso{

\code{\link{dba.contrast}}, \code{\link{dba.report}}

}


\examples{
data(tamoxifen_counts)

tamoxifen <- dba.analyze(tamoxifen)
tamoxifen

data(tamoxifen_counts)
tamoxifen <- dba.contrast(tamoxifen,categories=DBA_CONDITION,block=tamoxifen$masks$MCF7)
tamoxifen <- dba.analyze(tamoxifen,method=DBA_ALL_METHODS)
tamoxifen
}

