\name{dba.show}
\alias{dba.show}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List attributes of peaksets of contrasts associated with a DBA object
}
\description{
Returns attributes of peaksets and/or contrasts associated with a DBA object.
}
\usage{
dba.show(DBA, mask, attributes, bContrasts=FALSE, bDesign=FALSE, 
         th=DBA$config$th) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


  \item{DBA}{
DBA object
}
%\subsection{Optional/defaulted arguments}{}
  \item{mask}{
mask of peaksets for which to get attributes (used when obtaining peakset attributes, 
i.e. \code{bContrasts=FALSE}).
}
\item{attributes}{
attribute or vector of attributes to retrieve. 
Number of intervals is always shown. 
Used when obtaining peakset attributes, i.e. \code{bContrasts=FALSE}. Values:
\itemize{ 
\item \code{\link{DBA_ID}}
\item \code{\link{DBA_TISSUE}}
\item \code{\link{DBA_FACTOR}}
\item \code{\link{DBA_CONDITION}}
\item \code{\link{DBA_TREATMENT}}
\item \code{\link{DBA_REPLICATE}}
\item \code{\link{DBA_CONSENSUS}}
\item \code{\link{DBA_CALLER}}
\item \code{\link{DBA_CONTROL}}
\item \code{\link{DBA_READS}}
\item \code{\link{DBA_INTERVALS}}
\item \code{\link{DBA_FRIP}}
}
}

 \item{bContrasts}{
logical indicating whether peaksets or contrast attributes are to be retrieved.
\code{TRUE} retrieves a dataframe of contrast information instead of peakset attributes. 
If no contrasts are set, returns possible contrasts. See  \code{\link{dba.contrast}}.
}

 \item{bDesign}{
logical indicating whether the model design should be returned, if present. 
\code{bContrasts} must be \code{FALSE} for this parameter to be used.
}

  \item{th}{
if \code{bContrasts} is \code{TRUE}, then \code{th} is used as the threshold for determining how many
significant sites there are for each contrast.
Only relevant when obtaining contrast attributes (\code{bContrasts=TRUE}) 
and  \code{\link{dba.analyze}} has been run.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
MODE: Return attributes of peaksets associated with a DBA object:

\code{dba.show(DBA, mask, attributes)}

MODE: Return contrasts associated with a DBA object:

\code{dba.show(DBA,bContrasts=TRUE, th)}

MODE: Return design associated with a DBA object:

\code{dba.show(DBA,bDesign=TRUE)}
}
\value{
dataframe with peakset attributes. 

If \code{bContrasts == FALSE}, each row represents a peakset, and each column is an attributes, with the final column, Intervals, indicating how many sites there are in the peakset.

If \code{bContrasts == TRUE}, each row represent a contrast, with the following columns:
\item{\code{Group1}}{Label for first group of contrast}
\item{\code{Members1}}{Number of samples in first group of contrast}
\item{\code{Group2}}{Label for first group of contrast}
\item{\code{Members3}}{Number of samples in first group of contrast}

if \code{\link{dba.analyze}} has been successfully run, there there will be up to 
four more columns showing the number of significant differentially bound (DB) sites identified for 
\item{\code{DB.edgeR}}{Number of significantly differentially bound sites identified using edgeR}
\item{\code{DB.DESeq}}{Number of significantly differentially bound sites identified using DESeq}
\item{\code{DB.edgeR.block}}{Number of significantly differentially bound sites identified for blocking analysis using edgeR}
\item{\code{DB.DESeq.block}}{Number of significantly differentially bound sites identified for blocking analysis using DESeq}
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
\code{\link{dba}}, \code{\link{dba.peakset}}, \code{\link{dba.contrast}}
\code{\link{dba.analyze}},  \code{\link{DBA.config}}.
}

\examples{
data(tamoxifen_peaks)
dba.show(tamoxifen)
dba.show(tamoxifen,tamoxifen$masks$Responsive)
dba.show(tamoxifen,attributes=c(DBA_TISSUE,DBA_REPLICATE,DBA_CONDITION))

data(tamoxifen_analysis)
dba.show(tamoxifen,bContrasts=TRUE)

#alternatively:
data(tamoxifen_analysis)
tamoxifen
tamoxifen$config$th <- .01
tamoxifen
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
