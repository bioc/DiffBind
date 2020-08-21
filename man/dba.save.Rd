\name{dba.save}
\alias{dba.save}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
save DBA object
}
\description{
Writes out DBA object
}
\usage{
dba.save(DBA, file='DBA', dir='.', pre='dba_', ext='RData', 
         bRemoveAnalysis=FALSE, bRemoveBackground=FALSE,
         bCompress=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


  \item{DBA}{
DBA object
}
%\subsection{Optional/defaulted arguments}{}
  \item{file}{
main filename
}
  \item{dir}{
directory to save model in
}
  \item{pre}{
string to pre-pend to filename
}
 \item{ext}{
extensions to use
}
 \item{bRemoveAnalysis}{
if \code{TRUE}, will remove the global \code{DESeq2} and/or \code{edgeR}
analysis objects. 
The analysis results will be retained.
If the analysis objects are required after re-loading, they will
be automatically re-generated.
}
 \item{bRemoveBackground}{
if \code{TRUE}, will remove the global binned background
counts used for normalization.
Any normalization factors calculated using these counts
will be retained.
If the the normalization factors need to be re-re-calculated after re-loading, 
the binned background counts will
be automatically re-generated.
}
\item{bCompress}{
logical indicating saved DBA object should be compressed as much as possible.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
string containing full path and filename.
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
 \code{\link{dba.load}}
}

\examples{
\dontrun{
data(tamoxifen_peaks)
savefile <- dba.save(tamoxifen,'tamoxifenPeaks')
savefile
rm(tamoxifen)
tamoxifen <- dba.load('tamoxifenPeaks')
unlink(savefile)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
