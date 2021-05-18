\name{DBA tamoxifen resistance dataset}
\alias{tamoxifen}
\alias{tamoxifen_peaks}
\alias{tamoxifen_counts}
\alias{tamoxifen_analysis}
\alias{tamoxifen_greylist}
\alias{tamoxifen.greylist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tamoxifen resistance dataset used for DBA examples
}
\description{
Tamoxifen resistance dataset used for DBA examples
}
\usage{
data(tamoxifen_peaks)

data(tamoxifen_counts)

data(tamoxifen_analysis)

data(tamoxifen_greylist)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tamoxifen_peaks}{
load tamoxifen resistance dataset DBA object with peak (occupancy) data
}
\item{tamoxifen_counts}{
load tamoxifen resistance dataset DBA object with count (affinity) data.
Also includes background bins counts for background normalization.
}

\item{tamoxifen_analysis}{
load tamoxifen resistance dataset DBA object with count (affinity) data 
and \code{DESeq2}-based differential binding analysis results.
This analysis uses a blacklists, computed greylists,
 background normalization, and a two-factor design.
}

\item{tamoxifen_greylist}{
load greylist for tamoxifen dataset. 
Generated as shown in \code{dba.blacklist} example: \code{\link{dba.blacklist}}.
}
}
\details{
The tamoxifen resistance dataset is used for the DBA vignette 
and man page examples. 

Data used to create these objects can be downloaded at
\url{https://content.cruk.cam.ac.uk/bioinformatics/software/DiffBind/DiffBind_vignette_data.tar.gz}.
}
\value{
loads a DBA object named \code{tamoxifen} (or \code{tamoxifen.greylist}).
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
\note{
The script for generating these files (\code{GenerateDataFiles.R})
is included with the package in the \code{inst/extras} directory.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(tamoxifen_peaks)
tamoxifen
data(tamoxifen_counts)
plot(tamoxifen)
data(tamoxifen_analysis)
dba.plotMA(tamoxifen)
data(tamoxifen_greylist)
tamoxifen.greylist$master
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
