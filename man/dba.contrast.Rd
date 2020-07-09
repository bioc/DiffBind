\name{dba.contrast}
\alias{dba.contrast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set up contrasts for differential binding affinity analysis
}
\description{
Sets up contrasts for differential binding affinity analysis
}
\usage{
dba.contrast(DBA, design=missing(group1), contrast,
             group1, group2=!group1, name1, name2,
             minMembers=3, block, bNot=FALSE, bComplex=FALSE,
             categories=c(DBA_TISSUE,DBA_FACTOR,DBA_CONDITION,DBA_TREATMENT),
             bGetCoefficients=FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

%\subsection{Required arguments}{}


  \item{DBA}{
DBA object with count data
}
%\subsection{Optional/defaulted arguments}{}
\item{design}{
Either a logical value, or a character string containing a valid design formula. 

If a logical value is specified, \code{TRUE} indicates that a design should automatically be generated. If \code{contrast} is missing, contrasts will automatically be added and an appropriate design computed. If a \code{contrast} is specified, it must consist of a character vector of length three, containing a factor and two factor values. No \code{group}s can be specified. If set to \code{FALSE}, the contrast will be added between the \code{groups}, if specified; otherwise, if \code{group} is missing, all possible contrasts will be added.

If a design formula is specified, it must be composed from the following allowable factors:
   \itemize{ 
     \item Tissue
     \item Factor
     \item Condition
     \item Treatment
     \item Replicate
     \item Caller
  }
}
\item{contrast}{
If a \code{design} has been specified (previously or in the current call), the following contrasts forms may be indicated:

\itemize{
  \item Character vector of length three. The first element is a factor from the \code{design}. 
  The second and third elements are values for that factor associated with sample groups.
  \item List of length 1, containing a design matrix column name (as obtained using \code{bGetCoefficients}).
  \item List of length 2, containing two design matrix column names 
  (as obtained using \code{bGetCoefficients}), first the numerator and the second the denominator.
  \item Character vector of length one, containing a design matrix column name 
  (as obtained using \code{bGetCoefficients}).
  \item Numeric vector of the same length as the list of design matrix column names (as obtained using \code{bGetCoefficients}), with a weighting for each column.
  }
}
  \item{group1}{
mask of samples in first group (when adding a specific contrast).
See \code{\link{dba.mask}}.
Can not be used with an explicit \code{design}.
}
  \item{group2}{
mask of samples in second group (when adding a specific contrast). 
See \code{\link{dba.mask}}.
Can not be used with an explicit \code{design}.
}
  \item{name1}{
label for samples in first group (when adding a specific contrast). 
}
  \item{name2}{
label for samples in second group (when adding a specific contrast). 
}
  \item{minMembers}{
when automatically generating contrasts, minimum number of unique samples in a group. 
Must be at least 2, as replicates are strongly advised.
If you wish to do an analysis with no replicates, you can set the group1 and group2 parameters explicitly.
}
  \item{bNot}{
include contrasts consisting of a group and all other samples not in that group 
(indicated by a ! in the contrast name).
}

  \item{bComplex}{
include complex contrasts where groups include samples with the same values for multiple factors.
}

  \item{categories}{
when automatically generating contrasts, attribute or vector of attributes to base contrasts on:
\itemize{ 
\item DBA_ID
\item DBA_TISSUE
\item DBA_FACTOR
\item DBA_CONDITION
\item DBA_TREATMENT
\item DBA_REPLICATE
\item DBA_CALLER
}
}
\item{block}{
blocking attribute for multi-factor analysis.  This may be specified as either a value, a vector, or a list. 

If block is a value, the specified metadata field is used to derive the blocking factor. One of:
\itemize{ 
\item DBA_TISSUE
\item DBA_FACTOR
\item DBA_CONDITION
\item DBA_TREATMENT
\item DBA_REPLICATE
\item DBA_CALLER
}

If block is a vector, it can either be a mask (logical vector) or a vector of peakset numbers. In this case, the peaksets indicated in the blocking vector are all given the same factor value (true), while any peaksets not included in the vector take the alternative factor value (false).

If block is a list, it should be a list of vectors (either logical masks or vectors of peakset numbers), with each indicating a set of peaksets that should share the same value. Each peakset should appear at most once, and any peaksets not specified will be given an default value (other).

}

\item{bGetCoefficients}{
If \code{TRUE}, return the names of the columns (coeffcients) associated with the \code{design}.
These can be used to specify a \code{contrast}. 
If \code{bGetCoefficients=TRUE}, all other parameters 
(except \code{DBA} and \code{design}, if specified) will be ignored.
}
}
\details{

MODE: Set up a specific contrast using a design:

\code{dba.contrast(DBA, design, contrast)}

MODE: Set up all possible contrasts:

\code{dba.contrast(DBA, minMembers, categories)}

MODE: Set up a specific contrast without an explicit design:

\code{dba.contrast(DBA, design=FALSE, group1, group2, name1, name2, block)}

}
\value{
DBA object with contrast(s) set as \code{DBA$contrasts}. 

Contrast list can be retrieved using \code{dba.show(DBA, bContrasts=TRUE)}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Rory Stark
}
\note{
Contrasts will only be set up for peaksets where \code{DBA_CALLER == "counts"}.

Contrasts can be cleared by \code{DBA$contrasts <- NULL}.
}

\seealso{

\code{\link{dba.analyze}}

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Set up an explicit contrast
data(tamoxifen_counts)
tamoxifen <- dba.contrast(tamoxifen, contrast=c("Condition","Responsive","Resistant"))
tamoxifen
tamoxifen <- dba.analyze(tamoxifen)
tamoxifen

# Add another contrast
tamoxifen <- dba.contrast(tamoxifen, contrast=c("Tissue","MCF7","BT474"))
dba.show(tamoxifen,bDesign=TRUE)

# Change design
tamoxifen <- dba.contrast(tamoxifen,design="~Tissue + Condition")
tamoxifen <- dba.analyze(tamoxifen)
tamoxifen

# Automatically add all contrasts between sample groups 
# where at least THREE samples have the same factor value
data(tamoxifen_counts)
tamoxifen <- dba.contrast(tamoxifen)
tamoxifen

# Automatically add all contrasts between sample groups 
# where at least TWO samples have the same factor value
tamoxifen <- dba.contrast(tamoxifen, minMembers=2)
dba.show(tamoxifen,bContrasts=TRUE)

# Use of complex contrasts                          
data(tamoxifen_counts)
tamoxifen <- dba.contrast(tamoxifen, contrast=c("Tissue","MCF7","BT474"))
dba.show(tamoxifen,bDesign=TRUE)
tamoxifen <- dba.contrast(tamoxifen,design="~Tissue + Condition")
tamoxifen

dba.contrast(tamoxifen, bGetCoefficients=TRUE)
tamoxifen <- dba.contrast(tamoxifen,contrast="Tissue_MCF7_vs_BT474")
tamoxifen <- dba.contrast(tamoxifen,contrast=list("Tissue_MCF7_vs_BT474"))
tamoxifen <- dba.contrast(tamoxifen,contrast=c(0,1,0,0,0))
tamoxifen <- dba.contrast(tamoxifen,
                          contrast=list("Tissue_MCF7_vs_BT474","Tissue_T47D_vs_BT474"))
tamoxifen <- dba.contrast(tamoxifen,contrast=c(0,1,1,1,0))
tamoxifen <- dba.contrast(tamoxifen,contrast=c(0,0,0,0,1))
dba.show(tamoxifen,bContrasts=TRUE)
tamoxifen <- dba.analyze(tamoxifen)
tamoxifen
tamoxifen <- dba.contrast(tamoxifen, 
                          contrast=c("Condition","Responsive","Resistant"))
tamoxifen <- dba.analyze(tamoxifen)
dba.show(tamoxifen,bContrasts=TRUE)[7:8,]


# Explicit contrast, without design
data(tamoxifen_counts)
tamoxifen <- dba.contrast(tamoxifen, design=FALSE,
                          group1=tamoxifen$masks$Responsive, name1="Responsive",
                          group2=tamoxifen$masks$Resistant,  name2="Resistant",
                          block=DBA_TISSUE)
dba.show(tamoxifen, bContrasts=TRUE)
tamoxifen <- dba.analyze(tamoxifen)
dba.show(tamoxifen,bContrasts=TRUE)
dba.plotVenn(tamoxifen,contrast=1,method=c(DBA_DESEQ2,DBA_DESEQ2_BLOCK))
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
