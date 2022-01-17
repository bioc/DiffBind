\name{dba.normalize}
\alias{dba.normalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specify parameters for normalizing a dataset;
calculate library sizes and normalization factors.
}
\description{
Enables normalization of datasets using a variety of  methods,
including background, spike-in, and parallel factor normalization.
Alternatively, allows a user to specify library sizes 
and normalization factors directly, or retrieve computed ones.
}
\usage{
dba.normalize(DBA, method = DBA$config$AnalysisMethod,
              normalize = DBA_NORM_DEFAULT, library = DBA_LIBSIZE_DEFAULT, 
              background = FALSE, spikein = FALSE, offsets = FALSE,
              libFun=mean, bRetrieve=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DBA}{
DBA object that includes count data for a consensus peakset.
}
\item{method}{
Underlying method, or vector of methods, for which to normalize.

Supported methods:
\itemize{ 
  \item{\code{\link{DBA_EDGER}}} {use \code{edgeR} package for analysis}
  \item{\code{\link{DBA_DESEQ2}}} {use \code{DESeq2} package for analysis}
  \item{\code{\link{DBA_ALL_METHODS}}} {normalize for both 
 both \code{edgeR} and \code{DESeq2}}
}
}

  \item{normalize}{
Either user-supplied normalization factors in a numeric vector,
or a specification of a method to use to calculate normalization factors.


Methods can be specified using one of the following:
 \itemize{
   \item{\code{DBA_NORM_RLE} ("RLE")} 
   {RLE normalization (native to \code{DBA_DESEQ2},
   and available for \code{DBA_EDGER}).}
   \item{\code{DBA_NORM_TMM} ("TMM")} 
   {TMM normalization (native to \code{DBA_EDGER},
   and available for \code{DBA_DESEQ2}).}
   \item{\code{DBA_NORM_NATIVE} ("native")} 
   {Use native method based on \code{method}: 
   \code{DBA_NORM_RLE} for \code{DBA_DESEQ2} or
   \code{DBA_NORM_TMM} for \code{DBA_EDGER}.}
   \item{\code{DBA_NORM_LIB} ("lib")}     
   {Normalize by library size only.
   Library sizes can be specified using the \code{library} parameter.
   Normalization factors will be calculated to give each equal weight
   in a manner appropriate for the analysis \code{method}.
   See also the \code{libFun} parameter, which can be used to
   scale the normalization factors for \code{DESeq2.}}
   \item{\code{DBA_NORM_DEFAULT} ("default")} 
   {Default method: The "preferred" normalization approach
   depending on \code{method} and whether an explicit
   design is present. See \code{Details} below.}
   \item{\code{DBA_NORM_OFFSETS} ("offsets")} 
   {Indicates that offsets have been specified using the \code{offsets}
   parameter, and they should be used without alteration.}
   \item{\code{DBA_NORM_OFFSETS_ADJUST} ("adjust offsets")} 
   {Indicates that offsets have been specified using the \code{offsets}
   parameter, and they should be adjusted for library size and mean centering
   before being used in a \code{DBA_DESEQ2} analysis.}
 }
}

\item{library}{
Either user-supplied library sizes in a numeric vector,
or a specification of a method to use to calculate library sizes.


Library sizes can be based on one of the following:
 \itemize{
  \item{\code{DBA_LIBSIZE_FULL} ("full")} 
  {Use the full library size (total number of reads in BAM/SAM/BED file)}
  \item{\code{DBA_LIBSIZE_PEAKREADS} ("RiP")} 
  {Use the number of reads that overlap consensus peaks.}
  \item{\code{DBA_LIBSIZE_BACKGROUND} ("background")} 
  {Use the total number of reads aligned to the chromosomes for which there is
  at least one peak. This required a background bin calculation (see parameter
  \code{background}). These values are usually the same or similar to
  \code{DBA_LIBSIZE_FULL}.}
  \item{\code{DBA_LIBSIZE_DEFAULT} ("default")} 
  {Default method: The "preferred" library size
   depending on \code{method}, \code{background}, 
   and whether an explicit design is present. 
   See \code{Details} below.}
 }
}

\item{background}{
This parameter controls the option to use "background" bins, which
should not have differential enrichment between samples,
as the basis for normalizing (instead of using reads counts
overlapping consensus peaks). 
When enabled, the chromosomes for which there are peaks
in the consensus peakset are tiled into large bins and reads
overlapping these bins are counted.

If present, \code{background} can either be a logical value, a numeric value,
or a previously computed \code{$background} object.

If \code{background} is a logical value and set to \code{TRUE}, background bins will be
computed using the default bin size of 15000bp.
Setting this value to \code{FALSE} will prevent background mode
from being used in any default settings.

If \code{background} is a numeric value, it will be used as the bin size.

If \code{background} is a previously  computed
\code{$background} object, these counts will be used as the background.
A \code{$background} object can be obtained by calling 
\code{\link{dba.normalize}} with \code{bRetrieve=TRUE}
and \code{method=DBA_ALL_METHODS}.

After counting (or setting) background bins,
both the \code{normalize} and \code{library} parameters
will be used to determine how the
final normalization factors are calculated.

If \code{background} is missing, it will be set to \code{TRUE} if 
\code{library=DBA_LIBSIZE_BACKGROUND}, or if 
\code{library=DBA_LIBSIZE_DEFAULT} and certain conditions are met
(see \code{Details} below).

If \code{background} is not \code{FALSE}, then
the library size will be set to \code{library=DBA_LIBSIZE_BACKGROUND}
}

\item{spikein}{
Either a logical value, a character vector of chromosome names,
a \code{\link{GRanges}} object containing peaks for a parallel factor,
or a \code{$background} object containing previously computed 
spike-in read counts.

If \code{spikein} is a logical value set to \code{FALSE}, 
no spike-in normalization is performed.

If \code{spikein} is a logical value set to  \code{TRUE},
background normalization is performed using spike-in tracks.
There must be a spike-in track for each sample.
see \code{\link{dba}} and/or \code{\link{dba.peakset}} for
details on how to include a spike-in track with a sample
(eg. by including a \code{Spikein} column in the sample sheet.)
All chromosomes in the spike-in bam files will be used.

If \code{spikein} is a character vector of one or more chromosome names,
only reads on the named chromosome(s) will be used for background normalization.
If spike-in tracks are available, reads on chromosomes with these names in the
spike-in track will be counted.
If no spike-in tracks are available, reads on chromosomes with these names
in the main \code{bamReads} bam files will be counted.

If \code{spikein} is a \code{\link{GRanges}} object containing peaks
for a parallel factor, then background normalization is performed counting reads
in the spike-in tracks overlapping peaks in this object.

If \code{spikein} is a previously  computed
\code{$background} object, these counts will be used as the spikein background.
A \code{$background} object can be obtained by calling 
\code{\link{dba.normalize}} with \code{bRetrieve=TRUE}
and \code{method=DBA_ALL_METHODS}.

Note that if \code{spikein} is not \code{FALSE}, then
the library size will be set to \code{library=DBA_LIBSIZE_BACKGROUND}.
When \code{spikein}s are used, the default is \code{normalize=DBA_NORM_NATIVE}.
}
\item{offsets}{
This parameter controls the use of offsets (matrix of normalization factors)
instead of a single normalization factor for each sample. It can either 
be a logical value, a \code{matrix}, or a \code{\link{SummarizedExperiment}}.

If it is a logical value and set to \code{FALSE}, 
no offsets will be computed or used. A value of \code{TRUE}
indicates that an offset matrix should be computed using a \code{loess} fit.

Alternatively, user-calculated normalization offsets can be supplied
as a \code{matrix} or as a  \code{\link{SummarizedExperiment}} 
(containing an \code{assay} named "offsets").
In this case, the user may also set the 
\code{normalize} parameter to indicate whether the 
offsets should be applied as-is to a  \code{DESeq2} analysis 
(\code{DBA_NORM_OFFSETS}, default),
or if they should be adjusted for library size and mean centering
(\code{DBA_NORM_OFFSETS_ADJUST}).

}

\item{libFun}{
When \code{normalize=DBA_NORM_LIB}, normalization factors are
calculated by dividing the library sizes for each sample by
a common denominator, obtained by applying \code{libFun}
to the vector of library sizes.

For \code{method=DBA_EDGER}, the normalization factors are further 
adjusted so as to make all the effective library sizes (library sizes multiplied by normalization factors) the same, and adjusted to multiply to 1.
}

\item{bRetrieve}{
If set to \code{TRUE}, information about the current normalization will be 
returned. 
The only other relevant parameter in this case is the \code{method}.

If \code{method=DBA_DESEQ2} or \code{method=DBA_EDGER}, 
a record will be returned including normalization values 
for the appropriate analysis method. This record is a \code{list}
consists of the following elements:
 \itemize{
   \item{$norm.method}
   {A character string corresponding to the normalization method,
   generally one of the values that can be supplied as a value to 
   \code{normalize}.}
   \item{$norm.factors}
   {A vector containing the computed normalization factors.}
   \item{$lib.method}
    {A character string corresponding to the value of the method used to
    calculate the library size,
   generally one of the values that can be supplied as a value to 
   \code{library}.}
   \item{$lib.sizes}
   {A vector containing the computed library sizes.}
%   \item{$filter.value}
%   {A numeric value indicating the minimum values used to 
%   filter consensus peaks.} 
   \item{$background}
   {If the normalization if based on binned background reads,
   this field will be \code{TRUE}.}
   \item{$control.subtract}
   {If control reads were subtracted from the read counts,
   this field will be \code{TRUE}.}
 }
 
 If \code{method=DBA_ALL_METHODS}, the record be a list
 with one of the above records for each \code{method} for
 which normalization factors have been computed 
 (\code{$DESeq2} and \code{edgeR}).
 

 If \code{background} bins have been calculated,
 this will include an element called \code{$background}.
  This element can be passed in as the value to \code{background} or 
 \code{spikein} to re-use a previously computed set of reads.
 It contains three subfields:
 \itemize{
   \item{$background$binned}
   {a \code{SummarizedExperiment} object containing the binned counts.}
   \item{$background$bin.size}
   {a numeric value with the bin size used.}
      \item{$background$back.calc}
   {character string indicating how the background was calculated 
   (bins, spike-ins, or parallel factor).}
 }

 
 If \code{offsets} are available,
 this will include an element called \code{$offsets} with two
 subfields:
 \itemize{
   \item{$offsets$offsets}
   {a \code{matrix} or a \code{SummarizedExperiment} object 
   containing the offsets.}
   \item{offsets$offset.method}
   {a character string indicating the source of the offsets, either
   \code{"loess"} or \code{"user"}.}
 }
 
}

\item{...}{
Extra parameters to be passed to \code{limma::loessFit}
when computing offsets.
}
}
\details{
The default normalization parameters are as follows:
\itemize{
  \item{\code{normalize=DBA_NORM_LIB}}
  \item{\code{library=DBA_LIBSIZE_FULL}}
  \item{\code{background=FALSE}}
}

If \code{background=TRUE}, then the default becomes
\code{library=DBA_LIBSIZE_BACKGROUND}.

If \code{\link{dba.contrast}} has been
used to set up contrasts with \code{design=FALSE} (pre-3.0 mode),
then the defaults are:
\itemize{
  \item{\code{normalize=DBA_NORM_DEFAULT}}
  \item{\code{library=DBA_LIBSIZE_FULL}}
  \item{\code{background=FALSE}}
}
In this case, \code{normalize=DBA_NORM_LIB} will be set for 
\code{method=DBA_DESEQ2} for backwards compatibility.
}

\value{
Either a \code{DBA} object with normalization terms added,
or (if \code{bRetrieve=TRUE}), a record
or normalization details.
}

\author{
Rory Stark
}

\note{
The \code{csaw} package is used to compute  
\code{background} bins and \code{offsets} based on 
\code{limma::loessFit}.

See the \code{DiffBind} vignette for technical details of how this
is done, and the \code{csaw} vignette for details on 
background bins and loess offsets can be used to address 
different biases in ChIP-seq data.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dba.count}}, \code{\link{dba.analyze}}, \code{\link{dba.save}}
}

\examples{
# load DBA object with counts 
data(tamoxifen_counts)
tamoxifen <- dba.contrast(tamoxifen,design="~Tissue + Condition")

# default normalization: Full library sizes
tamoxifen <- dba.normalize(tamoxifen)
dba.normalize(tamoxifen, bRetrieve=TRUE)
dba.analyze(tamoxifen)

# RLE/TMM using Reads in Peaks
tamoxifen <- dba.normalize(tamoxifen, method=DBA_ALL_METHODS,
                           normalize=DBA_NORM_NATIVE, 
                           library=DBA_LIBSIZE_PEAKREADS)
dba.normalize(tamoxifen, method=DBA_DESEQ2, bRetrieve=TRUE)
dba.normalize(tamoxifen, method=DBA_EDGER, bRetrieve=TRUE)
tamoxifen <- dba.analyze(tamoxifen, method=DBA_ALL_METHODS)
dba.show(tamoxifen,bContrasts=TRUE)
dba.plotVenn(tamoxifen,contrast=1,method=DBA_ALL_METHODS,bDB=TRUE)

# TMM in Background using precomputed background
norm <- dba.normalize(tamoxifen,method=DBA_ALL_METHODS,bRetrieve=TRUE)
tamoxifen <- dba.normalize(tamoxifen, background=norm$background,
                           normalize="TMM", method=DBA_ALL_METHODS)
tamoxifen <- dba.analyze(tamoxifen)
dba.show(tamoxifen,bContrasts=TRUE)
dba.plotMA(tamoxifen)

# LOESS offsets
tamoxifen <- dba.normalize(tamoxifen, method=DBA_ALL_METHODS, offsets=TRUE)
tamoxifen <- dba.analyze(tamoxifen, method=DBA_ALL_METHODS)
dba.show(tamoxifen,bContrasts=TRUE)

par(mfrow=c(3,1))
dba.plotMA(tamoxifen,th=0,bNormalized=FALSE)
dba.plotMA(tamoxifen,method=DBA_DESEQ2)
dba.plotMA(tamoxifen,method=DBA_EDGER)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
