\name{DiffBind3}
\alias{DiffBind3}

\title{Differences between DiffBind 3.0 and earlier versions}


\description{Notes on the differences between 
DiffBind 3.0 and previous versions, and how run in a 
"backward compatible" manner.}

\section{Overview}{

Beginning with version 3.0, \code{\link{DiffBind}} introduces 
substantial updates and new features that may cause scripts
written for earlier versions to function differently 
(or not at all), as well as altering the results.
This page givens details on these changes,
and how to approximate results computed
with earlier version if desired.

The major change in version 3.0 is in how the data are modeled.
In previous versions, a separate model was derived for each contrast,
including data only for those samples present in the contrast.
Model design options were implicit and limited to either a single
factor, or a subset of two-factor "blocked" designs.
Starting in version 3.0, the default mode is to 
include all the data in a single model,
allowing for any allowable design formula and any set of
allowable contrasts.

Another change starting from version 3.0 is in how normalization is done.
There are more normalization options, and more explicit control over them.
The default normalization options have also changed, 
so reproducing a pre-3.0 analysis requires that normalization 
parameters to be specified.

It is recommended that existing analyses be re-run with the current software.
Existing scripts should execute (with the exception of two
normalization parameters which have been moved from \code{\link{dba.analyze}}
to the new interface function \code{\link{dba.normalize}}.)

See the \code{\link{DiffBind}} vignette for more information on 
processing and analyzing ChIP-seq (and ATAC-seq) experiments.

}

\section{Changes to Defaults}{
\itemize{
 \item{\bold{blacklist} is applied by default, if available, using automatic 
detection of reference genome.}
 \item{\bold{greylists} are generated from controls and applied by default.}
 \item{\bold{minimum read counts} are now 0 instead of being rounded up to 1
 (this is now controllable).}
 \item{\bold{centering peaks around summits} is now done by default using 401-bp 
 wide peaks (recommend to use `summits=100` for ATAC-seq).}
 \item{\bold{read counting} is now performed by `summarizeOverlaps()` by default, 
 with single-end/paired-end counting automatically detected.}
 \item{\bold{filtering} is performed by default; consensus peaks where no peak has
 and RPKM value of at least 1 in any sample are filtered.}
 \item{\bold{control read subtraction} is now turned off by default if a greylist 
 is present}
 \item{\bold{normalization} is based on full library sizes by default for both 
 `edgeR` and `DESeq2`analyses.}
 \item{\bold{score} is set to normalized values by default.}
 }
}

\section{Backward compatibility}{

Most existing \code{\link{DiffBind}} scripts and saved objects will
run correctly using version 3.0, but there may be differences
in the results. 

This section describes how to approximate earlier
results for existing scripts and objects.

\subsection{Running with saved DBA objects:}{
If a \code{DBA} object was created with an earlier version of 
\code{\link{DiffBind}}, and saved using the \code{\link{dba.save}}
function, and loaded using the \code{\link{dba.load}} function,
all settings should be preserved, such that running the analysis anew
will yield the same results.

In order to re-run the analysis using the post-version 3.0 settings,
the original script should be used to re-create the \code{DBA} object.
}

\subsection{Re-running DiffBind scripts:}{
By default, if you re-run a \code{\link{DiffBind}} script, it will use the
new defaults from version 3.0 and beyond.
In order to re-analyze an experiment in the pre-version 3.0 mode,
a number of defaults need to be changed.

When calling \code{\link{dba.count}}, the following defaults are changed:
\itemize{
  \item{\code{summits}:} {This parameter is now set by default. Setting 
     \code{summits=FALSE} will preempt re-centering each 
      peak interval around its point of highest pileup.}  
  \item{\code{filter}:}  {The new default for this parameter is \code{1}
  and is based on \code{RPKM} values; previously
       it was set to \code{filter=0} and was based on read counts.}
  \item{\code{minCount}:} {This is a new parameter representing a minimum read count value.
                   It now default to \code{0}; 
                   to get the previous behavior, set \code{minCount=1}. }
}

The easiest way to perform subsequent processing in a pre-version 3.0 
manner is to set a configuration option:

\code{DBA$config$design <- FALSE}

This will result in the appropriate defaults being set for
the new interface function, \code{\link{dba.normalize}} (which
does not need to be invoked explicitly.)
The pre-version 3.0 settings for \code{\link{dba.normalize}} parameters
are as follows:
\itemize{
  \item{\code{normalize}:} {\code{DBA_NORM_DEFAULT}}
  \item{\code{library}:} {\code{DBA_LIBSIZE_FULL}}
  \item{\code{background}:} {\code{FALSE}}
}
Note that two parameters that used to be available when calling 
\code{\link{dba.analyze}} have been moved:

\itemize{
  \item{\code{bSubControl}:} {now integrated into  \code{\link{dba.count}}. \code{FALSE} by default (unless a greylist has
  been added using \code{\link{dba.blacklist}}).}
  \item{\code{bFullLibrarySize}:} {now integrated into  
  \code{\link{dba.normalize}} as an option for  the \code{library} parameter. 
  \code{library=DBA_LIBSIZE_FULL} is equivalent to \code{bFullLibrarySize=TRUE},
  and \code{library=DBA_LIBSIZE_PEAKREADS} is equivalent to 
  \code{bFullLibrarySize=FALSE}.}
}
}
}
\author{
Rory Stark
}

\seealso{
The \code{\link{DiffBind}} vignette has been updated to show how to
analyze experiments using version 3.0.
}

