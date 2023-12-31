\name{DiffBindNews}
\title{DiffBind News}
\encoding{UTF-8}

\section{Version 3.9.4}{
  \itemize{
    \item{Rollover bugfixes}
    \item{add option to return ggplot in dba.plotVolcano}
  }
}

\section{Version 3.8}{

Mostly bug fixes
   \itemize{
     \item{Rollover bugfixes from Release branch}
     \item{Re-normalize after applying black/greylists}
     \item{Apply Rhtslib patch}
     \item{Update vignette}
     \item{Remove contact info for Gord}
  }
}

\section{Version 3.6}{

Mostly bug fixes, defaults for spike-in normalization.

\itemize{
 \item{Change default spike-in normalization to NATIVE (TLE/TMM)}
 \item{Change spike-in normalization to use reference library sizes}
 \item{Fix bug involving Called columns in reports}
 \item{Maintain FRiP calculations}
 \item{Improve dba.plotProfile sample labels}
 \item{Fix issue with package:parallel}
 \item{Fix error/warning in dba.blacklist relating to non-matching 
 chromosome names}
 \item{Update man page for dba.report to clarify how fold changes are 
 reported when bNormalized=FALSE.}
 \item{Add some new test conditions to GenerateDataFiles.R}
}
}

\section{Version 3.4}{
\itemize{
 \item{Fix bFlip issues}
 \item{Fix bug where mode not passed to summarizeOverlaps}
 \item{Add $inter.feature config parameter}
 \item{Re-compute FDR when fold change other than 0 is specified}
 \item{Remove most gc() calls for performance}
}
}


\section{Version 3.2}{
\itemize{
 \item{New type of plot: dba.plotProfile()}
 \item{Can mix single-end and paired-end bam files}
 \item{Various bug fixes}
}
}

\section{Version 3.0}{
\itemize{
\item{
This is a major release of DiffBind and the new version will be 3.0.

The main upgrade involves how the modelling and normalization are done. 
DiffBind now supports models and contrasts of arbitrary complexity 
using either/both DESeq2 or/and edgeR,
as well as a myriad of normalization options.

NB:

The previous methods for modelling are maintained for backward compatibility, 
however they are not the default. 
To repeat earlier analyses, dba.contrast() must be 
called explicitly with design=FALSE.
See ?DiffBind3 for more information.

The default mode for dba.count() is now to center around summits 
(resulting in 401bp intervals). 
To to avoid recentering around summits as was the previous default, 
set summits=FALSE (or to a more appropriate value).

Normalization options have been moved from dba.analyze() 
to the new interface function dba.normalize().
Any non-default normalization options must be specified using dba.normalize().
}
 \item{Summary of Changes:}
   \itemize{
       \item{dba.analyze():}
          \itemize{
             \item{Automatic mode can start at any point, including from a 
                   sample sheet, and continue default analysis}
              \item{Remove normalization options bSubControl, bFullLibrarySize,
                    filter, and filterFun from dba.analyze(), 
                    now set in dba.normalize().}
              \item{Add support to analyze using full model design formula.}
              \item{Update DESeq2 analysis method.}
              \item{Update edgeR analysis method.}
              \item{Moved edgeR bTagwise parameter to $config option}
              \item{Remove support for DESeq analysis method.}
              \item{Add ability to retrieve DESeq2 or edgeR object}
          }
          
        \item{dba.contrast():} {
           \itemize{
              \item{Add design parameter to set design formula}
              \item{Add contrast parameter to specify contrast using 
                    variety of methods}
              \item{Add reorderMeta parameter to set factor value order}
              \item{Add bGetCoefficients parameter to get design coefficient
                   names to use in contrasts}
           }
        }
        
        \item{NEW: dba.normalize():}
           \itemize{
              \item{Support TMM, RLE, and Library size noramlization for 
                    both DESeq2 and edgeR}
              \item{Support background bin normalization using csaw}
              \item{Support offsets and loess fit normalization}
              \item{Support spike-in normalization with combined or separate reference genomes}
              \item{Support parallel factor normalization}
              \item{bSubControl default depend on presence of Greylist}
           }
           
        \item{dba.count():}
           \itemize{
             \item{Change default to summits=250; 
                   to avoid recentering around summits,must set summits=FALSE}
             \item{Default for bUseSummarizeOverlaps in dba.count is now TRUE}
             \item{Automatically detect single-end/paired-end in dba.count}
             \item{Automatically index unindexed bams in dba.count and dba.normalize}
             \item{move bSubControl parameter}
             \item{Default score is now new score DBA_SCORE_NORMALIZED}
             \item{Add minCount parameter to dba.count(), default now 0 instead of 1}
             \item{Filtering peak by read count thresholds only available in dba.count()}
             \item{Fix bug in dba.count() with user-supplied peakset and summits=n}
           }
          
        \item{NEW: dba.blacklist():}
           \itemize{
             \item{Apply ENCODE blacklist}
             \item{Automatically detect reference genome for blacklist}
             \item{Apply Greylists}
             \item{Generate Greylists from controls using GreyListChIP package}
           }

        \item{Plotting changes:}
          \itemize{
            \item{Add loess fit line to dba.plotMA()}
            \item{Add ability in dba.plotMA() to plot aribitrary samples 
                  (without contrast).}
            \item{Add mask parameter to dba.plotBox() } 
            \item{Support negative scores, eg Fold changes in report-based objects,
                  to enable fold-change heatmaps.}
            \item{Removed bCorPlot as a parameter to dba(), dba.count(), 
                  and dba.analyze(). Use config.}
          }
          
        \item{dba.show() / print changes:}
           \itemize{
             \item{Updated dba.show() and print() to deal with designs and 
                   different contrast types}
             \item{Add ability to retrieve design formula in dba.show()}
             \item{Removed bUsePval parameter in dba.show()}
             \item{Added constant variable DBA_READS to access library sizes}
           }
        
        \item{Vignette and help pages:} 
           \itemize{
              \item{Replace multi-factor analysis section}
              \item{Add extensive normalization section}
              \item{Add blacklist/greylist section.}
              \item{Add pike-in and parallel normalization examples}
              \item{Add DiffBind3 help page and vignette section with 
                    information on backward compatibility.}
              \item{Update technical details sections}
              \item{General updates to all sections}
              \item{Add GenerateDataFiles.R to package}
           }
        \item{Various bugfixes and cosmetic changes.}
    }
 }
}


\section{Version 2.15.1}{
   \itemize{
      \item{Fix issue when extracting SummarizedExperiment from DBA}
    }
}


\section{Version 2.13.1}{
   \itemize{
      \item{Change dba.plotPCA to use proper loading}
    }
}

\section{Version 2.8.0}{
Minor features, documentation updates, bug fixes:

  \itemize{
    \item{Features}    
      \itemize{
        \item{dba.report: add precision option to dba.report}
        \item{dba.plot*: make plots use same precision as reports when thresholding}
        \item{dba: add dir option to dba}
      }
      
    \item{Documentation updates}
      \itemize{
        \item{Vignette: change default to bFullLibrarySize=TRUE in 
              description of DESeq2 analysis}
         \item{Vignette: update vignette to not change dir}
         \item{dba.report: clean up description of bCalled in man page}
         \item{dba.report: modify example inman page to be clearer}
         \item{dba: update man page to not change dir}
         \item{dba.save: dontrun example code for dba.save writing into LIB}
         \item{dba: dontrun example code for dba setting wd to LIB}
      }
      
    \item{Bug fixes}
      \itemize{
        \item{dba.report: Sort report by FDR instead of p_value}
        \item{dba.peakset: fix bug when adding consensus peaks with chromosomes 
              not present in some peaksets}
        \item{dba.count: fix bug when recentering passed-in peaks}
        \item{dba.report: fix bug when using filtering from dba.analyze}
        \item{dba.count: fix bug when passing in peaks using factors}
        \item{dba.count: fix bug caused by not registering one of the C
              routines correctly.}
      }
   }
}

\section{Version 2.6.0}{
Changes since last release include bug fixes and some minor feature refinements.

   \itemize{
      \item{Feature changes}
         \itemize{
            \item{Change sortFun parameter in dba.plotHeatmap to default sd, with
                  FALSE as option for no sorting}
             \item{Sort peaks when adding directly via dba.peakset}
             \item{don't add _ if no initString in dba.report}
             \item{Internal feature: alternate peak counts: pv.resetCounts}
         }
      \item{Bug Fixes}
         \itemize{
             \item{Fix bug in dba.plotHeatmap is all values in a row are zero}
             \item{Change authors in vignette to conform to new standard}
             \item{Fix single peak boundary conditions}
             \item{Subset config$fragmentSize when masking}
             \item{Update example peaks to match data objects}
             \item{Fix bug in dba.plotHeatmap if all values in a row are zero}
             \item{Bugfix when returning report as GRanges with only one site}
             \item{Bugfix when plotting venns of consensus peaksets}
             \item{Fix buffer overrun causing segfault on MacOS}
         }
   }
}

\section{Version 2.4.0}{

Changes since least release include bug fixes and one new feature:

  \itemize{
     \item{Feature: add new plot - dba.plotVolcano}
   }
}

\section{Version 2.2.0}{

Changes since least release include bug fixes and some new features:

  \itemize{
     \item{Feature: Control which principal components are plotted using components parameter in dba.plotPCA}
     \item{Feature: Control axis range using xrange and yrange parameters in dba.plotMA}
     \item{Feature: Filtering per-contrast using filter and filterFun parameters in dba.analyze)}
     \item{Feature: Flip which group in contrast shows gain/loss (sign of fold change) using bFlip parameter in dba.report}
     \item{Feature: Flip which group in contrast shows gain/loss (sign of fold change) using bFlip parameter in dba.plotMA}
   }
}

\section{Version 2.0.0}{
  DiffBind has had substantial internal changes since the last release.
  These are focused on improving memory efficiency for large experiments.
  \itemize{
     \item{Feature changes}
        \itemize{
           \item{Change default analysis method to DESeq2 for its more conservative normalization}
           \item{Designate DESeq method as obsolete (in favor of DESeq2); alter documentation and vignette accordingly.}
           \item{Change default FDR threshold to 0.05}
           \item{Add bNot parameter to dba.contrast to remove ! contrasts by default}
           \item{Remove bReturnPeaksets parameter from dba.plotVenn (does this by default)}
           \item{Change bCorPlot default to FALSE (no more automatic clustering heatmaps)}
        }
     \item{Internal changes}
       \itemize{
          \item{Bump version number to 2.0}
          \item{Update vignette}
          \item{Remove $allvectors and $vectors; replace with $merged (without score matrix) and $binding}
          \item{Upgrade peaksort to use inplace peakOrder}
          \item{Optimize peak merging memory usage}
          \item{Change PCA method from princomp() to prcomp()}
          \item{maxGap implemented in merge}
          \item{Include the beginnings of some unit tests using the testthat package.}
      }
     \item{Bug fixes}
       \itemize{
          \item{Fix bug in retrieving SummarizedExperiment}
          \item{Fix bug when no peaks}
          \item{Fix bugs in non-standard chromosome names and chromosome order}
          \item{Fix bugs in Called logical vectors}
          \item{Ensure loading sample sheet doesn't treat hash as comment char.}
          \item{Tildes in file paths now accepted.}
          \item{Spaces trimmed from entries in sample sheets (with warning).}
          \item{Functions added to importFrom stats to satisfy BiocCheck.}
       }
  }
}

\section{Version 1.16.0}{
\itemize{
   \item{Roll up bugfixes}
   \item{dba.plotHeatmap returns binding sites in row order}
   }
}


\section{Version 1.14.0}{
\itemize{
   \item{Add support for reading Excel-format sample sheets (.xls, .xlsx extensions}
   \item{Update DESeq2 reference in vignette; fix vignette samplesheet}
   \item{use vennPlot from systemPiper}
   \item{Fix Makevars to avoid gnu-specific extensions}
   \item{Replace 'require' with 'requireNamespace' to eliminate NOTEs regarding misuse of 'require'}
   \item{Remove non-ASCII characters from a couple of comments}
   \item{Change Gord's email address}
   \item{New: color vector lists for dba.plotHeatmap and colors for dba.plotPCA labels}
   \item{Fix: bug causing two plots when changing score in dba.plotHeatmap and dba.plotPCA}
 }
}

\section{Version 1.12.0}{
\itemize{
   \item{Mostly bug fixes!}
 }
}

\section{Version 1.10.0}{

Changes since last release:

\itemize{
   \item{Counting}
      \itemize{
         \item{New: option to compute summits}
         \item{New: option to center peaks with fixed width around summits}
         \item{New: scores for summits (height, position) and CPM for TMM values}
         \item{New: filter reads by mapping quality (mapQCth)}
         \item{New: support for PE bam data using summarizeOverlaps}
         \item{Remove: bCalledMask (now always TRUE)}
         \item{Change: insertLength to fragmentSize}
         \item{Add: fragmentSize can be a vector with a size for each sample}
         \item{Change: fragmentSize default is 125 bp}
      }
   \item{Plotting}
      \itemize{
        \item{Change: colors based on CRUK color scheme}
        \item{PCA plots}
           \itemize{
              \item{New: legend}
              \item{New: label parameter for adding text labels of points in 2D plot}
           }
         \item{Venn diagrams}
         \itemize{
            \item{New: plot overlaps of differentially bound sites by specifying contrasts, thresholds etc.}
            \item{New: able to return overlapping peaksets as GRanges directly}
            \item{New: able to generate new DBA object consisting of overlapping peaks}
            \item{New: labelAttributes for controlling default labels}
            \item{New: default main and sub titles}
         }
        \item{Heatmaps}
         \itemize{
            \item{Fix: don’t plot column vector for attributes where every sample has a different value}
         }
      }

   \item{General}
   \itemize{
     \item{New: add attribute value: DBA_ALL_ATTRIBUTES}
     \item{Change: SN  (signal/noise) to FRIP (fraction of reads in peaks)}
     \item{Change: “Down” to “Loss” and Up” to “Gain”}
     \item{Vignette}
       \itemize{
          \item{Change: vignette uses BiocStyles and dynamically generated figures}
           \item{Change: example data based on hg19 instead of hg18}
           \item{Change: example reads from bam files instead of bed files}
            \item{New: section on using DiffBind and ChIPQC together}
     }
     \item{New configuration defaults options (DBA$config):}
        \itemize{
          \item{Metadata name strings: ID, Tissue, Factor, Condition, Treatment, Caller}
          \item{th: significance threshold}
          \item{bUsePval}
          \item{fragmentSize}
          \item{mapQCth: filter reads by mapping quality}
          \item{fragments (for summarizeOverlaps)}
      }
   }
   \item{Bugs/Issues}
   \itemize{
      \item{Fix: bRemoveDuplicates had some unpredictable behaviour}
      \item{Fix: chrN_random were being counted against chrN}
      \item{Disable: tamoxifen_GEO.R doesn’t work after SRA changed format of archived data}
   }
}
}


\section{Version 1.8.0}{

Changes since last release:

\itemize{

     \item{Add support for DESeq2:}{
     \itemize{
      \item{New: Add DBA_DESEQ2, DBA_ALL_METHODS and DBA_ALL_BLOCK method constants}
      \item{Change: dba.analyze can analyze using DESeq2}
      \item{Change: all reporting and plotting functions support DESeq2 results}
      \item{Change: vignette includes comparison of edgeR, DESeq, and DESeq2}
     }
   }

     \item{Changes to counting using dba.count:}{
       \itemize{
        \item{Change: optimize built-in counting code to use much less memory and run faster}
        \item{Change: deprecate bLowMem, replaced with bUseSummarizeOverlaps}
        \item{New: add readFormat parameter to specify read file type (instead of using file suffix)}
       }
     }

   \item{New: generation of result-based DBA object using dba.report (makes it easier to work with differentially bound peaksets)}{}

    \item{Changes to defaults:}{
      \itemize{
        \item{Change: default score is now DBA_SCORE_TMM_MINUS_FULL instead of DBA_SCORE_TMM_MINUS_EFFECTIVE in dba.count}
       \item{Change: default value for bFullLibrarySize is now TRUE in dba.analyze}
        \item{New: add bCorPlot option to DBA$config to turn off CorPlot by default}
      }
    }

   \item{Various bugfixes, improved warnings, updated documentation}{}
  }
}

\section{Version 1.6.0}{
   \itemize{
      \item{New: Low memory counting of bam files using Rsamtools and summarizeOverlaps (bLowMem in dba.count)}
      \item{New: Ability to read in externally derived counts (e.g. from htSeq) (dba.count)}
      \item{Improved: Features to deal with filtering intervals based on read scores (dba.count)}
         \itemize{
            \item{Change parameter name: maxFilter -> filter}
            \item{Allow maxFilter to be a numerical vector to retrieve filtering rate}
            \item{Add parameter: filterFun to control filtering method}
        }
      \item{New: Support for SummarizedExperiment objects (dba and dba.report)}
          \itemize{
            \item{Add bSummarizedExperiment option to dba() to convert DBA object}
            \item{Add DataType = DBA_DATA_SUMMARIZED_REPORT option to dba.report() to return SummarizedExperiment}
         }
      \item{Documentation: Add section to vignette showing how to obtain full tamoxifen resistance dataset}
         \itemize{
            \item{Add section to vignette showing how to obtains full tamoxifen dataset}
            \item{Add script (tamoxifen_GEO.R) and sample sheet (tamoxifen_GEO.csv) to extras for full tamoxifen dataset}
            \item{Add examples to man page for dba.count to show filtering}
            \item{Add examples to man pages for dba and dba.report to show retrieval of SummarizedExperiment objects}
            \item{Update and cleanup vignette and man pages}
         }
      \item{Various bugfixes and improved warnings}
   }
}

\section{Version 1.4.0}{
   \itemize{
      \item{Plotting}
      \itemize{
          \item{dba.plotMA}
          \itemize{
             \item{Smooth plots now default}
             \item{Added fold parameter in addition to th (threshold)}
          }
          \item{dba.plotHeatmap}
          \itemize{
             \item{Side colorbars added}
                \item{Add support for specifying sample mask to include any subset of samples in a contrast plot, including samples that were not in the original contrast}
           }
           \item{dba.plotVenn}
           \itemize{
                \item{Changed plotter from limma to T. Girke's overLapper}
                \item{Added support for 4-way Venns (also in dba/overlap)}
	  }
	  \item{dba.plotPCA}
	  \itemize{
               \item{Add support for specifying sample mask to include any subset of samples in a contrast plot, including samples that were not in the original contrast}
            }
      }
      \item{Peaksets (dba and dba.peakset)}
      \itemize{
         \item{Peakset formats}
         \itemize{
            \item{narrowPeaks format supported}
            \item{Can override file format, score column, and score orientation defaults for supported peak callers}
         }
         \item{Consensus peaksets}
         \itemize{
            \item{Added ability to generate sets of consensus peaksets based on metadata attributes: for example create consensus peaksets for each tissue type and/or condition, or for all unique samples by taking the consensus of their replicate peaksets}
         }
      }
      \item{Read counting (dba.count)}
      \itemize{
         \item{Compute Signal-to-Noise ratio when counting}
         \item{Added bScaleControl to down-scale control reads by default}
         \item{Add option to specify a mask in peak parameter to limit which peaksets are used to for a consensus by overlap. Works with new consensus peakset options in dba.peakset}
         \item{Remove references to support for SAM files}
      }
      \item{Analysis (dba.analyze)}
      \itemize{
         \item{edgeR: updated calls to math change sin edgeR; updated vignette and references}
         \item{DESeq: updated to work with current DESeq; use pooled-CR dispersion estimation method for blocking analysis; update vignette}
      }
      \item{Various bug fixes; more informative warnings; update documentation including vignette, new examples and cross-referencing in man pages}
   }
 }


\section{version 1.2.0:(2012-03-30)}{
 \itemize{
   \item{GRanges is default class for peaksets and reports instead of RangedData, controlled by DataType parameter.}

   \item{Both analysis methods (edgeR and DESeq) use generalized linear models (GLMs) for two-group contrasts by default.}

   \item{Blocking factors (for two-factor analysis) can be specified flexibly such that arbitrary blocking factors can be used. }

   \item{Section added to vignette showing an analysis using a blocking factor.}

   \item{Added new metadata type, DBA_TREATMENT.}

   \item{New DBA_SCORE_ options for specifying scoring method, including TMM normalized counts,
and ability to change scoring method on the fly in dba.plotHeatmap and dba.plotPCA when plotting global binding matrix.}

   \item{bRemoveDuplicates parameter in dba.count allows duplicate reads to be discarded when computing counts}

   \item{More efficient use of memory when analyzing (controlled by bReduceObjects parameter in dba.analyze).}

   \item{various bugs fixed, man pages updated, and warning messages added.}
  }
}

\section{version 1.0.0}{
   \itemize{
     \item{Initial release}
  }
}

