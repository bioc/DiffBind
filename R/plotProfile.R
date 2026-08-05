pv.style          <- "point"
pv.nOfWindows     <- 100
pv.bin_size       <- 20
pv.distanceAround <- 1500
pv.distanceUp     <- 1000
pv.distanceDown   <- 1000


pv.plotProfile <- function(pv, mask, sites, maxSites=1000, labels,
                           scores="Score", absScores=TRUE, annotate=TRUE,
                           normalize=TRUE,merge=DBA_REPLICATE,
                           doPlot=TRUE, returnVal="profiles",
                           ...) {

  ## Profiles are computed natively (BAM coverage binned with
  ## EnrichedHeatmap::normalizeToMatrix) and rendered with
  ## EnrichedHeatmap/ComplexHeatmap, replacing the withdrawn
  ## profileplyr/soGGi backend. See pv.profiles() and pv.profileHeatmap().
  if (!requireNamespace("EnrichedHeatmap", quietly = TRUE)) {
    stop("Package 'EnrichedHeatmap' is required for dba.plotProfile().",
         call. = FALSE)
  }

  if(missing(sites)) {
    sites <- NULL
  }
  
  if(!is(pv,"SummarizedExperiment")) { # Need to compute profiles
    
    if(missing(mask)) {
      mask <- NULL
    }
    
    if(missing(labels)) {
      labels <- NULL
    }
    
    sitelabels <- NULL
    if(is(labels,"list")) {
      if(!is.null(labels$sites)) {
        sitelabels <- labels$sites
      }
      if(!is.null(labels$samples)) {
        labels <- labels$samples
      } else {
        labels <- NULL
      }
    }
    
    groups <- FALSE
    condition1 <- condition2 <- NULL
    
    # Check for default contrast 
    if(is.null(sites)) { 
      if(!is.null(pv$contrasts[[1]])) {
        if(pv$config$AnalysisMethod[1] == DBA_EDGER) {
          if(!is.null(pv$contrasts[[1]]$edgeR)) {
            sites <- 1
          }
        } else {
          if(!is.null(pv$contrasts[[1]]$DESeq2)) {
            sites <- 1
          }
        }
      }
    }
    
    # If sites is unspecified
    if(is.null(sites)) { 
      # Get all sites from binding matrix and scramble
      sites <- pv.peakMatrix_toGR(pv,pv$binding)
      sites <- sites[sample(1:length(sites)),]
      scores <- NULL
    }  else if(is(sites,"numeric")) {
      # Generate Report-based
      con <- sites
      sites <- dba.report(pv, contrast=con,
                          bDB=TRUE,bGain=TRUE,bLoss=TRUE,bAll=FALSE)
      #sites$peaks[[2]]$Score <- abs(sites$peaks[[2]]$Score)
      
      if(is.null(mask)) {
        if(!is.null(pv$contrasts[[con]]$group1)) {
          condition1 <- which(pv$contrasts[[con]]$group1)
          condition2 <- which(pv$contrasts[[con]]$group2)
          name1 <- pv$contrasts[[con]]$name1
          name2 <- pv$contrasts[[con]]$name2
          mask <- list(name1=condition1,name2=condition2)
          names(mask) <- c(name1,name2)
        }
      }
    }
    
    ## Which samples
    samplegroups <- NULL
    if(is.null(mask)) {
      mask <- 1:length(pv$peaks)
    } else if(is(mask,"logical")) {
      mask <- which(mask)
    } else if(is(mask,"list")) {
      samplegroups <- lapply(mask,function(x){if(is(x,"logical")){which(x)}else{x}})
      mask <- unlist(samplegroups)
    }
    
    samples <- pv$class["bamRead",mask]
    
    # Normalization factors
    normfacs <- rep(1, length(mask))
    if(is(normalize,"logical")) {
      if(normalize != FALSE) {
        normfacs <- dba.normalize(pv,bRetrieve=TRUE)$norm.factors[mask]
      }
    } else {
      if(length(normalize) == length(mask)) {
        normfacs <- normalize
      } else {
        stop("normalize must include one factor for each sample.")
      }
    }
    
    ## Sample labels
    if(is.null(labels)) {
      sampnames <- pv$class[DBA_ID,mask]
    } else if(is(labels[1],"character")) {
      sampnames  <- pv$class[DBA_ID,mask]
    } else if(length(labels) == 1) {
      if(labels > 0) {
        sampnames <- pv$class[labels,mask]
      } else {
        sampnames  <- pv$class[DBA_ID,mask]
      }
    } else {
      if(labels[1] > 0) {
        sampnames <-  apply(pv$class[labels,mask], 2, paste,collapse="_")
      } else {
        sampnames  <- pv$class[DBA_ID,mask]
      }
    }
    names(samples) <- sampnames
    
    ## Get sites
    if (is(sites,"vector")) { # specified sites
      # Get specified sites from binding matrix
      sites  <- pv.peakMatrix_toGR(pv,pv$binding[sites,])
      scores <- NULL
    } else if(is(sites, "GRanges")) { # Supplied sites
      # User supplied sites
      # just pass them through
    } else if(is(sites,"GRangesList")) { # Supplied groups of sites
      groups <- TRUE
      groupnames <- names(sites)
    } else if(is(sites,"DBA")) { # Report DBA object
      
      # check if report-based
      if(is.null(sites$resultObject)) {
        stop("sites is not a report-based DBA object.")
      } else if(sites$resultObject==FALSE) {
        stop("sites is not a report-based DBA object.")
      }
      
      # Get site groups from report DBA
      groups <- TRUE
      
      # Get each peakset, form GRangesList
      peaks <- sites$peaks
      for(i in 1:length(peaks)) {
        peaks[[i]] <- GRanges(peaks[[i]])
        #names(mcols(peaks[[i]])) <- "score"
      }
      peaks <- GRangesList(peaks)
      
      # Groupnames from report DBA $ class
      names(peaks) <- groupnames <- pv.groupNames(sites)
      sites <- peaks
    }
    
    if(is(sites,"GRanges")) {
      sites <- GRangesList(sites)
    }
    
    if(is(sites,"GRangesList")){ 
      if(!is.null(sitelabels)) {
        if(length(sitelabels) != length(sites)) {
          warning("Should be same number of site group labels as site groups",
                  .call=FALSE)
        } else {
          names(sites) <- groupnames <- sitelabels
        }
      }
    }
    
    # Limit sites
    for(i in 1:length(sites)) {
      sites[[i]] <- sites[[i]][1:min(maxSites,length(sites[[i]])),]
    }
    
    # save sites in bedfiles
    bedfiles <- NULL
    for(i in 1:length(sites)) {
      bedfile <- tempfile(as.character(Sys.getpid()))
      rtracklayer::export.bed(sites[[i]],con=bedfile)
      bedfiles <- c(bedfiles,bedfile)
    }
    
    # Generate mergelist from attributes if required
    if(!is.null(merge)) {
      merge <- pv.mergelist(pv, mask, merge, labels)
      if(!is.null(samplegroups)) {
        if(length(samplegroups) == length(merge)) {
          if(!is.null(names(samplegroups))) {
            names(merge) <- names(samplegroups)
          }
        }
      }
      sampnames <- names(merge)
    } 
    
    ## Generate profiles
    profiles <- pv.profiles(pv, samples, bedfiles,
                            mergelist = merge, normfacs = normfacs,
                            ...) 
    names(assays(profiles)) <- sampnames
    
    # Add group labels
    if(!groups) {
      groupnames <- names(sites)[1]
      if(is.null(groupnames)) {
        rowData(profiles)$sgGroup <- "Sites"
      } else {
        rowData(profiles)$sgGroup <- groupnames
      }
    } else {
      fns <- strsplit(bedfiles,"/")
      grlabels <- as.character(rowData(profiles)$sgGroup)
      for(i in 1:length(fns)) {
        fn <- fns[[i]]
        grp <- fn[length(fn)]
        grlabels[grlabels==grp] <- groupnames[i]
      }
      rowData(profiles)$"Binding Sites" <- factor(grlabels,
                                                  levels=unique(grlabels),
                                                  labels=groupnames)
      metadata(profiles)$params$rowGroupsInUse <- "Binding Sites"
    }
    
    # Add samplegroup labels
    if(!is.null(samplegroups)) {
      # convert to merged samples 
      if(!is.null(merge)) {
        if(!is.null(merge)) {
          conditions <- NULL
          sgroups <- samplegroups
          sampnum <- 1
          for(sgroup in 1:length(sgroups)) {
            sgroups[[sgroup]] <- sampnum:(sampnum+length(sgroups[[sgroup]])-1)
            sampnum <- sampnum + length(sgroups[[sgroup]])
          }
          for(tomerge in merge) {
            gnum <- which(unlist(lapply(sgroups,function(x){
              if(tomerge[1] %in% x){TRUE}else{FALSE}})))
            conditions <- c(conditions,names(sgroups)[gnum])
          }
          gnames <- NULL
          for(gname in names(samplegroups)) {
            gnames <- c(gnames,rep(gname, sum(conditions %in% gname)))
          }
          conditions <- gnames
        }
      }
      else {
        conditions <- NULL
        for(sgroup in 1:length(samplegroups)) {
          sname <- names(samplegroups)[sgroup]
          conditions <- c(conditions,rep(sname,length(samplegroups[[sgroup]])))
        }
      }
      # add samplegroup labels
      metadata(profiles) <- c(metadata(profiles),
                              list("Sample Group"=factor(conditions, 
                                                         levels=unique(conditions))))
    }
    
    # delete bedfiles
    unlink(bedfiles)
    
  } else {  # Profiles already defined -- just plotting
    profiles <- pv
    sampnames <- names(assays(profiles))
    gnames <- unique(rowData(profiles)$"Binding Sites")
    if(is.null(gnames)) {
      gnames <- unique(rowData(profiles)$sgGroup)
      groups <- FALSE
    } else {
      groups <- TRUE
    }
    groupnames <- as.character(gnames)
    
    # Normalization factors
    normfacs <- rep(1, length(sampnames))
    if(!is(normalize,"logical")) {
      if(length(normalize) == length(sampnames)) {
        normfacs <- normalize
        if(length(normfacs) != length(assays(profiles))) {
          warning("Wrong number of normalization factors, skipping.",call.=FALSE)
        } else {
          for(i in 1:length(normfacs)) {
            assay(profiles,i) <-
              assay(profiles,i) / normfacs[i]
          }
        }
      } else {
        stop("normalize must include one factor for each sample.")
      }
    }
  }
  
  
  # Scores
  profiles <- pv.siteScore(profiles, sites, scores, absScores)
  
  # Annotation
  profiles <- pv.annotate(pv, profiles, annotate)
  
  # Set up Groups
  if(groups) {
    profiles <- pv.groupProfiles(profiles)
  }
  
  # Generate plot
  if(doPlot) {
    if(returnVal=="HeatmapList") {
      return_ht_list <- TRUE
    } else {
      return_ht_list <- FALSE
    }
    doAnnotation <- annotate != FALSE
    profilehm <- pv.profileHeatmap(profiles, 
                                   samples_names=sampnames,
                                   group_names=groupnames,
                                   annotate=doAnnotation,
                                   return_ht_list,
                                   ...)
  }
  
  #return
  
  if(returnVal == "profiles") {
    return(profiles)
  }

  # "EnrichedHeatmap" returns the drawn heatmap; "HeatmapList" returns the
  # undrawn HeatmapList (see the return_ht_list flag set above).
  if(returnVal %in% c("EnrichedHeatmap", "HeatmapList")) {
    return(profilehm)
  }

}

pv.profiles <- function(pv, samples, sites, mergelist=NULL, normfacs=NULL,
                        style=pv.style, nOfWindows=pv.nOfWindows, 
                        bin_size=pv.bin_size, 
                        distanceAround=pv.distanceAround, 
                        distanceUp=pv.distanceUp, distanceDown=pv.distanceDown,
                        ...) {
  
  # check config
  if(!is.null(pv$config$pp.style)) {
    style <- pv$config$pp.style
  }
  if(!is.null(pv$config$pp.nOfWindows)) {
    nOfWindows <- pv$config$pp.nOfWindows
  }
  if(!is.null(pv$config$pp.bin_size)) {
    bin_size <- pv$config$pp.bin_size
  }
  if(!is.null(pv$config$pp.distanceAround)) {
    distanceAround <- pv$config$pp.distanceAround
  }
  if(!is.null(pv$config$pp.distanceUp)) {
    distanceUp <- pv$config$pp.distanceUp
  }
  if(!is.null(pv$config$pp.distanceDown)) {
    distanceDown <- pv$config$pp.distanceDown
  }
  
  args <- pv.sepProfilingArgs(list(...), remove=FALSE)
  if(length(args) > 0) {
    if(!is.null(args$nOfWindows)) {
      nOfWindows <- args$nOfWindows
    }
    if(!is.null(args$bin_size)) {
      bin_size <- args$bin_size
    }
    if(!is.null(args$distanceAround)) {
      distanceAround <- args$distanceAround
    }
    if(!is.null(args$distanceUp)) {
      distanceDown <- args$distanceDown
    }
  }
  
  # Check bam files: exist, PE/SE, BAI

  if(is.null(pv$config$singleEnd)) {
    bfile <- pv.BamFile(samples[1], bIndex=TRUE)
    pv$config$singleEnd	<- !suppressMessages(
      Rsamtools::testPairedEndBam(bfile))
  }
  paired <- !pv$config$singleEnd

  # Fragment length for single-end read extension
  fragment <- pv$config$fragmentSize
  if(is.null(fragment) || !is.numeric(fragment) || any(fragment <= 0)) {
    fragment <- 125
  }
  fragment <- as.integer(fragment[1])

  # Get MulticoreParam

  if(pv$config$RunParallel) {
    if(is.null(pv$config$cores)) {
      cores <- BiocParallel::multicoreWorkers()
    } else {
      cores <- pv$config$cores
    }
    param <- BiocParallel::MulticoreParam(workers=cores)
  } else {
    param <-  BiocParallel::SerialParam()
  }

  # Assemble the target sites from the per-group bed files, tracking which
  # group (bed file) each site came from so downstream group labelling works.
  sitegr <- NULL
  for(i in seq_along(sites)) {
    gr <- rtracklayer::import.bed(sites[i])
    gr$sgGroup <- basename(sites[i])
    sitegr <- if(is.null(sitegr)) gr else c(sitegr, gr)
  }
  if(is.null(sitegr$name)) {
    sitegr$name <- as.character(seq_along(sitegr))
  }

  # Profiling target: a point (peak centre) or the region body
  if(style == "point") {
    target <- GenomicRanges::resize(sitegr, width=1, fix="center")
  } else {
    target <- sitegr
  }

  # Compute one binned signal matrix per sample using native BAM coverage +
  # EnrichedHeatmap::normalizeToMatrix (replaces profileplyr/soGGi).
  message("Generating profiles...")
  oneMatrix <- function(bamfile) {
    reads <- pv.profileReads(bamfile, paired, fragment)
    if(style == "point") {
      # 'point' style: fixed-size bins extending distanceAround bp up- and
      # down-stream of each peak centre (default 1500 -> 3000bp / 150 bins).
      EnrichedHeatmap::normalizeToMatrix(reads, target,
                                         extend = distanceAround,
                                         w = bin_size, mean_mode = "coverage",
                                         background = 0, smooth = TRUE)
    } else {
      # 'percentOfRegion' style: the region body is divided into nOfWindows
      # bins, and the flanks extend distanceAround percent of the region width
      # on each side, binned at that same resolution.  For the documented
      # example (400bp peaks, nOfWindows=20, distanceAround=300): 20bp bins,
      # flanks of 1200bp, giving 60 + 20 + 60 = 140 bins across 2800bp.
      # target_ratio must describe the real geometry, otherwise the body is
      # sampled at a different resolution than the flanks.
      medw <- stats::median(GenomicRanges::width(target))
      ext  <- round(medw * distanceAround / 100)
      EnrichedHeatmap::normalizeToMatrix(reads, target,
                                         extend = ext,
                                         w = max(1, round(medw / nOfWindows)),
                                         k = nOfWindows,
                                         target_ratio = medw / (medw + 2*ext),
                                         mean_mode = "coverage",
                                         background = 0, smooth = TRUE)
    }
  }
  res <- tryCatch(
    suppressMessages(
      matrices <- BiocParallel::bplapply(samples, oneMatrix, BPPARAM = param)),
    error = function(x){stop("Error generating profiles: ",
                             conditionMessage(x), call.=FALSE)}
  )

  # Wrap each sample as a one-assay SummarizedExperiment so the existing
  # merge/normalize machinery (which operates via assay()) works unchanged.
  profiles <- lapply(matrices, function(m)
    SummarizedExperiment::SummarizedExperiment(assays = list(m),
                                               rowRanges = sitegr))

  # Normalize (divide each sample's signal by its normalization factor)
  if(!is.null(normfacs)) {
    if(length(normfacs) != length(profiles)) {
      warning("Wrong number of normalization factors, skipping.",call.=FALSE)
    } else {
      for(i in seq_along(normfacs)) {
        assay(profiles[[i]],1) <- assay(profiles[[i]],1) / normfacs[i]
      }
    }
  }

  # Merge replicate samples (means of the signal matrices)
  if(!is.null(mergelist)) {
    profiles <- pv.mergeProfiles(profiles, mergelist)
  }

  # Combine per-sample matrices into a single multi-assay object
  matlist <- lapply(profiles, function(x) assay(x,1))
  combined <- SummarizedExperiment::SummarizedExperiment(assays = matlist,
                                                         rowRanges = sitegr)
  rowData(combined)$name    <- sitegr$name
  rowData(combined)$sgGroup <- sitegr$sgGroup
  S4Vectors::metadata(combined)$params <-
    list(style = style, nOfWindows = nOfWindows, bin_size = bin_size,
         distanceAround = distanceAround,
         distanceUp = distanceUp, distanceDown = distanceDown)

  return(combined)
}

# Read a BAM into fragment-level GRanges: paired-end read pairs give the
# fragment span; single-end reads are extended to the fragment length.
pv.profileReads <- function(bamfile, paired, fragment) {
  if(paired) {
    ga <- GenomicAlignments::readGAlignmentPairs(bamfile)
    reads <- GenomicRanges::granges(ga)
  } else {
    ga <- GenomicAlignments::readGAlignments(bamfile)
    reads <- GenomicRanges::resize(GenomicRanges::granges(ga), width = fragment)
  }
  suppressWarnings(GenomicRanges::trim(reads))
}

pv.siteScore <- function(profiles, sites, scores, absScores) {

  if(!is.null(sites)) { # Add per-site metadata (e.g. fold-change scores)
    psites <- rowData(profiles)$name
    if(is(sites,"GRangesList")) {
      names(sites) <- NULL
    }
    sites  <- unlist(sites)
    ssites <- names(sites)

    matches <- match(psites, ssites)
    if(sum(is.na(matches)) > 0) {
      rowData(profiles) <- cbind(rowData(profiles), mcols(sites))
    } else {
      rowData(profiles) <- cbind(rowData(profiles),
                                 mcols(sites[matches,]))
    }
  }

  metadata(profiles)$params$mcolToOrderBy <- NULL
  if(!is.null(scores)) {
    if(!scores %in% names(rowData(profiles))) {
      message(scores," not a valid score column, using mean signal.")
    } else {
      metadata(profiles)$params$mcolToOrderBy <- "score"
      whichscore <- match(scores, names(rowData(profiles)))
      rowData(profiles)$score <- rowData(profiles)[[whichscore]]
      if(absScores) {
        rowData(profiles)$score <- abs(rowData(profiles)$score)
      }
    }
  }

  # Order sites (rows) by the chosen score, else by mean signal. Reordering
  # the SummarizedExperiment reorders rowRanges, rowData and every assay
  # matrix consistently (replaces profileplyr::orderBy()).
  if(!is.null(metadata(profiles)$params$mcolToOrderBy)) {
    ord <- order(rowData(profiles)$score, decreasing=TRUE)
  } else {
    ord <- order(rowMeans(assay(profiles,1)), decreasing=TRUE)
  }
  profiles <- profiles[ord, ]

  return(profiles)
}

pv.annotate <- function(pv, profiles, annotate) {

  ## Genomic-feature annotation of the profiled sites was previously provided
  ## by profileplyr::annotateRanges() (via ChIPseeker). That optional
  ## annotation strip is not yet re-implemented in the native backend, so the
  ## profiles are returned without a Features column; pv.profileHeatmap()
  ## detects this and omits the feature annotation. The core profile plot
  ## (per-sample heatmaps, site groups, sample-group colouring) is unaffected.
  if(!is.null(rowData(profiles)$Features)) {
    return(profiles)
  }
  if(!isFALSE(annotate)) {
    message("Site feature annotation is not available in this version; ",
            "plotting profiles without the feature annotation column.")
  }
  return(profiles)
}

pv.setAnno <- function(dataset) {
  Features <-  as.character(rowData(dataset)$annotation_short)
  Features[Features=="3p UTR"] <- "Gene Body"
  Features[Features=="5p UTR"] <- "Gene Body"
  Features[Features=="Downstream"] <- "Gene Body"
  Features[Features=="Exon"] <- "Gene Body"
  Features[Features=="Intron"] <- "Gene Body"
  Features[Features=="Distal Intergenic"] <- "Intergenic"
  Features <- factor(Features,levels=c("Promoter","Gene Body","Intergenic"))
  rowData(dataset)$Features <- Features
  
  Distances <- rowData(dataset)$distanceToTSS
  Distances[Distances < -100000] <- -1000000
  Distances[Distances > 100000]  <- 1000000
  rowData(dataset)$"TSS Distance" <- Distances
  
  return(dataset)
}

pv.groupColors <- c(crukBlue,crukCyan,crukMagenta, crukGrey, 
                    "green", "forestgreen")

pv.conditionColors <- list(c("white",crukMagenta), 
                           c("white",crukBlue), 
                           c("white",crukGrey), 
                           c("white",crukCyan))

pv.profileHeatmap <- function(profiles, samples_names, group_names,
                              annotate=FALSE,
                              ret_ht_list=FALSE, ...){
  message("Plotting...")

  mats  <- assays(profiles)
  nsamp <- length(mats)
  if(is.null(samples_names)) {
    samples_names <- names(mats)
  }
  if(is.null(samples_names)) {
    samples_names <- paste0("Sample", seq_len(nsamp))
  }

  # Row split by binding-site group (e.g. Gain / Loss).  When no analysis has
  # been done there is a single group (named "Sites"), which is still shown as
  # a one-colour bar with its own legend.  The legend is titled after whichever
  # rowData column supplies the grouping, as in the original plots.
  split_name <- "Binding Sites"
  row_split  <- rowData(profiles)$"Binding Sites"
  if(is.null(row_split)) {
    sg <- rowData(profiles)$sgGroup
    if(!is.null(sg)) {
      row_split  <- factor(sg, levels=unique(sg))
      split_name <- "sgGroup"
    }
  }
  # Gain sites are conventionally shown above Loss sites
  if(!is.null(row_split) && all(c("Gain","Loss") %in% levels(row_split))) {
    ord <- c("Gain","Loss", setdiff(levels(row_split), c("Gain","Loss")))
    row_split <- factor(row_split, levels=ord)
  }

  # Site-group colours (e.g. Gain=blue, Loss=cyan) for the left annotation bar
  # and the composite profile lines -- shared across all sample panels.
  group_cols <- NULL
  if(!is.null(row_split)) {
    glev <- levels(row_split)
    group_cols <- pv.groupColors[seq_along(glev)]
    names(group_cols) <- glev
  }

  # Per-sample colouring by contrast condition (heatmap body fill)
  conditions <- S4Vectors::metadata(profiles)$"Sample Group"
  if(!is.null(conditions)) {
    cond_int <- as.integer(factor(conditions, levels=unique(conditions)))
  } else {
    cond_int <- rep(1L, nsamp)
  }
  endcol <- vapply(pv.conditionColors, function(x) x[[2]], character(1))

  # Shared colour scale across panels so samples are comparable
  qmax <- suppressWarnings(
    stats::quantile(unlist(lapply(mats, as.vector)), 0.99, na.rm=TRUE))
  if(!is.finite(qmax) || qmax <= 0) {
    qmax <- max(vapply(mats, function(m) max(m, na.rm=TRUE), numeric(1)))
  }

  # Shared y-axis for the composite profile curves, so between-sample
  # differences in overall signal strength remain visible (rather than each
  # top plot being auto-scaled to its own maximum).
  if(!is.null(row_split)) {
    grp_idx <- split(seq_len(nrow(mats[[1]])), row_split)
  } else {
    grp_idx <- list(seq_len(nrow(mats[[1]])))
  }
  prof_max <- max(vapply(mats, function(m)
    max(vapply(grp_idx, function(idx)
      max(colMeans(m[idx, , drop = FALSE]), na.rm = TRUE), numeric(1))),
    numeric(1)))
  prof_ylim <- c(0, prof_max * 1.05)

  # use_raster is controllable (rasterization needs a working bitmap device)
  dots <- list(...)
  dots <- pv.sepProfilingArgs(dots, remove=TRUE)
  use_raster <- if(!is.null(dots$use_raster)) dots$use_raster else TRUE

  # Rasterizing keeps these heatmaps small, but ComplexHeatmap writes each
  # heatmap body to a temporary PNG using type="cairo", and cairo is not usable
  # in every R installation (on macOS it needs X11/libXrender, frequently
  # absent), which makes the raster step fail with an opaque error.  Choose a
  # bitmap type that works here, and only fall back to un-rasterized output if
  # none does.
  raster_param <- pv.rasterParam()
  if(use_raster && is.null(raster_param)) {
    use_raster   <- FALSE
    raster_param <- list()
  }

  # all_color_scales_equal (default TRUE): every panel shares the one colour
  # scale (qmax) so intensities are directly comparable between samples.  When
  # FALSE, each panel is scaled to its own 99th-percentile signal -- matching
  # the original profileplyr behaviour of independent per-sample colour scales.
  all_equal <- if(!is.null(dots$all_color_scales_equal))
                 isTRUE(dots$all_color_scales_equal) else TRUE

  # EnrichedHeatmap draws a dashed vertical line at the target position by
  # default; the original plots suppressed it (profileplyr's equivalent
  # matrices_pos_line was set FALSE).  Keep it off, but let it be turned on.
  pos_line <- if(!is.null(dots$pos_line)) isTRUE(dots$pos_line) else FALSE

  # Shrink the sample-name column titles if the figure is too narrow to hold
  # them side by side (overridable by passing column_title_gp).
  title_gp <- if(!is.null(dots$column_title_gp)) {
                dots$column_title_gp
              } else {
                grid::gpar(fontsize=pv.titleFontsize(samples_names, nsamp))
              }
  if(all_equal) {
    panel_qmax <- rep(qmax, nsamp)
  } else {
    panel_qmax <- vapply(mats, function(m) {
      q <- suppressWarnings(stats::quantile(as.vector(m), 0.99, na.rm=TRUE))
      if(!is.finite(q) || q <= 0) q <- max(m, na.rm=TRUE)
      q
    }, numeric(1))
  }

  pr <- metadata(profiles)$params
  if(!is.null(pr$style) && pr$style == "point") {
    ext <- if(!is.null(pr$distanceAround)) pr$distanceAround else pr$distanceUp
    axis_name <- c(paste0("-", ext), "0", as.character(ext))
  } else {
    # percentOfRegion: label the flanks by their width in bins, which is
    # distanceAround percent of the nOfWindows bins spanning the region body.
    nflank <- NULL
    if(!is.null(pr$nOfWindows) && !is.null(pr$distanceAround)) {
      nflank <- round(pr$nOfWindows * pr$distanceAround / 100)
    }
    axis_name <- if(!is.null(nflank) && nflank > 0) {
      c(paste0("-", nflank), "start", "end", as.character(nflank))
    } else {
      c("upstream", "start", "end", "downstream")
    }
  }

  # Composite profile lines (top of each panel) coloured by SITE GROUP
  # (Gain/Loss), the same in every panel -- not by sample condition.
  topGp <- if(!is.null(group_cols)) grid::gpar(col = unname(group_cols), lwd = 2)
           else grid::gpar(col = "black", lwd = 2)
  topAnno <- ComplexHeatmap::HeatmapAnnotation(
    enrich = EnrichedHeatmap::anno_enriched(gp = topGp, ylim = prof_ylim,
                                            pos_line = pos_line))

  # Left colour bar identifying the site groups, with its own legend
  # (drawn once, on the first panel).
  leftAnno <- NULL
  if(!is.null(row_split)) {
    # The annotation is keyed by split_name, so build the arguments by name.
    anno_arg <- list(row_split)
    names(anno_arg) <- split_name
    col_arg <- list(group_cols)
    names(col_arg) <- split_name
    leg_arg <- list(list(title = split_name))
    names(leg_arg) <- split_name
    leftAnno <- do.call(ComplexHeatmap::rowAnnotation,
                        c(anno_arg,
                          list(col = col_arg,
                               show_annotation_name = FALSE,
                               annotation_legend_param = leg_arg)))
  }

  htlist <- NULL
  legend_conds_seen <- integer(0)
  for(i in seq_len(nsamp)) {
    col2   <- endcol[[ ((cond_int[i]-1) %% length(endcol)) + 1 ]]
    colfun <- circlize::colorRamp2(c(0, panel_qmax[i]), c("white", col2))
    # Shared scale: draw one legend per condition.  Independent scales: every
    # panel needs its own legend (each has a different maximum).
    if(all_equal) {
      show_leg  <- !(cond_int[i] %in% legend_conds_seen)
      legend_conds_seen <- c(legend_conds_seen, cond_int[i])
      leg_title <- if(!is.null(conditions)) conditions[i] else "Signal"
    } else {
      show_leg  <- TRUE
      leg_title <- samples_names[i]
    }
    ehm <- EnrichedHeatmap::EnrichedHeatmap(
      mats[[i]], name = samples_names[i], column_title = samples_names[i],
      column_title_gp = title_gp,
      col = colfun, axis_name = axis_name, pos_line = pos_line,
      use_raster = use_raster, raster_quality = 5,
      raster_device_param = raster_param,
      top_annotation = topAnno,
      left_annotation = if(i == 1) leftAnno else NULL,
      show_heatmap_legend = show_leg,
      heatmap_legend_param = list(title = as.character(leg_title)))
    htlist <- if(is.null(htlist)) ehm else htlist + ehm
  }

  if(ret_ht_list) {
    return(htlist)
  }

  ComplexHeatmap::draw(htlist, split = row_split,
                       ht_gap = grid::unit(4, "mm"),
                       main_heatmap = 1)
  return(invisible(htlist))
}


## Find a bitmap device type usable for ComplexHeatmap's raster temporary
## images.  ComplexHeatmap asks for type="cairo", which is unavailable in R
## builds without a working cairo (notably macOS installations lacking X11), so
## probe once per session and remember the answer:
##   list()            -- cairo works; leave ComplexHeatmap's default alone
##   list(type=<type>) -- override with a type that does work
##   NULL              -- no usable type; the caller should not rasterize
## Sample names are drawn as column titles at a fixed point size, so several
## samples in a narrow figure (e.g. the 7-inch figures in the plotProfileDemo
## notebook) make the titles run into each other.  Scale the font down to what
## the width available per panel can hold, never going above ComplexHeatmap's
## own default, so wide figures are unaffected.  A character advance of
## 0.6 * fontsize is a deliberately conservative estimate for a proportional
## font -- erring small costs a little legibility, erring large overlaps.
pv.titleFontsize <- function(labels, nsamp, default=13.2, minimum=5) {
  width <- tryCatch(grDevices::dev.size("in")[1], error=function(e) NA_real_)
  if(!is.finite(width) || width <= 0) {
    return(default)
  }
  # leave room for the legends on the right and the site-group bar on the left
  avail <- (width - 1.75) / max(1L, nsamp)
  chars <- max(c(nchar(as.character(labels)), 1L))
  max(minimum, min(default, avail * 72 / (chars * 0.6)))
}

pv.rasterParam <- function() {
  cached <- getOption("DiffBind.rasterDeviceType", NULL)
  if(is.null(cached)) {
    # Some devices (quartz) only write the file once something has been drawn,
    # so the probe has to draw -- which makes it essential to confirm that the
    # probe's own device really opened.  An unusable type (cairo without X11)
    # only warns, leaving the caller's device current: drawing then lands on
    # that device and closing it destroys the caller's plot.  So compare
    # dev.cur() before and after, and never close a device we did not open.
    probe <- function(type) {
      file   <- tempfile(fileext=".png")
      before <- grDevices::dev.cur()
      ok     <- FALSE
      tryCatch({
        suppressWarnings(grDevices::png(file, width=64, height=64, type=type))
        if(!identical(grDevices::dev.cur(), before)) {
          graphics::par(mar=rep(0,4))  # default margins exceed a 64px canvas
          graphics::plot.new()
          grDevices::dev.off()
          ok <- file.exists(file) && file.info(file)$size > 0
        }
      }, error=function(e) NULL)
      # Restore the device that was current on entry.
      guard <- 0L
      while(!identical(grDevices::dev.cur(), before) &&
            grDevices::dev.cur() != 1L && guard < 8L) {
        try(grDevices::dev.off(), silent=TRUE)
        guard <- guard + 1L
      }
      unlink(file)
      isTRUE(ok)
    }
    cached <- "none"
    for(type in c("cairo","quartz","Xlib")) {
      if(probe(type)) {
        cached <- type
        break
      }
    }
    options(DiffBind.rasterDeviceType=cached)
  }
  if(cached == "cairo") {
    return(list())
  }
  if(cached == "none") {
    return(NULL)
  }
  list(type=cached)
}

pv.addArg <- function(addarg, param, val, args=NULL) {
  
  if(!is.null(args)) {
    if(param %in% names(args)) {
      return(addarg)
    }
  }
  
  if(is.null(addarg)) {
    addarg <- list(x=val)
    names(addarg) <- param
  } else {
    addarg <- pv.listadd(addarg, val)
    names(addarg)[length(addarg)] <- param
  }
  return(addarg)
}

pv.groupNames <- function(pv) {
  if(length(unique(pv$class[DBA_ID,])) == ncol(pv$class)) {
    return(pv$class[DBA_ID,])
  }
  
  labels <- NULL
  for(i in c(DBA_ID,DBA_TISSUE, DBA_FACTOR, 
             DBA_CONDITION, DBA_TREATMENT)) {
    if(length(unique(pv$class[i,])) > 1) {
      labels <- paste(labels,pv$class[i,],sep="_")
    }
  }
  
  labels <- sub('.', '',labels)
  
  return(labels)
}

pv.groupProfiles <- function(profiles) {
  return(profiles)
}

pv.mergelist <- function(pv, samples, merge, labels=NULL){
  
  if(!is(merge,"list")) { # Make list from excluded attributes
    attributes <- c(DBA_TISSUE,DBA_FACTOR,DBA_CONDITION,DBA_TREATMENT,DBA_REPLICATE)
    atts <- pv$class[attributes[! attributes %in% merge],samples]
    allmerge <- apply(atts,2,paste,collapse="")
    tomerge  <- unique(allmerge)
    if(is.null(labels)) {
      labels <- -merge
    }
    
    merge <- NULL
    for(val in tomerge) {
      merge <- pv.listadd(merge,which(allmerge %in% val))
    }
  }
  
  # Attach labels to mergelist
  
  names(merge) <- pv.makeMergeLabel(pv, merge, samples,labels)
  
  # check uniqueness of labels
  names(merge) <- make.unique(as.character(names(merge)))
  
  # return list of merge groups with label names
  return(merge)
}

pv.makeMergeLabel <- function(pv, merge, samples, labels=NULL) {
  
  if(is(labels[1],"character")) {
    if(length(labels) != length(merge)) {
      warning("Should be same number of sample labels as samples.",call.=FALSE)
    } else {
      return(labels)
    }
  }
  
  if(is.null(names(merge))) {
    
    atts <- c(DBA_TISSUE,DBA_FACTOR,DBA_CONDITION,DBA_TREATMENT,DBA_REPLICATE)
    
    #samples <- samples[unlist(merge)]
    
    if(is.null(labels)) {
      labels <- atts
    }
    
    if(labels[1] < 0) { # labels are negative, remove
      
      labels <-  atts[which(!atts %in% abs(labels))]
      
      # Strip common attributes as well
      if(length(labels) > 1) {
        uatts <- apply(pv$class[labels,samples],1,unique)
        for(i in length(labels):1) {
          if(length(uatts[[i]])==1) {
            labels <- labels[-i]
          }
        }
      }
    }
    
    res <- lapply(merge, pv.doMakeMergeLabel, pv$class, samples, labels)
    
    return(res)
  }
}

pv.doMakeMergeLabel <- function(tomerge, meta, samples, labels) {
  res <- NULL
  for(label in labels) {
    res <- c(res, unique(meta[label,samples[tomerge]]))
  }
  if(length(res) > 1) {
    res <- paste(res, collapse="_")
  }
  return(res)
}


pv.mergeProfiles <- function(profiles, mergelist, bMean=TRUE) {
  
  totalsamps <- length(profiles)
  mergelist <- mergelist[order(unlist(lapply(mergelist,min)))]
  mergesamps <- sort(unlist(mergelist))
  if(length(mergesamps) != totalsamps) {
    toadd <- which(!(1:totalsamps %in% mergesamps))
  } else toadd <- NULL
  
  resultlist <- NULL
  
  for(tomerge in mergelist) {
    
    if(length(tomerge) > 1) {
      
      adding <- TRUE
      while(adding){
        if(length(toadd) > 0) {
          if(toadd[1] < tomerge[1]) {
            resultlist <- pv.listadd(resultlist, profiles[[toadd[1]]])
            toadd <- toadd[-1]
          } else {
            adding <- FALSE
          }
        } else adding <- FALSE
      }
      
      profile <- profiles[[tomerge[1]]]
      
      merged <- assay(profile)
      for(i in 2:length(tomerge)) {
        merged <- merged + assay(profiles[[tomerge[i]]])
      }
      
      if(bMean) {
        merged <- merged / length(tomerge)
      }
      
      assay(profile) <- merged
      
    } else {
      profile <- profiles[[tomerge[1]]]
    }
    
    resultlist <- pv.listadd(resultlist, profile)
  }
  
  
  while(length(toadd) > 0) {
    resultlist <- pv.listadd(resultlist, profiles[[toadd[1]]])
    toadd <- toadd[-1]
  }
  
  return(resultlist)
  
}

pv.ProfilingArgs <- c("style","nOfWindows","bin_size",
                      "distanceAround","distanceUp","distanceDown")

pv.sepProfilingArgs <- function(arglist, remove=FALSE) {
  profiling <- which(names(arglist) %in% pv.ProfilingArgs)
  if(length(profiling) > 0) {
    proargs <- arglist[profiling]
    plotargs <- arglist[!profiling]
  } else {
    proargs <- NULL
    plotargs <- arglist
  }
  
  if(remove) {
    return(plotargs)
  } else {
    return(proargs)
  }
}

