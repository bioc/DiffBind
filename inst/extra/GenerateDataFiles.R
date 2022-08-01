##
## Script for generating data files included with DiffBind package
##
## To be run from working directory containing vignette data (including bam files)

# setwd("~/Work/DiffBind/Vignette/DiffBind_Vignette/")

library(DiffBind)
library(tools)
NUMCORES <- 18

GENERATE_KARYOTYPES <- FALSE
GENERATE_ANALYSIS   <- TRUE
GENERATE_SPIKEINS   <- TRUE


if(GENERATE_KARYOTYPES) {
  lib <- "/data/personal/stark01/rlibs"
  alllibs <- .libPaths()
  .libPaths(c(lib,alllibs))
  alllibs <- .libPaths()
  library(BSgenome)
  
  genomes <- c("BSgenome.Hsapiens.UCSC.hg19",
               "BSgenome.Hsapiens.UCSC.hg38",
               "BSgenome.Hsapiens.NCBI.GRCh38",
               "BSgenome.Mmusculus.UCSC.mm9",
               "BSgenome.Mmusculus.UCSC.mm10",
               "BSgenome.Celegans.UCSC.ce10",
               "BSgenome.Celegans.UCSC.ce11",
               "BSgenome.Dmelanogaster.UCSC.dm3",
               "BSgenome.Dmelanogaster.UCSC.dm6")
  
  dba.ktypes <- NULL
  for(genome in genomes) {
    bsgenes <- BSgenome::getBSgenome(genome)
    if(is.null(dba.ktypes)) {
      dba.ktypes <- list(seqinfo(bsgenes))
    } else {
      dba.ktypes <- DiffBind:::pv.listadd(dba.ktypes,seqinfo(bsgenes))
    }
  }
  
  installed <- BSgenome::installed.genomes()
  othergenomes <- BSgenome::available.genomes(splitNameParts = TRUE)
  othergenomes <- othergenomes[!othergenomes[,5],1]
  othergenomes <- othergenomes[-match(genomes,othergenomes)]
  othergenomes <- othergenomes[othergenomes %in% installed]
  for(genome in othergenomes) {
    bsgenes <- BSgenome::getBSgenome(genome)
    dba.ktypes <- DiffBind:::pv.listadd(dba.ktypes,seqinfo(bsgenes))
  }
  names(dba.ktypes) <- c(genomes, othergenomes)
  save(dba.ktypes,file="ktypes.rda")
}

if(GENERATE_ANALYSIS) {
  
  load("ktypes.rda")
  
  ## Load sample sheet to generate peak data
  tamoxifen <- dba(sampleSheet = "tamoxifen.csv")
  tamoxifen$config$RunParallel <- FALSE
  tamoxifen$config$cores <- NULL
  config <- tamoxifen$config
  tam.bl <- tamoxifen
  
  ## Generate greylist
  tamoxifen$config$RunParallel <- TRUE
  tamoxifen$config$cores <- NUMCORES
  ktype <- dba.ktypes$BSgenome.Hsapiens.UCSC.hg19["chr18"]
  tamoxifen$config$greylist.pval = .999
  tamoxifen <- dba.blacklist(tamoxifen,
                             blacklist=DBA_BLACKLIST_HG19,greylist=ktype, 
                             cores=NUMCORES)
  tamoxifen.greylist <- dba.blacklist(tamoxifen, Retrieve=DBA_GREYLIST)
  
  save(tamoxifen.greylist,file="tamoxifen_greylist.rda")
  
  tam.bl$blacklist <- tamoxifen$blacklist
  tam.bl$greylist  <- tamoxifen$greylist
  tamoxifen <- tam.bl
  tamoxifen$config <- config
  save(tamoxifen,file="tamoxifen_peaks.rda")
  
  ## TEST: Counting using separate consensus peaks and different summits
  tamoxifen$config$RunParallel <- TRUE
  tamoxifen$config$cores <- NUMCORES
  tamoxifen_consensus <- dba.peakset(tamoxifen,
                                     consensus = c(DBA_TISSUE, DBA_CONDITION),
                                     minOverlap = 0.66)
  
  tamoxifen_consensus <- dba(tamoxifen_consensus,
                             mask = tamoxifen_consensus$masks$Consensus,
                             minOverlap = 1)
  
  consensus_peaks <- dba.peakset(tamoxifen_consensus, bRetrieve = TRUE)
  
  tam1 <- dba.count(tamoxifen, summits=TRUE, peaks=consensus_peaks)
  tam2 <- dba.count(tam1, summits=100)
  tam3 <- dba.count(tam1, summits=1000)
  if(nrow(tam3$binding) >= nrow(tam2$binding)) {
    stop("Problem counts/centering peaks")
  }
  gc()
  
  ## Generate count data with background normalization
  tamoxifen$config$RunParallel <- TRUE
  tamoxifen$config$cores <- NUMCORES
  tamoxifen <- dba.count(tamoxifen)
  tam <- dba.normalize(tamoxifen, background=TRUE)
  tamoxifen$norm$background <- tam$norm$background
  tam  <- dba.contrast(tamoxifen, reorderMeta = list(Condition="Responsive"))
  tamoxifen$meta <- tam$meta
  tamoxifen$config <- config
  save(tamoxifen,file="tamoxifen_counts.rda")
  
  ## Generate analysis data with blacklists/greylists/background normalization
  tamoxifen <- tam.bl
  tamoxifen$config$RunParallel <- TRUE
  tamoxifen$config$cores <- NUMCORES
  tamoxifen <- dba.count(tamoxifen)
  tamoxifen$norm$background <- tam$norm$background
  tamoxifen <- dba.contrast(tamoxifen,design="~Tissue + Condition",
                            reorderMeta = list(Condition="Responsive"))
  tamoxifen <- dba.contrast(tamoxifen)
  tamoxifen <- dba.analyze(tamoxifen)
  tamoxifen$config <- config
  save(tamoxifen,file="tamoxifen_analysis.rda")
}

if(GENERATE_SPIKEINS) {
  ## spikein and parallel factors ##
  load("tamoxifen_peaks.rda")
  source("GenerateSpikein.R")
}

## Compress files
resaveRdaFiles(".", compress="auto",compression_level = 9)

# cp tamoxifen_*.rda DiffBind/data/
# cp ktypes.rda spikes.rda parallelFactor.rda DiffBind/inst/extra