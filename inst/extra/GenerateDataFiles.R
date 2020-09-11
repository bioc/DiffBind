##
## Script for generating data files included with DiffBind package
##
## To be run from working directory containing vignette data (including bam files)

library(DiffBind)
library(tools)
NUMCORES <- 18

setwd("~/Work/DiffBind/Vignette/DiffBind_Vignette/")

## Greylists: Gather karyotypes for internally supported reference genomes

genomes <- c("BSgenome.Hsapiens.UCSC.hg19",
             "BSgenome.Hsapiens.UCSC.hg38",
             "BSgenome.Hsapiens.NCBI.GRCh38",
             "BSgenome.Mmusculus.UCSC.mm9",
             "BSgenome.Mmusculus.UCSC.mm10",
             "BSgenome.Celegans.UCSC.ce10",
             "BSgenome.Celegans.UCSC.ce11",
             "BSgenome.Dmelanogaster.UCSC.dm3",
             "BSgenome.Dmelanogaster.UCSC.dm6")

library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Hsapiens.NCBI.GRCh38")
library("BSgenome.Mmusculus.UCSC.mm9")
library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Celegans.UCSC.ce10")
library("BSgenome.Celegans.UCSC.ce11")
library("BSgenome.Dmelanogaster.UCSC.dm3")
library("BSgenome.Dmelanogaster.UCSC.dm6")

dba.ktypes <- list(seqinfo(BSgenome.Hsapiens.UCSC.hg19))
dba.ktypes <- DiffBind:::pv.listadd(dba.ktypes,seqinfo(BSgenome.Hsapiens.UCSC.hg38))
dba.ktypes <- DiffBind:::pv.listadd(dba.ktypes,seqinfo(BSgenome.Hsapiens.NCBI.GRCh38))
dba.ktypes <- DiffBind:::pv.listadd(dba.ktypes,seqinfo(BSgenome.Mmusculus.UCSC.mm9))
dba.ktypes <- DiffBind:::pv.listadd(dba.ktypes,seqinfo(BSgenome.Mmusculus.UCSC.mm10))
dba.ktypes <- DiffBind:::pv.listadd(dba.ktypes,seqinfo(BSgenome.Celegans.UCSC.ce10))
dba.ktypes <- DiffBind:::pv.listadd(dba.ktypes,seqinfo(BSgenome.Celegans.UCSC.ce11))
dba.ktypes <- DiffBind:::pv.listadd(dba.ktypes,seqinfo(BSgenome.Dmelanogaster.UCSC.dm3))
dba.ktypes <- DiffBind:::pv.listadd(dba.ktypes,seqinfo(BSgenome.Dmelanogaster.UCSC.dm6))
names(dba.ktypes) <- genomes

save(dba.ktypes,file="ktypes.rda")

## Load sample sheet to generate peak data
tamoxifen <- dba(sampleSheet = "tamoxifen.csv")
tamoxifen$config$RunParallel <- FALSE
config <- tamoxifen$config
save(tamoxifen,file="tamoxifen_peaks.rda")

## Generate greylist
tamoxifen$config$RunParallel <- TRUE
tamoxifen$config$cores <- NUMCORES
ktype <- dba.ktypes$BSgenome.Hsapiens.UCSC.hg19["chr18"]
tamoxifen$config$greylist.pval = .999
tamoxifen <- dba.blacklist(tamoxifen,
                           blacklist=DBA_BLACKLIST_HG19,greylist=ktype, 
                           cores=3)
tamoxifen.greylist <- dba.blacklist(tamoxifen, Retrieve=DBA_GREYLIST)

save(tamoxifen.greylist,file="tamoxifen_greylist.rda")

tam.bl <- tamoxifen

## Generate count data with background normalization
load("tamoxifen_peaks.rda")
tamoxifen$config$RunParallel <- TRUE
tamoxifen$config$cores <- NUMCORES
tamoxifen <- dba.count(tamoxifen)
tam <- dba.normalize(tamoxifen, background=TRUE)
tamoxifen$norm$background <- tam$norm$background
tamoxifen$config <- config
save(tamoxifen,file="tamoxifen_counts.rda")

## Genereate analysis data with blacklists/greylists/background normalization
tamoxifen <- tam.bl
tamoxifen$config$RunParallel <- TRUE
tamoxifen$config$cores <- NUMCORES
tamoxifen <- dba.count(tamoxifen)
tam <- dba.normalize(tamoxifen)
tamoxifen <- dba.contrast(tamoxifen,design="~Tissue + Condition")
tamoxifen <- dba.contrast(tamoxifen)
tamoxifen <- dba.analyze(tamoxifen)
tamoxifen$config <- config
save(tamoxifen,file="tamoxifen_analysis.rda")

## Compress files
resaveRdaFiles(".", compress="auto",compression_level = 9)

# cp tamoxifen_*.rda DiffBind/data/
# cp ktypes.rda DiffBind/install/extra