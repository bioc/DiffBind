library(DiffBind)

crwd <- getwd()
setwd("holding/BrundleData/inst/extdata/")

############# Drosophila spike-ins #############
samples <- read.csv("samplesheet/samplesheet_SLX8047_hs.csv")
dmsamples <- read.csv("samplesheet/samplesheet_SLX8047_dm.csv")

samples$Spikein <- dmsamples$bamReads
spikes <- dba(sampleSheet = samples)
spikes$config$doBlacklist <- spikes$config$doGreylist <- FALSE
spikes <- dba.count(spikes)
spikes <- dba.normalize(spikes, spikein = TRUE)
spikes.spikeins <- spikes$norm$background
spikes$norm <- NULL
spikes <- dba.normalize(spikes, normalize="RLE", background=TRUE)
spikes.background <- spikes$norm$background
spikes$config <- tamoxifen$config
spikes$config$doBlacklist <- spikes$config$doGreylist <- FALSE
spikes$config$cores <- NULL
save(spikes,spikes.background,spikes.spikeins,file="spikes.rda")

############# parallel factor with CTCF #############
library(Brundle)
data(dbaExperiment,package="Brundle")
pfac <- dba(dbaExperiment, minOverlap=1)
pfac$class["Spikein",] <- pfac$class["bamRead",]
ctcf <- GRanges(jg.controlPeakset[,1:3])
parallelFactor.peaks <- ctcf
pfac <- dba.normalize(pfac, normalize="RLE", background=TRUE)
parallelFactor.background <- pfac$norm$background
pfac$norm <- NULL
pfac <- dba.normalize(pfac,spikein = ctcf)
parallelFactor.ctcf <- pfac$norm$background
pfac <- dba.contrast(pfac,reorderMeta = list(Condition="Fulvestrant"))
pfac$design <- NULL
pfac$contrasts <- NULL
pfac$config <- tamoxifen$config
parallelFactor <- pfac
parallelFactor$config$doBlacklist <-
  parallelFactor$config$doGreylist <- FALSE
parallelFactor$config$cores <- NULL
save(parallelFactor, parallelFactor.peaks,
     parallelFactor.background, parallelFactor.ctcf,
     file="parallelFactor.rda")

setwd(crwd)
system("cp holding/BrundleData/inst/extdata/*.rda .")



