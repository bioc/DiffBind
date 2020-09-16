library(DiffBind)

crwd <- getwd()
setwd("holding/BrundleData/inst/extdata/")

samples <- read.csv("samplesheet/samplesheet_SLX8047_hs.csv")
dmsamples <- read.csv("samplesheet/samplesheet_SLX8047_dm.csv")

samples$Spikein <- dmsamples$bamReads
spikes <- dba(sampleSheet = samples)
spikes <- dba.count(spikes)

par(mfrow=c(3,2))

dba.plotMA(spikes,contrast=list(Fulvestrant=spikes$masks$Fulvestrant),
           bNormalized=FALSE,sub="RAW" )

spikes <- dba.contrast(spikes,contrast=c("Condition","Fulvestrant","none"),
                       reorderMeta = list(Condition="none"))

spikes$norm <- NULL
spikes <- dba.normalize(spikes, normalize=DBA_NORM_LIB, library="full")
spikes <- dba.analyze(spikes)
dba.plotMA(spikes,bNormalized=TRUE, sub="LIB full")

spikes$norm <- NULL
spikes <- dba.normalize(spikes, normalize="RLE")
spikes <- dba.analyze(spikes)
dba.plotMA(spikes,bNormalized=TRUE, sub="RLE RiP")

spikes$norm <- NULL
spikes <- dba.normalize(spikes, normalize="RLE", background=TRUE)
spikes.background <- spikes$norm$background
spikes <- dba.analyze(spikes)
dba.plotMA(spikes,bNormalized=TRUE, sub="RLE BG")

## drosophila spike-ins
peaks <- read.table("peaks/drosophila/SLX-8047.D705_D506.C81G5ANXX.s_1.r_1.fq_peaks.narrowPeak")
dmchrs <- unique(peaks[,1])

spikes$norm <- NULL
spikes <- dba.normalize(spikes, spikein = TRUE)
spikes.spikeins <- spikes$norm$background

spikes$norm$background <- spikes.background
save(spikes,spikes.background,spikes.spikeins,file="spikes.rda")
spikes$norm$background <- spikes.spikeins

spikes <- dba.analyze(spikes)
dba.plotMA(spikes,bNormalized=TRUE, sub="LIB spikein")

spikes <- dba.normalize(spikes, normalize=DBA_NORM_RLE, spikein = TRUE)
spikes <- dba.analyze(spikes)
dba.plotMA(spikes,bNormalized=TRUE, sub="RLE spikein")

############# parallel factor with CTCF #############
library(Brundle)
data(dbaExperiment,package="Brundle")
pfac <- dba(dbaExperiment, minOverlap =1)
pfac$class["Spikein",] <- pfac$class["bamRead",]
ctcf <- GRanges(jg.controlPeakset[,1:3])
parallelFactor.peaks <- ctcf

par(mfrow=c(3,2))
dba.plotMA(pfac,contrast=list(Fulvestrant=pfac$masks$Fulvestrant),
           bNormalized=FALSE,sub="RAW" )

pfac <- dba.contrast(pfac, contrast=c("Condition","Fulvestrant","none"),
                     reorderMeta = list(Condition="none"))
pfac <- dba.normalize(pfac)
pfac <- dba.analyze(pfac)
dba.plotMA(pfac,bNormalized=TRUE, sub="LIB full")

pfac$norm <- NULL
pfac <- dba.normalize(pfac, normalize="RLE")
pfac <- dba.analyze(pfac)
dba.plotMA(pfac,bNormalized=TRUE, sub="RLE RiP")

pfac$norm <- NULL
pfac <- dba.normalize(pfac, normalize="RLE", background=TRUE)
parallelFactor.background <- pfac$norm$background
pfac <- dba.analyze(pfac)
dba.plotMA(pfac,bNormalized=TRUE, sub="RLE BG")

pfac$norm <- NULL
pfac <- dba.normalize(pfac,spikein = ctcf)
parallelFactor.ctcf <- pfac$norm$background
pfac <- dba.analyze(pfac)
dba.plotMA(pfac,bNormalized=TRUE, sub="LIB CTCF")
parallelFactor <- pfac
save(parallelFactor, parallelFactor.peaks,
     parallelFactor.background, parallelFactor.ctcf,
     file="parallelFactor.rda")

pfac <- dba.normalize(pfac,spikein = ctcf, norm=DBA_NORM_RLE)
pfac <- dba.analyze(pfac)
dba.plotMA(pfac,bNormalized=TRUE, sub="RLE CTCF")

setwd(crwd)
system("cp holding/BrundleData/inst/extdata/*.rda .")



