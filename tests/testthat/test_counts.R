test_that("generating DBA from counts is same as reads",{
  
  data(tamoxifen_counts)
  tamoxifen <- dba.count(tamoxifen,peaks=NULL, score=DBA_SCORE_READS)
  data <- dba.peakset(tamoxifen, bRetrieve=TRUE)
  peaks <- data[,0]
  counts <- mcols(data)
  meta <- dba.show(tamoxifen)
  
  tam <- NULL
  for(i in 1:ncol(counts)) {
    tam <- dba.peakset(tam, peaks=peaks,
                       sampID=meta$ID[i],tissue=meta$Tissue[i],
                       factor=meta$Factor[i], condition=meta$Condition[i],
                       treatment=meta$Treatment[i],replicate=meta$Replicate[i],
                       peak.caller=meta$Caller[i], counts=counts[,i])
  }
  
  tam <- dba.normalize(tam, library="RiP")
  tam <- dba.analyze(tam)
  tam.res <- dba.show(tam,bContrasts = TRUE)
  
  tamoxifen <- dba.normalize(tamoxifen, library="RiP")
  tamoxifen <- dba.analyze(tamoxifen)
  tamoxifen.res <- dba.show(tamoxifen,bContrasts = TRUE)
  
  expect_equal(tam.res$DB.DESeq2,tamoxifen.res$DB.DESeq2)
  
  samples <- tamoxifen$samples
  data <- dba.peakset(tamoxifen, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
  for(i in 1:ncol(counts)) {
    countFile=sprintf("%s.txt",meta$ID[i])
    write.table(data[,c(1:3,3+i)],countFile,
                quote=FALSE,row.names=FALSE,col.names=FALSE)
    samples$Counts[i] <- countFile
  }
  samples$bamReads <- NA
  samples$bamControl <- NA
  samples$Peaks <- NA
  
  tam2 <- dba(sampleSheet = samples)
  
  tam2 <- dba.normalize(tam2, library="RiP")
  tam2 <- dba.analyze(tam2)
  tam2.res <- dba.show(tam2,bContrasts = TRUE)
  
  expect_equal(tam2.res$DB.DESeq2,tamoxifen.res$DB.DESeq2)
  
  for(i in 1:ncol(counts)) {
    countFile=sprintf("%s.txt",meta$ID[i])
    unlink(countFile)
  }

})