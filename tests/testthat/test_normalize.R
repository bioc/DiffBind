test_that("normalizing data",{
  
  if(Sys.info()["sysname"] != "Windows") {
    data(tamoxifen_counts)
    tamoxifen <- dba.contrast(tamoxifen)
    tamoxifen <- dba.contrast(tamoxifen,design="~Tissue + Condition")
    
    for (nf in c(DBA_NORM_DEFAULT, DBA_NORM_LIB, DBA_NORM_TMM, DBA_NORM_RLE)) {
      tam <- dba.normalize(tamoxifen,method = DBA_ALL_METHODS,
                           library = DBA_LIBSIZE_PEAKREADS, normalize = nf)
      set.seed(3793)
      tam <- dba.analyze(tam,method=DBA_ALL_METHODS)
      rep1.e <- dba.report(tam,method=DBA_EDGER,bCounts=TRUE, 
                           bCalled=TRUE, bCalledDetail=TRUE)
      rep1.d <- dba.report(tam,method=DBA_DESEQ2,bCounts=TRUE, 
                           bCalled=TRUE, bCalledDetail=TRUE)
      
      expect_true(length(rep1.e) > 500)
      expect_true(length(rep1.d) > 500)
      
      norm <- dba.normalize(tam, bRetrieve=TRUE, method=DBA_ALL_METHODS)
      tam$norm <- NULL
      tam <- dba.normalize(tam,method = DBA_EDGER,
                           library = norm$edgeR$lib.sizes, 
                           normalize = norm$edgeR$norm.factors,
                           background=FALSE)
      tam <- dba.normalize(tam,method = DBA_DESEQ2,
                           library = norm$DESeq2$lib.sizes, 
                           normalize = norm$DESeq2$norm.factors,
                           background=FALSE)
      set.seed(3793)
      tam <- dba.analyze(tam,method=DBA_ALL_METHODS)
      rep2.e <- dba.report(tam,method=DBA_EDGER,bCounts=TRUE, 
                           bCalled=TRUE, bCalledDetail=TRUE)
      rep2.d <- dba.report(tam,method=DBA_DESEQ2,bCounts=TRUE, 
                           bCalled=TRUE, bCalledDetail=TRUE)
      
      expect_equal(sum(rep1.e != rep2.e),0)
      expect_equal(sum(rep1.d != rep2.d),0)
      
      if(nf==DBA_NORM_LIB) {
        tam <- dba.normalize(tamoxifen,method = DBA_ALL_METHODS,
                             library = DBA_LIBSIZE_FULL, 
                             normalize = nf, offsets=FALSE)
        set.seed(3793)
        tam <- dba.analyze(tam,method=DBA_ALL_METHODS)
        rep1.e <- dba.report(tam,method=DBA_EDGER,bCounts=TRUE, 
                             bCalled=TRUE, bCalledDetail=TRUE)
        rep1.d <- dba.report(tam,method=DBA_DESEQ2,bCounts=TRUE, 
                             bCalled=TRUE, bCalledDetail=TRUE)
        
        expect_true(length(rep1.e) > 500)
        expect_true(length(rep1.d) > 500)
        
        norm <- dba.normalize(tam, bRetrieve=TRUE, method=DBA_ALL_METHODS)
        tam$norm <- NULL
        tam <- dba.normalize(tam,method = DBA_EDGER,
                             library = norm$edgeR$lib.sizes, 
                             normalize = norm$edgeR$norm.factors,
                             background=FALSE)
        tam <- dba.normalize(tam,method = DBA_DESEQ2,
                             library = norm$DESeq2$lib.sizes, 
                             normalize = norm$DESeq2$norm.factors,
                             background=FALSE)
        set.seed(3793)
        tam <- dba.analyze(tam,method=DBA_ALL_METHODS)
        rep2.e <- dba.report(tam,method=DBA_EDGER,bCounts=TRUE, 
                             bCalled=TRUE, bCalledDetail=TRUE)
        rep2.d <- dba.report(tam,method=DBA_DESEQ2,bCounts=TRUE, 
                             bCalled=TRUE, bCalledDetail=TRUE)
        
        expect_equal(sum(rep1.e != rep2.e),0)
        expect_equal(sum(rep1.d != rep2.d),0)
        
        tam <- dba.normalize(tamoxifen,method = DBA_ALL_METHODS,
                             library = DBA_LIBSIZE_PEAKREADS, 
                             normalize = nf, offsets=TRUE)
        set.seed(3793)
        tam <- dba.analyze(tam,method=DBA_ALL_METHODS)
        rep1.e <- dba.report(tam,method=DBA_EDGER,bCounts=TRUE, 
                             bCalled=TRUE, bCalledDetail=TRUE)
        rep1.d <- dba.report(tam,method=DBA_DESEQ2,bCounts=TRUE, 
                             bCalled=TRUE, bCalledDetail=TRUE)
        
        expect_true(length(rep1.e) > 500)
        expect_true(length(rep1.d) > 500)
        
        norm <- dba.normalize(tam, bRetrieve=TRUE, method=DBA_ALL_METHODS)
        tam$norm <- NULL
        tam <- dba.normalize(tam,method = DBA_EDGER,
                             library = norm$edgeR$lib.sizes, 
                             offsets   = norm$offsets$offsets,
                             normalize = DBA_NORM_OFFSETS_ADJUST,
                             background=FALSE)
        tam <- dba.normalize(tam,method = DBA_DESEQ2,
                             library = norm$DESeq2$lib.sizes, 
                             normalize = DBA_NORM_OFFSETS_ADJUST,
                             offsets   = norm$offsets$offsets,
                             background=FALSE)
        set.seed(3793)
        tam <- dba.analyze(tam,method=DBA_ALL_METHODS)
        rep2.e <- dba.report(tam,method=DBA_EDGER,bCounts=TRUE, 
                             bCalled=TRUE, bCalledDetail=TRUE)
        rep2.d <- dba.report(tam,method=DBA_DESEQ2,bCounts=TRUE, 
                             bCalled=TRUE, bCalledDetail=TRUE)
        
        expect_equal(sum(rep1.e != rep2.e),0)
        expect_equal(sum(rep1.d != rep2.d),0)
      }
    }
  }
})