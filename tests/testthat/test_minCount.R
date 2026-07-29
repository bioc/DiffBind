## Regression tests for minCount preservation.
##
## Background: pv.vectors() rebuilds a DBA object from its parts and is called
## by many DiffBind operations (dba(), dba.count(), dba.peakset() retrieval,
## masking, analysis modelling). A copy-paste typo (pv$minCount <- pv$minCount,
## reading from the freshly-nulled object instead of the saved local) caused it
## to drop $minCount. Most callers happened to re-stamp it, so the defect was
## latent -- but it could silently reset a user-set minCount.
##
## These tests exercise pv.vectors() directly (the exact function fixed), using
## the bundled counted dataset so they run under R CMD check without BAM files.
## The count-path aspect (dba.count(minCount=) with summit re-centering) needs
## BAMs and is verified separately outside the package test suite.

test_that("pv.vectors() preserves a user-set minCount", {
  data(tamoxifen_counts)

  tamoxifen$minCount <- 7L
  rebuilt <- DiffBind:::pv.vectors(tamoxifen, minOverlap = 1, bAllSame = TRUE)
  expect_identical(rebuilt$minCount, 7L)
})

test_that("pv.vectors() round-trips minCount = 0 (default) unchanged", {
  data(tamoxifen_counts)

  tamoxifen$minCount <- 0
  rebuilt <- DiffBind:::pv.vectors(tamoxifen, minOverlap = 1, bAllSame = TRUE)
  expect_identical(rebuilt$minCount, 0)
})
