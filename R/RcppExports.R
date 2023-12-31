# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

mergePeaks <- function(data, maxGap) {
    .Call('_DiffBind_mergePeaks', PACKAGE = 'DiffBind', data, maxGap)
}

mergeScores <- function(sMerged, sScore, sPeaks, abs = NULL) {
    .Call('_DiffBind_mergeScores', PACKAGE = 'DiffBind', sMerged, sScore, sPeaks, abs)
}

peakOrder <- function(schrom, sleft, sright) {
    .Call('_DiffBind_peakOrder', PACKAGE = 'DiffBind', schrom, sleft, sright)
}

