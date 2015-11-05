context("Loading Sample Sheets")

collectWarnings <- function(fun,data) {
  msgs <- NULL
  wHandler <- function(w) {
    msgs <<- c(msgs,w$message)
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(do.call(fun,data),warning=wHandler)
  return(list("value"=val,"warnings"=msgs))
}

test_that("spaces are trimmed from sample sheets, basic sanity",{
  data <- data.frame(c('a','b','c '),c('x',' y','z'))
  ans <- data.frame(c('a','b','c'),c('x','y','z'))
  colnames(data) <- c('alpha','bravo')
  colnames(ans) <- c('alpha','bravo')
  res <- suppressWarnings(stripSpaces(data))
  expect_equal(sum(res$alpha==ans$alpha),3)
  expect_equal(sum(res$bravo==ans$bravo),3)
})

test_that("warning is generated correctly",{
  data <- data.frame(c('a','b','c '),c('x','y','z'))
  colnames(data) <- c('alpha','bravo')
  expect_warning(stripSpaces(data),
                 'Removed white space from c in column alpha (row 3)',
                 fixed=TRUE)
})

test_that("warning multiple rows is generated correctly",{
  data <- data.frame(c('a','b ','c '),c('x','y','z'))
  colnames(data) <- c('alpha','bravo')
  expect_warning(stripSpaces(data),
                 'Removed white space from b,c in column alpha (rows 2,3)',
                 fixed=TRUE)
})

test_that("warning multiple columns generates multiple warnings",{
  data <- data.frame(c('a','b ','c '),c('x ','y','z'))
  colnames(data) <- c('alpha','bravo')
  res <- collectWarnings(stripSpaces,list(data))
  msgs <- res$warnings
  expect_equal(2,length(msgs))
  expect_equal(msgs[1],"Removed white space from b,c in column alpha (rows 2,3)",fixed=TRUE)
  expect_equal(msgs[2],"Removed white space from x in column bravo (row 1)",fixed=TRUE)
})

test_that("loading csv sample sheet with spaces generates warnings",{
  wd <- getwd()
  setwd(system.file('extra','testdata',package='DiffBind'))
  res <- collectWarnings(dba,list('sampleSheet'='test_sampleSheet_Spaces.csv','bCorPlot'=FALSE))
  expect_equal(2,length(res$warnings))
  expect_equal(res$warnings[1],"Removed white space from treated in column Condition (row 3)",fixed=TRUE)
  expect_equal(res$warnings[2],"Removed white space from bed_with_header.bed.gz,bed_without_header.bed.gz,raw_with_header.raw.gz in column Peaks (rows 1,2,3)",fixed=TRUE)
  expect_equal(res$value$samples$Condition[3],"treated")
  expect_equal(res$value$samples$Condition[4],"treated")
  expect_equal(res$value$samples$Peaks[1],"bed_with_header.bed.gz")
  expect_equal(res$value$samples$Peaks[2],"bed_without_header.bed.gz")
  expect_equal(res$value$samples$Peaks[3],"raw_with_header.raw.gz")
  expect_equal(res$value$samples$Peaks[4],"raw_without_header.raw.gz")
  setwd(wd)
})

test_that("loading xlsx sample sheet with spaces generates warnings",{
  if (require(XLConnect)) {
    wd <- getwd()
    setwd(system.file('extra','testdata',package='DiffBind'))
    res <- collectWarnings(dba,list('sampleSheet'='test_sampleSheet_Spaces.csv','bCorPlot'=FALSE))
    expect_equal(2,length(res$warnings))
    expect_equal(res$warnings[1],"Removed white space from treated in column Condition (row 3)",fixed=TRUE)
    expect_equal(res$warnings[2],"Removed white space from bed_with_header.bed.gz,bed_without_header.bed.gz,raw_with_header.raw.gz in column Peaks (rows 1,2,3)",fixed=TRUE)
    expect_equal(res$value$samples$Condition[3],"treated")
    expect_equal(res$value$samples$Condition[4],"treated")
    expect_equal(res$value$samples$Peaks[1],"bed_with_header.bed.gz")
    expect_equal(res$value$samples$Peaks[2],"bed_without_header.bed.gz")
    expect_equal(res$value$samples$Peaks[3],"raw_with_header.raw.gz")
    expect_equal(res$value$samples$Peaks[4],"raw_without_header.raw.gz")
    setwd(wd)
  }
})

test_that("sample sheet loading does not interpret hash as a comment character",{
  wd <- getwd()
  setwd(system.file('extra','testdata',package='DiffBind'))
  fn = 'test_sampleSheet_comment_char.csv'
  dobj = dba(sampleSheet=fn)
  expect_is(dobj,"DBA")
  expect_equal(dobj$samples$SampleID[2],"bravo#")
  expect_equal(dobj$samples$SampleID[3],"char#lie")
  setwd(wd)
})
