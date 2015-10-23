context("Interval Sorting")

is_sorted <- function(peaks) {
  for (i in 1:(nrow(peaks)-1)) {
    if (peaks[i,1] > peaks[i+1,1] ||
        (peaks[i,1] == peaks[i+1,1] &&
         (peaks[i,2] > peaks[i+1,2] ||
          (peaks[i,2] == peaks[i+1,2] && peaks[i,3] > peaks[i+1,3])))) {
      return(FALSE)
    }
  }
  return(TRUE)
}
                
test_that("basic sanity of sorting",{
  x <- data.frame(c(4,3,2,1),c(5,3,2,7),c(7,5,2,7))
  y <- cpp_peakOrder(x)
  z <- c(4,3,2,1)
  expect_equal(y,z)
})

test_that("sort large set of peaks",{
  chromNames = seq(1,95)
  n <- 10000
  chroms <- sample(chromNames,n,replace=TRUE)
  starts <- sample.int(200000,n,replace=TRUE)
  widths <- sample.int(1000,n,replace=TRUE)
  peaks <- data.frame(chroms,
                      starts,
                      starts + widths)
  names(peaks) <- c('chrom','left','right')
  ord <- cpp_peakOrder(peaks)
  sorted <- peaks[ord,]
  expect_true(is_sorted(sorted))
})

test_that("sort equal chromosomes",{
  x = data.frame(c(6,6,6,6,6),c(6,3,9,1,3),c(3,8,3,9,9))
  y = cpp_peakOrder(x)
  z = c(4,2,5,1,3)
  expect_equal(y,z)
})

test_that("sort equal left endpoints",{
  x = data.frame(c(6,6,6,6,6),c(9,9,9,9,9),c(3,8,1,7,9))
  y = cpp_peakOrder(x)
  z = c(3,1,4,2,5)
  expect_equal(y,z)
})
