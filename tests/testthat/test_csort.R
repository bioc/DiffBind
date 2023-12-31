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
  y <- peakOrder(x[,1],x[,2],x[,3])
  z <- c(4,3,2,1)
  expect_equal(y,z)
})

test_that("sort large set of peaks",{
  chromNames = seq(1,95)
  n <- 1000
  chroms <- sample(chromNames,n,replace=TRUE)
  starts <- sample.int(200000,n,replace=TRUE)
  widths <- sample.int(1000,n,replace=TRUE)
  peaks <- data.frame(chroms,
                      starts,
                      starts + widths)
  names(peaks) <- c('chrom','left','right')
  ord <- peakOrder(peaks[,1],peaks[,2],peaks[,3])
  sorted <- peaks[ord,]
  expect_true(is_sorted(sorted))
})

test_that("sort equal chromosomes",{
  x = data.frame(c(6,6,6,6,6),c(6,3,9,1,3),c(3,8,3,9,9))
  y = peakOrder(x[,1],x[,2],x[,3])
  z = c(4,2,5,1,3)
  expect_equal(y,z)
})

test_that("sort equal left endpoints",{
  x = data.frame(c(6,6,6,6,6),c(9,9,9,9,9),c(3,8,1,7,9))
  y = peakOrder(x[,1],x[,2],x[,3])
  z = c(3,1,4,2,5)
  expect_equal(y,z)
})

#test_that("sort non-numeric chromosomes",{
#  x = data.frame(c("chr1","chr2","chr1","chr2","chr3"),c(9,9,9,9,9),c(3,8,1,7,9))
#  y = peakOrder(x[,1],x[,2],x[,3])
#  z = c(3,1,4,2,5)
#  expect_equal(y,z)
#})

test_that("test incomplete data",{
  x = data.frame(c(6,6,6,6,6),c(9,9,9,9,9))
  expect_error(peakOrder(x[,1],x[,2],x[,3]),"undefined columns selected")
})

test_that("test no data",{
  x = data.frame()
  expect_error(peakOrder(x[,1],x[,2],x[,3]),"undefined columns selected")
})

test_that("test no data",{
  x = NA
  expect_error(peakOrder(x[,1],x[,2],x[,3]),"incorrect number of dimensions")
})

test_that("test no data",{
  x = NA
  expect_error(peakOrder(x[,1],x[,2],x[,3]),"incorrect number of dimensions")
})
