context("Peak Merging")

test_that("basic sanity of merging",{
  x <- data.frame(c(1,2,3,4),c(5,3,2,7),c(7,5,9,8))
  names(x) <- c('chr','left','right')
  y <- mergePeaks(x,-1)
  expect_equal(x,y)
})

test_that("merge overlapping intervals",{
  x <- data.frame(c(1,2,2,4),c(5,3,5,7),c(7,8,9,8))
  names(x) <- c('chr','left','right')
  y <- mergePeaks(x,-1)
  z <- data.frame(c(1,2,4),c(5,3,7),c(7,9,8))
  names(z) <- c('chr','left','right')
  expect_equal(z,y)
})

test_that("merge last intervals",{
  x <- data.frame(c(1,2,4,4),c(5,3,5,7),c(7,8,9,8))
  names(x) <- c('chr','left','right')
  y <- mergePeaks(x,-1)
  z <- data.frame(c(1,2,4),c(5,3,5),c(7,8,9))
  names(z) <- c('chr','left','right')
  expect_equal(z,y)
})

test_that("merge multiple intervals",{
  x <- data.frame(c(1,2,2,2,4),c(5,3,5,7,9),c(7,8,9,10,20))
  names(x) <- c('chr','left','right')
  y <- mergePeaks(x,-1)
  z <- data.frame(c(1,2,4),c(5,3,9),c(7,10,20))
  names(z) <- c('chr','left','right')
  expect_equal(z,y)
})

test_that("don't merge touching intervals",{
  x <- data.frame(c(1,2,2,2,4),c(5,3,5,7,9),c(7,5,9,10,20))
  names(x) <- c('chr','left','right')
  y <- mergePeaks(x,-1)
  z <- data.frame(c(1,2,2,4),c(5,3,5,9),c(7,5,10,20))
  names(z) <- c('chr','left','right')
  expect_equal(z,y)
})

test_that("merge 1bp intervals",{
  x <- data.frame(c(1,2,2,2,4),c(5,3,5,7,9),c(7,6,9,10,20))
  names(x) <- c('chr','left','right')
  y <- mergePeaks(x,-1)
  z <- data.frame(c(1,2,4),c(5,3,9),c(7,10,20))
  names(z) <- c('chr','left','right')
  expect_equal(z,y)
})

test_that("merge with maxgap less than -1 (require bigger overlap)",{
  x <- data.frame(c(1,2,2,2,4),c(5,3,5,7,9),c(7,7,11,12,20))
  names(x) <- c('chr','left','right')
  y <- mergePeaks(x,-3)
  z <- data.frame(c(1,2,2,4),c(5,3,5,9),c(7,7,12,20))
  names(z) <- c('chr','left','right')
  expect_equal(z,y)
})

test_that("merge with positive maxgap",{
  x <- data.frame(c(1,2,2,2,2),c(5,3,13,23,34),c(7,11,20,30,40))
  names(x) <- c('chr','left','right')
  y <- mergePeaks(x,3)
  z <- data.frame(c(1,2,2),c(5,3,34),c(7,30,40))
  names(z) <- c('chr','left','right')
  expect_equal(z,y)
})

test_that("merge scores basic sanity",{
  x <- data.frame(c(1,2,2,2,4),c(20,100,150,200,40),c(40,120,170,220,60))
  s <- c(20,30,40,50,60)
  y <- data.frame(c(1,2,2,4),c(20,100,155,40),c(30,120,165,60),c(21,29,100,100))
  z <- mergeScores(x,s,y)
  expect_equal(z$score,c(21,30,100,50,100))
  expect_equal(z$included,c(1,1,1,0,1))
})

test_that("merge scores with multiple peaks in y matching one in x",{
  x <- data.frame(c(1,2,2,2,4),c(20,100,150,200,40),c(40,120,170,220,60))
  s <- c(20,30,40,50,60)
  y <- data.frame(c(1,2,2,2,2,4),c(20,100,111,112,155,40),c(30,110,120,120,165,60),c(21,40,39,38,100,100))
  z <- mergeScores(x,s,y)
  expect_equal(z$score,c(21,40,100,50,100))
  expect_equal(z$included,c(1,1,1,0,1))
})

test_that("merge scores with multiple incoming peak sets",{
  x <- data.frame(c(1,2,2,2,4),c(20,100,150,200,40),c(40,120,170,220,60))
  s <- c(0,0,0,0,0)
  names(x) <- c('chr','left','right')
  y <- data.frame(c(2,2,2,2,4),c(100,111,112,155,40),c(110,120,120,165,60),c(40,39,38,100,100))
  names(y) <- c('chr','left','right','score')
  z <- data.frame(c(1),c(28),c(30),c(9))
  names(z) <- c('chr','left','right','score')
  a <- mergeScores(x,s,y)
  s <- a$score
  b <- mergeScores(x,s,z)
  expect_equal(b$score,c(9,40,100,0,100))
  expect_equal(a$included,c(0,1,1,0,1))
  expect_equal(b$included,c(1,0,0,0,0))
})

test_that("merge scores with multiple incoming peak sets and abs flag set to false",{
  x <- data.frame(c(1,2,2,2,4),c(20,100,150,200,40),c(40,120,170,220,60))
  s <- c(0,0,0,0,0)
  names(x) <- c('chr','left','right')
  y <- data.frame(c(2,2,2,2,4),c(100,111,112,155,40),c(110,120,120,165,60),c(40,39,38,100,100))
  names(y) <- c('chr','left','right','score')
  z <- data.frame(c(1),c(28),c(30),c(9))
  names(z) <- c('chr','left','right','score')
  a <- mergeScores(x,s,y,FALSE)
  s <- a$score
  b <- mergeScores(x,s,z,FALSE)
  expect_equal(b$score,c(9,40,100,0,100))
  expect_equal(a$included,c(0,1,1,0,1))
  expect_equal(b$included,c(1,0,0,0,0))
})

test_that("merge scores with multiple incoming peak sets and abs flag set to true",{
  x <- data.frame(c(1,2,2,2,4),c(20,100,150,200,40),c(40,120,170,220,60))
  s <- c(0,0,0,0,0)
  names(x) <- c('chr','left','right')
  y <- data.frame(c(2,2,2,2,4),c(100,111,112,155,40),c(110,120,120,165,60),c(40,39,38,100,100))
  names(y) <- c('chr','left','right','score')
  z <- data.frame(c(1),c(28),c(30),c(9))
  names(z) <- c('chr','left','right','score')
  a <- mergeScores(x,s,y,TRUE)
  s <- a$score
  b <- mergeScores(x,s,z,TRUE)
  expect_equal(b$score,c(9,40,100,0,100))
  expect_equal(a$included,c(0,1,1,0,1))
  expect_equal(b$included,c(1,0,0,0,0))
})

test_that("merge scores with negative scores and abs flag set to false",{
  x <- data.frame(c(1,2,2,2,4),c(20,100,150,200,40),c(40,120,170,220,60))
  s <- c(0,0,0,0,0)
  names(x) <- c('chr','left','right')
  y <- data.frame(c(2,2,2,2,4),c(100,111,112,155,40),c(110,120,120,165,60),c(40,39,38,100,100))
  names(y) <- c('chr','left','right','score')
  z <- data.frame(c(1),c(28),c(30),c(-9))
  names(z) <- c('chr','left','right','score')
  a <- mergeScores(x,s,y,FALSE)
  s <- a$score
  b <- mergeScores(x,s,z,FALSE)
  expect_equal(b$score,c(0,40,100,0,100))
  expect_equal(a$included,c(0,1,1,0,1))
  expect_equal(b$included,c(1,0,0,0,0))
})

test_that("merge scores with negative scores and abs flag set to true",{
  x <- data.frame(c(1,2,2,2,4),c(20,100,150,200,40),c(40,120,170,220,60))
  s <- c(0,0,0,0,0)
  names(x) <- c('chr','left','right')
  y <- data.frame(c(2,2,2,2,4),c(100,111,112,155,40),c(110,120,120,165,60),c(-40,39,38,100,100))
  names(y) <- c('chr','left','right','score')
  z <- data.frame(c(1,2),c(28,100),c(30,110),c(-9,40))
  names(z) <- c('chr','left','right','score')
  a <- mergeScores(x,s,y,TRUE)
  s <- a$score
  b <- mergeScores(x,s,z,TRUE)
  expect_equal(a$score,c(0,-40,100,0,100))
  expect_equal(b$score,c(-9,40,100,0,100))
  expect_equal(a$included,c(0,1,1,0,1))
  expect_equal(b$included,c(1,1,0,0,0))
})

test_that("merge scores with double negative scores and abs flag set to true",{
  x <- data.frame(c(1,2,2,2,4),c(20,100,150,200,40),c(40,120,170,220,60))
  s <- c(0,0,0,0,0)
  names(x) <- c('chr','left','right')
  y <- data.frame(c(2,2,2,2,4),c(100,111,112,155,40),c(110,120,120,165,60),c(-40,39,38,100,100))
  names(y) <- c('chr','left','right','score')
  z <- data.frame(c(1,2),c(28,100),c(30,110),c(-9,-40))
  names(z) <- c('chr','left','right','score')
  a <- mergeScores(x,s,y,TRUE)
  s <- a$score
  b <- mergeScores(x,s,z,TRUE)
  expect_equal(a$score,c(0,-40,100,0,100))
  expect_equal(b$score,c(-9,-40,100,0,100))
  expect_equal(a$included,c(0,1,1,0,1))
  expect_equal(b$included,c(1,1,0,0,0))
})



#test_that("merge scores with non-numeric arguments",{
#  x <- data.frame(c("chr1","chr2","chr2","chr2","chr4"),c(20,100,150,200,40),c(40,120,170,220,60))
#  s <- c(20,30,40,50,60)
#  y <- data.frame(c('chr1','chr2','chr2','chr4'),c(20,100,155,40),c(30,120,165,60),c(21,29,100,100))
#  z <- mergeScores(x,s,y)
#  expect_equal(z$score,c(21,30,100,50,100))
#  expect_equal(z$included,c(1,1,1,0,1))
#})

test_that("merge scores with integers",{
  aa <- as.integer(c(1,2,2,2,5))
  bb <- as.integer(c(20,100,150,200,40))
  cc <- as.integer(c(40,120,170,220,60))
  x <- data.frame(aa,bb,cc)
  s <- as.integer(c(20,30,40,50,60))
  y <- data.frame(c(1,2,2,5),c(20,100,155,40),c(30,120,165,60),c(21,29,100,100))
  z <- mergeScores(x,s,y)
  expect_equal(z$score,c(21,30,100,50,100))
  expect_equal(z$included,c(1,1,1,0,1))
})

test_that("merge scores with factor",{
  aa <- as.factor(c(1,2,2,2,5))
  bb <- as.integer(c(20,100,150,200,40))
  cc <- as.integer(c(40,120,170,220,60))
  x <- data.frame(aa,bb,cc)
  s <- as.integer(c(20,30,40,50,60))
  y <- data.frame(as.factor(c(1,2,2,5)),c(20,100,155,40),c(30,120,165,60),c(21,29,100,100))
  z <- mergeScores(x,s,y)
  expect_equal(z$score,c(21,30,100,50,100))
  expect_equal(z$included,c(1,1,1,0,1))
})

