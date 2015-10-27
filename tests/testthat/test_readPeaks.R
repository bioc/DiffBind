context("Reading Peaks")

test_that("bed with header is detected correctly",{
  fn <- system.file('extra','testdata','bed_with_header.bed.gz',package='DiffBind')
  expect_true(hasHeader(fn))
})

test_that("bed without header is detected correctly",{
  fn <- system.file('extra','testdata','bed_without_header.bed.gz',package='DiffBind')
  expect_false(hasHeader(fn))
})

test_that("raw with header is detected correctly",{
  fn <- system.file('extra','testdata','raw_with_header.raw.gz',package='DiffBind')
  expect_true(hasHeader(fn))
})

test_that("raw without header is detected correctly",{
  fn <- system.file('extra','testdata','raw_without_header.raw.gz',package='DiffBind')
  expect_false(hasHeader(fn))
})

test_that("bed with header is read correctly",{
  fn <- system.file('extra','testdata','bed_with_header.bed.gz',package='DiffBind')
  data <- pv.readbed(fn,checkHeader=TRUE)
  expect_equal(nrow(data),10)
  expect_equal(data[1,2],97113)
  expect_equal(data[1,3],100326)
})

test_that("bed without header is read correctly",{
  fn <- system.file('extra','testdata','bed_without_header.bed.gz',package='DiffBind')
  data <- pv.readbed(fn,checkHeader=TRUE)
  expect_equal(nrow(data),10)
  expect_equal(data[1,2],97113)
  expect_equal(data[1,3],100326)
})

test_that("raw with header is read correctly",{
  fn <- system.file('extra','testdata','raw_with_header.raw.gz',package='DiffBind')
  data <- pv.readbed(fn,checkHeader=TRUE)
  expect_equal(nrow(data),10)
  expect_equal(data[1,2],97113)
  expect_equal(data[1,3],100326)
})

test_that("raw without header is read correctly",{
  fn <- system.file('extra','testdata','raw_without_header.raw.gz',package='DiffBind')
  data <- pv.readbed(fn,checkHeader=TRUE)
  expect_equal(nrow(data),10)
  expect_equal(data[1,2],97113)
  expect_equal(data[1,3],100326)
})

