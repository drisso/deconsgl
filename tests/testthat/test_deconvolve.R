context("Test deconvolve function.")

BiocParallel::register(BiocParallel::SerialParam())

test_that("coefs for random covariates are 0", {
  set.seed(13124)
  # No graph information
  Y <- matrix(rnorm(30), ncol=3, nrow=10)
  X <- matrix(rnorm(30), ncol=3, nrow=10)
  rownames(X) <- rownames(Y) <- paste0("gene", 1:10)
  colnames(X) <- paste0("celltype", 1:3)
  colnames(Y) <- paste0("sample", 1:3)

  groups <- c(1, 2, 2)

  res1 <- deconvolve(Y, X, purity = 0.5)
  res2 <- deconvolve(Y, X, groups = groups, purity = 0.5)

  expect_true(all(res1$absolute == 0))
  expect_true(all(res2$absolute == 0))
})

test_that("coefs for important covariates are not 0", {
  set.seed(13124)
  # No graph information
  X <- matrix(c(rep(1, 10), rnorm(30)), ncol=4, nrow=10)
  y <- 10 + 2 * X[,2] + rnorm(10, sd=.2)
  Y <- matrix(y, ncol=1)
  rownames(X) <- rownames(Y) <- paste0("gene", 1:10)
  colnames(X) <- paste0("celltype", 1:4)
  colnames(Y) <- paste0("sample", 1)

  groups <- c(1, 2, 2, 2)

  res1 <- deconvolve(Y, X, purity = 0.5)
  res2 <- deconvolve(Y, X, groups = groups, purity = 0.5)

  expect_true(res1$absolute[2,1] > 0)
  expect_true(res2$absolute[2,1] > 0)

  expect_true(all(res1$absolute >= 0))
  expect_true(all(res2$absolute >= 0))

})
