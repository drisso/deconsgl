#' Deconvolution of tumor samples
#'
#' Given a matrix of tumor samples and a matrix of cell type gene signatures,
#' this function returns a vector of cell type composition.
#'
#' @details Optionally, the results can be standardize to sum to one, to give
#'   the proportion of (non-tumor) cell types in each sample. However, often the
#'   unstandardized coefficients are more informative as they are inversely
#'   related to the tumor purity of each sample.
#'
#' @details The method is based on the sparse group lasso (SGL) as implemented
#'   in the lsgl package (negative coefficients are forced to 0).
#'
#' @details Note that both X and Y need to have rownames. The common rownames
#'   will be used as the genes used for the deconvolution.
#'
#' @author Davide Risso
#'
#' @param Y a matrix of log gene expression to deconvolute. Genes in row,
#'   samples in columns.
#' @param X a matrix with the cell type gene signatures (log expression) to be
#'   used for deconvolution. Genes in row and cell types in columns.
#' @param purity numeric vector of length the number of samples in Y. It
#'   specifies the purity of each sample.
#' @param groups vector specifying the groups of cell types
#'   (a la group lasso).
#' @param intercept logical. Whether an intercept should be added to X.
#' @param ... additional parameters to be passed to or \code{\link[lsgl]{fit}}.
#'   See its man pages for details.
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#'
#' @details The deconvolution uses the purity of the samples as prior
#'   information to penalize the coefficients (the more purity the more likely
#'   the coefficients to be zero).
#'
#' @return A list with two components: the estimated proportion of cell types
#'   with and without standardization to sum to one.
#' @export
#'
#' @importFrom lsgl cv fit lsgl.algorithm.config
#' @importFrom sglOptim best_model
#' @import BiocParallel
#'
#' @examples
#'
#' X <- matrix(rnorm(30), ncol=3, nrow=10)
#' Y <- matrix(ncol=2, nrow=10)
#' Y[,1] <- 10 + 2 * X[,2] + rnorm(10, sd=.2)
#' Y[,2] <- X[,1] + 3 * X[,3] + rnorm(10, sd=.2)
#' rownames(X) <- rownames(Y) <- paste0("gene", 1:10)
#' colnames(Y) <- paste0("sample", 1:2)
#' colnames(X) <- paste0("celltype", 1:3)
#' res <- deconvolve(Y, X, purity = 0.5,
#'                   BPPARAM = BiocParallel::SerialParam())
#'
deconvolve <- function(Y, X,
                       purity = 1,
                       groups = NULL,
                       intercept = TRUE, ...,
                       BPPARAM = BiocParallel::bpparam()) {

  if(any(is.null(rownames(X)) | is.null(rownames(Y)))) {
    stop("X and Y need to have rownames.")
  }

  if(!(all(rownames(X) %in% rownames(Y)))) {
    stop("Not all genes in the signature are present in the expression matrix.",
         "Keeping only the genes present in both.")
  }
  genes <- intersect(rownames(X), rownames(Y))

  X <- X[genes, , drop = FALSE]
  Y <- Y[genes, , drop = FALSE]

  penalty <- purity

  if(length(penalty) == 1) {
    penalty <- rep(penalty, NCOL(Y))
  } else if(length(penalty) != NCOL(Y)) {
    stop("`penalty` must be either a scalar or a vector of length the number of samples.")
  }

  res <- bplapply(seq_len(NCOL(Y)), function(j) {
      fit_cv <- lsgl::cv(X, Y[,j], lambda = penalty[j], grouping = groups,
                         ..., intercept = intercept,
                         algorithm.config = lsgl.algorithm.config(verbose=FALSE))
      res <- lsgl::fit(X, Y[,j], lambda = penalty[j], grouping = groups,
                       ..., intercept = intercept,
                       algorithm.config = lsgl.algorithm.config(verbose=FALSE))
      return(res$beta[[best_model(fit_cv)]])
  }, BPPARAM = BPPARAM)
  res <- sapply(res, as.vector)
  res[res<0] <- 0

  names <- colnames(X)
  if(intercept) {
    names <- c("Intercept", names)
  }

  rownames(res) <- names

  colnames(res) <- colnames(Y)

  if(intercept) {
    absolute <- res[-1, , drop = FALSE]
  } else {
    absolute <- res
  }

  retval <- list(coefficients = res,
                 absolute = absolute,
                 relative = scale(absolute, scale = pmax(colSums(absolute), 1),
                                  center = FALSE))
  return(retval)
}
