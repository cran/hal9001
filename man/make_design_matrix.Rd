% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{make_design_matrix}
\alias{make_design_matrix}
\title{Build HAL Design Matrix}
\usage{
make_design_matrix(X, blist, p_reserve = 0.5)
}
\arguments{
\item{X}{Matrix of covariates containing observed data in the columns.}

\item{blist}{List of basis functions with which to build HAL design matrix.}

\item{p_reserve}{Sparse matrix pre-allocation proportion. Default value is 0.5.
If one expects a dense HAL design matrix, it is useful to set p_reserve to a higher value.}
}
\value{
A \code{dgCMatrix} sparse matrix of indicator basis functions
corresponding to the design matrix in a zero-order highly adaptive lasso.
}
\description{
Make a HAL design matrix based on original design matrix X and a list of
basis functions in argument blist
}
\examples{
\donttest{
gendata <- function(n) {
  W1 <- runif(n, -3, 3)
  W2 <- rnorm(n)
  W3 <- runif(n)
  W4 <- rnorm(n)
  g0 <- plogis(0.5 * (-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4))
  A <- rbinom(n, 1, g0)
  Q0 <- plogis(0.15 * (2 * A + 2 * A * W1 + 6 * A * W3 * W4 - 3))
  Y <- rbinom(n, 1, Q0)
  data.frame(A, W1, W2, W3, W4, Y)
}
set.seed(1234)
data <- gendata(100)
covars <- setdiff(names(data), "Y")
X <- as.matrix(data[, covars, drop = FALSE])
basis_list <- enumerate_basis(X)
x_basis <- make_design_matrix(X, basis_list)
}

}
