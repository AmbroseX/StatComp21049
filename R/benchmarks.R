#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to comp
#' @examples
#' \dontrun{
#' data(data); attach(data)
#' tm1 <- microbenchmark::microbenchmark(
#' rnR = gibbsR(100,10),
#' ...
#' }
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm rgamma
#' @useDynLib StatComp21049
NULL
