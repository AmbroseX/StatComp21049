#' @title A function that calculates matrix eigenvalues and eigenvectors
#' @description Returns a list of eigenvalues and eigenvectors for a matrix
#' @param X a matrix with T*N size
#' @return a list with eigenvalues and eigenvectors
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(200 * 100), 200, 100)
#' ev <- Get_Matrix_Eigens(X)
#' eigenvalues <- ev$values
#' eigenvectors <- ev$vectors
#' }
Get_Matrix_Eigens <- function(X) {
  CovX <- t(X) %*% t(t(X))
  ev <- eigen(CovX)
  return(ev)
}
