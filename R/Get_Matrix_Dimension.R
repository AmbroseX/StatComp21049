#' @title A function that Calculates the dimension of the matrix
#' @description Returns the dimension of the matrix X
#' @param X a matrix with  T*N size
#' @return a number, the dimension of the matrix X
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(200 * 100), 200, 100)
#' ev <- Get_Matrix_Dimension(X)
#' }
#' @export
Get_Matrix_Dimension <- function(X) {
  ev <- Get_Matrix_Eigens(X)
  eigenvalues <- ev$values
  PR <- sum(eigenvalues)^2 / (sum(eigenvalues^2))
  return(PR)
}
