#' @title Embedding the data after PCA dimensionality reduction
#' @description Embedding the data after PCA dimensionality reduction
#' @param X a matrix with T*N size
#' @param K_mode_dimension is the dimensionality of PCA dimensionality reduction
#' @param K_delay is the delay time
#' @return a matrix
#' @examples
#' \dontrun{
#' data("curve_data", "time_Elapsed", "w")
#' ## you can use your data to Proceed to the steps below
#' rec <- Downsample(time_Elapsed, curve_data, 16)
#' time_downsampled <- rec$time_downsampled
#' data_downsampled <- rec$data_downsampled
#' ev <- Get_Matrix_Eigens(data_downsampled)
#' eigenvalues <- ev$values
#' eigenvectors <- ev$vectors
#' PCAmode <- 5
#' K_delay <- 12
#' Xeigenvector <- eigenvectors[, 1:PCAmode]
#' X_PCA <- data_downsampled %*% Xeigenvector
#' test_CurveFilter_embeddingdelay <- Embedding_delay(X_PCA, PCAmode, K_delay)
#' X <- test_CurveFilter_embeddingdelay %*% w
#' ## Plot the projection
#' plot(X[, 3], X[, 4], type = "l")
#' plot(X[, 6], X[, 7], type = "l")
#' }
#' @export
Embedding_delay <- function(X, K_mode_dimension = 5, K_delay = 12) {
  num_Frame <- dim(X)[1]
  Ybase <- matrix(0, num_Frame - K_delay + 1, K_mode_dimension * K_delay)
  j <- 1
  for (k in 1:K_delay) {
    Ybase[, j:(j + K_mode_dimension - 1)] <- X[k:(k + num_Frame - K_delay), 1:K_mode_dimension]
    j <- j + K_mode_dimension
  }
  return(Ybase)
}
