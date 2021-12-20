#' @title A function to downsample the raw data
#' @description Downsample the time sequence and data
#' @param timeseq a vector with T*1 size
#' the unit is (s)
#' @param dataseq a matrix with T*N size
#' @param dft downsample to dfz (Hz)
#' @return a list,contain $time_downsampled, $data_downsampled
#' @examples
#' \dontrun{
#' data("curve_data", "time_Elapsed")
#' rec <- time_downsample(time_Elapsed, curve_data, 16)
#' time_downsampled <- rec$time_downsampled
#' data_downsampled <- rec$data_downsampled
#' }
#' @export
Downsample <- function(timeseq, dataseq, dft = 16) {
  sizedata <- dim(dataseq)
  m <- sizedata[1]
  n <- sizedata[2]


  dt <- diff(timeseq)
  meanft <- mean(1. / dt)
  rfolde <- round(meanft / dft)


  if (rfolde > 1) {
    T <- ceiling(m / rfolde)
    time_downsampled <- numeric(T)
    data_downsampled <- matrix(0, T, n)
    for (i in 1:T) {
      time_downsampled[i] <- timeseq[4 * i - 3]
      for (j in 1:n) {
        data_downsampled[i, j] <- dataseq[4 * i - 3, j]
      }
    }
  } else {
    time_downsampled <- timeseq
    data_downsampled <- dataseq
  }
  rec <- list(time_downsampled = time_downsampled, data_downsampled = data_downsampled)
  return(rec)
}
