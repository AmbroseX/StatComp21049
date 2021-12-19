#' @title A function that Plot the number of dimensions from a segment to a segment
#' @description Plot the dimensions from a segment to a segment
#' @param X a matrix with  T*N size
#' @param Sstart the Start segment
#' @param Send the end segment, 1< Sstart<= Send < 100
#' @return a number, the dimension of the matrix X
#' @examples
#' \dontrun{
#' X <- Plot_Matrix_Dimension()
#' }
#' @importFrom utils data
#' @export
Plot_Matrix_Dimension <- function(X = 100, Sstart=1,Send = 100){
  if(X == 100){
    X <- StatComp21049::curve_data
    #data(curve_data)
    # X <- curve_data
  }
  n <- Send-Sstart+1

  Pdim <- numeric(n)
  flag <- 1
  for(i in Sstart:Send){
    Pdim[flag] <- Get_Matrix_Dimension(X[ ,Sstart:i])
    flag <- flag+1
  }
  plot(c(Sstart:Send), Pdim, xlab = "Segments", ylab = "Dimension")
}
