#' @title A Gibbs sampler using R
#' @description A Gibbs sampler using R
#' @param N the length of the chain
#' @param thin the number of Binomial(thin, y)
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' rnR <- gibbsR(100,10)
#' plot(rnR[1,],type='l')
#' plot(rnR[2,],type='l')
#' }
#' @export
#' @importFrom stats rbinom rbeta
gibbsR <- function(N, thin){
  X <- matrix(0,2,N);
  X[1,1] <- 1
  X[2,1] <- 0.05
  a <- 1
  b <- 2
  X1 <- 0
  Y1 <- 0
  for(i in 2:N){
    Y1 <- X[2,i-1]
    X[1,i] <- rbinom(1,thin,Y1)
    X1 <- X[1,i]
    X[2,i] <- rbeta(1,X1+a,thin-X1+b)
  }
  return(X);
}
