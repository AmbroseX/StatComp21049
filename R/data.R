#' Example data of the curve fillter for a worm
#'
#' Each line represents the curvature of 100 segments of the body at a moment
#' Column i represents the curvature of the body of the ith segment over time.
#'
#' @format A Matrix with 43,288 rows times 100 columns.
#' @usage data(curve_data)
#' @seealso
#' [time_Elapsed],[w]
#' @source \url{http://www.wenlab.org/}
"curve_data"


#' Example data of the time series for a worm
#'
#' The unit is seconds.
#'
#' @format A Matrix with 43288*100.
#' @usage data(time_Elapsed)
#' @seealso
#' [curve_data],[w]
#' @source \url{http://www.wenlab.org/}
"time_Elapsed"

#' Example Projection matrix for a worm
#'
#' @format A Matrix with 60*7
#' @usage data(w)
#' @seealso
#' [curve_data],[time_Elapsed]
#' @source \url{http://www.wenlab.org/}
"w"








