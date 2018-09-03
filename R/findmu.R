#' Determine equidistant means of trinormal ROC data simulation
#'
#' A function that computes the equidistant means \code{muy} and \code{muz} for
#' a specific \code{mux}. The VUS as well as the set of standard errors are
#' given as arguments to the function.
#'
#' @details Defaults are: VUS = 1/6, standard errors for all three classes equal
#'   1. The searching algorithm is stepwise increasing the differences
#'   \code{muy-mux} and \code{muz-mux} according to the variable \code{step}.
#'   The algorithm stops when the computed VUS exceeds the preferred VUS. The
#'   according parameters \code{mux, muy, muz} are returned with the computed
#'   VUS.
#'
#'   Remark: The bigger \code{VUS} and the smaller \code{step} is chosen, the
#'   longer the computation lasts.
#'
#' @param mux  The numeric mean of the healthy class. Default is zero.
#' @param sdx,sdy,sdz The numeric standard errors of the healthy, intermediate
#'   and diseased class, for which the according means have to be determined
#'   given a specifiv VUS.
#' @param VUS The Volume Under the Surface. A numeric value between 1/6 and 1. Default
#' is 1/6.
#' @param step A numeric indicating the step size each iteration takes in order to
#' find the closest set of means. Default set to 0.001.
#' @return A data frame with the following components:
#'   \item{mux}{The initial mean of the healthy class}
#'   \item{muy}{The mean of the intermediate class computed for the specified \code{VUS}.}
#'   \item{muz}{The mean of the diseased class computed for the specified \code{VUS}.}
#'   \item{VUS}{The VUS computed for \code{mux}, \code{muy} and \code{muz}.}
#' @export
#' @examples
#' # find equidistant means with mux=2.7 and VUS = 0.45:
#' findmu(mux = 2.7, VUS = 0.45)
#' # specify standard errors:
#' findmu(mux = 2.7, sdx = 1.1, sdy = 1.3, sdz = 1.5, VUS = 0.45)


findmu <- function(mux=0, sdx=1, sdy=1, sdz=1, VUS = 1/6, step=0.001) {

  # check if sd's are appropriatly set:
  if (sdx <= 0 | sdy <= 0 | sdz <=0)
    stop("Standard errors must be positive.")
  # check if VUS is appropriatly set:
  if (!missing(VUS) & (VUS < 1/6 | VUS > 1))
    stop("VUS must lie inside [1/6, 1].")

  dvus <- function(x) {
    A <- sdy/sdx
    B <- (mux-muy)/sdx
    C <- sdy/sdz
    D <- (muz-muy)/sdz
    pnorm(A*x - B)*pnorm(-C*x+D)*dnorm(x)
  }

  tempVus <- 1/6
  muy <- 0
  muz <- 0

  while(tempVus < VUS & tempVus < 0.999) {

    muy <- muy + step
    muz <- muz + 2*step
    tempVus <- integrate( dvus, -Inf, Inf)$value

  }

  return(data.frame(Par=c("mux", "muy", "muz", "VUS"),
                    Coeff = c(mux, muy, muz, tempVus)))
}

