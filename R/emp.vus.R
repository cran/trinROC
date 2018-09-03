#' Empirical VUS calculation
#'
#' This function computes the empirical Volume Under the Surface (VUS)
#' of three-class ROC data.
#'
#' @details This function computes the empirical VUS of three-class ROC data
#'   using the \code{expand.grid} function. It has been shown to be faster than
#'   computation using the \code{merge} function (\code{VUS.merge()}) or direct
#'   geometrical imlementation. The measurements can be input as seperate vectors
#'   \code{x, y, z} or as a data frame \code{dat}.
#' @param x,y,z Numeric vectors contaning the measurements from the healthy,
#'   intermediate and diseased class.
#' @param dat A data frame of the following structure: The first column
#'   represents a factor with three levels, containing the true class membership
#'   of each measurement. The levels are ordered according to the convention of
#'   higher values for more severe disease status. The second column contains
#'   all measurements obtained from Classifier.
#' @param old.version A logical to switch computation method to the old version,
#'   which is up to 50\% faster in computation (at N=50).
#' @return It returns the numeric \code{VUS} of the data.
#' @references Scurfield, B. K. (1996). Multiple-event forced-choice tasks in
#'   the theory of signal detectability. \emph{Journal of Mathematical
#'   Psychology} \bold{40.3}, 253–269.
#' @references Nakas CT and Yiannoutsos CT (2004) Ordered multiple-class roc
#'   analysis with continuous measurements. \emph{Statistics in Medicine} \bold{23}(22):
#'   3437–3449.
#' @export
#' @examples
#' data(krebs)
#' x1 <- with(krebs, cancer[trueClass=="healthy", 4])
#' y1 <- with(krebs, cancer[trueClass=="intermediate", 4])
#' z1 <- with(krebs, cancer[trueClass=="diseased", 4])
#'
#' emp.vus(x1, y1, z1)
#' # Alternatively:
#' emp.vus(dat = krebs[,c(1,4)])


emp.vus <- function(x, y, z, dat = NULL, old.version = TRUE) {

  # if data comes in a data.frame, unpack it:
  ## Important: levels symbolize the correctly ordered classes
  if (!is.null(dat)) {
    if (class(dat) != "data.frame" || class(dat[,1]) != "factor" | ncol(dat) <= 1)
      stop("Data should be organized as a data frame with the group index factor at
           the first column and marker measurements at the second column.")

    data.temp <- split(dat[,2], dat[,1], drop=FALSE)
    x <- data.temp[[1]]
    y <- data.temp[[2]]
    z <- data.temp[[3]]
  }

  # check for NA's:
  x.1 <- as.numeric(na.omit(x))
  y.1 <- as.numeric(na.omit(y))
  z.1 <- as.numeric(na.omit(z))

  if (!old.version) {
  dat1 <- expand.grid(x=x.1, y=y.1, z=z.1)

  sumup.grid <- function(x,y,z) {
    sum.vus <- 1 * (x<y & y<z) +
             0.5 * ((x<y & y==z ) | (x==y & y<z)) +
             1/6 * (x==y & y==z)
    VUS <- mean(sum.vus)
    return(VUS)
  }

  return(VUS = do.call(sumup.grid, dat1))
  } else {

    dat1 <- expand.grid(x=x.1, y=y.1, z=z.1, KEEP.OUT.ATTRS = FALSE)
    x.lt.y <- dat1$x<dat1$y
    y.lt.z <- dat1$y<dat1$z
    x.eq.y <- dat1$x==dat1$y
    y.eq.z <- dat1$y==dat1$z

    sum.vus <- 1 * (x.lt.y & y.lt.z) +
             0.5 * ((x.lt.y & y.eq.z ) | (x.eq.y & y.lt.z)) +
             1/6 * (x.eq.y & y.eq.z)
    return(VUS = mean(sum.vus))
  }
}
