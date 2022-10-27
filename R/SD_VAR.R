# #' Package standard error
# #'
# #' A function used to set the MLE or the sample variance case.
# #' @param x A vector to compute the standard error.
# #' @param ... Further arguments.
# #' @export

SD <- function(x) {
  n <- length(x)
  if (getOption("trinROC.MLE")) {
    sd(x)*sqrt((n-1)/n)
  } else {
    sd(x)
  }
}

# #' Package variance
# #'
# #' A function used to set the MLE or the sample variance case
# #' @param x A vector to compute the variance.
# #' @param ... Further arguments.
# #' @export

VAR <- function(x) {
  n <- length(x)
  if (getOption("trinROC.MLE")) {
    var(x)*(n-1)/n
  } else {
    var(x)
  }
}


# #' Package covariance
# #'
# #' A function used to set the MLE or the sample covariance case
# #' @param x A vector to compute the covariance.
# #' @param y A vector of same length to compute the covariance.
# #' @param ... Further arguments.
# #' @export

COV <- function(x, y=NULL) {
  n <- length(x)
  if (getOption("trinROC.MLE")) {
    cov(x,y)*(n)/(n-1)
  } else {
    cov(x,y)
  }
}
