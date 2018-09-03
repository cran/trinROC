#' Synthetic small data set to investigate three-class ROC data.
#'
#' A dataset containing randomly generated measurements from three diagnostic classes as they may
#' arise in a cancer investigation. For illustration, this dataset has been
#' chosen to be smaller than the data set \code{\link{cancer}}.
#'
#' @format A data frame with 50 rows and 5 variables (4 classifiers):
#' \describe{
#' \item{trueClass}{A factor, indicating the class membership of the
#' individuals.}
#' \item{Fac1, Fac2, Fac3, Fac4}{Measurements obtained from the patients that underwent the
#' clinical study.}}
#' @keywords datasets
"krebs"
