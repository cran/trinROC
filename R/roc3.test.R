#' Statistical test function for computing multiple tests on three-class ROC data
#'
#' A statistical test function that assesses three-class ROC data with the
#' trinormal based ROC test, the trinormal VUS test and the Bootstrap test.
#'
#' @details For the preliminary assessment of a classifier, different
#'   statistical tests have been proposed in the literature. This function can
#'   be used for either comparison of single classifiers to a null hypothesis of
#'   being not better than a random allocation function or comparison of two
#'   classifiers under the null hypothesis of having equal discriminatory power.
#'   Depending on the specification of the user, (s)he can apply the trinormal
#'   based ROC test (LINK), the test developed by Xiong et. al. or the Bootstrap
#'   test or any combination of these tests. More information of the specific
#'   tests can be obtained by calling \code{?functionname}. If more than two
#'   markers are present, a pairwise comparison between each marker is realized.
#' @param dat  A data frame of the following structure: The first column
#'   represents a factor with three levels, containing the true class membership
#'   of each measurement. The levels are ordered according to the convention of
#'   higher values for more severe disease status.
#'
#' @param type A character, specifying which tests are applied to \code{dat}.
#'   \code{"ROC"} implies the trinormal based ROC test, \code{"VUS"} the trinormal
#'   VUS test and \code{"Bootstrap"} the Bootstrap test.
#' @param paired A logical indicating whether data arose from a paired setting.
#' If data is paired, each class must have equal sample size for both classifiers.
#' @param conf.level confidence level of the interval. A numeric value between (0,1)
#'   yielding the significance level \eqn{\alpha=1-\code{conf.level}}.
#' @param n.boot An integer incicating the number of Bootstrap replicates sampled
#'  to obtain the variance of the VUS. Default is 1000.
#' @param p.adjust A logical, indicating whether a FDR adjustment
#'    should be applied to the p-values. Default is \code{FALSE}.
#' @section Note: If \code{type = "Bootstrap"}, the Bootstrap test is evaluated. This
#'  may take some time, especially with sample sizes > 100.
#' @return A list with components:
#'   \item{Overview}{a data frame with number of columns according to number of
#'   markers. Rows contain the following information about the makers:
#'   \enumerate{
#'     \item Index according to smallest VUS
#'     \item VUS
#'     \item P-values of statistical test specified by \code{type}
#'     \item Number of NA's
#'     }}
#'   \item{O.orig}{the unsorted \code{Overview.}}
#'   \item{P.values}{a list, containing the upper triangular matrices of the optionally adjusted
#'   p-values of the statistical tests chosen by \code{type}.}
#'   \item{Test.Values}{a list, containing the upper triangular matrices of the
#'   test values of the statistical tests chosen by \code{type}.}
#' @examples
#' data(krebs)
#' roc3.test(krebs, type = c("ROC", "VUS"), paired = TRUE)[c("Overview","P.values")]
#' @export

roc3.test <- function(dat, type = c("ROC","VUS","Bootstrap"),
                      paired = FALSE, conf.level = 0.95, n.boot = 1000,
                      p.adjust = FALSE) {
#                      alternative = c("two.sided", "less", "greater"))

  # checking arguments:
  type <- match.arg(type, several.ok = TRUE)

  if (paired & ncol(dat) <= 2)
    stop("There exists no paired method for single classifier assessment")

  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level <= 0 || conf.level >= 1))
    stop("'conf.level' must be a single number between 0 and 1")

  if ( !inherits(dat,"data.frame") | !inherits(dat[,1],"factor") | ncol(dat) <= 1)
      stop("Data should be organized as a data frame with the group index factor at
         the first column and marker measurements at the second and third column.")
  # check for classes of column vectors of dat, should all be numeric:
  if (any(sapply(1 : (ncol(dat)-1), function(i) class(dat[, i+1])!="numeric")) ) {
    not.num <- sapply(1 : (ncol(dat)-1), function(i) class(dat[, i+1])!="numeric")
    dat[, not.num==FALSE] <- as.numeric(dat[, not.num==FALSE])
    warning("Some measurements were not numeric. Forced to numeric.")
  }

  # count nr of classes of dat:
  nclass     <- ncol(dat)-1
  classnames <- names(dat)[-1]
  clns  <- substr(classnames, 1, 4)
  # indicater vector of the three applied tests:
  counttype  <- c(any(type=="ROC"), any(type=="VUS"),
                  any(type=="Bootstrap"))

  # preparing output:
  Overview <- data.frame(matrix(NA, (4+length(type)), nclass))
  p.values <- sapply(1:sum(counttype), function(i) {
    temp <- type[counttype]
    paste("p.value",temp[i], "test") })

  row.names(Overview) <- c("Charts","Emp. VUS", "Trin. VUS", p.values,"Nr. of NA's")
  colnames(Overview) <- classnames

  # fill in empirical VUS:
  Overview[2,] <- sapply(1:nclass, function(i) emp.vus(dat=dat[, c(1,(i+1))]) )
  Overview[3,] <- sapply(1:nclass, function(i) trinVUS.test(dat=dat[, c(1,(i+1))])$estimate[1])
  Overview[(4+sum(counttype)),] <- sapply(1:nclass, function(i) sum(is.na(dat[,(i+1)])))

  # list of p-value comparisons and test values:
  pv.comp <- list()
  test.comp <- list()

  # computing the tests:
  if (any(type == "ROC")) {
    # trinROCMat  <- matrix(list(NA), nclass,nclass)
    # trinROCMat  <- rep(list(list()), nclass*nclass)
    if ( nclass == 1) {
      temptrinROC <- suppressWarnings(
                   trinROC.test(dat = dat, conf.level = conf.level))
      # if "roc", third position is taken:
      Overview[4, 1] <- temptrinROC$p.value

    } else {
      trinROCPval <- data.frame(matrix(NA, nclass,nclass))
      row.names(trinROCPval) <- classnames
      colnames(trinROCPval)  <- classnames
      trinROCTest <- trinROCPval
      # fill p.values for single tests.
      Overview[4,] <- sapply(1:nclass, function(i) { suppressWarnings(
        trinROC.test(dat = dat[, c(1, (i+1))], conf.level=conf.level)$p.value)
      })

      # fill results from pairwise evaluation of classifiers
      for (i in 1:(nclass-1)) {
        for (j in (i+1):nclass) {
          #trinROCMat[[i+nclass*(j-1)]] <- trinROC.test(dat = dat[, c(1,(i+1),(j+1))],
          #                        conf.level = conf.level, paired = paired)
          #names(trinROCMat)[i+nclass*(j-1)] <- paste0(clns[i],"/",clns[j])
          temptrinROC <- suppressWarnings(
                        trinROC.test(dat = dat[, c(1,(i+1),(j+1))],
                        conf.level = conf.level, paired = paired) )
          trinROCPval[i,j] <- temptrinROC$p.value
          trinROCTest[i,j] <- temptrinROC$statistic
        }}

      # should p adjustment be applied:
      if (p.adjust) {
        Pval <- p.adjust(na.omit(as.vector(as.matrix(trinROCPval))), method = "fdr")
        trinROCPval[upper.tri(trinROCPval, diag=F)] <- Pval
      }

      pv.comp$trinROC   <- trinROCPval
      test.comp$trinROC <- trinROCTest
    }
  }

  if (any(type == "VUS")) {
    if ( nclass == 1) {
      temptrinVUS <- suppressWarnings(
                    trinVUS.test(dat = dat, conf.level = conf.level))
      Overview[(3+sum(counttype[1:2])),1] <- temptrinVUS$p.value

    } else {
      trinVUSPval <- data.frame(matrix(NA, nclass,nclass))
      row.names(trinVUSPval) <- classnames
      colnames(trinVUSPval)  <- classnames
      trinVUSTest <- trinVUSPval
      # fill p.values for single tests.
      Overview[(3+sum(counttype[1:2])),] <- sapply(1:nclass, function(i) {
        suppressWarnings(
          trinVUS.test(dat = dat[, c(1, (i+1))], conf.level=conf.level)$p.value)
      })

      # fill results from pairwise evaluation of classifiers
      for (i in 1:(nclass-1)) {
        for (j in (i+1):nclass) {
          temptrinVUS <- suppressWarnings(
                          trinVUS.test(dat = dat[, c(1,(i+1),(j+1))],
                          conf.level = conf.level, paired = paired) )
          trinVUSPval[i,j] <- temptrinVUS$p.value
          trinVUSTest[i,j] <- temptrinVUS$statistic
        }}

      # should p adjustment be applied:
      if (p.adjust) {
        Pval <- p.adjust(na.omit(as.vector(as.matrix(trinVUSPval))), method = "fdr")
        trinVUSPval[upper.tri(trinVUSPval, diag=F)] <- Pval
      }

      pv.comp$trinVUS   <- trinVUSPval
      test.comp$trinVUS <- trinVUSTest
    }
  }

  if (any(type == "Bootstrap")) {
    if ( nclass == 1) {
      tempBoot <- suppressWarnings(
                    boot.test(dat = dat, conf.level = conf.level) )
      Overview[(3+sum(counttype)),1] <- tempBoot$p.value

    } else {
      BootPval <- data.frame(matrix(NA, nclass,nclass))
      row.names(BootPval) <- classnames
      colnames(BootPval)  <- classnames
      BootTest <- BootPval
      # fill p.values for single tests.
      Overview[(3+sum(counttype)),] <- sapply(1:nclass, function(i) {
        suppressWarnings(
          boot.test(dat = dat[, c(1, (i+1))], conf.level=conf.level)$p.value )
      })

      # fill results from pairwise evaluation of classifiers
      for (i in 1:(nclass-1)) {
        for (j in (i+1):nclass) {
          #print(paste0("Step: ", i*nclass-1)j))
          tempBoot <- suppressWarnings (
                              boot.test(dat = dat[, c(1,(i+1),(j+1))],
                              conf.level = conf.level, paired = paired) )
          BootPval[i,j] <- tempBoot$p.value
          BootTest[i,j] <- tempBoot$statistic
        }}

      # should p adjustment be applied:
      if (p.adjust) {
        Pval <- p.adjust(na.omit(as.vector(as.matrix(BootPval))), method = "fdr")
        BootPval[upper.tri(BootPval, diag=F)] <- Pval
      }

      pv.comp$Boot   <- BootPval
      test.comp$Boot <- BootTest
    }
  }

  OverviewSorted <- Overview[order(Overview[2,], decreasing = T)]
  OverviewSorted[1, ] <- 1:nclass
  OverviewSorted <- round(OverviewSorted, 4)


  res <- list(Overview = OverviewSorted, O.orig = Overview[-1,], P.values = pv.comp,
              Test.Values = test.comp)
  return(res)

}
