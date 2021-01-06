#' Trinormal based ROC test
#'
#' A statistical test function to assess three-class ROC data. It is possible to
#' investigate a single classifier or make a comparison of two independent /
#' correlated classifiers.
#'
#' @details The trinormal ROC model is a parametric model in three-class ROC
#' analysis. It is based on normality in each of the trhee classes D_-
#' (healthy), D_0 (intermediate) and D_+ (diseased) with denoted distributions
#' \eqn{N(\mu_-,\sigma_-^2)}, \eqn{N(\mu_0,\sigma_0^2)} and
#' \eqn{N(\mu_+,\sigma_+^2)}. A classifier of a trinormal ROC model classifies
#' individuals into one of the three ordered classes based on two cut-off points
#' \eqn{c_- < c_+}. We define \eqn{t_-=F_-(c_-)} and \eqn{t_+
#' =1-F_+(c_+)=G_+(c_+)}. Now, the ROC surface can be written as
#'
#' \deqn{ROCs(t_-,t_+) = \Phi \left(\frac{\Phi^{-1} (1-t_+) +d}{c} \right) -
#' \Phi \left(\frac{\Phi^{-1} (t_-)+b}{a} \right)}
#'
#' whith parameters a, b, c and c given by \eqn{a =
#' \frac{\hat{\sigma}_0}{\hat{\sigma}_-}, b = \frac{ \hat{\mu}_- -
#' \hat{\mu}_0}{\hat{\sigma}_-}, c = \frac{\hat{\sigma}_0}{\hat{\sigma}_+}, d =
#' \frac{ \hat{\mu}_+ - \hat{\mu}_0}{\hat{\sigma}_+} }. It is a surface in the
#' unit cube that plots the probability of a measurement to get assigned to the
#' intermediate class as the two thresholds \eqn{c_-,c_+} are varying.
#'
#' Based on the reference standard, the trinormal based ROC test can be used
#' to assess the discriminatory power of such classifiers. It distinguishes
#' between single classifier assessment, where a classifier is compared to some
#' hypothetical distributions in the classes, and comparison between two
#' classifiers. The latter case tests for equality between the parameters a, b,
#' c and d of the ROC curves. The data can arise in a unpaired or paired
#' setting. If \code{paired} is \code{TRUE}, a correlation is introduced which
#' has to be taken into account. Therefore the sets of the two classifiers have
#' to have classwise equal size. The data can be input as the data frame
#' \code{dat} or as single vectors \code{x1, y1, z1, ...}.
#'
#' As the Chi-squared test is by definition a one-sided test, the variable
#' \code{alternative} cannot be specified in this test. For this 'goodness of
#' fit' test, we assume the parameters \eqn{a_1, \dots , d_1} and \eqn{a_2, \dots , d_2} to have a
#' pairwise equivalent normal distribution (in large sample sets).
# Thus, when
# the realized Chi-squared value is way out on the right tail of its
# distribution, it indicates a poor fit, and if it is far enough, relative to
# some pre-specified threshold, we might conclude that it is so poor that we
# don't believe the data are from that reference distribution. If we were to
# use the Chi-squared test as a two-sided test, we would also be worried if the
# statistic were too far into the left side of the chi-squared distribution.
# This would mean that we are worried the fit might be too good. This is simply
# not something we are typically worried about.
#'
#' @param dat  a data frame of the following structure: The first column
#'   represents a factor with three levels, containing the true class membership
#'   of each measurement. The levels are ordered according to the convention of
#'   higher values for more severe disease status. The second column contains
#'   all measurements obtained from Classifier 1 (in the case of single marker
#'   assessment). In the case of comparison of two markers, column three
#'   contains the measurementss from the Classifier.
#' @param x1,y1,z1  (non-empty) numeric vectors of data from the healthy,
#'   intermediate and diseased class from Classifier 1.
#' @param x2,y2,z2  numeric vectors of data from the healthy, intermediate and
#'   diseased class from Classifier 2.
#' @param paired a logical indicating whether data arose from a paired setting.
#'   If \code{TRUE}, each class must have equal sample size for both
#'   classifiers.
#' @importFrom stats complete.cases cor cov sd dnorm integrate na.omit nlm
#' pchisq pnorm qchisq qnorm var p.adjust shapiro.test
#' @importFrom grDevices cm.colors
#' @importFrom rgl open3d surface3d grid3d axes3d rgl.snapshot
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot  aes coord_flip facet_grid geom_boxplot
#'   geom_density geom_histogram geom_jitter labs scale_colour_manual
#'   scale_fill_manual stat_boxplot aes_string
#' @param conf.level confidence level of the interval. A numeric value between (0,1)
#'   yielding the significance level \eqn{\alpha=1-}\code{conf.level}.
#'
#' @return A list of class \code{"htest"} containing the following components:
#'   \item{statistic}{the value of the chi-squared statistic.}
#'   \item{parameter}{the degrees of freedom for the chi-squared statistic.}
#'   \item{p.value}{the p-value for the test.}
#'   \item{conf.int}{a confidence interval for the test.}
#'   \item{estimate}{a data frame containing the estimated VUS and parameters
#'   a, b, c and d from Classifier 1 and Classifier 2 (if specified).}
#'   \item{null.value}{a character expressing the null hypothesis.}
#'   \item{alternative}{a character string describing the alternative hypothesis.}
#'   \item{method}{a character string indicating what
#'   type of trinormal based ROC test was performed.}
#'   \item{data.name}{a character string giving the names of the data.}
#'   \item{CovMat}{the covariance matrix of the chi-squared statistic.}
#'   \item{Summary}{a data frame representing the number of NA's as well as
#'   the means and the standard deviations per class.}
#'
#' @references Noll, S., Furrer, R., Reiser, B. and Nakas, C. T. (2019).
#'   Inference in ROC surface analysis via a trinormal model-based testing approach.
#'   \emph{Stat}, \bold{8}(1), e249.
#' @seealso \code{\link{trinVUS.test}}, \code{\link{boot.test}}.
#'
#' @export
#' @examples
#' data(cancer)
#' data(krebs)
#'
#' # investigate a single marker:
#' trinROC.test(dat = cancer[,c(1,3)])
#' trinROC.test(dat = krebs[,c(1,5)])
#'
#' # result is equal to:
#' x1 <- with(cancer, cancer[trueClass=="healthy", 3])
#' y1 <- with(cancer, cancer[trueClass=="intermediate", 3])
#' z1 <- with(cancer, cancer[trueClass=="diseased", 3])
#' trinROC.test(x1, y1, z1)
#'
#' # comparison of marker 2 and 6:
#' trinROC.test(dat = cancer[,c(1,3,5)], paired = TRUE)
#' trinROC.test(dat = cancer[,c(1,3,5)], paired = FALSE)
#'
#' # result is equal to:
#' x2 <- with(cancer, cancer[trueClass=="healthy", 5])
#' y2 <- with(cancer, cancer[trueClass=="intermediate", 5])
#' z2 <- with(cancer, cancer[trueClass=="diseased", 5])
#' trinROC.test(x1, y1, z1, x2, y2, z2, paired = TRUE)



trinROC.test <- function(x1, y1, z1, x2 = 0, y2 = 0, z2 = 0, dat = NULL,
                        paired = FALSE, conf.level = 0.95) {

  alternative <- "two.sided" # Normally 3 poss. for test objects. chi-squar. has only 1,
  # namely "less", here "two.sided" is only used intern for correct print output.
  # to understand how output of a htest object is generated: stat:::print.htest

  # check if confidence level is appropriatly set:
  if (!missing(conf.level) & (length(conf.level) != 1 | !is.finite(conf.level) |
                               conf.level <= 0 | conf.level >= 1))
    stop("'conf.level' must be a single number between 0 and 1")


  # if data comes in a data.frame, unpack it:
  ## Important: levels symbolize the correctly ordered classes
  if (!is.null(dat)) {
    if ( !inherits(dat,"data.frame") | !inherits(dat[,1],"factor") | ncol(dat) <= 1)
      stop("Data should be organized as a data frame with the group index factor at
           the first column and marker measurements at the second and third column.")
    if (any(sapply(1 : (ncol(dat)-1), function(i) class(dat[, i+1])!="numeric")) ) {
      for (i in 1 : (ncol(dat)-1)) {
        if (class(dat[, i+1])!="numeric") dat[,(i+1)] <- as.numeric(dat[,(i+1)]) }
      warning("Some measurements were not numeric. Forced to numeric.")
    }

    data.temp <- split(dat[,2], dat[,1], drop=FALSE)
    x1 <- data.temp[[1]]
    y1 <- data.temp[[2]]
    z1 <- data.temp[[3]]
    if (ncol(dat) > 2) {
      data.temp <- split(dat[,3], dat[,1], drop=FALSE)
      x2 <- data.temp[[1]]
      y2 <- data.temp[[2]]
      z2 <- data.temp[[3]]
      }
  }

  # check for NA's:
  if (any(is.na(c(x1,x2)))) {
    warning("there are NA's in the healthy classes which are omitted.")
    if  (paired) {
      delx <- complete.cases(cbind(x1,x2))
      x.1 <- x1[delx]; x.2 <- x2[delx] }
    else {
      x.1 <- as.numeric(na.omit(x1)); x.2 <- as.numeric(na.omit(x2)) }
  } else { x.1 <- x1; x.2 <- x2 }

  if (any(is.na(c(y1,y2)))) {
    warning("there are NA's in the intermediate classes which are omitted.")
    if  (paired) {
      dely <- complete.cases(cbind(y1,y2))
      y.1 <- y1[dely]; y.2 <- y2[dely]
    } else {
      y.1 <- as.numeric(na.omit(y1)); y.2 <- as.numeric(na.omit(y2)) }
  } else { y.1 <- y1; y.2 <- y2 }

  if (any(is.na(c(z1,z2)))) {
    warning("there are NA's in the diseased classes which are omitted.")
    if (paired) {
      delz <- complete.cases(cbind(z1,z2))
      z.1 <- z1[delz]; z.2 <- z2[delz]
    } else {
    z.1 <- as.numeric(na.omit(z1)); z.2 <- as.numeric(na.omit(z2)) }
  } else { z.1 <- z1; z.2 <- z2 }


  # check if we are in single curve assessment or comparison of two classifiers:
  if (length(x.2)==1 && length(y.2)==1 && length(z.2)==1 ) {
    twocurves <- FALSE
    method <- "Trinormal based ROC test for single classifier assessment"
    dname <- if (!is.null(dat)) { paste(c(levels(dat[,1]), "of",
                                  names(dat)[2]), collapse = "  ")
      } else { paste(deparse(substitute(x1)), deparse(substitute(y1)),
                    "and", deparse(substitute(z1))) }
  } else {
  # in a case of two paired classifiers, check if data has equal size:
    if ((paired) &&
       (length(x.1)!=length(x.2) || length(y.1)!=length(y.2) ||
        length(z.1)!= length(z.2)) )
      stop('The two test sets do not have equal size.')

    twocurves <- TRUE
    method <- "Trinormal based ROC test for comparison of two independent classifiers"
    dname <- if (!is.null(dat)) { paste( c(levels(dat[,1]), "of",names(dat)[2],
              " and ",levels(dat[,1]),"of",names(dat)[3]), collapse = " ")
    } else { paste(deparse(substitute(x1)), deparse(substitute(y1)),
                   deparse(substitute(z1)), "and", deparse(substitute(x2)),
                   deparse(substitute(y2)), deparse(substitute(z2)) ) }
  }

  # compute means and standard errors of Classifier 1:
  mux1 <- mean(x.1); sdx1 <- SD(x.1) # healthy ind.
  muy1 <- mean(y.1); sdy1 <- SD(y.1) # early diseased ind.
  muz1 <- mean(z.1); sdz1 <- SD(z.1) # diseased ind.

  # sample size:
  nh1 <- length(x.1); nh2 <- length(x.2)
  n01 <- length(y.1); n02 <- length(y.2)
  nd1 <- length(z.1); nd2 <- length(z.2)
  if (nh1 < 2 || n01 < 2 || nd1 < 2 )
    stop("not enough observations (you need at least two observations in each class).")

  summary.dat <- data.frame(n = c(nh1,n01,nd1), mu=c(mux1,muy1,muz1),
                            sd=c(sdx1,sdy1,sdz1),
                            row.names=c("healthy", "intermediate", "diseased"))

  # check variability of the data:
  if (any(c(sdx1,sdy1,sdz1) < 10 * .Machine$double.eps * abs(c(mux1,muy1,muz1)) ) )
    stop("data of Classifier 1 is essentially constant")

  # check ordering of classes:
  # if (mux1 > muy1 || muy1 > muz1) {
  #   stop(paste(c("convention of ordered classes in Classifier 1 is violated:\n",
  #        "Mean x1 = ", round(mux1,3), ", Mean y1 = ", round(muy1,3),
  #        ", Mean z1 = ", round(muz1,3)), collapse = "") ) }
  mux2 <- 0; sdx2 <- 0
  muy2 <- 0; sdy2 <- 0
  muz2 <- 0; sdz2 <- 0
  if (twocurves) {

    # compute means and standard errors of Classifier 2:
    mux2 <- mean(x.2); sdx2 <- SD(x.2)
    muy2 <- mean(y.2); sdy2 <- SD(y.2)
    muz2 <- mean(z.2); sdz2 <- SD(z.2)
    if (any(c(sdx2,sdy2,sdz2) < 10 * .Machine$double.eps * abs(c(mux2,muy2,muz2)) ) )
      stop("data of Classifier 2 is essentially constant")
    summary.dat2 <-  data.frame(n = c(nh2,n02,nd2), mu=c(mux2,muy2,muz2),
                              sd=c(sdx2,sdy2,sdz2),
                              row.names=c("healthy", "intermediate", "diseased"))
    summary.dat <- list(Classifier1 = summary.dat, Classifier2 = summary.dat2)

    # if (mux2 > muy2 || muy2 > muz2)
    #   stop(paste(c("convention of ordered classes in Classifier 1 is violated:\n",
    #                "Mean x2 = ", round(mux2,3), ", Mean y2 = ", round(muy2,3),
    #                ", Mean z2 = ", round(muz2,3)), collapse = "") )
    }

  # compute correlation:
  rho.h <- 0; rho.0 <- 0; rho.d <- 0
  if (paired) {
    rho.h <- COV(x.1,x.2)/(sdx1*sdx2)  # correlation of healthy ind
    rho.0 <- COV(y.1,y.2)/(sdy1*sdy2)  # correlation of early disease ind
    rho.d <- COV(z.1,z.2)/(sdz1*sdz2)  # correlation of diseased ind
    method <- "Trinormal based ROC test for comparison of paired ROC surfaces" }

  # define the parameters A-D:
  A1 <- sdy1/sdx1
  B1 <- (mux1-muy1)/sdx1
  C1 <- sdy1/sdz1
  D1 <- (muz1-muy1)/sdz1
  A2 <- 1; B2 <- 0; C2 <- 1; D2 <- 0 # H0 when singe marker assessment
  # if (var.equal) {
  #    A2 <- A1; C2 <- C1 # consider different variances as H0
  #    method <- "Extended Metz-Kronman test for single curve assessment, adjusted sd's" }
  if (twocurves) {
    A2 <- sdy2/sdx2
    B2 <- (mux2-muy2)/sdx2
    C2 <- sdy2/sdz2
    D2 <- (muz2-muy2)/sdz2
  }

  # compute variances of A1-A2,..,D1-D2 with the DELTA METHOD:
  varA <- (A1^2/2)*(1/n01 + 1/nh1) +
          twocurves*( (A2^2/2)*(1/n02 + 1/nh2)) -
          paired*(A1*A2*rho.0^2/n01 + A1*A2*rho.h^2/nh1)
  varB <- A1^2/n01 + 1/nh1 + B1^2/(2*nh1) +
          twocurves*( A2^2/n02 + 1/nh2 + B2^2/(2*nh2)) -
          paired*2*(A1*A2*rho.0/n01 + rho.h/nh1 + B1*B2*rho.h^2/(2*nh1) )
  varC <- (C1^2/2)*(1/nd1 + 1/n01) +
          twocurves*( (C2^2/2)*(1/nd2 + 1/n02) ) -
          paired*(C1*C2*rho.0^2/n01 + C1*C2*rho.d^2/nd1)
  varD <- C1^2/n01 + 1/nd1 + D1^2/(2*nd1) +
          twocurves*( C2^2/n02 + 1/nd2 + D2^2/(2*nd2)) -
          paired*2*(C1*C2*rho.0/n01 + rho.d/nd1 + D1*D2*rho.d^2/(2*nd1) )

  # compute covariances of A1-A2,B1-B2,C1-C2,D1-D2 with DELTA METHOD:
  covAB <- A1*B1/(2*nh1) +
          twocurves*( A2*B2/(2*nh2)) -
          paired*(A1*B2*rho.h^2/(2*nh1) + A2*B1*rho.h^2/(2*nh1))
  covAC <- A1*C1/(2*n01) +
          twocurves*( A2*C2/(2*n02)) -
          paired*(A1*C2*rho.0^2/(2*n01) + A2*C1*rho.0^2/(2*n01))
  covAD <- 0
  covBC <- 0
  covBD <- A1*C1/n01 + twocurves*( A2*C2/n02) -
            paired*(A1*C2*rho.0/n01 + A2*C1*rho.0/n01)
  covCD <- C1*D1/(2*nd1) + twocurves*( C2*D2/(2*nd2)) -
            paired*(C2*D1*rho.d^2/(2*nd1) + C1*D2*rho.d^2/(2*nd1))
  covariance <- matrix(c(varA,covAB,covAC,covAD,
                         covAB,varB,covBC,covBD,
                         covAC,covBC,varC,covCD,
                         covAD,covBD,covCD,varD),4,4)

  # compute Chi-squared test statistic:
  if (!inherits(try(solve(covariance),silent=TRUE),"matrix")) {
    chi2 <- 0
    print("matrix is not invertible. Conclude complete alikeness.")
  } else {

  test <- t(c(A1-A2, B1-B2, C1-C2, D1-D2))%*% solve(covariance) %*%
            c(A1-A2, B1-B2, C1-C2, D1-D2)
  chi2 <- as.numeric(test)
  }
  # computation of VUS1 and VUS2:
  vus.1 <- 0; vus.2 <- 1/6

  dvus1 <- function(x) {  pnorm(A1*x - B1)*pnorm(-C1*x+D1)*dnorm(x) }
  vus.1 <- integrate( dvus1, -Inf, Inf)$value

  if (twocurves) {
    dvus2 <- function(x) {  pnorm(A2*x - B2)*pnorm(-C2*x+D2)*dnorm(x) }
    vus.2 <- integrate( dvus2, -Inf, Inf)$value }

  # P value of the Chi square test:
  df <- 4
  names(df)   <- "df"
  names(chi2) <- "Chi-Squared test"
  mekro.p     <- pchisq(chi2, df, lower.tail = FALSE)
  cint        <- c(0, qchisq(conf.level, df))
  attr(cint, "conf.level") <- conf.level

  param       <- as.data.frame(matrix( c(vus.1,A1,B1,C1,D1),nrow = 1))
  colnames(param) <- c("VUS", "a","b","c","d")
  rownames(param) <- "Classifier1:"
  if (twocurves) {
    param[2,] <- c(vus.2,A2,B2,C2,D2)
    rownames(param)[2] <- "Classifier2:" }

  # naming estimates:
  if (!is.null(dat)) {
    rownames(param)[1]   <- names(dat)[2]
    if (class(summary.dat)=="list") names(summary.dat)[1] <- names(dat)[2]
    if (ncol(dat) > 2) {
      rownames(param)[2] <- names(dat)[3]
      names(summary.dat)[2] <- names(dat)[3] }
  }

  null.value <- 0
  names(null.value) <- "a1-a2, b1-b1, c1-c2 and d1-d2"
  rval <- list(statistic = chi2, parameter = df, p.value = unname(mekro.p),
           conf.int = NULL, estimate =  param,
           null.value = null.value, alternative = alternative,
           method = method, data.name = dname,
           CovMat = covariance, Summary = summary.dat)

  class(rval) <- "htest" #hypothesis test element
  return(rval)

}
