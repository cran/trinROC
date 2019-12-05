#' Trinormal VUS test
#'
#' A statistical test function to assess three-class ROC data. It can be used
#' for assessment of a single classifier or comparison of two independent /
#' correlated classifiers, using the statistical test developed by Xiong et al.
#' (2007).
#'
#' @details
#   The trinormal ROC model is a parametric model in three-class ROC
#   analysis. It is based on normality in each of the trhee classes D_-
#   (healthy), D_0 (intermediate) and D_+ (diseased) with denoted distributions
#   \eqn{N(\mu_-,\sigma_-^2)}, \eqn{N(\mu_0,\sigma_0^2)} and
#   \eqn{N(\mu_+,\sigma_+^2)}. A classifier of a trinormal ROC model classifies
#   individuals into one of the three ordered classes based on two cut-off
#   points \eqn{c_- < c_+}. We define \eqn{t_-=F_-(c_-)} and \eqn{t_+
#   =1-F_+(c_+)=G_+(c_+)}. Now, the ROC surface can be written as
#
#   \deqn{ROCs(t_-,t_+) = \Phi \left(\frac{\Phi^{-1} (1-t_+) +d}{c} \right) -
#   \Phi \left(\frac{\Phi^{-1} (t_-)+b}{a} \right)}{ROCs(t_-,t_+) = \Phi((\Phi^{-1} (1-t_+) +d)/c ) -
#   \Phi( (\Phi^{-1} (t_-)+b)/a )}
#
#   whith parameters a, b, c and d given by \eqn{a =
#   \frac{\hat{\sigma}_0}{\hat{\sigma}_-}}{a =
#   \hat{\sigma}_0/\hat{\sigma}_-}, \eqn{b = \frac{ \hat{\mu}_- -
#   \hat{\mu}_0}{\hat{\sigma}_-}}{b = ( \hat{\mu}_- -
#   \hat{\mu}_0)/\hat{\sigma}_-}, \eqn{c = \frac{\hat{\sigma}_0}{\hat{\sigma}_+}}{c = \hat{\sigma}_0/\hat{\sigma}_+}, \eqn{d
#   = \frac{ \hat{\mu}_+ - \hat{\mu}_0}{\hat{\sigma}_+} }{d
#   = (\hat{\mu}_+ - \hat{\mu}_0)/\hat{\sigma}_+}. It is a surface in
#   the unit cube that plots the probability of a measurement to get assigned
#   to the intermediate class as the two thresholds \eqn{c_-,c_+} are varying.
#'   Based on the reference standard, this trinormal VUS test assesses the
#'   discriminatory power of classifiers by comparing the volumes under the ROC
#'   surfaces (VUS). It distinguishes between single classifier assessment,
#'   where a classifier is compared to the chance plane with VUS=1/6, and
#'   comparison between two classifiers. The latter case tests the equality
#'   between VUS_1 and VUS_2. The data can arise in a unpaired or paired
#'   setting. If \code{paired} is \code{TRUE}, a correlation is introduced which
#'   has to be taken into account. Therefore the sets of the two classifiers
#'   have to have classwise equal size. The data can be input as the data
#'   frame \code{dat} or as single vectors \code{x1, y1, z1, ...}.
#'
#' @param dat  a data frame of the following structure: The first column
#'   represents a factor with three levels, containing the true class membership
#'   of each measurement. The levels are ordered according to the convention of
#'   higher values for more severe disease status. The second column contains
#'   all measurements obtained from Classifier 1 (in the case of single marker
#'   assessment). In the case of comparison of two markers, column three
#'   contains the measurementss from the Classifier.
#' @param x1,y1,z1  non-empty numeric vectors of data from the healthy,
#'   intermediate and diseased class from Classifier 1.
#' @param x2,y2,z2  numeric vectors of data from the healthy, intermediate and
#'   diseased class from Classifier 2, only needed in a comparison of two
#'   classifiers.
#' @param paired logical; indicating whether data arose from a paired setting.
#'   If \code{TRUE}, each class must have equal sample size for both
#'   classifiers.
#' @param conf.level confidence level of the interval. A numeric value between (0,1)
#'   yielding the significance level \eqn{\alpha=1-\code{conf.level}}.
#' @param alternative character string specifying the alternative hypothesis,
#'   must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}. You can specify
#'   just the initial letter. For two sided test, notice \eqn{H0: Z = (VUS_1-VUS_2) /
#'   (Var(VUS_1)+Var(VUS_2)-2Cov(VUS_1,VUS_2))^{0.5}}.
#' @return A list of class \code{"htest"} containing the following components:
#'   \item{statistic}{the value of the Z-statistic.}
#'   \item{p.value}{the p-value for the test.}
#'   \item{conf.int}{a confidence interval for the test.}
#'   \item{estimate}{a data frame containing the
#'   estimated VUS from Classifier 1 and Classifier 2 (if specified).}
#'   \item{null.value}{a character expressing the null hypothesis.}
#'   \item{alternative}{a character string describing the alternative  hypothesis.}
#'   \item{method}{a character string indicating what type of trinormal VUS test was performed.}
#'   \item{data.name}{a character string giving the names of the data.} \item{Summary}{a data frame representing the
#'   number of NA's as well as the means and the standard deviations per class.}
#'   \item{Sigma}{the covariance matrix of the VUS.}
#'
#' @seealso \code{\link{trinROC.test}}, \code{\link{boot.test}}.
#'
#' @export
#' @references Xiong, C., Van Belle, G.  Miller J. P., Morris, J. C. (2006). Measuring and estimating
#'   diagnostic accuracy when there are three ordinal diagnostic groups.
#'   \emph{Statistics in Medicine}, \bold{25}(7), 1251–1273.
#'
#' @references Xiong, C., van Belle, G.,  Miller, J. P.,  Yan, Y.,  Gao, F., Yu, K., and Morris, J. C. (2007). A parametric comparison
#'   of diagnostic accuracy with three ordinal diagnostic groups.
#'   \emph{Biometrical Journal}, \bold{49}(5), 682–693. \url{http://doi.org/10.1002/bimj.200610359}.
#' @examples
#' data(cancer)
#' data(krebs)
#'
#' # investigate a single marker:
#' trinVUS.test(dat = cancer[,c(1,3)])
#' trinVUS.test(dat = krebs[,c(1,5)])
#'
#' # result is equal to:
#' x1 <- with(cancer, cancer[trueClass=="healthy", 3])
#' y1 <- with(cancer, cancer[trueClass=="intermediate", 3])
#' z1 <- with(cancer, cancer[trueClass=="diseased", 3])
#' trinVUS.test(x1, y1, z1)
#'
#' # comparison of marker 2 and 6:
#' trinVUS.test(dat = cancer[,c(1,3,5)], paired = TRUE)
#' trinVUS.test(dat = cancer[,c(1,3,5)], paired = FALSE)
#'
#' # result is equal to:
#' x2 <- with(cancer, cancer[trueClass=="healthy", 5])
#' y2 <- with(cancer, cancer[trueClass=="intermediate", 5])
#' z2 <- with(cancer, cancer[trueClass=="diseased", 5])
#' trinVUS.test(x1, y1, z1, x2, y2, z2, paired = TRUE)


trinVUS.test <- function(x1, y1, z1, x2 = 0, y2 = 0, z2 = 0, dat = NULL,
                       paired = FALSE, conf.level = 0.95,
                       alternative = c("two.sided", "less", "greater")) {

  alternative <- match.arg(alternative)
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
    method <- "Trinormal VUS test for single classifier assessment"
    dname <- if (!is.null(dat)) {
      paste(c(levels(dat[,1]), "of", names(dat)[2]), collapse = " ")
    } else { paste(deparse(substitute(x1)), deparse(substitute(y1)),
                   "and", deparse(substitute(z1))) }
  } else {
    # in a case of two paired classifiers, check if data has equal size:
    if (paired &&
        (length(x.1)!=length(x.2) || length(y.1)!=length(y.2) ||
         length(z.1)!= length(z.2)) )
      stop('The two test sets do not have equal size.')

    twocurves <- TRUE
    method <- "Trinormal VUS test for comparison of two independent classifiers"
    dname <- if (!is.null(dat)) { paste( c(levels(dat[,1]), "of",names(dat)[2],
                   " and ",levels(dat[,1]),"of",names(dat)[3]), collapse = " ")
    } else { paste(deparse(substitute(x1)), deparse(substitute(y1)),
                   deparse(substitute(z1)), "and", deparse(substitute(x2)),
                   deparse(substitute(y2)), deparse(substitute(z2)) ) }
  }

  # compute means and standard errors of Classifier 1:
  mux1<-mean(x.1); sdx1<-SD(x.1) # healthy ind.
  muy1<-mean(y.1); sdy1<-SD(y.1) # early diseased ind.
  muz1<-mean(z.1); sdz1<-SD(z.1) # diseased ind.

  # sample size:
  nh1 <- length(x.1); nh2 <- length(x.2)
  n01 <- length(y.1); n02 <- length(y.2)
  nd1 <- length(z.1); nd2 <- length(z.2)
  if (nh1 < 2 || n01 < 2 || nd1 < 2)
    stop("not enough observations (you need at least two observations in each class).")

  summary.dat <- data.frame(n = c(nh1,n01,nd1), mu=c(mux1,muy1,muz1),
                            sd=c(sdx1,sdy1,sdz1),
                            row.names=c("healthy", "intermediate", "diseased"))

  # check variability of the data:
  if (any(c(sdx1,sdy1,sdz1) < 10 * .Machine$double.eps * abs(c(mux1,muy1,muz1)) ) )
    stop("data of Classifier 1 is essentially constant")

  mux2<-0; sdx2<-0
  muy2<-0; sdy2<-0
  muz2<-0; sdz2<-0
  if (twocurves) {
    # compute means and standard errors of Classifier 2:
    mux2<-mean(x.2); sdx2<-SD(x.2)
    muy2<-mean(y.2); sdy2<-SD(y.2)
    muz2<-mean(z.2); sdz2<-SD(z.2)
      # check variability of the data:
    if (any(c(sdx2,sdy2,sdz2) < 10 * .Machine$double.eps * abs(c(mux2,muy2,muz2)) ) )
      stop("data of Classifier 2 is essentially constant")
    summary.dat2<- data.frame(n = c(nh2,n02,nd2), mu=c(mux2,muy2,muz2),
                              sd=c(sdx2,sdy2,sdz2),
                              row.names=c("healthy", "intermediate", "diseased"))
    summary.dat <- list(Classifier1 = summary.dat, Classifier2 = summary.dat2)
  }

  # compute correlation:
  rho.h<-0; rho.0<-0; rho.d<-0
  if (paired) {
    rho.h <- COV(x.1,x.2)/(sdx1*sdx2)  # correlation of healthy ind
    rho.0 <- COV(y.1,y.2)/(sdy1*sdy2)  # correlation of early disease ind
    rho.d <- COV(z.1,z.2)/(sdz1*sdz2)  # correlation of diseased ind
    method <- "Trinormal VUS test for comparison of paired ROC data" }

  # define the parameters A-D:
  A1 <- sdy1/sdx1
  B1 <- (mux1-muy1)/sdx1
  C1 <- sdy1/sdz1
  D1 <- (muz1-muy1)/sdz1
  A2 <- 1; B2<-0; C2<-1; D2<-0 # H0 when singe marker assessment
  if (twocurves) {
    A2 <- sdy2/sdx2
    B2 <- (mux2-muy2)/sdx2
    C2 <- sdy2/sdz2
    D2 <- (muz2-muy2)/sdz2
  }

  # (co)variances of parameters for Var(VUS1):
  var.A1 <- (A1^2/2)*(1/n01 + 1/nh1)
  var.B1 <- A1^2/n01 + 1/nh1 + B1^2/(2*nh1)
  var.C1 <- (C1^2/2)*(1/nd1 + 1/n01)
  var.D1 <- C1^2/n01 + 1/nd1 + D1^2/(2*nd1)
  cov.A1B1 <- A1*B1/(2*nh1)
  cov.A1C1 <- A1*C1/(2*n01)
  cov.B1D1 <- A1*C1/n01
  cov.C1D1 <- C1*D1/(2*nd1)

  if (twocurves) {
    # (co)variances of parameters for Var(VUS2):
    var.A2 <- (A2^2/2)*(1/n02 + 1/nh2)
    var.B2 <- A2^2/n02 + 1/nh2 + B2^2/(2*nh2)
    var.C2 <- (C2^2/2)*(1/nd2 + 1/n02)
    var.D2 <- C2^2/n02 + 1/nd2 + D2^2/(2*nd2)
    cov.A2B2 <- A2*B2/(2*nh2)
    cov.A2C2 <- A2*C2/(2*n02)
    cov.B2D2 <- A2*C2/n02
    cov.C2D2 <- C2*D2/(2*nd2)

    if (paired) {
      # covariances of parameters for Cov(VUS1,VUS2):
      cov.A1A2 <- A1*A2*rho.0^2/(2*n01) + A1*A2*rho.h^2/(2*nh1)
      cov.B1B2 <- A1*A2*rho.0/n01 + rho.h/nh1 + B1*B2*rho.h^2/(2*nh1)
      cov.C1C2 <- C1*C2*rho.0^2/(2*n01) + C1*C2*rho.d^2/(2*nd1)
      cov.D1D2 <- C1*C2*rho.0/n01 + rho.d/nd1 + D1*D2*rho.d^2/(2*nd1)
      cov.A1B2 <- A1*B2*rho.h^2/(2*nh1)
      cov.A2B1 <- A2*B1*rho.h^2/(2*nh1)
      cov.A1C2 <- A1*C2*rho.0^2/(2*n01)
      cov.A2C1 <- A2*C1*rho.0^2/(2*n01)
      cov.B1D2 <- A1*C2*rho.0/n01
      cov.B2D1 <- A2*C1*rho.0/n01
      cov.C1D2 <- C1*D2*rho.d^2/(2*nd1)
      cov.C2D1 <- C2*D1*rho.d^2/(2*nd1)

      } }

  # Define partial derivatives of A,B,C,D of the VUS:
  Der.A1 <- function(x) { dnorm(A1*x - B1)*pnorm(-C1*x+D1)*dnorm(x)*x }
  DA1    <- integrate( Der.A1, -Inf, Inf)$value
  Der.B1 <- function(x) { - dnorm(A1*x - B1)*pnorm(-C1*x+D1)*dnorm(x) }
  DB1    <- integrate( Der.B1, -Inf, Inf)$value
  Der.C1 <- function(x) { - pnorm(A1*x - B1)*dnorm(-C1*x+D1)*dnorm(x)*x }
  DC1    <- integrate( Der.C1, -Inf, Inf)$value
  Der.D1 <- function(x) { pnorm(A1*x - B1)*dnorm(-C1*x+D1)*dnorm(x) }
  DD1    <- integrate( Der.D1, -Inf, Inf)$value

  if(twocurves) {
  Der.A2 <- function(x) { dnorm(A2*x - B2)*pnorm(-C2*x+D2)*dnorm(x)*x }
  DA2    <- integrate( Der.A2, -Inf, Inf)$value
  Der.B2 <- function(x) { - dnorm(A2*x - B2)*pnorm(-C2*x+D2)*dnorm(x) }
  DB2    <- integrate( Der.B2, -Inf, Inf)$value
  Der.C2 <- function(x) { - pnorm(A2*x - B2)*dnorm(-C2*x+D2)*dnorm(x)*x }
  DC2    <- integrate( Der.C2, -Inf, Inf)$value
  Der.D2 <- function(x) { pnorm(A2*x - B2)*dnorm(-C2*x+D2)*dnorm(x) }
  DD2    <- integrate( Der.D2, -Inf, Inf)$value
  }

  # computation of Var(VUS1) and Var(VUS2):
  Var.vus1 <- 0; Var.vus2 <- 0; cov.vus1vus2 <- 0

  Var.vus1 <- DA1^2*var.A1 + DB1^2*var.B1 + DC1^2*var.C1 + DD1^2*var.D1 +
              2*(DA1*DB1*cov.A1B1 + DA1*DC1*cov.A1C1 + DB1*DD1*cov.B1D1 +
                   DC1*DD1*cov.C1D1)
  if (twocurves) {
  Var.vus2 <- DA2^2*var.A2 + DB2^2*var.B2 + DC2^2*var.C2 + DD2^2*var.D2 +
    2*(DA2*DB2*cov.A2B2 + DA2*DC2*cov.A2C2 + DB2*DD2*cov.B2D2 + DC2*DD2*cov.C2D2)
  if (paired) {
  cov.vus1vus2 <- DA1*DA2*cov.A1A2 + DA1*DB2*cov.A1B2 + DA1*DC2*cov.A1C2 +
                  DB1*DA2*cov.A2B1 + DB1*DB2*cov.B1B2 + DB1*DD2*cov.B1D2 +
                  DC1*DA2*cov.A2C1 + DC1*DC2*cov.C1C2 + DC1*DD2*cov.C1D2 +
                  DD1*DB2*cov.B2D1 + DD1*DC2*cov.C2D1 + DD1*DD2*cov.D1D2
  } }


  # computation of VUS1 and VUS2:
  vus.1 <- 0; vus.2 <- 1/6

  dvus1 <- function(x) {  pnorm(A1*x - B1)*pnorm(-C1*x+D1)*dnorm(x) }
  vus.1 <- integrate( dvus1, -Inf, Inf)$value
  if (twocurves) {
  dvus2 <- function(x) {  pnorm(A2*x - B2)*pnorm(-C2*x+D2)*dnorm(x) }
  vus.2 <- integrate( dvus2, -Inf, Inf)$value }

  # Compute test value and p-value:
  x.test  <- (vus.1-vus.2) / (Var.vus1 + Var.vus2 - 2*cov.vus1vus2)^0.5
  names(x.test) <- "Z-stat"

  if (alternative == "two.sided") {
    p.value <- 2 * pnorm(-abs(x.test), lower.tail = TRUE)
    alpha <- 1 - conf.level
    cint <- qnorm(1 - alpha/2) * (Var.vus1 + Var.vus2 - 2*cov.vus1vus2)^0.5
    cint <- (vus.1-vus.2) + c(-cint, cint)
  } else if (alternative == "less") {
    p.value <- pnorm(x.test, lower.tail = TRUE)
    cint <- c(-Inf, (vus.1-vus.2) + qnorm(conf.level) *
                (Var.vus1 + Var.vus2 - 2*cov.vus1vus2)^0.5)
  } else if (alternative == "greater") {
    p.value <- pnorm(x.test, lower.tail = FALSE)
    cint <- c((vus.1-vus.2) - qnorm(conf.level) *
                (Var.vus1 + Var.vus2 - 2*cov.vus1vus2)^0.5, Inf)
  }
  attr(cint, "conf.level") <- conf.level

  if (is.na(p.value))
    {p.value<-1; print("Denominator is zero. Conclude complete alikeness.")}

  estimate <- vus.1
  names(estimate) <- "VUS of Classifier 1"
  if (twocurves) {
    estimate <- c(vus.1, vus.2)
    names(estimate) <- c("VUS of Classifier 1", "VUS of Classifier 2")}

  sigma <- diag(c(Var.vus1, Var.vus2))
  sigma[1,2] <- sigma[2,1] <- cov.vus1vus2

  # naming estimates:
  if (!is.null(dat)) {
    names(estimate)[1] <- paste("VUS of", names(dat)[2])
    if (class(summary.dat)=="list") names(summary.dat)[1] <- names(dat)[2]
    if (ncol(dat) > 2) {
      names(estimate)[2] <- paste("VUS of", names(dat)[3])
      names(summary.dat)[2] <- names(dat)[3] }
  }

  null.value <- 0
  names(null.value) <- "Difference in VUS"
  rval <- list(statistic = x.test, p.value = unname(p.value),
               conf.int = NULL, estimate = estimate,
               null.value=null.value, alternative = alternative,
               method = method, data.name = dname,
               Summary = summary.dat, Sigma = sigma)
  class(rval) <- "htest" #hypothesis test element
  return(rval)
}


