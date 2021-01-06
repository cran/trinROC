#' Bootstrap test for three-class ROC data
#'
#' A statistical test function to assess three-class ROC data. It can be used
#' for assessment of a single classifier or comparison of two independent /
#' correlated classifiers, using the Bootstrap test.
#'
#' @details Based on the reference standard, the Bootstrap test assesses the
#'   discriminatory power of classifiers by comparing the volumes under the ROC
#'   surfaces (VUS). It distinguishes between single classifier assessment,
#'   where a classifier is compared to the chance plane with VUS=1/6, and
#'   comparison between two classifiers. The latter case tests the equality
#'   between VUS_1 and VUS_2. The data can arise in a unpaired or paired
#'   setting. If \code{paired} is \code{TRUE}, a correlation is introduced which
#'   has to be taken into account. Therefore the sets of the two classifiers
#'   have to have classwise equal size. The data can be input as the data
#'   frame \code{dat} or as single vectors \code{x1, y1, z1, ...}. The
#'   implemented methods to evaluate the \code{VUS} and \code{var(VUS),
#'   cov(vus.1,vus.2)} are based on the empirical model assumptions and
#'   resampling techniques. This means, there are no underlying distributions
#'   assumed in any of the classes.
#'
#' @param dat  A data frame of the following structure: The first column
#'   represents a factor with three levels, containing the true class membership
#'   of each measurement. The levels are ordered according to the convention of
#'   higher values for more severe disease status. The second column contains
#'   all measurements obtained from Classifier 1 (in the case of single marker
#'   assessment). In the case of comparison of two markers, column three
#'   contains the measurementss from the Classifier.
#' @param x1,y1,z1  Non-empty numeric vectors of data from the healthy,
#'   intermediate and diseased class from Classifier 1.
#' @param x2,y2,z2  Numeric vectors of data from the healthy, intermediate and
#'   diseased class from Classifier 2, only needed in a comparison of two
#'   classifiers.
#' @param paired A logical indicating whether data arose from a paired setting.
#'   If \code{TRUE}, each class must have equal sample size for both
#'   classifiers.
#' @param n.boot An integer incicating the number of bootstrap replicates
#'   sampled to obtain the variance of the VUS. Default is 1000.
#' @param conf.level confidence level of the interval. A numeric value between (0,1)
#'   yielding the significance level \eqn{\alpha=1-\code{conf.level}}.
#' @param alternative character string specifying the alternative hypothesis,
#'   must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}. You can specify
#'   just the initial letter. For two sided test, notice \eqn{H0: Z = (VUS_1-VUS_2) /
#'   (Var(VUS_1)+Var(VUS_2)-2Cov(VUS_1,VUS_2))^{0.5}}.
#' @return A list of class \code{"htest"} containing the following components:
#'   \item{statistic}{the value of the Z-statistic.}
#'   \item{p.value}{the p-value for the test.}
# #'   \item{conf.int}{a confidence interval for the test.}
#'   \item{estimate}{a data frame containing the
#'   estimated parameters from Classifier 1 and Classifier 2 (if specified).}
#'   \item{null.value}{a character expressing the null hypothesis.}
#'   \item{alternative}{a character string describing the alternative
#'   hypothesis.}
#'   \item{method}{a character string indicating what type of extended
#'   Metz--Kronman test was performed.}
#'   \item{data.name}{a character string giving the names of the data.}
#'   \item{Summary}{A data frame representing the number of NA's as well as the
#'   means and the standard deviations per class.}
#'   \item{Sigma}{The covariance matrix of the VUS.}
#'
#' @seealso \code{\link{trinROC.test}}, \code{\link{trinVUS.test}}.
#' @export
#' @references Nakas, C. T. and C. T. Yiannoutsos (2004). Ordered multiple-class
#'  ROC analysis with continuous measurements. \emph{Statistics in
#'  Medicine}, \bold{23}(22), 3437–3449.
#' @examples
#' data(cancer)
#' data(krebs)
#'
#' # investigate a single marker:
#' boot.test(dat = krebs[,c(1,2)], n.boot=500)
#'
#' # result is equal to:
#' x1 <- with(krebs, krebs[trueClass=="healthy", 2])
#' y1 <- with(krebs, krebs[trueClass=="intermediate", 2])
#' z1 <- with(krebs, krebs[trueClass=="diseased", 2])
#' \donttest{boot.test(x1, y1, z1, n.boot=500) }
#'
#' # comparison of marker 2 and 6:
#' \donttest{boot.test(dat = krebs[,c(1,2,5)], paired = TRUE) }
#'
#' # result is equal to:
#' x2 <- with(krebs, krebs[trueClass=="healthy", 5])
#' y2 <- with(krebs, krebs[trueClass=="intermediate", 5])
#' z2 <- with(krebs, krebs[trueClass=="diseased", 5])
#' \donttest{boot.test(x1, y1, z1, x2, y2, z2, paired = TRUE) }


boot.test <- function(x1, y1, z1, x2 = 0, y2 = 0, z2 = 0, dat = NULL,
                      paired = FALSE, n.boot = 1000, conf.level = 0.95,
                      alternative = c("two.sided", "less", "greater")) {

  alternative <- match.arg(alternative)
  # check if confidence level is appropriatly set:
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level <= 0 || conf.level >= 1))
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
  if (length(x.2)==1 & length(y.2)==1 & length(z.2)==1 ) {
    twocurves <- FALSE
    method <- "Bootstrap test for single classifier assessment"
    dname <- if (!is.null(dat)) { paste(c(levels(dat[,1]), "of",
                                          names(dat)[2]), collapse = " ")
    } else { paste(deparse(substitute(x1)), deparse(substitute(y1)),
                   "and", deparse(substitute(z1))) }
  } else {
    # in a case of two paired classifiers, check if data has equal size:
    if ((paired) &
        (length(x.1)!=length(x.2) | length(y.1)!=length(y.2) |
         length(z.1)!= length(z.2)) )
      stop('The two test sets do not have equal size.')

    twocurves <- TRUE
    method <- "Bootstrap test for comparison of two independent classifiers"
    dname <- if (!is.null(dat)) { paste( c(levels(dat[,1]), "of",names(dat)[2],
                                           " and ",levels(dat[,1]),"of",names(dat)[3]), collapse = " ")
    } else { paste(deparse(substitute(x1)), deparse(substitute(y1)),
                   deparse(substitute(z1)), "and", deparse(substitute(x2)),
                   deparse(substitute(y2)), deparse(substitute(z2)) ) }
  }


  nh1 <- length(x.1)
  n01 <- length(y.1)
  nd1 <- length(z.1)
  if (twocurves) {
    nh2 <- length(x.2)
    n02 <- length(y.2)
    nd2 <- length(z.2) }

  mux1<-mean(x.1); sdx1<-SD(x.1) # healthy ind.
  muy1<-mean(y.1); sdy1<-SD(y.1) # early diseased ind.
  muz1<-mean(z.1); sdz1<-SD(z.1) # diseased ind.

  summary.dat <- data.frame(n = c(nh1,n01,nd1), mu=c(mux1,muy1,muz1),
                            sd=c(sdx1,sdy1,sdz1), row.names=c("D-", "D0", "D+"))
  # check variability of the data:
  if (any(c(sdx1,sdy1,sdz1) < 10 * .Machine$double.eps * abs(c(mux1,muy1,muz1)) ) )
    stop("data of Classifier 1 is essentially constant")

  if (twocurves) {
    # compute means and standard errors of Classifier 2:
    mux2<-mean(x.2); sdx2<-SD(x.2)
    muy2<-mean(y.2); sdy2<-SD(y.2)
    muz2<-mean(z.2); sdz2<-SD(z.2)
    if (any(c(sdx2,sdy2,sdz2) < 10 * .Machine$double.eps * abs(c(mux2,muy2,muz2)) ) )
      stop("data of Classifier 2 is essentially constant")
    summary.dat2<- data.frame(n = c(nh2,n02,nd2), mu=c(mux2,muy2,muz2),
                              sd=c(sdx2,sdy2,sdz2), row.names=c("D-", "D0", "D+"))
    summary.dat <- list(Classifier1 = summary.dat, Classifier2 = summary.dat2)
  } else{
      mux2<-0; sdx2<-0
      muy2<-0; sdy2<-0
      muz2<-0; sdz2<-0
  }


    emp.vus.INT <- function(x, y, z) {
        dat1 <- expand.grid(x=x, y=y, z=z, KEEP.OUT.ATTRS = FALSE)
        x.lt.y <- dat1$x<dat1$y
        y.lt.z <- dat1$y<dat1$z
        x.eq.y <- dat1$x==dat1$y
        y.eq.z <- dat1$y==dat1$z

        sum.vus <- (x.lt.y & y.lt.z) +
            0.5 * ((x.lt.y & y.eq.z ) | (x.eq.y & y.lt.z)) +
            1/6 * (x.eq.y & y.eq.z)
        return(VUS = mean(sum.vus))
    }


  # compute empirical VUS1 and VUS2:
  vus.1 <- emp.vus.INT(x.1, y.1, z.1)
  vus.2 <- 1/6 # under H0
  if (twocurves) {
      vus.2 <- emp.vus.INT(x.2, y.2, z.2)
  }

  boot.sample <- function (x, seed) {
    if (!missing(seed) & !is.null(seed))
      set.seed(seed)
    return(sample(x, replace = TRUE))
  }

  boot.sample2 <- function (x, seed=1){
    n <- length(x)
    if (!missing(seed) | !is.null(seed))
      set.seed(seed)
    id <- sample(1:n, n, replace = TRUE)
    res <- x[id]
    return(list(res = res, id = id))
  }

  # new we compute the variance Var(VUS) of the Bootstrap test:
  # for paired data:
  if (paired) {
    boot.vus <- sapply(1:n.boot, function(i) {
      new.x  <- boot.sample2(x.1, i)
      new.x1 <- new.x$res
      new.x2 <- x.2[new.x$id]
      new.y  <- boot.sample2(y.1, (i+n.boot+172))
      new.y1 <- new.y$res
      new.y2 <- y.2[new.y$id]
      new.z  <- boot.sample2(z.1, (i+2*n.boot+172))
      new.z1 <- new.z$res
      new.z2 <- z.2[new.z$id]
      c(emp.vus.INT(new.x1, new.y1, new.z1), emp.vus.INT(new.x2, new.y2, new.z2))
    })
    boot.vus1 <- boot.vus[1,]
    boot.vus2 <- boot.vus[2,]
  }

  # for unpaired data:
  if (twocurves & !paired) {
    boot.vus1 <- sapply(1:n.boot, function(i) {
      new.x <- boot.sample(x.1, i)
      new.y <- boot.sample(y.1, (i+n.boot+172))
      new.z <- boot.sample(z.1, (i+2*n.boot+172))
      emp.vus.INT(new.x, new.y, new.z)
    })
    boot.vus2 <- sapply(1:n.boot, function(i) {
      new.x <- boot.sample(x.2, (i+3*n.boot+172))
      new.y <- boot.sample(y.2, (i+4*n.boot+172))
      new.z <- boot.sample(z.2, (i+5*n.boot+172))
      emp.vus.INT(new.x, new.y, new.z)
    })
  }

  # for single marker:
  if (!twocurves) {
    boot.vus1 <- sapply(1:n.boot, function(i) {
      new.x <- boot.sample(x.1, i)
      new.y <- boot.sample(y.1, (i+n.boot+172))
      new.z <- boot.sample(z.1, (i+2*n.boot+172))
      emp.vus.INT(new.x, new.y, new.z)
    })
  }

  # compute variances and covariance:
  var.vus1 <- VAR(boot.vus1)
  var.vus2 <- 0; cov.vus1vus2 <- 0 # default values
  if (twocurves) {
    var.vus2 <- VAR(boot.vus2) }
  if (paired) {
    cov.vus1vus2 <- COV(boot.vus1, boot.vus2) }

  # Compute test value and p-value:
  z.value  <- (vus.1-vus.2) / (var.vus1 + var.vus2 - 2*cov.vus1vus2)^0.5
  names(z.value) <- "Z-stat"

  if (alternative == "two.sided") {
    p.value <- 2 * pnorm(-abs(z.value), lower.tail = TRUE)
    alpha <- 1 - conf.level
    cint <- qnorm(1 - alpha/2) * (var.vus1 + var.vus2 - 2*cov.vus1vus2)^0.5
    cint <- (vus.1-vus.2) + c(-cint, cint)
  } else if (alternative == "less") {
    p.value <- pnorm(z.value, lower.tail = TRUE)
    cint <- c(-Inf, (vus.1-vus.2) + qnorm(conf.level) *
                (var.vus1 + var.vus2 - 2*cov.vus1vus2)^0.5)
  } else if (alternative == "greater") {
    p.value <- pnorm(z.value, lower.tail = FALSE)
    cint <- c((vus.1-vus.2) - qnorm(conf.level) *
                (var.vus1 + var.vus2 - 2*cov.vus1vus2)^0.5, Inf)
  }
  attr(cint, "conf.level") <- conf.level

  if (is.na(p.value))
  {p.value <- 1; print("Denominator is zero. Conclude complete alikeness.")}

  estimate <- vus.1
  names(estimate) <- "VUS of Classifier 1"
  if (twocurves) {
    estimate <- c(vus.1, vus.2)
    names(estimate) <- c("VUS of Classifier 1", "VUS of Classifier 2")}

  sigma <- diag(c(var.vus1, var.vus2))
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
  rval <- list(statistic = z.value, p.value = unname(p.value),
               conf.int = NULL, estimate =  estimate,
               null.value = null.value, alternative = alternative,
               method = method, data.name = dname,
               Summary = summary.dat, Sigma = sigma)
  class(rval) <- "htest" #hypothesis test element

  return(rval)
}



