#' Box-Cox transformation on three-class ROC data
#'
#' A transformation function for three-class ROC data in order to obtain normally
#' distributed classes.
#'
#' @details A Box-Cox transformation computing
#'
#' \deqn{X^{(\lambda)} = \left\{ \begin{array}{ll} (X^\lambda -1)/\lambda, 	&
#' \mbox{if } \; \lambda \neq 0,\\ \log(X),	& \mbox{else } \; \lambda = 0,
#' \end{array} \right.}{X^{(\lambda)} = log(X) if \lambda = 0 and X^{(\lambda)}
#' = (X^\lambda -1)/\lambda otherwise}
#'
#' with optimal \eqn{\lambda} estimated from the likelihood kernel function,
# \eqn{l_P(\lambda) = C - \frac{n}{2} \log({\sigma}_{(\lambda)}^2) +
# (\lambda-1)\sum_{i=1}^n \log(x_i)}{C -   n/2 log(\sigma^2_{(\lambda)}) +
# (\lambda-1)\sum_i log(x_i)},
#' as formally described in the supplementary
#' material in Bantis et al. (2017). If the data include any nonpositive
#' observations, a shifting parameter \code{lambda2} can be included in the
#' transformation given by:
#'
#' \deqn{X^{(\lambda)} = \left\{ \begin{array}{ll} ((X+\lambda_2)^\lambda -1)/\lambda, &	\mbox{if }
#' \, \lambda \neq 0,\\ \log(X+\lambda_2),	& \mbox{else } \; \lambda = 0. \end{array}
#' \right.\\
#' }{X^{(\lambda)} = log(X+\lambda_2), if \lambda = 0
#'  and X^{(\lambda)} =  ((X+\lambda_2)^\lambda -1)/\lambda
#'  otherwise.
#' }
#'
#' @param x,y,z  vectors containing the data of the three classes "healthy",
#'  "intermediate" and "diseased" to be transformed. In two-class ROC analysis only.
#' @param lambda vector of possible lambdas the log-likelihood function is evaluated.
#' @param lambda2 numeric shifting parameter. For the implemented Box-Cox
#'   transformation positive measurements in \code{x, y, z} are required.
#'   \code{lambda2} is used to shift these measurements.
#'@param eps numeric; indicating the bandwith around zero, where \code{lambda}
#'  is treated to be zero and the data is log-transformed.
#'@param verbose logical; indicating whether output should be displayed (default) or
#'  not.
# #'@section Warning:
# #'  There is no guarantee that the data behaves similar to
# #'  normal distributed data after transforming.
#' @references Bantis LE, Nakas CT, Reiser B, Myall D and Dalrymple-Alford JC
#'   (2015) Construction of joint confidence regions for the optimal true class
#'   fractions of receiver operating characteristic (roc) surfaces and
#'   manifolds. \emph{Statistical Methods in Medical Research} \bold{26}(3): 1429–1442.
#' @references Box, G. E. P. and Cox, D. R.  (1964). An analysis of
#'   transformations (with discussion). \emph{Journal of the Royal Statistical Society,
#'   Series B}, \bold{26}, 211–252.
#' @return A list with  components:
#'   \item{xbc, ybc, zbc}{The transformed vectors.}
#'   \item{lambda}{estimated optimal parameter.}
#'   \item{shapiro.p.value}{p-values obtained from \code{shapiro.test()} of
#'   the original and transformed data.}
#' @seealso
#'    \code{\link{shapiro.test}} and \code{\link[MASS:boxcox]{boxcox}} from the package \code{MASS}.
#' @export
#' @examples
#' data(cancer)
#' x1 <- with(cancer, cancer[trueClass=="healthy", 9])
#' y1 <- with(cancer, cancer[trueClass=="intermediate", 9])
#' z1 <- with(cancer, cancer[trueClass=="diseased", 9])
#'
#' boxcoxROC(x1, y1, z1)

boxcoxROC <- function(x, y, z, lambda = seq(-2., 2., 0.05),
                      lambda2 = NULL, eps = 0.02, verbose = TRUE) {

  if (any(is.na(c(x,y,z))))
    stop("Data includes NA's. Please omit first.")

  if (!is.null(lambda2)) {
    x <- x + lambda2
    y <- y + lambda2
    z <- z + lambda2
  } else {
    lambda2 <- 0
  }

  if (any(min(c(x,y,z)) <= 0))
    stop(paste("Data must be strictly positive. Specify shifting parameter lambda2 at least ",
               abs(min(x,y,z))))

  # length vectors:
  n <- length(x)
  m <- length(y)
  l <- length(z)


  roxlik <- function(h) {
    xh <- ((x^h)-1)/h
    yh <- ((y^h)-1)/h
    zh <- ((z^h)-1)/h
    loglik<- -(-n/2*log(sum((xh-sum(xh)/n)^2)/n) -
                m/2*log(sum((yh-sum(yh)/m)^2)/m) -
                l/2*log(sum((zh-sum(zh)/l)^2)/l) +
                (h-1)*(sum(log(x))+sum(log(y))+sum(log(z))) )

    return(loglik)
  }

  # using a non-linear minimization function to get optimal lambda:
  # a<-seq(-2,2, length.out = 100)
  # plot(a,roxlik(a), type="l")
  # opt <- optimize(roxlik, c(-2,2), tol=0.001)$minimum

  #parto <- nlm(roxlik, p=-.5) # p=1 is starting value of newton-type minimization
  #lambda<- parto$estimate
  ll <- sapply(lambda, roxlik)
  lambda <- lambda[which.min(ll)]

  # transforming data with optimal lambda:
  if (abs(lambda) <= eps) {
    xbc <- log(x)
    ybc <- log(y)
    zbc <- log(z)
  } else {
    xbc <- ((x^lambda)-1)/lambda
    ybc <- ((y^lambda)-1)/lambda
    zbc <- ((z^lambda)-1)/lambda
  }

  shap1 <- c(shapiro.test(x)$p, shapiro.test(y)$p, shapiro.test(z)$p)
  shap2 <- c(shapiro.test(xbc)$p, shapiro.test(ybc)$p, shapiro.test(zbc)$p)
  shap  <- rbind(shap1, shap2)
  colnames(shap) <- c("x", "y", "z")
  rownames(shap) <- c("Original data", "Box-Cox transformed")

  # output:
  if (verbose) {
    cat("---------------------------------------------------------------------",
        "\n", sep = " ")
    cat(" Optimal lambda       = ", lambda, "\n", sep = "")
    cat(" Shift param. lambda2 = ", lambda2, "\n\n", sep = "")
    cat(" Shapiro p-values for original data: ", "\n", sep = "")
    cat(" x = ",shap[1,1], ", y = ", shap[1,2], ", z = ", shap[1,3], "\n\n", sep = "")
    cat(" Shapiro p-values for Box-Cox transformed data: ", "\n", sep = "")
    cat(" x = ",shap[2,1], ", y = ", shap[2,2], ", z = ", shap[2,3], "\n", sep = "")
    cat("---------------------------------------------------------------------",
        "\n", sep = " ") }

  rval <- list(xbc = xbc,ybc = ybc,zbc = zbc,
               lambda = lambda, shapiro.p.value = shap)
  invisible(rval)
}
