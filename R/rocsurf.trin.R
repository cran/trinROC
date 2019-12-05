#' Trinormal ROC surface plot
#'
#' Function for computation of the trinormal ROC surface.
#'
#'@details This function takes three-class ROC data and computes the three
#'  dimentional surface using the R-package \code{rgl}. The ROC surface is
#'  defined as
#'
#'  \deqn{z = ROCs(t_-,t_+) = F_0(c_+) - F_0(c_-)=F_0(G_+^{-1}(t_+) )
#'  -F_0(F_-^{-1}(t_-) ),}
#'
#'  where \eqn{c_-, c_+} are the two cut-off points and \eqn{F_-, F_0, F_+} the
#'  cdf of the three classes with \eqn{G = 1-F}.
#'@param x,y,z  Vectors containing the data of the three classes "healthy",
#'  "intermediate" and "diseased".
#'@param p An integer for the precision of the surface. \code{p} gives the
#'  number of gridpoints per axis.
#'@param plot logical. If TRUE (default), the VUS is plotted using \code{surface3d} from the package \code{rgl}.
#'@param saveVUS A logical whether to save a PNG of the VUS in your current
#'  working directory (default is \code{FALSE}).
#'@return A list with the following components:
#'\item{t1, t2}{The vectors \eqn{t_-=F_-^{(c_-)}} and \eqn{t_+=F_+^{(c_+)}}}
#'   \item{zVUS}{The matix containing the surface values.}
#'   \item{x, y, z}{The original data.}
#'
#' @export
#' @references Xiong, C., G. Van Belle, et al. (2006). Measuring and estimating
#'   diagnostic accuracy when there are three ordinal diagnostic groups.
#'   \emph{Statistics in Medicine} 25(7), 1251–1273.
#' @examples
#' data(cancer)
#' x1 <- with(cancer, cancer[trueClass=="healthy", 8])
#' y1 <- with(cancer, cancer[trueClass=="intermediate", 8])
#' z1 <- with(cancer, cancer[trueClass=="diseased", 8])
#'
#' rocsurf.trin(x1, y1, z1)

rocsurf.trin <- function(x, y, z, p=300, plot = TRUE, saveVUS = FALSE) {

  # compute estimated means and variances:
  mu_1 <- mean(x) ; sigma_1 <- VAR(x)
  mu_2 <- mean(y) ; sigma_2 <- VAR(y)
  mu_3 <- mean(z) ; sigma_3 <- VAR(z)

  # probabilities t1 = F_-(c1), t2 = G_+(c2) for cut-offs c1 < c2:
  t1 <- seq(0,1, length.out = p)
  t2 <- seq(0,1, length.out = p)
  # set matrix of values for ROC surface:
  Z <- matrix(0, nrow = p, ncol = p)

  # compute the parameters a,b,c,d from xiong et al. 2006:
  a <- sigma_2/sigma_1
  b <- (mu_1 - mu_2)/sigma_1
  c <- sigma_2/sigma_3
  d <- (mu_3 - mu_2)/sigma_3

  for (i in 1:p) {
    j <- 1 # set counter

    repeat{
      Z[i,j] <- pnorm( (qnorm(1-t2[j]) + d)/c ) - pnorm( (qnorm(t1[i])+b)/a )
      j <- j+1
      # only defined in t_+ < 1-Phi(F_+^-1(t_-))
      if (j == p+1 || t2[j] >= 1-pnorm(qnorm(t1[i], mean=mu_1,sd=sigma_1),
                                       mean=mu_3,sd=sigma_3) )   break
    }
  }

  # delete some 0-entries of z:
  for (i in 1:p) {
    for (j in 1:p) {
      # set first condition
      if (Z[i,j] == 0) {
        # border value of i:
        if (i == 1 ) {
          if (Z[i,j-1]==0 || is.na(Z[i,j-1])) Z[i,j] <- NA
          # border value of j
        } else if (j == 1 ) {
          if (Z[i-1,j]==0 || is.na(Z[i-1,j])) Z[i,j] <- NA
          # set second condition for inner points of z:
        } else
          if ( (Z[i-1,j]==0 || is.na(Z[i-1,j])) &&
               (Z[i,j-1]==0 || is.na(Z[i,j-1])) &&
               (Z[i-1,j-1]==0 || is.na(Z[i-1,j-1])) )  Z[i,j] <- NA
      }
    }
  }

  if (plot) {
  # draw plot:
    colorlut <- cm.colors(50)
    col <- colorlut[ cut(Z, 50, labels = FALSE) ]

    # fix viewpoint of the visualisation:
    userMatrix <- matrix(0, 4, 4)
    userMatrix[1,] <- c(-0.6, .6, 0, 0)
    userMatrix[2,] <- c(-0, 0, 1, 0)
    userMatrix[3,] <- c(0.65, .65, 0, 0)
    userMatrix[4,] <- c(0, 0, 0, 1)

    open3d( userMatrix = userMatrix, windowRect = c(0,0,550,600))
    surface3d(t1, t2, Z, color = col,
              shade = 0.75, smooth=F, shininess = 100 )
    grid3d(c("x", "y", "z"), n =10)
    axes3d(c("x+", "y+", "z+"), labels=T, color = "darkgray")

    if (saveVUS == T) rgl.snapshot("trinVUS.png")
  }
  rval <- list(t1 = t1, t2 = t2, zVUS = Z,
               x = x, y = y, z = z)
  invisible(rval)
}

