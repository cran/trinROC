#' Empirical ROC surface plot
#'
#' Function for computation of the empirical ROC surface.
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
#'@param plot logical. If TRUE (default), the VUS is plotted using \code{surface3d} from the package \code{rgl}.
#'@param saveVUS A logical whether to save a PNG of the VUS in your current
#'  working directory (default is \code{FALSE}).
#'@return A list with the following components:
#'   \item{t1, t2, zVUS}{The matice containing the surface values.}
#'   \item{x, y, z}{The original data.}
#'
#' @export
#' @seealso \code{\link[rgl:surface3d]{surface3d}}.
#' @examples
#' data(cancer)
#' x1 <- with(cancer, cancer[trueClass=="healthy", 9])
#' y1 <- with(cancer, cancer[trueClass=="intermediate", 9])
#' z1 <- with(cancer, cancer[trueClass=="diseased", 9])
#'
#' rocsurf.emp(x1, y1, z1)

rocsurf.emp <- function(x,y,z, plot=TRUE, saveVUS = FALSE) {

  # lengths of the class vectors:
  nh <- length(x) # healthy
  n0 <- length(y) # intermediate
  nd <- length(z) # diseased

  #construct vector of possible cut-offs:
  con <- c(x,y,z)
  # add an additional minimum value for the case c_- being below min(con)
  con <- c(con, min(con)-1)
  socon <- unique(sort(con))
  total_cutoff <- length(socon)

  # construct matrices that form will carry all point triples of the
  # ROC surface. Entries that wont be filled are left with NA entries
  t1 <- matrix(NA, total_cutoff, total_cutoff)
  t2 <- t1
  Z <- t1
  # martices are needed since we have two cutoffs (2dimensions) to cover
  # all possible combinations of cutofs.

  for (j in 1:(total_cutoff)){
    for (i in 1:j) {
      # fills columwise the rates in (upper right triangle, since c_-<c_+):
      # socon[i] represents c_-, socon[j] represents c_+
      t1[i,j] <- sum(x <= socon[i])/nh # F_-(c_-)
      Z[i,j] <- sum(y > socon[i] & y <= socon[j])/n0 # F_0(c_+)-F_0(c_-)
      t2[i,j] <- sum(z > socon[j])/nd # G_+(c_+)
    }
  }

  if (plot) {
  # draw plot:
    colorlut <- cm.colors(50)
    col <- colorlut[ cut(Z[,ncol(Z):1], 50, labels = FALSE) ]

    # fix viewpoint of the visualisation:
    userMatrix <- matrix(0, 4, 4)
    userMatrix[1,] <- c(-0.6, .6, 0, 0)
    userMatrix[2,] <- c(-0, 0, 1, 0)
    userMatrix[3,] <- c(0.65, .65, 0, 0)
    userMatrix[4,] <- c(0, 0, 0, 1)

    open3d( userMatrix = userMatrix, windowRect = c(0,0,550,600))
    surface3d(t1[,ncol(t1)], t2[,ncol(t2):1][1,], Z[,ncol(Z):1], color = col,
              shade = 0.75, smooth = FALSE, shininess = 100 )
    grid3d(c("x", "y", "z"), n =10)
    axes3d(c("x+", "y+", "z+"), labels = TRUE, color = "darkgray")

    if (saveVUS == TRUE) rgl.snapshot("empVUS.png")
  }

  rval <- list(t1 = t1, t2 = t2, zVUS = Z,
               x = x, y = y, z = z)
  invisible(rval)
}

# res <- ROCs(x1,y1,z1,p=300)
# (p <- plot_ly(x=res$t1, y=res$t2, z=res$z, colors="RdBu") %>% add_surface() )
