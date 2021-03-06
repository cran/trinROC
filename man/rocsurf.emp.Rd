% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rocsurf.emp.R
\name{rocsurf.emp}
\alias{rocsurf.emp}
\title{Empirical ROC surface plot}
\usage{
rocsurf.emp(x, y, z, plot = TRUE, saveVUS = FALSE)
}
\arguments{
\item{x, y, z}{Vectors containing the data of the three classes "healthy",
"intermediate" and "diseased".}

\item{plot}{logical. If TRUE (default), the VUS is plotted using \code{surface3d} from the package \code{rgl}.}

\item{saveVUS}{A logical whether to save a PNG of the VUS in your current
working directory (default is \code{FALSE}).}
}
\value{
A list with the following components:
  \item{t1, t2, zVUS}{The matice containing the surface values.}
  \item{x, y, z}{The original data.}
}
\description{
Function for computation of the empirical ROC surface.
}
\details{
This function takes three-class ROC data and computes the three
 dimentional surface using the R-package \code{rgl}. The ROC surface is
 defined as

 \deqn{z = ROCs(t_-,t_+) = F_0(c_+) - F_0(c_-)=F_0(G_+^{-1}(t_+) )
 -F_0(F_-^{-1}(t_-) ),}

 where \eqn{c_-, c_+} are the two cut-off points and \eqn{F_-, F_0, F_+} the
 cdf of the three classes with \eqn{G = 1-F}.
}
\examples{
data(cancer)
x1 <- with(cancer, cancer[trueClass=="healthy", 9])
y1 <- with(cancer, cancer[trueClass=="intermediate", 9])
z1 <- with(cancer, cancer[trueClass=="diseased", 9])

rocsurf.emp(x1, y1, z1)
}
\seealso{
\code{\link[rgl:surface3d]{surface3d}}.
}
