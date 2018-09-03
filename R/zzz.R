.onLoad <- function(libname, pkgname) {
  op <- options()

  op.trinROC <- list(
    trinROC.MLE = TRUE
  )
  toset <- !(names(op.trinROC) %in% names(op))
  if(any(toset)) options(op.trinROC[toset])

  invisible()
}

#
# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("")
#
# }
