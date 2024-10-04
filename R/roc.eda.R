#' Exploratory data analysis for a three-class ROC marker
#'
#' A function that investigates data that arose from a single marker and
#' contains the reference standard of the three classes "healthy",
#' "intermediate" and "diseased".
#'
#'@details For the preliminary assessment of a classifier, exporatory
#'  data analysis  (EDA) on the markers is necessary. This function assesses
#'  measurements from a single marker and computes the VUS, statistical tests
#'  and returns a summary table as well as some plots of the data.
#'@param x,y,z numeric vectors contaning the measurements from the healthy,
#'  intermediate and diseased class.
#'@param dat  a data frame of the following structure: The first column
#'  represents a factor with three levels, containing the true class membership
#'  of each measurement. The levels are ordered according to the convention of
#'  higher values for more severe disease status.
#'@param type a character, specifying if the \code{empirical} VUS and tests or
#'  the \code{trinormal} VUS and tests are computed.
#'@param plotVUS a logical whether to evaluate and plot the VUS (default is
#'  \code{FALSE}). Note: To save a png \code{plotVUS} needs to be \code{TRUE} too.
#'@param saveVUS a logical whether to save a PNG of the VUS in your current
#'  working directory (default is \code{FALSE}).
#'@param sep.dens a logical indicating if the density plots should be plotted
#'  on separate x-axes (\code{TRUE}) or on a common axis (\code{FALSE}, is
#'  default).
#'@param scatter a logical indicating if the measurements per class plot should
#'  be plotted as a boxplot (default) or as a scatterplot (\code{scatter =
#'  TRUE}).
#' @param conf.level A numeric value between 0 and 1 yielding the significance
#'   level \eqn{\alpha=1-\code{conf.level}}.
#'@param n.boot an integer incicating the number of bootstrap replicates sampled
#'  to obtain the variance of the VUS. Default is 1000.
#'@param verbose a logical, indicating whether output should be displayed or
#'  not. Default is \code{TRUE}.
#'@param alternative a character string specifying the alternative hypothesis,
#'  must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.
#'@param class.labels a vector of character strings of length 3 specifying the labels
#'  of the three classes used for the axis labeling in the boxplot.
#'@return A list with class "htest" containing the following components:
#'  \item{statistic}{The value of the test(s).}
#'  \item{p.value}{The p-value for the test(s).}
#'  \item{VUS}{the VUS computed with the specific method defined in
#'  \code{type}.}
#'  \item{dat.summary}{A data frame displaying size, mean and standard deviation
#'  of the three classes.}
#'  \item{alternative}{The alternative hypothesis.}
#'  \item{type}{a character containing the the method used for the exploratory
#'  data analysis.}
#'  \item{data.name}{a character containing the name of the data.}
#'  \item{xVUS, yVUS, zVUS}{(if \code{plotVUS = TRUE}) numeric vectors and
#'  matrices computed by \code{rocsurf.emp} or \code{rocsurf.trin}, used for
#'  displaying the surface with package \code{rgl}.}
#'  \item{histROC}{a \code{ggplot2} object, displaying the historgrams and
#'  densities of the three classes.}
#'  \item{meas.overview}{A ggplot2 object, displaying the boxplots (if
#'  \code{scatter = FALSE}) or scatter plots of the three classes (if
#'  \code{scatter = TRUE}).}
#'
#' @section Warning:
#' If \code{type = "empirical"}, computation may take a while, as \code{roc.eda} calls
#' the function \code{boot.test()}.
#' @seealso \code{\link{trinROC.test}}, \code{\link{trinVUS.test}} for trinormal
#'   data investigation, \code{\link{boot.test}} for empirical data analysis.
#'   \code{\link{rocsurf.emp}}, \code{\link{rocsurf.trin}} for the surface plot.
#' @export
#' @examples
#' data(krebs)
#'
#' # empirical EDA:
#' x <- with(krebs, krebs[trueClass=="healthy", 5])
#' y <- with(krebs, krebs[trueClass=="intermediate", 5])
#' z <- with(krebs, krebs[trueClass=="diseased", 5])
#' roc.eda(x, y, z, type = "e", plotVUS = FALSE)
#'
#' # equal data input via:
#' # roc.eda(dat = krebs[,c(1,5)], type = "e", plotVUS = FALSE)
#'
#' data(cancer)
#' # trinormal EDA:
#' roc.eda(dat = cancer[,c(1,10)], type = "trin", plotVUS = FALSE)
#' # trinormal EDA with different plots:
#' roc.eda(dat = cancer[,c(1,5)], type = "t", sep.dens = TRUE, scatter = TRUE)

#@param transform.data A logical, if TRUE a box-cox transformation is performed
#  on the measurements. This option only affects the results if \code{type =
#  "trinormal"}.


# RF:
# - changes from [0,1] -> (0,1)
# - class to inherits

roc.eda <- function(x, y, z, dat = NULL, type = c("empirical", "trinormal"),
                      plotVUS = FALSE, saveVUS = FALSE, sep.dens = FALSE,
                      scatter = FALSE, conf.level = 0.95, n.boot = 1000,
                      verbose = TRUE, alternative = c("two.sided", "less", "greater"),
                      class.labels = c("healthy", "intermediate", "diseased")) {

  # checking arguments:
  type <- match.arg(type)
  Marker <- "Classifier"

  if (!missing(conf.level) & (length(conf.level) != 1 | !is.finite(conf.level) |
                               conf.level <= 0 | conf.level >= 1))
    stop("'conf.level' must be a single number between 0 and 1.")

  alternative <- match.arg(alternative)

  if (!is.null(dat)) {
    if ( !inherits(dat,"data.frame") | !inherits(dat[,1],"factor") | ncol(dat) <= 1)
      stop("Data should be organized as a data frame with the group index factor at
           the first column and marker measurements at the second and third column.")
    # check for classes of column vectors of dat, should all be numeric:
    if (any(sapply(1 : (ncol(dat)-1), function(i) !inherits(dat[, i+1],"numeric"))) ) {
      not.num <- sapply(1 : (ncol(dat)-1), function(i) inherits(dat[, i+1],"numeric"))
      dat[, not.num==FALSE] <- as.numeric(dat[, not.num==FALSE])
      warning("Some measurements were not numeric. Forced to numeric.")
    }

  if (!length(class.labels)==3 | !inherits(class.labels,"character")){
    class.labels <- c("healthy", "intermediate", "diseased")
    warning("'class.labels' was not a vector of character strings of length 3.
            Forced to default.")
  }
    if (length(class.labels)!= length(unique(class.labels))){
      class.labels <- c("healthy", "intermediate", "diseased")
      warning("'class.labels' was not a vector of 3 different character strings.
              Forced to default.")
    }

    data.temp <- split(dat[,2], dat[,1], drop=FALSE)
    x <- data.temp[[1]]
    y <- data.temp[[2]]
    z <- data.temp[[3]]
    Marker <- names(dat)[2]
  }

  # check for NA's:
  if (any(is.na(x))) {
    warning("there are NA's in the healthy classes which are omitted.")
    x <- as.numeric(na.omit(x)) }

  if (any(is.na(y))) {
    warning("there are NA's in the intermediate classes which are omitted.")
    y <- as.numeric(na.omit(y)) }

  if (any(is.na(z))) {
    warning("there are NA's in the diseased classes which are omitted.")
    z <- as.numeric(na.omit(z)) }

    # computing the tests:
    if (any(type == "trinormal")) {

      method <- "Data overview of trinormal ROC Classifier"
      method2<- "Applied tests: Trinormal based ROC and VUS test"

      # if (transform.data) {
      #   transformed <- boxcoxROC(x, y, z)
      #   x <- transformed[[1]]
      #   y <- transformed[[2]]
      #   z <- transformed[[3]]
      #   lambda <- transformed[[4]]
      # }

      temptrinROC <- suppressWarnings(
                   trinROC.test(x, y, z, conf.level = conf.level))
      temptrinVUS <- suppressWarnings(
                    trinVUS.test(x, y, z, conf.level = conf.level,
                                 alternative = alternative))

      VUS       <- temptrinVUS$estimate[[1]]
      statistic <- c(temptrinROC$statistic, temptrinVUS$statistic)
      p.value   <- c(temptrinROC$p.value, temptrinVUS$p.value)
      conf.int  <- c(temptrinVUS$conf.int.estimate)
      names(VUS)      <- "trinormal VUS: "
      names(statistic)<- c("ROC test statistic: ", "VUS test statistic: ")
      names(p.value)  <- c("ROC p.value: ", "VUS p.value: ")
      summary   <- temptrinROC$Summary

      # compute 3dim surface:
      surf <- list(t1 = NULL, t2 = NULL, zVUS = NULL)

      if (plotVUS) {
        (surf <- rocsurf.trin(x ,y , z, saveVUS = saveVUS))
      }
  }

  if (any(type == "empirical")) {

    method <- "Data overview of empirical ROC Classifier"
    method2<- "Applied test: Bootstrap test"

    # compute the statistical summary:
    tempBoot <- suppressWarnings(
              boot.test(x, y, z, conf.level = conf.level, n.boot = n.boot,
                        alternative = alternative) )

    VUS      <- tempBoot$estimate[1]
    statistic<- tempBoot$statistic
    p.value  <- tempBoot$p.value
    conf.int <- tempBoot$conf.int.estimate
    names(VUS)      <- "empirical VUS: "
    names(statistic)<- "Boot statistic: "
    names(p.value)  <- "Boot p.value: "
    summary  <- tempBoot$Summary

    # compute 3dim surface:
    surf <- list(t1 = NULL, t2 = NULL, zVUS = NULL)

    if (plotVUS) {
      (surf <- rocsurf.emp(x ,y , z, saveVUS = saveVUS))
    }
  }

    ## compute 2dim plots:
    # construct data frame for plotting:
    trueClass <- factor(c(rep(class.labels[1], length(x)),
                        rep(class.labels[2], length(y)),
                        rep(class.labels[3], length(z))),
                        levels = class.labels )
    value <- c(x,y,z)
    data <- data.frame(trueClass = trueClass, value = value)

    if (sep.dens == FALSE) {
      # common x-axe histograms & densities:
      histROC <- ggplot(data, aes(x = value, colour=trueClass, fill=trueClass)) +
        geom_histogram(aes(y =after_stat(density)), binwidth=(max(data$value)-min(data$value))/15,
                       position = "dodge", alpha=0.7, show.legend = FALSE) +
        scale_colour_manual(values=c("#79AB67", "#6EA3D0", "#D68898"), guide = "none") +
        scale_fill_manual(values=c("#79AB67", "#6EA3D0", "#D68898"), name = "Class") +
        labs(y="Count",x=paste(Marker, "measurements")) +
        #facet_grid(. ~ trueClass, scales = "free") +
        geom_density(aes(col=trueClass), show.legend = FALSE, alpha=0.2)

    # separate x-axe histograms & densities:
    } else {
      histROC <- ggplot(data, aes(x = value, colour=trueClass, fill=trueClass)) +
        geom_histogram(aes(y =after_stat(density)),binwidth=(max(data$value)-min(data$value))/30,
                       show.legend = FALSE) +
        scale_colour_manual(values=c("#79AB67", "#6EA3D0", "#D68898"), guide = "none") +
        scale_fill_manual(values=c("#79AB67", "#6EA3D0", "#D68898"), name="Class") +
        labs(y="Density",x=paste(Marker, "measurements")) +
        facet_grid(. ~ trueClass, scales = "free") +
        geom_density(col=2, show.legend = FALSE, fill = NA)
    }

    if (scatter == FALSE) {
    # boxplots:
    meas.overview <- ggplot(data, aes(y=value, x=trueClass, fill = trueClass)) +
      stat_boxplot(geom ='errorbar') + geom_boxplot() +
      scale_fill_manual(values=c("#79AB67", "#6EA3D0", "#D68898"), guide = "none") +
      coord_flip() + labs(x="", y=paste(Marker, "measurements"))

    } else {
    # scattter plot of data:
    meas.overview <- ggplot(data, aes(y=value, x=trueClass, color = trueClass)) +
      geom_jitter(width = 0.25) +
      scale_colour_manual(values=c("#79AB67", "#6EA3D0", "#D68898"), guide = "none") +
        coord_flip() + labs(x="", y=paste(Marker, "measurements"))
    }

    grid.arrange(histROC, meas.overview, ncol = 2)


  # prepare output:
  dname <- if (!is.null(dat)) levels(dat[,1])
  else c(deparse(substitute(x)), deparse(substitute(y)), deparse(substitute(z)))

  null.value <- 0
  names(null.value) <- "Difference in VUS"


  if (verbose) {
    cat("\n", method, "\n", sep = " ")
    cat("---------------------------------------------------------------------",
        "\n", sep = " ")
    cat("\n", method2, "\n", sep = " ")
    cat(" Significance level: ", 1-conf.level, "\n", sep = "")
    cat(" Alternative hypothesis: ", alternative, "\n", sep = "")
    cat("---------------------------------------------------------------------",
        "\n", sep = " ")
    cat(" data: ", dname[1], ", ", dname[2], " and ", dname[3], "\n\n", sep = "")
    cat(" ", rbind(names(statistic), round(statistic,3), c(", ",", ") ,
                   names(p.value), round(p.value,5))[,1], "\n ", sep = "")
    if (type == "trinormal") {
      cat("", rbind(names(statistic), round(statistic,3), c(", ",", ") ,
                   names(p.value), round(p.value,5))[,2], "\n", sep = "") }

    cat("\n", " ", names(VUS), round(VUS,3), ", 95% CI: ",sep="")
    cat(round(conf.int,3), "\n", sep = " ")
    if (type == "trinormal") {
      cat("\n", " Parameters: ", "\n", sep = "")
      cat(" ", names(temptrinROC$estimate[-1]),"\n",sep="\t")
      cat(" ", as.numeric(round(temptrinROC$estimate[-1],4)),"\n",sep="\t") }
    cat("---------------------------------------------------------------------",
        "\n", sep = " ")
  }

  rval <- list(statistic = statistic, p.value = p.value,
               VUS=VUS, conf.int=conf.int, dat.summary=summary, alternative = alternative,
               type = type, data.name = dname, method = method,
               xVUS = surf$t1, yVUS = surf$t2, zVUS = surf$zVUS,
               histROC = histROC, meas.overview = meas.overview,
               x = x, y = y, z = z)

  invisible(rval)

}
