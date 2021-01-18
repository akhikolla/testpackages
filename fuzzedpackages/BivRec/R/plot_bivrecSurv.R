########################    plot.bivrecSurv     ########################
#' Plot Bivariate Alternating Recurrent Series
#'
#' @description
#' This function plots bivariate recurrent event gap times from a \verb{bivrecSurv} object with an option to create separate plots based on categorical covariates.
#'
#' @import graphics
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix
#'
#' @param x A \verb{bivrecSurv} object.
#' @param y Either empty or NULL.
#' @param by A vector or data frame with categorical variables. See details.
#' @param main Optional string with plot title. Default is no title (will go to default when using \verb{by} argument).
#' @param xlab Optional string with label for horizontal axis. Default is "Time".
#' @param ylab Optional string with label for vertical axis. Default is "Individual".
#' @param type Optional vector of strings to label Type I and Type II gap times. Default is c("Type I", "Type II").
#' @param ... Additional arguments to be passed to graphical methods if needed.
#'
#' @export
#'
#' @details
#' Argument \verb{by} must be a vector or data frame with one or several categorical variables (up to 6 categories each).
#' Plots of the bivariate alternating recurrent outcome will be created by category for each variable.
#' When the \verb{by} argument is used it overrides the \verb{main} argument and sets it to its default \verb{main=""}.
#' To avoid errors make sure the vectors used for the \verb{bivrecSurv} object have the same length as the categorical variables and no missing values.
#'
#' @examples
#' \dontrun{
#'# Simulate bivariate alternating recurrent event data
#' library(BivRec)
#' set.seed(28)
#' bivrec_data <- simBivRec(nsize=100, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63,
#'                set=1.1)
#' plot(x = with(bivrec_data,bivrecSurv(id, epi, xij, yij, d1, d2)), main="Example",
#'      type = c("In Hospital", "Out of Hospital"))
#'
#' #Present the data by subgroups
#' #Note that the covariate a2 in the function will be dropped because it is continuous
#' attach(bivrec_data)
#' plot(x = bivrecSurv(id, epi, xij, yij, d1, d2), by = data.frame(a1, a2),
#'      type = c("In Hospital", "Out of Hospital"))
#' detach(bivrec_data)
#' }

plot.bivrecSurv <- function(x, y=NULL, by, type = NULL, main = NULL,
                            xlab = NULL, ylab = NULL, ...){

  if (!inherits(x, "bivrecSurv")) stop("Must be a bivrecSurv object.")
  object <- x

  #check arguments for labels
  if (missing(type)) {type=c("Type I","Type II")}
  if (missing(xlab)) {xlab="Time"}
  if (missing(ylab)) {ylab="Individual"}
  if (missing(main)) {main=""}
  args = c(main, xlab, ylab, type)

  if (missing(by)) {
    ##EXTRACT VECTORS FOR PLOTTING FUNCTION
  parameters <- object$data4Creg[,-(5:7)]
  colnames(parameters) <- c("id", "episode", "xij", "yij", "ci")
  ctimes <- object$data4Lreg$ctime
  nsubject <- object$data4Lreg$n

  basicplot(parameters=parameters, ctimes=ctimes,
            nsubject=nsubject, temp = NULL, args = args,
            c=0.95, cm=1.5, byp=FALSE)

  } else {

    if (is.data.frame(by)) {
      nrows_c <-  c(nrow(object$data4Creg), nrow(na.omit(by)))
      if (nrows_c[1]!= nrows_c[2]) {
        stop("Non-conformable arguments. Variables to create bivrecSurv object and categorical variables to split plots must have the same length and no missingness.")
      }
    } else {
      if (is.vector(by)) {
        nrows_c <-  c(nrow(object$data4Creg), length(na.omit(by)))
        if (nrows_c[1]!= nrows_c[2]) {
          stop("Non-conformable arguments. Vector variables to create bivrecSurv object and categorical variables to plot by must have the same length.")
        }
      } else {stop("Parameter by must be a vector or data frame.")}
      }

    #colnames(df) <- c("id", "episode", "xij", "yij", "ci")
    df = as.data.frame(cbind(object$data4Creg[,-(5:7)], na.omit(by)))
    plotBy(df, args)

  }





}
