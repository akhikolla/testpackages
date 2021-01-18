#' @title Partial association analysis between ordinal responses after adjusting for a set of covariates
#'
#' @description This function is mainly designed for conducting the partial association analysis.
#' It provides two ways of using:
#'
#' 1. A user-friendly way: only need "responses", "adjustments", and
#' "data". All the rest of the argument will be setted as default (see Arguments for details of
#' default).
#'
#' 2. An advanced way: user can input a list of fitted models by "fitted.models", then the
#' "responses" and "adjustments" are not necessary. Supported class of cumulative link models in
#' \code{\link[ordinal]{clm}}, \code{\link[stats]{glm}}, \code{\link[rms]{lrm}},
#' \code{\link[rms]{orm}}, \code{\link[MASS]{polr}}, \code{\link[VGAM]{vglm}}, .
#'
#' It generates an object that has partial association matrix, marginal association,
#' and some attributes: "arguments" saves c(association, method, resids.type). "responses"
#' contains the names of response variables. The attribute "adjustments" contains the names
#' of covariates. The "summary" function of "PAsso" class of object provides marginal association '
#' matrix for comparison purpose.
#'
#' @param responses A string vector that specifies response variables. It requires to be equal
#' or greater than two variables in the data frame.
#' @param adjustments A string vector specifies covariates/confounders that need to
#' be adjusted.
#' @param data A data.frame including responses and adjustments.
#' @param uni.model A character string specifying the universal model setting for all
#' responses. Default \code{"logit"} refers to cumulative logit model. \code{"probit"}
#' refers to cumulative probit model. \code{"acat"} fits an adjacent categories regression model.
#'
#' @param models A string vector contains default link functions of fitting models with
#' respect to each response variable. If \code{"models"} is missing or has any one of the model
#' unspecified, \code{"uni.model"} is used to specify same models for all responses automatically.
#' But, this argument has higher priority than the \code{"uni.model"} as long as the length of
#' \code{"models"} equals to the number of \code{"responses"}.
#'
#' @param method A string argument to specify correlation coefficient method.
#' Three choices \code{c("kendall", "pearson", "wolfsigma")}. The default is
#' \code{"kendall"}
#' @param resids.type A character string specifying which type of residuals to generate
#' Current options are \code{"latent"} and \code{"uniform"}. Default is \code{"latent"}.
#' \describe{
#'   \item{\code{surrogate}}{surrogate residuals (Liu and Zhang, 2017);}
#'   \item{\code{sign}}{sign-based residuals (Li and Shepherd, 2010, 2012);}
#'   \item{\code{general}}{generalized residuals (Franses and Paap, 2001);}
#'   \item{\code{deviance}}{deviance residuals (-2*loglik).}
#' }
#'
#' Although \code{"sign"}, \code{"general"}, and \code{"deviance"} are provided in
#' this package, these residuals are problematic for partial association analysis
#' between ordinal response (more discussions see Liu, Dungang, Li, Shaobo, Yu, Yan,
#' and Moustaki, Irini.(2020))
#'
#' @param jitter A character string specifying how to generate surrogate residuals.
#' Current options are \code{"latent"} and \code{"uniform"}. Default is \code{"latent"}.
#' \describe{
#'   \item{\code{latent}}{surrogate residuals.}
#'   \item{\code{uniform}}{sign-based residuals.}
#' }
#'
#' @param jitter.uniform.scale A character string specifying the scale on which to perform
#' the jittering whenever \code{jitter = "uniform"}. More details: \code{PAsso::residuals}.
#'
#' @param fitted.models A list contains all the models (S3 objects) you want to
#' assess for the partial association between ordinal responses after adjusting
#' for a set of covariates covariates. All of these models should be applied to the
#' same dataset, having same covariates, same sample size etc. The models in this
#' list can be an object of class \code{\link[ordinal]{clm}},
#' \code{\link[stats]{glm}}, \code{\link[rms]{lrm}}, \code{\link[rms]{orm}},
#' \code{\link[MASS]{polr}}, \code{\link[VGAM]{vglm}}.
#'
#' @param n_draws A number to specify draws of surrogate residuls
#' such that the partial correlation coefficients are calculated repeatedly. The final
#' correlation coefficients are the average of all partial correlation coefficients.
#' It is the \code{"nsim"} argument in \code{"residuals()"} function.
#' @param association An default argument to specify the partial association. Leave this
#' further development of package such that other association analyses can be embedded.
#' @param ... Additional optional arguments.
#'
#' @return An object of class \code{"PAsso"} is a list containing at least the following
#' components. It contains the partial correlation matrix and multiple repeats if
#' \code{n_draws} > 1. This object has "arguments"
#' attribute saved as c(association, method, resids.type), "responses" attribute, and
#' "adjustments" attribute.
#' The list contains:
#' \describe{
#'   \item{\code{corr}}{The estimated correlation matrix(average of \code{rep_MatCorr})
#'   of partial association after adjusting confounders;}
#'   \item{\code{rep_corr}}{The replications of estimated correlation matrix;}
#'   \item{\code{rep_SRs}}{The replications of surrogate residuals if partial association is applied;}
#'   \item{\code{fitted.models}}{The list stores all fitted.models;}
#'   \item{\code{data}}{The data.frame of original dataset;}
#'   \item{\code{mods_n}}{The sample size of each fitted model;}
#'   \item{\code{cor_func}}{The correlation function after assign different method;}
#'   \item{\code{Marg_corr}}{The marginal association matrix.}
#' }
#' @references
#' Liu, Dungang, Li, Shaobo, Yu, Yan, and Moustaki, Irini. Assessing partial association between
#' ordinal variables: quantification, visualization, and hypothesis testing, \emph{Journal of the
#' American Statistical Association}, Revision under review.
#'
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach. \emph{Journal of the American Statistical Association}.
#' \url{http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20}
#'
#' Li, Chun, and Bryan E. Shepherd. "Test of association between two ordinal variables while
#' adjusting for covariates." \emph{Journal of the American Statistical Association} 105, no.
#' 490 (2010): 612-620. \url{https://doi.org/10.1198/jasa.2010.tm09386}
#'
#' Li, Chun, and Bryan E. Shepherd. "A new residual for ordinal outcomes." \emph{Biometrika}
#' 99, no. 2 (2012): 473-480. \url{https://doi.org/10.1093/biomet/asr073}
#'
#' Franses, Philip Hans, and Richard Paap. Quantitative models in marketing research.
#' Cambridge University Press, 2001.
#' \url{https://pdfs.semanticscholar.org/dad0/820f287a8cf5a4e8039549e35fc111fd86e5.pdf}
#'
#' @importFrom ggplot2 aes_string geom_abline geom_boxplot geom_point
#'
#' @importFrom ggplot2 geom_smooth ggplot ggtitle guides labs xlab ylab
#'
#' @importFrom stats .checkMFClasses lowess median model.frame model.matrix
#'
#' @importFrom stats model.response nobs pbinom pcauchy plogis pnorm ppoints
#'
#' @importFrom stats predict qcauchy qlogis qnorm qqline qqplot qqnorm quantile
#'
#' @importFrom stats qunif runif
#'
#' @importFrom stats filter lag
#'
#' @importFrom MASS polr
#'
#' @importFrom copBasic wolfCOP
#'
#' @importFrom pcaPP cor.fk
#'
#' @importFrom methods slot
#'
#' @importFrom VGAM vglm acat summary
#'
#'
#' @export
#'
#' @examples
#'
#' ###########################################################
#' # User-friendly way of using
#' ###########################################################
#' library(MASS)
#'
#' # Import ANES2016 data in "PAsso"
#' data(ANES2016)
#'
#' # User-friendly way of using: Parial association analysis
#' PAsso_1 <- PAsso(responses = c("PreVote.num", "PID"),
#'                 adjustments = c("income.num", "age", "edu.year"),
#'                 data = ANES2016,
#'                 method = c("kendall"))
#'
#' print(PAsso_1, digits = 4)
#' summary(PAsso_1, digits = 4)
#'
#' ###########################################################
#' # Advanced way of using
#' ###########################################################
#'
#' fit.vote<- glm(PreVote.num ~ income.num+ age + edu.year, data = ANES2016,
#'                family = binomial(link = "probit"))
#' fit.PID<- polr(as.factor(PID) ~ income.num+age+edu.year, data = ANES2016,
#'                method="probit", Hess = TRUE)
#'
#' PAsso_adv1 <- PAsso(fitted.models=list(fit.vote, fit.PID),
#'                     method = c("kendall"),
#'                     resids.type = "surrogate")
#'
#' print(PAsso_adv1, digits = 4)
#' summary(PAsso_adv1, digits = 4)
#'
#' @useDynLib PAsso
PAsso <- function(responses, adjustments, data,
                  uni.model = c("probit", "logit", "acat"),
                  models = NULL,
                  method = c("kendall", "pearson", "wolfsigma"),
                  resids.type = c("surrogate", "sign", "general", "deviance"),
                  jitter = c("latent", "uniform"),
                  jitter.uniform.scale = c("probability", "response"),
                  fitted.models = NULL,
                  n_draws = 20,
                  association = "partial", ...){

  # TEST HEADER:
  # data=boot_data
  # responses = attr(object, "responses")
  # adjustments = attr(object, "adjustments")
  # models = attr(object, "models")
  # association = arguments[1]
  # method = arguments[2]
  # resids.type = arguments[3]

  # responses = c("PreVote.num", "PID")
  # models = c("probit", "acat")
  # adjustments <- c("income.num", "age", "edu.year")
  # association = "partial"; method = "kendall";
  # n_draws = 30; data = ANES2016;
  # resids.type = "surrogate"; jitter = "latent"
  # jitter.uniform.scale = "response"
  # models = c("probit", "probit")

  n_resp <- ifelse(missing(responses), length(fitted.models), length(responses))

  # Match arguments -------------------------------------------------
  call <- match.call()

  # Give default models for each response
  uni.model <- match.arg(uni.model)
  if (missing(models)) {
    models <- rep(uni.model, n_resp)
  }

  association <- match.arg(association)
  method <- match.arg(method)
  resids.type <- match.arg(resids.type)
  jitter <- match.arg(jitter)
  jitter.uniform.scale <- match.arg(jitter.uniform.scale)

  # Initialize the cor_func for different methods
  cor_func <- switch(method,
                     # kendall = function(X) cor(X, method = "kendall"),
                     kendall = function(X) cor.fk(X),
                     pearson = function(X) cor(X),
                     wolfsigma = function(X)
                       t(wolfCOP(para = data.frame(X), as.sample = TRUE)))



  # Main Body ---------------------------------------------------------------
  # Partial Association based on surrogate residuals
  # Pull out ingredients for cooking partial association --------------------------------------------
  if (!missing(fitted.models) & is.list(fitted.models)) { # If fitted.models is imported as a list, then done!

    conf_temp0 <- lapply(fitted.models, function(mod) colnames(model.frame(mod)[1,]))

    responses <- lapply(conf_temp0, `[[`, 1)
    n_responses <- length(responses)

    mods_ys <- sapply(fitted.models, function(mod) as.numeric(model.frame(mod)[,1]))
    colnames(mods_ys) <- c(responses)
    # Save numeric response for calculating marginal corr

    conf_temp <- lapply(conf_temp0, function(mod) mod[-1])


    # responses <- sapply(fitted.models, function(mod) colnames(model.frame(mod))[1])
    # n_responses <- length(responses)
    #
    # mods_ys <- sapply(fitted.models, function(mod) model.frame(mod)[,1])
    # colnames(mods_ys) <- c(responses)
    # conf_temp <- lapply(fitted.models, function(mod) colnames(model.frame(mod)[1,-1]))

    if (length(unique.default(conf_temp)) != 1L) { # If data are not same, stop!
      stop("The imported fitted.models have different confounders!")
    } else { # If confounders are same, use them!
      adjustments <- conf_temp[[1]]; n_adjustments <- length(adjustments)
    }

    data <- data.frame(mods_ys, getCovariates(fitted.models[[1]])) # Make sure data is data.frame to avoid issue!
    colnames(data) <- c(responses, adjustments)

    # Obtain models(links) from the "fitted.models", and update "models"
    distName <- c()
    for (i in 1:n_responses) {
      distName[i] <- getDistributionName(fitted.models[[i]])

      # Add "ord" as "ordinal response model"
      # Use "ord" as an extend of "glm", "clm", "lrm", "orm", "polr", "vglm"
      # class object such that "residuals.ord" can recognize them all.
      # Here is a test for other model objects: polr, vglm etc.
      class(fitted.models[[i]]) <- c("ord", class(fitted.models[[i]]))
    }
    models <- sapply(X = distName,
                     FUN = function(x) switch(x,
                                              "logis" = "logit",
                                              "norm" = "probit",
                                              "gumbel" = "loglog",
                                              "Gumbel" = "cloglog",
                                              "cauchy" = "cauchit"))

  } else { # Start to deal with "responses" and "adjustments"
    if (length(responses) != length(models)) {
      mess_str <- "Not all models of the responses are specified.
      Use same model as specified in 'uni.model' for all responses!
      If it is not you want, quit and specify again."
      message(paste(strwrap(mess_str), collapse = "\n"))
      models <- rep(uni.model, n_resp)
    }

    n_responses <- length(responses)
    n_adjustments <- length(adjustments)

    # mods_ys <- as.numeric(sapply(models, function(mod) model.frame(mod)[,1]))
    mods_ys <- apply(X = data[,responses], MARGIN = 2, FUN = as.numeric)
    # Save numeric response for calculating marginal corr

    formulaAll <-
      sapply(responses,
             function(X) paste(X, paste(adjustments, collapse = "+"), sep = "~"))

    for (i in responses) { # Change all response variables to factor type to avoid error!
      data[,i] <- as.factor(data[,i])
    }

    # Combine Models together -------------------------------
    # Use two different functions (polr, glm) if response has differen levels!
    # NEEDED FEATURE: Test for other packages.
    fitted.models <- list()
    for (i in 1:n_responses) {

      if (models[i] %in% "acat") { # Conduct Adjacent Categories regression model.
        data[,responses[i]] <- as.numeric(data[,responses[i]])
        fitted_temp <- do.call("vglm",
                               list(formula = as.formula(formulaAll[i]),
                                    data = quote(data),
                                    family = quote(acat(reverse = TRUE, parallel = TRUE))))

        # need to change the response back to factor to avoid issue!
        data[,responses[i]] <- as.factor(data[,responses[i]])

      } else if (length(unique(data[,responses[i]])) > 2) {
        # If response has more than 2 levels, use "polr" or "VGAM::acat", otherwise "glm".

        temp_model <- ifelse(models[i]=="logit", "logistic", models[i])
        # Change link from "logit" to "logistic" since "polr" only accept "logistic" link

        if (temp_model %in% c("logit", "logistic", "probit")) {
          fitted_temp <- do.call("polr", list(formula = as.formula(formulaAll[i]),
                                              Hess = TRUE, # Need this to draw coefficients and t-values as output
                                              method = temp_model, data = quote(data)))
        }

      } else {
        fitted_temp <- do.call("glm",
                               list(formula = as.formula(formulaAll[i]),
                                    family = quote(binomial(link = models[i])),
                                    data = quote(data)))

        # In order to reregister residuals for the glm class, I add "ord" class as "ordinal response model"
        # Use "ord" as an extend of "glm" class object such that "residuals.ord" can recognize it
        # This class of object can be extended to other model objects: polr, vglm etc.
        class(fitted_temp) <- c("ord", class(fitted_temp))
      }
      fitted.models[[i]] <- fitted_temp
      # FIXED: formula now is shown as what it is (by "do.call" and "quote")!
    }

    # Transfer the data again to obtain numeric responses!!!
    data <- data.frame(mods_ys, data[,adjustments]) # Make sure data is data.frame to avoid issue!
    colnames(data) <- c(responses, adjustments)
  }

  mods_n <- sapply(fitted.models, nobs)
  # mods_ys_names <- responses # sapply(fitted.models, function(mod) colnames(model.frame(mod))[1])

  if (length(unique.default(mods_n)) != 1L) { # if sample sizes of all models object are not same, stop!
    stop("Stop due to different data in the models!")
  }


  # Simulate surrogate residuals for calculating the correlation -----------------------------------
  MatCorr <- matrix(1, nrow = n_responses, ncol = n_responses, dimnames = list(responses, responses))

  if (n_draws==1) { # Use just one replication of SR to calcualte partial correlation!
    rep_SRs <- array(data = NA, dim = c(mods_n[1],n_draws,n_responses)) # Keep array for consistency.

    for (i in 1:n_responses) {
      # print(class(fitted.models[[i]]))
      if (isS4(fitted.models[[i]]) & inherits(fitted.models[[i]], "vglm")) {
        use_func <- "residualsAcat"
      } else {
        use_func <- "residuals"
      }
      temp_resids <- do.call(what = use_func,
                             args = list(object = fitted.models[[i]],
                                         type = resids.type, jitter = jitter,
                                         jitter.uniform.scale = jitter.uniform.scale,
                                         nsim = n_draws)
                             )
      rep_SRs[,1,i] <- temp_resids
    }

    dimnames(rep_SRs)[[1]] <- c(1:mods_n[1])
    dimnames(rep_SRs)[[2]] <- seq(n_draws)
    dimnames(rep_SRs)[[3]] <- responses

    if (method == "wolfsigma") { # wolfsigma only return one value!
      MatCorr[upper.tri(MatCorr)] <- MatCorr[lower.tri(MatCorr)] <- cor_func(rep_SRs[,1,])
    } else {
      MatCorr <- cor_func(rep_SRs[,1,])
    }
    rep_MatCorr <- MatCorr # Save Replication again!

  } else { # Repeat n_draws times to get SRs and calcualte average partial correlation!
    rep_SRs <- array(data = NA, dim = c(mods_n[1],n_draws,n_responses))
    for (i in 1:n_responses) {
      # print(class(fitted.models[[i]]))
      if (isS4(fitted.models[[i]]) & inherits(fitted.models[[i]], "vglm")) {
        use_func <- "residualsAcat"
      } else {
        use_func <- "residuals"
      }
      temp_resids <- do.call(what = use_func,
                             args = list(object = fitted.models[[i]],
                                         type = resids.type, jitter = jitter,
                                         jitter.uniform.scale = jitter.uniform.scale,
                                         nsim = n_draws)
                             )
      rep_SRs[,,i] <- attr(temp_resids, "draws")
    }

    dimnames(rep_SRs)[[1]] <- c(1:mods_n[1])
    dimnames(rep_SRs)[[2]] <- seq(n_draws)
    dimnames(rep_SRs)[[3]] <- responses

    rep_MatCorr_temp <- apply(X = rep_SRs, MARGIN = 2, cor_func) # For each copy, calculate correlation
    # For "wolfsigma", below code has bug!
    MatCorr <- matrix(apply(X = rep_MatCorr_temp, MARGIN = 1, mean),
                      nrow = n_responses, ncol = n_responses, dimnames = list(responses, responses))
    # ABOVE: Take mean as the estimate of partial correlation matrix!

    rep_MatCorr <- array(rep_MatCorr_temp, dim = c(n_responses, n_responses, n_draws),
                         dimnames = list(responses, responses))
  }

  # Add marginal association!
  # Marg_corr <- cor(mods_ys, method = method)
  Marg_corr <- cor_func(mods_ys)

  PartialAsso <- list(corr=MatCorr, rep_corr=rep_MatCorr, rep_SRs=rep_SRs,
                      fitted.models=fitted.models, data=data, mods_n=mods_n,
                      cor_func=cor_func, Marg_corr=Marg_corr)
  attr(PartialAsso, "arguments") <- c(association, method, resids.type,
                                      jitter, jitter.uniform.scale)
  attr(PartialAsso, "responses") <- responses
  attr(PartialAsso, "adjustments") <- adjustments
  attr(PartialAsso, "models") <- models

  class(PartialAsso) <- c("PAsso", class(PartialAsso))
  return(PartialAsso)

} ## end of function



#' @title Print partial association matrix
#' @param x A PAsso object for printing out results.
#'
#' @param digits A default number to specify decimal digit values.
#' @param ... Additional optional arguments.
#'
#' @name print
#' @method print PAsso
#'
#' @return Print partial association matrix of a PAsso object
#'
#' @export
#' @examples
#' # See PAsso for the example.
#'
print.PAsso <- function(x, digits = max(2, getOption("digits")-2), ...) {
  cat("-------------------------------------------- \n")
  cat("The partial correlation coefficient matrix: \n")

  # x$corr[lower.tri(x$corr)] <- NA

  # print.default(format(x$corr, digits = max(2, (digits))),
                # print.gap = 2, na.print = "",
                # quote = FALSE, ...)

  temp <- format(round(x$corr, digits=max(2, (digits))),
                 digits = max(2, (digits)), ...)
  temp[lower.tri(temp)] <- NA

  print.default(temp,
                print.gap = 2, na.print = "",
                quote = FALSE, ...)
}

#' @title Summary of partial association analysis
#' @description This function summarizes the partial association analysis by
#' providing partial association matrix, marginal association matrix, and a
#' matrix of the models' coefficients. The partial correlation coefficient
#' matrix displays the partial association between each pair of responses
#' after adjusting the covariates. While the marginal coefficient matrix displays association
#' before the adjustment.
#'
#' @param object A PAsso object to draw the summarized results, which includes partial association
#' matrix and a matrix of the models' coefficients.
#'
#' @param digits A default number to specify decimal digit values.
#' @param ... Additional optional arguments.
#'
#' @name summary
#' @method summary PAsso
#'
#' @return For a PAsso object, print its partial association matrix, marginal association
#' matrix, and a matrix of the models' coefficients.
#'
#' @export
#' @examples
#' # See PAsso for the example.
#'
summary.PAsso <- function(object, digits = max(3L, getOption("digits")-2L), ...) {
  cat("-------------------------------------------- \n")
  cat("The partial correlation coefficient matrix: \n\n")
  # print(signif(object$corr, ...))

  temp <- format(round(object$corr, digits=max(2, (digits))),
                 digits = max(2, (digits)), ...)
  temp[lower.tri(temp)] <- NA

  print.default(temp,
                print.gap = 2, na.print = "",
                quote = FALSE, ...)

  cat("-------------------------------------------- \n")
  cat("The marginal correlation coefficient matrix: \n\n")
  Marg_temp <- format(round(object$Marg_corr, digits=max(2, (digits))),
                      digits = max(2, (digits)), ...)
  Marg_temp[lower.tri(Marg_temp)] <- NA

  print.default(Marg_temp,
                print.gap = 2, na.print = "",
                quote = FALSE, ...)

  cat("\n--------------------------------------------\n")
  cat("--------------------------------------------\n")
  cat("\nThe coefficients of fitted models are: \n\n", sep = "")
  # print(object$fitted.models)

  # object=PAsso_1; digits=max(3L, getOption("digits")-2L)

  responses <- attr(object, "responses")
  adjustments <- attr(object, "adjustments")
  n_resp <- length(responses)
  n_adju <- length(adjustments)
  n_samp <- nobs(object$fitted.models[[1]])

  coefs_table <- matrix(NA, nrow = n_adju*3, ncol = n_resp)
  colnames(coefs_table) <- responses
  rownames_coefs_table <- rep(NA, n_adju*2)
  rownames_coefs_table[seq(1, n_adju*3, by = 3)] <- adjustments
  rownames_coefs_table[seq(2, n_adju*3, by = 3)] <- rep("Std. Error", n_adju)
  rownames_coefs_table[seq(3, n_adju*3, by = 3)] <- rep("---", n_adju)

  rownames(coefs_table) <- rownames_coefs_table

  for (i in 1:n_resp) {
    # i <- 1

    if (inherits(object$fitted.models[[i]], "polr")) {
      # Obtain results from the summary output
      sumry <- summary(object$fitted.models[[i]])$coefficients

      # Obtain coefficients and standard error
      coefs_se <- sumry[1:(n_adju), ]
      # Obtain p-value
      temp_p <- 2*pt(-abs(coefs_se[,3]), df = n_samp-1)
      `Pr` <- format.pval(temp_p, digits = max(1L, min(5L, digits - 1L)),
                  eps = .Machine$double.eps, ...)
    } else if (inherits(object$fitted.models[[i]], "glm")) {
      # Obtain results from the summary output
      sumry <- summary(object$fitted.models[[i]])$coefficients

      # Obtain coefficients and standard error
      coefs_se <- sumry[-1,1:3]
      # Obtain p-value
      temp_p <- sumry[-1,4]
      `Pr` <- format.pval(temp_p, digits = max(1L, min(5L, digits - 1L)),
                             eps = .Machine$double.eps, ...)
    } else if (inherits(object$fitted.models[[i]], "vglm")) {
      # Obtain results from the summary output of "vglm"

      sumry_temp <- summary(object$fitted.models[[i]])
      sumry <- slot(sumry_temp,"coef3")

      ncat_model <- ncat(object$fitted.models[[i]])
      # Obtain coefficients and standard error
      coefs_se <- sumry[-c(1:(ncat_model-1)),1:3]
      # Obtain p-value
      temp_p <- sumry[-c(1:(ncat_model-1)),4]
      `Pr` <- format.pval(temp_p, digits = max(1L, min(5L, digits - 1L)),
                          eps = .Machine$double.eps, ...)
    }
    Signif <- symnum(temp_p, corr = FALSE, na = FALSE,
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " "))

    coefs_se <- format(round(coefs_se, digits = max(2, (digits))),
                       digits = max(2, (digits)), ...)
    coefs_se <- cbind(coefs_se, `Pr`, Signif)

    coefs_table[seq(1, n_adju*3, by = 3),i] <- paste(coefs_se[,1], coefs_se[,5], sep = "")
    coefs_table[seq(2, n_adju*3, by = 3),i] <- coefs_se[,2]
  }
  print.default(coefs_table,
                print.gap = 2, na.print = "",
                quote = FALSE)

  # Add stars of significance level --------------------------------------------------
  if ((w <- getOption("width")) <
      nchar(sleg <- attr(Signif, "legend"))) {
    sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
  }
  cat("Signif. codes:  ", sleg, sep = "",
      fill = w + 4 + max(nchar(sleg, "bytes") - nchar(sleg)))
}



