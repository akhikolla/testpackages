#' @title Dose-Escalation Trials guided by Bayesian Logistic Regression Model
#'
#' @description \code{blrm_trial} facilitates the conduct of dose
#'     escalation studies guided by Bayesian Logistic Regression
#'     Models (BLRM). While the \code{blrm_exnex} only fits the BLRM
#'     model to data, the \code{blrm_trial} function standardizes the
#'     specification of the entire trial design and provides various
#'     standardized functions for trial data accrual and derivation of
#'     model summaries needed for dose-escalation decisions.
#'
#' @param data dose-toxicity data available at design stage of trial
#'
#' @param dose_info specificaion of the dose levels as
#'     planned for the ongoing trial arms.
#'
#' @param drug_info specification of drugs used in trial arms.
#'
#' @param simplified_prior logical (defaults to \code{FALSE})
#'     indicating whether a simplified prior should be employed based
#'     on the \code{reference_p_dlt} values provided in
#'     \code{drug_info}. \strong{Warning:} The simplified prior will
#'     change between releases. Please read instructions below in the
#'     respective section for the simplified prior.
#'
#' @param EXNEX_comp logical (default to \code{TRUE}) indicating
#'     whether EXchangeable-NonEXchangeable priors should be employed
#'     for all component parameters
#'
#' @param EX_prob_comp_hist prior weight (\eqn{[0,1]}, default to 1)
#'     on exchangeability for the component parameters in groups
#'     representing historical data
#'
#' @param EX_prob_comp_new prior weight (\eqn{[0,1]}, default to 0.8)
#'     on exchangeability for the component parameters in groups
#'     representing new or concurrent data
#'
#' @param EXNEX_inter logical (default to \code{FALSE}) indicating
#'     whether EXchangeable-NonEXchangeable priors should be employed
#'     for all interaction parameters
#'
#' @param EX_prob_inter prior weight (\eqn{[0,1]}, defaults to 0.8) on
#'     exchangeability for the interaction parameters
#'
#' @param formula_generator formula generation function (see for
#'     example \code{blrm_formula_linear} or 
#'     \code{blrm_formula_saturating}). The formula generator
#'     defines the employed interaction model.
#'
#' @param interval_prob defines the interval probabilities reported in
#'     the standard outputs. Defaults to \code{c(0, 0.16, 0.33, 1)}.
#'
#' @param interval_max_mass named vector defining for each interval of
#'     the \code{interval_prob} vector a maxmimal admissable
#'     probability mass for a given dose level. Whenever the posterior
#'     probability mass in a given interval exceeds the threshold,
#'     then the Escalation With Overdose Control (EWOC) criterion is
#'     considered to be not fullfilled. Dose levels not fullfilling
#'     EWOC are ineligible for the next cohort of patients. The
#'     default restricts the overdose probability to less than 0.25.
#'
#' @param ... Additional arguments are forwarded to \code{blrm_exnex},
#'     i.e. for the purpose of prior specification.
#' @param x \code{blrm_trial} object to print
#'
#' @details \code{blrm_trial} constructs an object of class
#'     \code{blrm_trial} which stores the compelte information about
#'     the study design of a dose-escalation trial. The study design
#'     is defined through the data sets (see sections below for a
#'     definition of the columns):
#'
#' \describe{
#'
#' \item{data (historical data)}{The \code{data} argument defines
#' available dose-toxicity data at the design stage of the
#' trial. Together with the prior of model (without any data) this
#' defines the prior used for the trial conduct.}
#'
#' \item{dose_info}{Definition of the pre-specified dose levels
#' explored in the ongoing trial arms. Thus, all dose-toxcitiy trial
#' data added to the object is expected correspond to one of the dose
#' levels in the pre-defined set of dose_info.}
#'
#' \item{drug_info}{Determines the drugs used in the trial, their
#' units, reference dose level and optionally defines the expected
#' probability for a toxicity at the reference dose.}
#'
#' }
#'
#' Once the \code{blrm_trial} object is setup the complete trial
#' design is specified and the model is fitted to the given
#' \code{data}. This allows evaluation of the pre-specified dose
#' levels of the trial design wrt. to safety, i.e. whether the
#' starting dose of the trial fullfills the escalate with overdose
#' criterion (EWOC) condition.
#'
#' The \code{blrm_trial} trial can also be constructed in a 2-step
#' process which allows for a more convenient specification of the
#' prior since meta data like number of drugs and the like can be
#' used. See the example section for details.
#'
#' After setup of the initial \code{blrm_trial} object additional data
#' is added through the use of the \code{update} method which has a
#' \code{add_data} argument intended to add data from the ongoing
#' trial. The \code{summary} function finally allows to extract
#' various model summaries. In particular, the EWOC criterion can be
#' calculated for the pre-defined dose levels of the trial.
#'
#' @section Simplified prior:
#'
#' As a convenience for the user, a simplified prior can be specifed
#' whenever the \code{reference_p_dlt} column is present in the
#' \code{drug_info} data set. However, the user is \strong{warned}
#' that the simplified prior will change in future releases of the
#' package and thus \strong{we strongly discourage the use of the
#' simplified prior for setting up trial designs}. The functionality
#' is intended to provide the user a quick start and as a starting
#' point. The actually instantiated prior can be seen as demonstrated
#' below in the examples.
#'
#' @section Input data:
#'
#' The data given to the \code{data} argument of \code{blrm_trial} is
#' considered as the available at design stage of the trial. The
#' collected input data thus does not necessarily need to have the
#' same dose levels as the pre-specified dose_info for the
#' ongoing trial(s). It's data columns must include, but are not
#' limited to:
#'
#' \describe{
#'   \item{group_id}{study}
#'
#'   \item{stratum_id}{optional, only required for differential discounting of groups}
#'
#'   \item{num_patients}{number of patients}
#'
#'   \item{num_toxicities}{number of toxicities}
#'
#'   \item{drug_A}{Columns for the dose of each treatment component,
#'   with column names matching the \code{drug_name} values specified
#'   in the \code{drug_info} argument of \code{blrm_trial}}
#'
#' }
#'
#' @section Drug info data:
#'
#' The drug information data-set defines drug properties. The fields
#' included are:
#'
#' \describe{
#'   \item{drug_name}{name of drug which is also used as column name for the dose}
#'   \item{dose_ref}{reference dose}
#'   \item{dose_unit}{units used for drug amounts}
#'   \item{reference_p_dlt}{optional; if provided, allows setup of a simplified prior}
#' }
#'
#' @section Dose info data:
#'
#' The \code{drug_info} data-set pre-specifies the dose levels of the
#' ongoing trial. Thus, all data added to the \code{blrm_trial}
#' through the \code{update} command must be consistent with the
#' pre-defined dose levels as no other than those pre-specified ones
#' can be explored in an ongoing trial.
#'
#' \describe{
#'
#' \item{dose_id}{optional column which assigns a unique id to each
#' group_id/dose combination. If not specified the column is
#' internally generated.}
#'
#' \item{group_id}{study}
#'
#' \item{drug_A}{Columns for the dose of each treatment component,
#'   with column names matching the \code{drug_name} values specified
#'   in the \code{drug_info} argument of \code{blrm_trial}}
#'
#' }
#'
#' @return The function returns an object of class \code{blrm_trial}.
#'
#' @family blrm_trial combo2 example
#'
#' @template ref-babb_ewoc
#' @template ref-mac
#'
#' @template start-example
#' @examples
#'
#' # construct initial blrm_trial object from built-in example datasets
#' combo2_trial_setup <- blrm_trial(
#'   data = hist_combo2,
#'   dose_info = dose_info_combo2,
#'   drug_info = drug_info_combo2,
#'   simplified_prior = TRUE
#' )
#'
#' # extract blrm_call to see setup of the prior as passed to blrm_exnex
#' summary(combo2_trial_setup, "blrm_exnex_call")
#'
#' # Warning: The simplified prior will change between releases!
#' # please refer to the combo2_trial example for a complete
#' # example. You can obtain this example with
#' # ?example-combo2_trial
#' # or by running
#' # example_model("combo2_trial")
#'
#' @template stop-example
#'
#' @export
blrm_trial <- function(
  data,
  dose_info,
  drug_info,

  simplified_prior = FALSE,
  EXNEX_comp = TRUE,
  EX_prob_comp_hist = 1,
  EX_prob_comp_new = 0.8,
  EXNEX_inter = FALSE,
  EX_prob_inter = 0.8,

  formula_generator = blrm_formula_linear,

  interval_prob = c(0, 0.16, 0.33, 1),
  interval_max_mass = c(prob_underdose = 1, prob_target=1, prob_overdose=0.25),

  ...
)
{
  trial <- list()

  # Prevent factor nonsense - enforce tibbles
  assert_tibble(data)
  assert_tibble(dose_info)
  assert_tibble(drug_info)

  assert_names(colnames(drug_info), must.include = c("drug_name", "dose_ref", "dose_unit"))

  # Assign drug names and reference doses
  assert_that(nrow(drug_info) > 0, msg="At least one drug must be defined")

  trial$ref_doses <- drug_info$dose_ref
  names(trial$ref_doses) <- drug_info$drug_name

  trial$dose_units <- drug_info$dose_unit
  trial$drug_info <- drug_info

  # Check data at t = 0

  # If no stratum is assigned, assign a single one
  if (!has_name(data, "stratum_id")) {
    data[["stratum_id"]] <- "all"
  }

  # If no historic data state is assigned, assign all given data to be historic
  # (= cohort time 0)
  if (!has_name(data, "cohort_time")) {
    data[["cohort_time"]] <- 0
  }

  assert_that(has_name(data, "stratum_id"))
  assert_that(has_name(data, "group_id"))
  assert_that(has_name(data, "num_patients"))
  assert_that(has_name(data, "num_toxicities"))

  # Data must have columns that correspond to drug name
  assert_that(has_name(data, trial$drug_info$drug_name), msg="Data must use same drug names as specified in drug_info")

  if (!has_name(data, "dose_id")) {
    data$dose_id <- integer(1)*NA ## Make sure dose_id is an integer column
  }

  # Check dose_info
  if (!has_name(dose_info, "stratum_id")) {
    message("No stratum defined - assigning all groups to single stratum \"all\"")
    dose_info[["stratum_id"]] <- "all"
  }

  assert_that(has_name(dose_info, "stratum_id"))
  assert_that(has_name(dose_info, "group_id"))
  assert_that(!has_name(dose_info, "num_patients"))
  assert_that(!has_name(dose_info, "num_toxicities"))

  # dose_info must have columns that correspond to drug name
  # Note: all() is needed for assertthat 0.2.0 compatibility
  assert_that(all(has_name(dose_info, trial$drug_info$drug_name)), msg="dose_info must use same drug names as specified in drug_info")

  if (!has_name(dose_info, "dose_id")) {
    # warning("dose_info do not contain \"dose_id\" - adding \"dose_id\" column")
    dose_info <- tibble::rowid_to_column(dose_info, "dose_id")
  }

  # Define helper functions
  merge_to_factor_levels <- function(df_1, df_2, variable_name)
  {
    if (is.factor(df_1[[variable_name]][0]) && is.factor(df_2[[variable_name]][0])){
      assert_that(
        nlevels(df_1[[variable_name]][0]) == nlevels(df_2[[variable_name]][0]) &&
        all(levels(df_1[[variable_name]][0]) == levels(df_2[[variable_name]][0])),
        msg = paste0("Factor levels in variable \"", variable_name, "\" inconsistent.")
      )
      all_factor_levels <- levels(df_1[[variable_name]][0])
    } else {
      all_level_df <- unique(
        bind_rows(
          select(df_1, variable_name),
          select(df_2, variable_name)
        )
      )
      all_level_strings <- all_level_df[[variable_name]]
      all_factor_levels <- factor(all_level_strings, levels=all_level_strings)
    }

    all_factor_levels
  }

  apply_factor_levels_to_variable <- function(df, variable_name, factor_levels)
  {
    df[[variable_name]] <- factor(df[[variable_name]], levels=factor_levels)

    df
  }

  # Apply consistent group_id and stratum_id factor levels to all data
  for (variable_name in c("group_id", "stratum_id"))
  {
    all_factor_levels <- merge_to_factor_levels(data, dose_info, variable_name)

    data <- apply_factor_levels_to_variable(data, variable_name, all_factor_levels)
    dose_info <- apply_factor_levels_to_variable(dose_info, variable_name, all_factor_levels)
  }



  # Assign data & dose_info
  trial$dose_info <- dose_info
  trial$data <- data
  trial$group_to_stratum_mapping <- .create_group_to_stratum_mapping(
    dose_info = dose_info,
    data = data,
    drug_names = trial$drug_info$drug_name
  )

  # Save group and stratum factor levels
  trial$group_id_factor <- trial$data[["group_id"]][0]
  trial$stratum_id_factor <- trial$data[["stratum_id"]][0]

  # Create formula for individual components
  assert_class(formula_generator, "function")

  trial$formula <- formula_generator(ref_doses = trial$ref_doses)
  assert_class(trial$formula, "blrm_formula")

  # Compute dimensionality of data and problem
  trial$num_components <- trial$formula$num_components
  trial$num_interaction_terms <- trial$formula$num_interaction_terms
  trial$num_strata <- nlevels(trial$stratum_id_factor)
  trial$num_groups <- nlevels(trial$group_id_factor)
  trial$num_groups_hist <- length(unique(trial$data[["group_id"]]))
  trial$num_groups_current <- length(unique(trial$dose_info[["group_id"]]))

  # Sanity check thresholds for EWOC
  assert_numeric(interval_max_mass, lower=0, upper=1, names="strict", any.missing=FALSE, finite=TRUE)
  assert_numeric(interval_prob, lower=0, upper=1, sorted=TRUE, any.missing=FALSE, finite=TRUE)

  assert_that(all(interval_max_mass >= 0) && all(interval_max_mass <= 1), msg="Maximum interval probability interval_max_mass must be 0 <= p <= 1")
  assert_that(all(interval_prob >= 0) && all(interval_prob <= 1) && all(diff(interval_prob) > 0), msg="interval_prob must be an ascending vector of probabilities")


  trial$interval_names <- NULL
  if (!is.null(names(interval_max_mass))) {
    trial$interval_names <- names(interval_max_mass)
  }
  assert_that((length(interval_max_mass) == length(interval_prob) - 1) || (length(interval_max_mass) == 0 && length(interval_prob) == 0), msg="Inconsistent interval number and maximal interval mass. The number of intervals should be length(interval_prob) - 1, or both interval_prob and interval_max_mass be empty.")
  trial$interval_max_mass <- interval_max_mass

  trial$interval_prob <- interval_prob

  trial <- structure(trial, class="blrm_trial")

  if (simplified_prior) {
    # Create simplified prior
    assert_that("reference_p_dlt" %in% colnames(trial$drug_info), msg='"reference_p_dlt" column is required for specifying simplified prior.')
    # Check reference dose toxicity probabilities
    assert_that(all(trial$drug_info$reference_p_dlt >= 0) && all(trial$drug_info$reference_p_dlt <= 1), msg="Reference dose DLT probability must be >= 0 and <= 1.")


    trial <- .blrm_trial_compute_simple_prior(
      object = trial,
      EXNEX_comp = EXNEX_comp,
      EX_prob_comp_hist = EX_prob_comp_hist,
      EX_prob_comp_new = EX_prob_comp_new,
      EXNEX_inter = EXNEX_inter,
      EX_prob_inter = EX_prob_inter,
      ...
    )

  } else {
    dot.args <- list(...)

    if (length(dot.args) > 0) {
      trial <- do.call(".blrm_trial_compute_prior", c(list(object=trial), dot.args))
    } else {
      message("Please configure blrm_exnex using the update() function.")
    }
  }

  trial
}

#'
#' @describeIn blrm_trial print function.
#' @method print blrm_trial
#' @export
print.blrm_trial <- function(x, ...)
{
  .assert_is_blrm_trial(x)
  pretty_formula <- gsub("~", "~\n", x$formula$blrm_formula)
  pretty_formula <- gsub("\\|", "|\n", pretty_formula)

  cat("Bayesian Logistic Regression Model trial with \n\n")

  cat("Number of observations   :", nrow(x$data), "\n")
  cat("Number of historic groups:", x$num_groups_hist, "\n")
  cat("Number of trial arms     :", x$num_groups_current, "\n")
  cat("Generated formula        :\n", pretty_formula, "\n\n")
  if (has_name(x, "blrmfit")) {
    print(x$blrmfit)
  } else {
    cat("BLRM prior undefined!\n")
  }
}



# Internal functions ------------------------------------------------------

.blrm_trial_compute_ewoc <- function(trial, trial_est) {
  .assert_is_blrm_trial(trial)

  trial_est_ewoc <- trial_est

  if (length(trial$interval_max_mass) > 0) {
    if (!is.null(trial$interval_names)) {
      cols <- colnames(trial_est_ewoc)
      cols[(length(cols)-length(trial$interval_names)+1):length(cols)] <- trial$interval_names
      colnames(trial_est_ewoc) <- cols
      interval_max_mass <- trial$interval_max_mass
    } else {
      interval_max_mass <- trial$interval_max_mass
      names(interval_max_mass) <- paste0("I", seq(1, length(interval_max_mass)))
    }
    trial_est_intervals <- trial_est_ewoc[, (length(colnames(trial_est_ewoc))-length(interval_max_mass)+1):length(colnames(trial_est_ewoc)) ]

    interval_max <- as_tibble(as.list(interval_max_mass))
    interval_max <- interval_max[rep(1, times=nrow(trial_est_ewoc)),]

    trial_est_ewoc$ewoc_ok <- rowSums( trial_est_intervals > interval_max ) == 0
  } else {
    trial_est_ewoc$ewoc_ok <- TRUE
  }

  trial_est_ewoc
}

.blrm_trial_predict <- function(trial, newdata, ...)
{
  .assert_is_blrm_trial_and_prior_is_set(trial)
  assert_tibble(newdata)

  if (length(trial$interval_prob) > 0)
  {
    trial_est <- summary(trial$blrmfit,
                         newdata = newdata,
                         interval_prob = trial$interval_prob,
                         ...
    )
  } else {
    trial_est <- summary(trial$blrmfit,
                         newdata = newdata,
                         ...
    )
  }


  trial_est_ewoc <- .blrm_trial_compute_ewoc(trial, trial_est)

  common_column_names <- intersect(colnames(trial_est_ewoc), colnames(newdata))
  newdata_prediction <- newdata

  if (length(common_column_names) > 0)
  {
    newdata_prediction <- select(newdata_prediction, -common_column_names)
  }

  newdata_prediction <- cbind(newdata_prediction, trial_est_ewoc)

  #
  ## Enrichment cohort predictive distribution - compute probability of P(>= 2 DLTs in cohort of size 10)
  #enrichment_data <- mutate(newdata, num_toxicities = 0, num_patients = 10)
  #enrichment_cohort_DLTs <- posterior_predict(blrmfit, newdata=enrichment_data)
  #
  # Enrichment is okay if P(>= 2 DLTs in cohort of size 10) <= 30%
  #enrichment_is_okay <- colMeans(enrichment_cohort_DLTs >= 2) <= 0.3
  #newdata_prediction <- mutate(newdata_prediction, enrichment_ok = enrichment_is_okay)

  as_tibble(newdata_prediction)
}



.blrm_trial_compute_simple_prior <- function(object,
  EXNEX_comp,
  EX_prob_comp_hist,
  EX_prob_comp_new,
  EXNEX_inter,
  EX_prob_inter, # TODO: Make EX_prob_inter for hist. data = 1 and new data .8?
  ...
)
{
  .assert_is_blrm_trial(object)
  trial <- object

  group_id <- 0 ## Suppress group_id binding warning

  ## Should not be used with more than one stratum
  assert_that(trial$num_strata == 1, msg="Simplified prior can only be set up with a single stratum. Omit the reference_p_dlt argument and set the prior up yourself using update().")

  ## This is primarily a help to get started with blrm_trial and should not be used in production - warn the user:
  warning("Simplified prior CAN and WILL change with releases. NOT recommended to use in production. Instantiating a simplified prior - run summary(trial, \"blrm_exnex_call\") to inspect arguments. ")

  ## Prior mean and sd on log mu_{alpha_i}, log mu_{beta_i}
  ref_p_dlt  <- trial$drug_info$reference_p_dlt
  names(ref_p_dlt)  <- trial$drug_info$drug_name
  prior_EX_mu_mean_comp  <- substitute(cbind(logit(refs), rep(0, trial[["num_components"]])), list(refs=ref_p_dlt))
  prior_EX_mu_sd_comp    <- substitute(matrix(c(2, 1), nrow = trial[["num_components"]], ncol = 2, byrow=TRUE))

  ## Prior mean and sd on tau_{alpha_{s,i}}, tau{beta_{s,i}}
  prior_EX_tau_mean_comp <- substitute(matrix(log(c(0.5, 0.25)), nrow = trial[["num_components"]], ncol = 2, byrow=TRUE))
  prior_EX_tau_sd_comp <- substitute(matrix(log(4)/1.96, nrow = trial[["num_components"]], ncol = 2))

  # Prior mean and sd on mu_{eta}
  prior_EX_mu_mean_inter  <- substitute(rep(0,   trial[["num_interaction_terms"]]))
  prior_EX_mu_sd_inter    <- substitute(rep(0.5, trial[["num_interaction_terms"]]))
  prior_EX_tau_mean_inter <- substitute(matrix(log(0.5)   , nrow = trial[["num_strata"]], ncol = trial[["num_interaction_terms"]]))
  prior_EX_tau_sd_inter   <- substitute(matrix(log(4)/1.96, nrow = trial[["num_strata"]], ncol = trial[["num_interaction_terms"]]))

  prior_is_EXNEX_comp <- rep(EXNEX_comp, trial[["num_components"]])

  ## By default, historical data have
  prior_EX_prob_comp <- matrix(NA, nrow = 0, ncol = trial[["num_components"]])
  ## Go in order of the factor levels
  for (group_name in levels(trial$group_id_factor)) {
    if (nrow(filter(trial$data, group_id == group_name)) > 0 &&
        nrow(filter(trial$dose_info, group_id == group_name)) == 0 ) {
      ## There is exclusively historic data
      prior_EX_prob_comp <- rbind(prior_EX_prob_comp, matrix(EX_prob_comp_hist, nrow = 1, ncol = trial[["num_components"]]))
    } else {
      ## There is no historic data for this group, or there is also new data expected
      prior_EX_prob_comp <- rbind(prior_EX_prob_comp, matrix(EX_prob_comp_new, nrow = 1, ncol = trial[["num_components"]]))
    }
  }

  prior_is_EXNEX_inter <- rep(EXNEX_inter, trial[["num_interaction_terms"]])
  prior_EX_prob_inter <- matrix(EX_prob_inter, nrow = trial[["num_groups"]], ncol = trial[["num_interaction_terms"]])

  prior_tau_dist <- 1
  ## Compute prior and return resulting trial object
  update(
    object = trial,

    prior_EX_mu_mean_comp = prior_EX_mu_mean_comp,
    prior_EX_mu_sd_comp = prior_EX_mu_sd_comp,

    prior_EX_tau_mean_comp = prior_EX_tau_mean_comp,
    prior_EX_tau_sd_comp = prior_EX_tau_sd_comp,

    prior_EX_mu_mean_inter = prior_EX_mu_mean_inter,
    prior_EX_mu_sd_inter = prior_EX_mu_sd_inter,

    prior_EX_tau_mean_inter = prior_EX_tau_mean_inter,
    prior_EX_tau_sd_inter = prior_EX_tau_sd_inter,

    prior_is_EXNEX_comp = prior_is_EXNEX_comp,
    prior_EX_prob_comp = prior_EX_prob_comp,

    prior_is_EXNEX_inter = prior_is_EXNEX_inter,
    prior_EX_prob_inter = prior_EX_prob_inter,

    prior_tau_dist = prior_tau_dist,

    ...
  )
}


.blrm_trial_compute_prior <- function(object, ...) {
  .assert_is_blrm_trial(object)
  trial <- object

  dot.args <- list(...) # Evaluate arguments here
  ref_doses <- trial[["ref_doses"]]
  formula <- as.formula(trial[["formula"]][["blrm_formula"]])
  data <- .blrm_trial_merge_data(trial, bind_rows(trial[["group_to_stratum_mapping"]], trial[["data"]]))
  ## SW: the as.name is something I picked up here:
  ## https://github.com/WinVector/wrapr/blob/master/extras/MacrosInR.md
  ## (just delete this comment if you fine with this change)
  dot.args <- modifyList(list(formula=formula, data=as.name("data")), dot.args)

  trial$blrmfit <- do.call("blrm_exnex", dot.args)

  trial$update_blrmfit <- function(trial, ...)
  {
    .assert_is_blrm_trial_and_prior_is_set(trial)

    args <- list(...)
    if (has_name(args, "data")) {
      args[["data"]] <- .blrm_trial_sanitize_data(trial=trial, data=args[["data"]])
      trial$data <- args[["data"]]
      args[["data"]] <- bind_rows(trial$group_to_stratum_mapping, args[["data"]])
    }

    if (has_name(args, "add_data")) {
      args[["add_data"]] <- .blrm_trial_sanitize_data(trial=trial, data=args[["add_data"]])

      if (!has_name(args[["add_data"]], "cohort_time")) {
        next_index <- max(select(trial$data, "cohort_time"))
        if (is.na(next_index)) {
          next_index <- 1
        } else {
          next_index <- next_index + 1
        }

        args[["add_data"]]$cohort_time <- next_index
      }

      trial$data <- bind_rows(trial$data, args[["add_data"]])
      args[["data"]] <- bind_rows(trial$group_to_stratum_mapping, trial$data)
      args[["add_data"]] <- NULL
    }

    if (has_name(args, "data")) {
      args[["data"]] <- .blrm_trial_merge_data(trial, args[["data"]])
    }

    trial$blrmfit <- do.call("update", c(list(trial$blrmfit), args))

    trial
  }

  trial
}

.blrm_trial_merge_data <- function(trial, data)
{
  .assert_is_blrm_trial(trial)
  assert_tibble(data)

  ## Silence R CRAN check notes
  num_patients <- num_toxicities <- NULL

  ## Summarize data for performance
  summarized_data <- ungroup(summarize(group_by_at(data, c("group_id", "stratum_id", names(trial$ref_doses))), num_patients = sum(num_patients), num_toxicities = sum(num_toxicities)))

  ## Return sorted data for improved reproducibility
  arrange_at(summarized_data, c("group_id", "stratum_id", names(trial$ref_doses)))
}

.blrm_trial_sanitize_data <- function(trial, data)
{
    .assert_is_blrm_trial(trial)
    assert_tibble(data)

    ## make R CMD check happy
    dose_id <- NULL

  if (!has_name(data, "stratum_id")) {
    if (nlevels(trial$stratum_id_factor) == 1) {
      data$stratum_id <- levels(trial$stratum_id_factor)[1]
      message("stratum_id not given, but only one stratum defined. Assigning first stratum.")
    } else {
      stop("stratum_id not defined even though there are multiple strata!")
    }
  }

  assert_that(has_name(data, "stratum_id"))
  assert_that(has_name(data, "group_id"))

  if (!is.factor(data$stratum_id[0])) {
    data$stratum_id <- factor(data$stratum_id, levels=levels(trial$stratum_id_factor))
    assert_that(!any(is.na(data$stratum_id)), msg="stratum_id level inconsistent or NA")
  }

  assert_factor(data$stratum_id[0], levels=levels(trial$stratum_id_factor), n.levels=nlevels(trial$stratum_id_factor))

  if (!is.factor(data$group_id[0])) {
    data$group_id <- factor(data$group_id, levels=levels(trial$group_id_factor))
    assert_that(!any(is.na(data$group_id)), msg="group_id level inconsistent or NA")
  }
  assert_factor(data$group_id[0], levels=levels(trial$group_id_factor), n.levels=nlevels(trial$group_id_factor))

  assert_that(has_name(data, "num_patients"))
  assert_that(has_name(data, "num_toxicities"))

  ## Data must have columns that correspond to drug name
  ## Note: all() is needed for assertthat 0.2.0 compatibility
  assert_that(all(has_name(data, colnames(trial$ref_doses))))

  ## Check that dose_id is consistent
  if(has_name(data, "dose_id")) {
    data_consistent_with_dose_info <- inner_join(data, trial$dose_info, by=c("dose_id", "group_id", "stratum_id", trial$drug_info$drug_name))
    assert_that(nrow(filter(data, !is.na(dose_id))) == nrow(data_consistent_with_dose_info),
                msg="dose_id inconsistent with dose combinations")

    if (any(is.na(data$dose_id))) {
      warning("dose_id NA was provided - cannot check consistency of new data with pre-specified dose_info")
    }
  } else {
    ## If no dose_id is provided, resolve it
    data_consistent_with_dose_info <- left_join(data, trial$dose_info, by=c("group_id", "stratum_id", trial$drug_info$drug_name))
    if(any(is.na(data_consistent_with_dose_info$dose_id))) {
      warning("Data that was provided does not correspond to a pre-specified dose!")
    }
    data <- data_consistent_with_dose_info
  }

  data
}


.create_group_to_stratum_mapping <- function(
  dose_info,
  data,
  drug_names
)
{
  assert_tibble(dose_info)
  assert_tibble(data)

  ## Create dummy data to map group_ids to stratum_ids in the model
  group_to_stratum_mapping <- unique(select(bind_rows(dose_info, data), one_of("group_id", "stratum_id")))

  ## Add group to stratum mapping so groups are mapped correctly even
  ## if they do not contain data (not yet)
  group_to_stratum_mapping <- bind_rows(data[0,], group_to_stratum_mapping)
  group_to_stratum_mapping <- mutate_at(group_to_stratum_mapping, vars("num_patients", "num_toxicities"), function(x)(0))

  group_to_stratum_mapping <- mutate_at(group_to_stratum_mapping, vars(drug_names), function(x)(1))

  group_to_stratum_mapping
}

.assert_is_blrm_trial_and_prior_is_set <- function(object)
{
  .assert_is_blrm_trial(object)
  assert_that(!is.null(object$blrmfit), msg="Prior must be specified. Use update() with arguments for blrm_exnex to configure the prior.")
}

.assert_is_blrm_trial <- function(object)
{
  assert_class(object, "blrm_trial")
}


