
context("blrm_trial tests")

set.seed(123144)

eps <- 1E-4


suppressPackageStartupMessages(library(dplyr))

## check that a blrm_trial object is produced
check_basic  <- function(example){
  expect_class(blrm_trial(example$histdata, example$dose_info, example$drug_info),
              "blrm_trial")
}

test_that("blrm_trial runs for a basic single-agent example",
          check_basic(examples$single_agent))

test_that("blrm_trial runs for a basic combo2 example",
          check_basic(examples$combo2))

test_that("blrm_trial runs for a basic combo3 example",
          check_basic(examples$combo3))

test_that("blrm_trial runs for a basic single-agent example with multiple strata",
          check_basic(examples$single_drug_with_strata))


# Check printing works ----------------------------------------------------
check_print_trial  <- function(example) {
  print(blrm_trial(example$histdata, example$dose_info, example$drug_info))
}

check_print_trial_with_prior  <- function(example) {
  print(blrm_trial(
    example$histdata,
    example$dose_info,
    example$drug_info,
    simplified_prior = TRUE
  ))
}

test_that("blrm_trial prints without prior",
          check_print_trial(examples$single_agent))

test_that("blrm_trial prints with prior",
          check_print_trial_with_prior(examples$single_agent))


# Check summary works for data & dose_info --------------------------------
check_trial_summary  <- function(example) {
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info)

    expect_list(summary(trial, summarize="dimensionality"))
    expect_tibble(summary(trial, summarize="data"))
    expect_tibble(summary(trial, summarize="dose_info"))
  })
}

check_trial_with_prior_summary  <- function(example) {
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info, simplified_prior = TRUE)

    expect_list(summary(trial, summarize="dimensionality"))
    expect_tibble(summary(trial, summarize="data"))
    expect_tibble(summary(trial, summarize="dose_info"))

    expect_data_frame(summary(trial, summarize="blrmfit"))
    expect_class(summary(trial, summarize="blrm_exnex_call"), "call")
    expect_tibble(summary(trial, summarize="data_prediction"))
    expect_tibble(summary(trial, summarize="dose_prediction"))
  })
}

test_that("blrm_trial summary without prior",
          check_trial_summary(examples$single_agent))

test_that("blrm_trial summary with prior",
          check_trial_with_prior_summary(examples$single_agent))


## Check for consistent behavior on group_id -------------------------------

## check that factor levels can be predefined, and the order of the levels is respected
check_factor_accepted <- function(example) {
  with(example, {
    # reverse the levels of the factor; separate levels for the two df's
    levels <- c(paste0("hist_", c("B", "A")), paste0("cur_", c("B", "A")))

    histdata$group_id <- factor(histdata$group_id, levels)
    dose_info$group_id <- factor(dose_info$group_id, levels)

    ## blrm_trial should accept group_id factors with consistent levels
    trial <- blrm_trial(histdata, dose_info, drug_info)

    ## make sure blrm_trial has not switched the order of factor levels
    expect_equal(levels(trial$group_id_factor), levels)
  })
}

## check that if factor levels are used, they must be consistent (error otherwise)
check_factor_not_inconsistent  <- function(example) {
  with(example, {
    # reverse the levels of the factor; separate levels for the two df's
    levels <- c(paste0("hist_", c("B", "A")), paste0("cur_", c("B", "A")))

    # Separate levels should error

    histdata$group_id <- factor(histdata$group_id, paste0("hist_", c("B", "A")))
    dose_info$group_id <- factor(dose_info$group_id, paste0("cur_", c("B", "A")))
    expect_error(blrm_trial(histdata, dose_info, drug_info), ".*inconsistent.*")
  })
}

## check ordering of factor levels
check_factor_group_id_reorder <- function(example) {
  with(example, {
    # reverse the levels of the factor; separate levels for the two df's
    levels <- c(paste0("hist_", c("B", "A")), paste0("cur_", c("B", "A")))
    levels_permuted <- c(paste0("hist_", c("A", "B")), paste0("cur_", c("A", "B")))
    # Separate levels should error
    histdata$group_id <- factor(histdata$group_id, levels)
    dose_info$group_id <- factor(dose_info$group_id, levels_permuted)
    expect_error(blrm_trial(histdata, dose_info, drug_info), ".*inconsistent.*")
  })
}


## check behavior when both are characters
check_character_group_id  <- function(example) {
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info)
    levs <- unique(c(histdata$group_id, dose_info$group_id))
    expect_equal(levels(trial$group_id_factor), levs)

    # reorder histdata and dose_info
    histdata <- histdata %>% arrange(desc(group_id))
    dose_info <- dose_info %>% arrange(desc(group_id))

    trial <- blrm_trial(histdata, dose_info, drug_info)
    levs <- unique(c(histdata$group_id, dose_info$group_id))
    expect_equal(levels(trial$group_id_factor), levs)
  })
}

## check behavior when histdata contains unobserved factor levels
check_factor_group_id_unobserved  <- function(example){
  with(example, {
    # check when there are unobserved levels in histdata$group_id
    levels <- c(paste0("hist_", c("A", "B", "C")), paste0("cur_", c("A", "B")))

    histdata$group_id <- factor(histdata$group_id, levels)
    dose_info$group_id <- factor(dose_info$group_id, levels)
    trial <- blrm_trial(histdata, dose_info, drug_info)

    # make sure blrm_trial did not remove the unobserved levels
    expect_equal(levels(trial$group_id_factor), levels)
  })
}


test_that("blrm_trial accepts group_id factor levels",
          check_factor_accepted(examples$combo2))

test_that("blrm_trial errors if inconsistent factor levels for group_id are input in data / dose_info",
          check_factor_not_inconsistent(examples$combo2))

test_that("blrm_trial does not reorder group_id factor levels from input group_id's",
          check_factor_group_id_reorder(examples$combo2))

test_that("blrm_trial handles group_id's as expected when they are input as characters",
          check_character_group_id(examples$combo2))

test_that("blrm_trial does not drop unobserved factor levels from input group_id's",
          check_factor_group_id_unobserved(examples$combo3))


## Check for consistent behavior on stratum_id -----------------------------


## check that factor levels can be predefined, and the order of the levels is respected
check_stratum_factor_accepted <- function(example) {
  with(example, {
    # reverse the levels of the factor; separate levels for the two df's
    levels <- c(paste0("strat_", c("B", "A")))

    histdata$stratum_id <- factor(histdata$stratum_id, levels)
    dose_info$stratum_id <- factor(dose_info$stratum_id, levels)

    ## blrm_trial should accept stratum_id factors with consistent levels
    trial <- blrm_trial(histdata, dose_info, drug_info)

    ## make sure blrm_trial has not switched the order of factor levels
    expect_equal(levels(trial$stratum_id_factor), levels)
  })
}

## check that if factor levels are used, they must be consistent (error otherwise)
check_stratum_factor_not_inconsistent  <- function(example) {
  with(example, {
    # reverse the levels of the factor; separate levels for the two df's
    levels <- c(paste0("strat_", c("A", "B", "C", "D")))
    levels_reversed <- c(paste0("strat_", c("A", "B")))

    # Separate levels should error

    histdata$stratum_id <- factor(histdata$stratum_id, levels)
    dose_info$stratum_id <- factor(dose_info$stratum_id, levels_reversed)
    expect_error(blrm_trial(histdata, dose_info, drug_info), ".*inconsistent.*")
  })
}

## check ordering of factor levels
check_stratum_factor_group_id_reorder <- function(example) {
  with(example, {
    # reverse the levels of the factor; separate levels for the two df's
    levels <- c(paste0("strat_", c("B", "A")))
    levels_permuted <- c(paste0("strat_", c("A", "B")))

    # Separate levels should error
    histdata$stratum_id <- factor(histdata$stratum_id, levels)
    dose_info$stratum_id <- factor(dose_info$stratum_id, levels_permuted)

    expect_error(blrm_trial(histdata, dose_info, drug_info), ".*inconsistent.*")
  })
}


## check behavior when both are characters
check_character_stratum_id  <- function(example){
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info)
    levs <- unique(c(histdata$stratum_id, dose_info$stratum_id))
    expect_equal(levels(trial$stratum_id_factor), levs)

    # reorder histdata and dose_info
    histdata <- histdata %>% arrange(desc(stratum_id))
    dose_info <- dose_info %>% arrange(desc(stratum_id))

    trial <- blrm_trial(histdata, dose_info, drug_info)
    levs <- unique(c(histdata$stratum_id, dose_info$stratum_id))
    expect_equal(levels(trial$stratum_id_factor), levs)
  })
}

## check behavior when histdata contains unobserved factor levels
check_factor_stratum_id_unobserved  <- function(example){
  with(example, {

    # check when there are unobserved levels in histdata$group_id
    levels <- c(paste0("strat_", c("A", "B", "C")))

    histdata$stratum_id <- factor(histdata$stratum_id, levels)
    dose_info$stratum_id <- factor(dose_info$stratum_id, levels)
    trial <- blrm_trial(histdata, dose_info, drug_info)

    # make sure blrm_trial did not remove the unobserved levels
    expect_equal(levels(trial$stratum_id_factor), levels)
  })
}


test_that("blrm_trial accepts stratum_id factor levels",
          check_stratum_factor_accepted(examples$single_drug_with_strata))

test_that("blrm_trial errors if inconsistent factor levels for stratum_id are input in data / dose_info",
          check_stratum_factor_not_inconsistent(examples$single_drug_with_strata))

test_that("blrm_trial does not reorder stratum_id factor levels from input stratum_ids",
          check_stratum_factor_group_id_reorder(examples$single_drug_with_strata))

test_that("blrm_trial handles stratum_ids as expected when they are input as characters",
          check_character_stratum_id(examples$single_drug_with_strata))

test_that("blrm_trial does not drop unobserved factor levels from input stratum_ids",
          check_factor_stratum_id_unobserved(examples$single_drug_with_strata))

# Check bind_rows_0  ------------------------------------------------------
check_bind_rows_0 <- function() {
  example_A <- tibble(
    group_id = "hist_A",
    drug1 = 1
  )

  example_A_full <- tibble(
    group_id = "hist_A",
    drug1 = 1,
    drug2 = 0
  )

  example_B <- tibble(
    group_id = "hist_B",
    drug2 = 100 * (1:2)
  )

  example_B_full <- tibble(
    group_id = "hist_B",
    drug1 = 0,
    drug2 = 100 * (1:2)
  )

  ## Check that bind_rows_0 creates non-existing columns and sets their value to 0
  example_bound <- bind_rows_0(example_A, example_B, example_A)

  expect_equal(example_bound, rbind(example_A_full, example_B_full, example_A_full))

  ## Check for error in case of NAs
}




check_bind_rows_0_error_on_NA <- function() {
  example_A <- tibble(
    group_id = "hist_A",
    drug1 = NA
  )

  ## Check for error in case of NAs present
  expect_error(bind_rows_0(example_A))
}

test_that("bind_rows_0 sets non-existing columns to 0 and merges tibbles correctly",
          check_bind_rows_0())

test_that("bind_rows_0 errors in case some provided elements are NA",
          check_bind_rows_0_error_on_NA())


# Test prior_summary ------------------------------------------------------

check_blrm_trial_prior_summary <- function(example)
{
  with(example, {
    # Test on trial that has no prior specified
    trial_no_prior <- blrm_trial(histdata, dose_info, drug_info)
    expect_error(prior_summary(trial_no_prior))

    trial_with_prior <- blrm_trial(histdata, dose_info, drug_info, simplified_prior = TRUE)
    ps <- prior_summary(trial_with_prior)
    expect_class(ps, "prior_summary.blrm_trial")
    expect_class(ps$prior_summary.blrmfit, "prior_summary.blrmfit")
  })
}

test_that("prior_summary calls blrmfita after setting up simplified prior",
          check_blrm_trial_prior_summary(examples$single_agent))


# Test an example with update, adding data --------------------------------
check_blrm_trial_update <- function(example, formula_generator)
{
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info, simplified_prior = TRUE, formula_generator = formula_generator)

    new_data <- filter(summary(trial, "dose_info"), dose_id == 1)
    new_data$num_patients <- 6
    new_data$num_toxicities <- 1

    trial_with_new_data <- update(trial, add_data=new_data)

    new_trial_data <- summary(trial_with_new_data, "data", prob = c(0.95, 0.5))

    # Data should have increased by one row
    expect_equal(nrow(new_trial_data), nrow(histdata) + nrow(new_data))

    # This row should be the original row (plus potentially additional information such as cohort_time)
    expect_equal(new_data, select(new_trial_data[nrow(new_trial_data),], colnames(new_data)))

    trial_with_replaced_data <- update(trial, data=new_trial_data)
    expect_equal(new_trial_data, summary(trial_with_replaced_data, "data"))
    
    # Test data aggregation works correctly
    trial_with_new_data_2 <- update(trial_with_new_data, add_data=new_data)
    new_trial_data_2 <- summary(trial_with_new_data_2, "data", prob = c(0.95, 0.5))
    
    # Data should have increased by one row
    expect_equal(nrow(new_trial_data_2), nrow(new_trial_data) + nrow(new_data))   
    
    # Check for internal consistency with blrmfit data
    blrmfit_data <- trial_with_new_data$blrmfit$data
    blrmfit_data_2 <- trial_with_new_data_2$blrmfit$data
    internal_data <- filter(blrmfit_data, group_id == new_data$group_id, drug1 == new_data$drug1)
    internal_trial_data <- filter(new_trial_data, group_id == new_data$group_id, drug1 == new_data$drug1)
    expect_equal(internal_trial_data$num_patients, new_data$num_patients)
    expect_equal(internal_trial_data$num_toxicities, new_data$num_toxicities)
    expect_equal(internal_data$num_patients, new_data$num_patients)
    expect_equal(internal_data$num_toxicities, new_data$num_toxicities)
    
    internal_data_2 <- filter(blrmfit_data_2, group_id == new_data$group_id, drug1 == new_data$drug1)
    internal_trial_data_2 <- filter(new_trial_data_2, group_id == new_data$group_id, drug1 == new_data$drug1)
    expect_equal(sum(internal_trial_data_2$num_patients), 2*new_data$num_patients)
    expect_equal(sum(internal_trial_data_2$num_toxicities), 2*new_data$num_toxicities)
    expect_equal(internal_data_2$num_patients, 2*new_data$num_patients)
    expect_equal(internal_data_2$num_toxicities, 2*new_data$num_toxicities)
    
    # Number of rows should not increase if same data is added twice
    expect_equal(nrow(blrmfit_data_2), nrow(blrmfit_data))   
  })
}

test_that("update function for blrm_trial adds data with add_data= and replaces data with data= as expected",
          check_blrm_trial_update(examples$single_agent, blrm_formula_linear))



# Test EWOC criterion -----------------------------------------------------
check_ewoc_criterion_defaults <- function(example) {
  with(example, {
    ## create basic blrm trial
    trial <- blrm_trial(histdata, dose_info, drug_info)


    # Test EWOC criterion is computed correctly
    interval_pred <- summary(trial, "data")

    # Check for reasonable default underdose probability behavior
    interval_pred_underdose <- interval_pred
    interval_pred_underdose[["prob_underdose"]] <- 1
    interval_pred_underdose[["prob_target"]] <- 0
    interval_pred_underdose[["prob_overdose"]] <- 0
    interval_pred_underdose <- OncoBayes2:::.blrm_trial_compute_ewoc(trial, interval_pred_underdose)
    expect_true(all(interval_pred_underdose$ewoc_ok))

    # Check for reasonable default target probability behavior
    interval_pred_target <- interval_pred
    interval_pred_target[["prob_underdose"]] <- 0
    interval_pred_target[["prob_target"]] <- 1
    interval_pred_target[["prob_overdose"]] <- 0
    interval_pred_target <- OncoBayes2:::.blrm_trial_compute_ewoc(trial, interval_pred_target)
    expect_true(all(interval_pred_target$ewoc_ok))

    # Check for reasonable default overdose probability behavior
    interval_pred_target <- interval_pred
    interval_pred_target[["prob_underdose"]] <- 0
    interval_pred_target[["prob_target"]] <- 0
    interval_pred_target[["prob_overdose"]] <- 1
    interval_pred_target <- OncoBayes2:::.blrm_trial_compute_ewoc(trial, interval_pred_target)
    expect_true(!any(interval_pred_target$ewoc_ok))
  })
}

check_ewoc_criterion_intervals <- function(example) {
  with(example, {
    for (named_interval in c(TRUE))
    {
      ## create basic blrm trial
      interval_prob <- c(0, 0.16, 0.33, 1)
      if (named_interval) {
        interval_max_mass <- c(prob_underdose = 1, prob_target=1, prob_overdose=0.25)
      } else {
        interval_max_mass <- c(1, 1, 0.25)
      }

      trial <- blrm_trial(
        histdata, dose_info, drug_info,
        interval_prob = interval_prob,
        interval_max_mass = interval_max_mass
      )

      # Test EWOC criterion is computed correctly
      interval_pred <- summary(trial, "data")

      # Check EWOC OK if barely not overdosing
      barely_not_overdosing <- interval_max_mass[[3]] * 0.95
      interval_pred_ok <- interval_pred
      interval_pred_ok[["prob_underdose"]] <- 1 - barely_not_overdosing
      interval_pred_ok[["prob_target"]] <- 0
      interval_pred_ok[["prob_overdose"]] <- barely_not_overdosing
      interval_pred_ok <- OncoBayes2:::.blrm_trial_compute_ewoc(trial, interval_pred_ok)
      expect_true(all(interval_pred_ok$ewoc_ok))

      # Check EWOC OK if barely overdosing
      barely_overdosing <- interval_max_mass[[3]] * 1.05
      interval_pred_not_ok <- interval_pred
      interval_pred_not_ok[["prob_underdose"]] <- 1 - barely_overdosing
      interval_pred_not_ok[["prob_target"]] <- 0
      interval_pred_not_ok[["prob_overdose"]] <- barely_overdosing
      interval_pred_not_ok <- OncoBayes2:::.blrm_trial_compute_ewoc(trial, interval_pred_not_ok)
      expect_true(!any(interval_pred_not_ok$ewoc_ok))
    }
  })
}

check_ewoc_criterion_empty_intervals <- function(example) {
  with(example, {
    ## create blrm trial with empty intervals
    interval_prob <- c(1)
    interval_max_mass <- numeric(0)

    trial <- blrm_trial(
      histdata, dose_info, drug_info,
      interval_prob = interval_prob,
      interval_max_mass = interval_max_mass
    )

    # Test EWOC criterion is computed correctly
    interval_pred <- summary(trial, "data")
    expect_true(all(interval_pred$ewoc_ok))

    interval_prob <- numeric(0)
    trial <- blrm_trial(
      histdata, dose_info, drug_info,
      interval_prob = interval_prob,
      interval_max_mass = interval_max_mass
    )

    # Test EWOC criterion is computed correctly
    interval_pred2 <- summary(trial, "data")
    expect_true(all(interval_pred2$ewoc_ok))
  })
}

test_that("EWOC defaults for 100% underdosing / target / overdosing probability are sane",
          check_ewoc_criterion_defaults(examples$single_agent))

test_that("EWOC interval limits are respected",
          check_ewoc_criterion_intervals(examples$single_agent))

test_that("EWOC intervals can be empty",
          check_ewoc_criterion_empty_intervals(examples$single_agent))


# Test update function with inconsistent dose_id, group_id, stratum_id --------------------
check_update_group_stratum_id_char <- function(example) {
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info, simplified_prior = TRUE)

    new_data <- filter(summary(trial, "dose_info"), dose_id == 1)
    new_data$num_patients <- 6
    new_data$num_toxicities <- 1

    # Providing group_id and stratum_id as characters should not error and resolve correctly
    new_data_char_factors <- new_data
    new_data_char_factors$group_id <- as.character(new_data_char_factors$group_id)
    new_data_char_factors$stratum_id <- as.character(new_data_char_factors$stratum_id)

    trial_with_new_data <- update(trial, add_data=new_data_char_factors)

    new_data_in_trial <- summary(trial_with_new_data, summarize="data")
    expect_equal(new_data$group_id, new_data_in_trial[nrow(new_data_in_trial), ]$group_id)
    expect_equal(new_data$stratum_id, new_data_in_trial[nrow(new_data_in_trial), ]$stratum_id)

    # Same for replacing the data
    new_data_in_trial_char_factors <- new_data_in_trial
    new_data_in_trial_char_factors$group_id <- as.character(new_data_in_trial_char_factors$group_id)
    new_data_in_trial_char_factors$stratum_id <- as.character(new_data_in_trial_char_factors$stratum_id)

    trial_with_new_data_replaced <- update(trial_with_new_data, data=new_data_in_trial_char_factors)
    new_data_in_trial_replaced <- summary(trial_with_new_data_replaced, summarize="data")
    expect_equal(new_data_in_trial, new_data_in_trial_replaced)
  })
}

check_update_group_stratum_id_consistency <- function(example) {
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info, simplified_prior = TRUE)

    new_data <- filter(summary(trial, "dose_info"), dose_id == 1)
    new_data$num_patients <- 6
    new_data$num_toxicities <- 1

    # Providing nonsense strings for group_id and stratum_id should error
    new_data_nonsense_group <- new_data
    new_data_nonsense_group$group_id <- "this group is nonsense"

    trial_copy <- trial
    expect_error(update(trial_copy, add_data=new_data_nonsense_group))
    trial_copy <- trial
    expect_error(update(trial_copy, data=new_data_nonsense_group))

    new_data_nonsense_stratum <- new_data
    new_data_nonsense_stratum$stratum_id <- "this stratum is nonsense"
    trial_copy <- trial
    expect_error(update(trial_copy, add_data=new_data_nonsense_stratum))
    trial_copy <- trial
    expect_error(update(trial_copy, data=new_data_nonsense_stratum))

    # Providing nonsense factors for group_id and stratum_id should error
    new_data_nonsense_group <- new_data
    new_data_nonsense_group$group_id <- factor(
      as.character(new_data$group_id), levels=c(levels(new_data$group_id), "nonsense_level")
    )
    trial_copy <- trial
    expect_error(update(trial_copy, add_data=new_data_nonsense_group))
    trial_copy <- trial
    expect_error(update(trial_copy, data=new_data_nonsense_group))

    new_data_nonsense_stratum <- new_data
    new_data_nonsense_stratum$stratum_id <- factor(
      as.character(new_data$stratum_id), levels=c(levels(new_data$stratum_id), "nonsense_level")
    )
    trial_copy <- trial
    expect_error(update(trial_copy, add_data=new_data_nonsense_stratum))
    trial_copy <- trial
    expect_error(update(trial_copy, data=new_data_nonsense_stratum))
  })
}

test_that("update() resolves string factors for group_id, stratum_id correctly",
          check_update_group_stratum_id_char(examples$single_agent))

test_that("update() errors on inconsistent factors for group_id, stratum_id",
          check_update_group_stratum_id_consistency(examples$single_agent))



# Test dose_id and dose_info for consistency, if provided ---------------------------------

check_update_resolves_doses_to_dose_id <- function(example) {
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info, simplified_prior = TRUE)

    new_data <- filter(summary(trial, "dose_info"), dose_id == 3)
    new_data$num_patients <- 6
    new_data$num_toxicities <- 1

    new_data_without_dose <- new_data
    new_data_without_dose$dose_id <- NULL

    ## Providing existing doses should resolve the dose_id
    trial_with_new_data <- update(trial, add_data=new_data_without_dose)

    new_data_in_trial <- summary(trial_with_new_data, summarize="data")

    expect_equal(select(new_data_in_trial[nrow(new_data_in_trial),], colnames(new_data)), new_data)

    ## Providing a non-predefined dose should create a warning
    new_data_not_predefined <- new_data_without_dose
    new_data_not_predefined$drug1 <- 1234
    expect_warning(update(trial, add_data=new_data_not_predefined))
  })
}

check_update_error_on_inconsistent_dose_id <- function(example) {
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info, simplified_prior = TRUE)

    new_data <- filter(summary(trial, "dose_info"), dose_id == 3)
    new_data$num_patients <- 6
    new_data$num_toxicities <- 1

    new_data_with_inconsistent_dose <- new_data
    new_data_with_inconsistent_dose$dose_id <- 1

    expect_error(update(trial, add_data=new_data_with_inconsistent_dose))

    new_data_with_NA_dose <- new_data
    new_data_with_NA_dose$drug1 <- 1234
    new_data_with_NA_dose$dose_id <- NA
    expect_warning(update(trial, add_data=new_data_with_NA_dose))
  })
}

test_that("update resolves dose combos to dose_id",
          check_update_resolves_doses_to_dose_id(examples$single_agent))

test_that("update prevents inconsistent dose combo to dose_id mapping and warns on NA dose_id",
          check_update_error_on_inconsistent_dose_id(examples$single_agent))


# Test the three ways of specifying a prior -----------------------------------------------------
check_simplified_prior <- function(example, formula_generator) {
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info, simplified_prior = TRUE, formula_generator = formula_generator)
    expect_tibble(summary(trial, "data_prediction"))
  })
}

check_full_prior_with_update <- function(example, formula_generator) {
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info, formula_generator = formula_generator)
    
    dims <- summary(trial, "dimensionality")
    
    prior_EX_mu_mean_comp  <- cbind(logit(0.1), rep(0, dims[["num_components"]]))
    prior_EX_mu_sd_comp    <- matrix(c(2, 1), nrow = dims[["num_components"]], ncol = 2, byrow=TRUE)
    
    ## Prior mean and sd on tau_{alpha_{s,i}}, tau{beta_{s,i}}
    prior_EX_tau_mean_comp <- matrix(log(c(0.5, 0.25)), nrow = dims[["num_components"]], ncol = 2, byrow=TRUE)
    prior_EX_tau_sd_comp <- matrix(log(4)/1.96, nrow = dims[["num_components"]], ncol = 2)
    
    # Prior mean and sd on mu_{eta}
    prior_EX_mu_mean_inter  <- rep(0,   dims[["num_interaction_terms"]])
    prior_EX_mu_sd_inter    <- rep(0.5, dims[["num_interaction_terms"]])
    prior_EX_tau_mean_inter <- matrix(log(0.5)   , nrow = dims[["num_strata"]], ncol = dims[["num_interaction_terms"]])
    prior_EX_tau_sd_inter   <- matrix(log(4)/1.96, nrow = dims[["num_strata"]], ncol = dims[["num_interaction_terms"]])
    
    prior_is_EXNEX_comp <- rep(FALSE, dims[["num_components"]])
    prior_EX_prob_comp <- matrix(1, nrow = dims[["num_groups"]], ncol = dims[["num_components"]])
    prior_is_EXNEX_inter <- rep(FALSE, dims[["num_interaction_terms"]])
    prior_EX_prob_inter <- matrix(1, nrow = dims[["num_groups"]], ncol = dims[["num_interaction_terms"]])
    
    prior_tau_dist <- 1
    
    trial <- update(trial, 
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
      
      prior_tau_dist = prior_tau_dist
    )
    
    expect_tibble(summary(trial, "data_prediction"))
  })
}

check_full_prior_direct <- function(formula_generator) {
  with(examples$single_agent, {
    dims <- list(
      num_components = 1,
      num_interaction_terms = 0,
      num_groups = 4,
      num_strata = 1)
    
    prior_EX_mu_mean_comp  <- cbind(logit(0.1), rep(0, dims[["num_components"]]))
    prior_EX_mu_sd_comp    <- matrix(c(2, 1), nrow = dims[["num_components"]], ncol = 2, byrow=TRUE)
    
    ## Prior mean and sd on tau_{alpha_{s,i}}, tau{beta_{s,i}}
    prior_EX_tau_mean_comp <- matrix(log(c(0.5, 0.25)), nrow = dims[["num_components"]], ncol = 2, byrow=TRUE)
    prior_EX_tau_sd_comp <- matrix(log(4)/1.96, nrow = dims[["num_components"]], ncol = 2)
    
    # Prior mean and sd on mu_{eta}
    prior_EX_mu_mean_inter  <- rep(0,   dims[["num_interaction_terms"]])
    prior_EX_mu_sd_inter    <- rep(0.5, dims[["num_interaction_terms"]])
    prior_EX_tau_mean_inter <- matrix(log(0.5)   , nrow = dims[["num_strata"]], ncol = dims[["num_interaction_terms"]])
    prior_EX_tau_sd_inter   <- matrix(log(4)/1.96, nrow = dims[["num_strata"]], ncol = dims[["num_interaction_terms"]])
    
    prior_is_EXNEX_comp <- rep(FALSE, dims[["num_components"]])
    prior_EX_prob_comp <- matrix(1, nrow = dims[["num_groups"]], ncol = dims[["num_components"]])
    prior_is_EXNEX_inter <- rep(FALSE, dims[["num_interaction_terms"]])
    prior_EX_prob_inter <- matrix(1, nrow = dims[["num_groups"]], ncol = dims[["num_interaction_terms"]])
    
    prior_tau_dist <- 1    

    trial <- blrm_trial(histdata, dose_info, drug_info, 
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
      formula_generator = formula_generator
    )
    
    expect_tibble(summary(trial, "data_prediction"))
  })
}


test_that("simplified prior specification enables prediction",
          check_simplified_prior(examples$single_agent, blrm_formula_linear))


test_that("Full prior specification with update() enables prediction",
          check_full_prior_with_update(examples$single_agent, blrm_formula_linear))

test_that("Full prior specification with blrm_trial() enables prediction",
          check_full_prior_direct(blrm_formula_linear))


# Test with saturating model ------------------------------------------------------------------

test_that("simplified prior specification enables prediction (single-agent, saturating interaction model)",
          check_simplified_prior(examples$single_agent, blrm_formula_saturating))

test_that("simplified prior specification enables prediction (combo2, saturating interaction model)",
          check_simplified_prior(examples$combo2, blrm_formula_saturating))

test_that("simplified prior specification enables prediction (combo3, saturating interaction model)",
          check_simplified_prior(examples$combo3, blrm_formula_saturating))


test_that("Full prior specification with update() enables prediction (saturating interaction model)",
          check_full_prior_with_update(examples$single_agent, blrm_formula_saturating))

test_that("Full prior specification with blrm_trial() enables prediction (saturating interaction model)",
          check_full_prior_direct(blrm_formula_saturating))
test_that("update function for blrm_trial adds data with add_data= and replaces data with data= as expected (saturating interaction model)",
          check_blrm_trial_update(examples$single_agent, blrm_formula_saturating))


# Test sorting in .blrm_trial_merge_data ------------------------------------------------------
set_prior <- function(trial, mu_sd_inter = 0.5) 
{
  dims <- summary(trial, "dimensionality")

  update(
    trial,
    # Prior mean and sd on mu_{alpha_i}, mu_{beta_i}
    prior_EX_mu_mean_comp  = matrix(c(logit(0.10), 0), nrow = dims$num_components, ncol = 2, TRUE),
    prior_EX_mu_sd_comp    = matrix(c(2, 1), nrow = dims$num_components, ncol = 2, TRUE),
    
    # Prior mean and sd on tau_{alpha_{s,i}}, tau{beta_{s,i}} 
    prior_EX_tau_mean_comp = do.call("abind", c(
      replicate(dims$num_strata, matrix(c(0, 0), nrow = dims$num_components, ncol = 2, TRUE), simplify=FALSE),
      list(along = 0))
    ),
    prior_EX_tau_sd_comp = do.call("abind", c(
      replicate(dims$num_strata, matrix(0, nrow = dims$num_components, ncol = 2, TRUE), simplify=FALSE),
      list(along = 0))
    ),
    
    # Prior mean and sd on mu_{eta}
    prior_EX_mu_mean_inter  = rep(0,   dims$num_interaction_terms),
    prior_EX_mu_sd_inter    = rep(mu_sd_inter, dims$num_interaction_terms), 
    prior_EX_tau_mean_inter = matrix(0, nrow = dims$num_strata, ncol = dims$num_interaction_terms),
    prior_EX_tau_sd_inter   = matrix(0, nrow = dims$num_strata, ncol = dims$num_interaction_terms),
    
    prior_is_EXNEX_comp = rep(FALSE, dims$num_components),
    prior_EX_prob_comp = matrix(1, nrow = dims$num_groups, ncol = dims$num_components),
    
    prior_is_EXNEX_inter = rep(FALSE, dims$num_interaction_terms),
    prior_EX_prob_inter = matrix(1, nrow = dims$num_groups, ncol = dims$num_interaction_terms),
    
    prior_tau_dist = 0
  )
}

set_prior_EXNEX <- function(trial, mu_sd_inter = 0.5) 
{
  dims <- summary(trial, "dimensionality")
  
  update(
    trial,
    # Prior mean and sd on mu_{alpha_i}, mu_{beta_i}
    prior_EX_mu_mean_comp  = matrix(c(logit(0.10), 0), nrow = dims$num_components, ncol = 2, TRUE),
    prior_EX_mu_sd_comp    = matrix(c(2, 1), nrow = dims$num_components, ncol = 2, TRUE),
    
    # Prior mean and sd on tau_{alpha_{s,i}}, tau{beta_{s,i}} 
    prior_EX_tau_mean_comp = do.call("abind", c(
      replicate(dims$num_strata, matrix(c(0, 0), nrow = dims$num_components, ncol = 2, TRUE), simplify=FALSE),
      list(along = 0))
    ),
    prior_EX_tau_sd_comp = do.call("abind", c(
      replicate(dims$num_strata, matrix(0, nrow = dims$num_components, ncol = 2, TRUE), simplify=FALSE),
      list(along = 0))
    ),
    
    # Prior mean and sd on mu_{eta}
    prior_EX_mu_mean_inter  = rep(0,   dims$num_interaction_terms),
    prior_EX_mu_sd_inter    = rep(mu_sd_inter, dims$num_interaction_terms), 
    prior_EX_tau_mean_inter = matrix(0, nrow = dims$num_strata, ncol = dims$num_interaction_terms),
    prior_EX_tau_sd_inter   = matrix(0, nrow = dims$num_strata, ncol = dims$num_interaction_terms),
    
    prior_is_EXNEX_comp = rep(TRUE, dims$num_components),
    prior_EX_prob_comp = matrix(0.8, nrow = dims$num_groups, ncol = dims$num_components),
    
    prior_is_EXNEX_inter = rep(TRUE, dims$num_interaction_terms),
    prior_EX_prob_inter = matrix(0.8, nrow = dims$num_groups, ncol = dims$num_interaction_terms),
    
    prior_tau_dist = 0
  )
}

check_data_sorting <- function(example) {
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info)
    trial <- set_prior(trial)
    
    new_data <- filter(summary(trial, "dose_info"))
    new_data <- mutate(new_data, num_patients = 2*dose_id)
    new_data$num_toxicities <- 1
    new_data <- arrange(new_data, dose_id)
    
    rev_new_data <- arrange(new_data, desc(dose_id))
    
    trial_2 <- update(trial, add_data = new_data)
    trial_2_rev <- update(trial, add_data = rev_new_data)
    
    expect_equal(trial_2$blrmfit$data, trial_2_rev$blrmfit$data)
  })
}

test_that(".blrm_trial_merge_data sorts data",
          check_data_sorting(examples$single_agent))
test_that(".blrm_trial_merge_data sorts data",
          check_data_sorting(examples$single_drug_with_strata))


# Test for corner cases -----------------------------------------------------------------------
check_trial_with_EXNEX_prior <- function(example) {
  with(example, {
    trial <- blrm_trial(histdata, dose_info, drug_info)
    print(trial)
    trial <- set_prior_EXNEX(trial)
  })
}

test_that("Try EXNEX with one group, but multiple components / interactions",
          check_trial_with_EXNEX_prior(examples$multi_drug_single_group))
