#' Update data of a BLRM analysis
#'
#' Adds data rows to a \code{\link{blrm_exnex}} or
#' \code{\link{blrm_trial}} analysis object.
#'
#' @param object blrmfit analysis object
#' @param add_data additional data added to analysis data of \code{object}
#' @param ... passed to default \code{update} command
#'
#' The data in \code{add_data} will be combined with data in
#' \code{object} using \code{bind_rows}. The indices for groups and
#' stratums (if defined) are matched between \code{add_data} and the
#' data of the analysis \code{object}.
#'
#' Note that the \code{add_data} argument must be named explicitly as
#' demonstrated in the example.
#'
#' @template start-example
#' @examples
#' example_model("single_agent")
#'
#' library(tibble)
#' new_cohort <- tibble(group_id="trial_A", drug_A=50, num_patients=4, num_toxicities=1)
#'
#' ## this would fail, since add_data argument must be named
#' ## new_blrmfit <- update(blrmfit, new_cohort)
#' new_blrmfit <- update(blrmfit, add_data=new_cohort)
#'
#' @template stop-example
#' @method update blrmfit
#' @export
update.blrmfit <- function(object, ..., add_data) {
    dots  <- list(...)
    if(!missing(add_data)) {
        if("data" %in% names(dots)) {
            current_data <- dots$data
            dots$data <- NULL
        } else {
            current_data <- object$data
        }

        ## ensure that group_id and stratum_id (if present) as used in
        ## object are part of the data in add_data
        labels <- all.vars(delete.response(terms(object$f, rhs=length(object$f)[2])))
        assert_set_equal(intersect(names(add_data), labels), labels, .var.name="grouping and/or stratum columns")
        has_stratum_fct  <- length(labels) == 2

        ## check add_data for consistency with existing data and
        ## matched against definitions of object
        strata_group_fct <- .get_strata_group_fct(object, add_data)

        ## convert the index columns into matched factors
        ## only do this if the object$data is a factor
        if (has_stratum_fct) {
            if(is.factor(current_data[[labels[1]]]))
                add_data[[labels[1]]] <- strata_group_fct$strata_fct
            if(is.factor(current_data[[labels[2]]]))
                add_data[[labels[2]]] <- strata_group_fct$group_fct
        } else {
            if(is.factor(current_data[[labels[1]]]))
                add_data[[labels[1]]] <- strata_group_fct$group_fct
        }
        .combined_data <- bind_rows(current_data, add_data)
        upd_call  <- do.call(update.default, c(list(object=object, data=.combined_data), dots, list(evaluate=FALSE)))
        return(eval(upd_call, c(list(.combined_data=.combined_data), dots), parent.frame()))
    }
    upd_call  <- do.call(update.default, c(list(object=object), dots, list(evaluate=FALSE)))
    return(eval(upd_call, dots, parent.frame()))
}

#' Update data and/or prior of a BLRM trial
#'
#' * Adds data rows to a \code{\link{blrm_trial}} object (add_data argument)
#' * Replaces data of a \code{\link{blrm_trial}} object (data argument)
#' * Sets the prior of a \code{\link{blrm_trial}} object (... argument will be
#'   passed to blrm_exnex)
#'
#' @param object blrm_trial object
#' @param ... passed to default \code{update} command of \code{blrm_exnex}
#'
#' @template start-example
#' @examples
#'
#' # the combo2_trial example demonstrates the use of add_data of
#' # update.blrmfit
#' example_model("combo2_trial")
#'
#' @template stop-example
#' @method update blrm_trial
#' @export
update.blrm_trial <- function(object, ...)
{
    .assert_is_blrm_trial(object)

    if (has_name(object, "update_blrmfit")) {
        # If blrmfit exists, update it
        object <- object$update_blrmfit(object, ...)

        return(object)
    } else {
        # Otherwise, create it
        object <- .blrm_trial_compute_prior(object, ...)

        return(object)
    }
}

