
#' @title Derive synthetic micro datasets for a given geography.
#' @description Derive synthetic micro datasets for each sub-geography of a given set of geographic 
#' macro data constraining tabulations. See Details... By default, micro dataset generation is run 
#' in parallel with load balancing. Macro data is assumed to have been pulled from the US Census API 
#' via the \code{acs} package.
#' 
#' @section Details:
#' In the absence of true micro level datasets for a given geographic area, synthetic datasets
#' can be used. This function uses conditional and marginal probability distributions (at the 
#' aggregate level) to generate synthetic micro population datasets, which are built one constraint 
#' at a time. Taking as input the macro level data (class \code{"macroACS"}), this function builds 
#' synthetic micro datasets for each lower level geographical area within the area of study.
#' 
#' In simplest terms, the goal is to generate a joint probability distribution for an attribute 
#' vector; and, to create synthetic individuals from this distribution. However, note that information
#' for the full joint distribution is typically not available, so we construct it as a product of 
#' conditional and marginal probabilities. This is done one attribute at a time; where it is assumed 
#' that there is some sort of continuum of attribute dependence. That is, some attributes are more 
#' important (eg. gender, age) in 'determining' others (eg. educational attainment, marital status, 
#' etc). These more important attributes need to be assigned first, whereas less important attributes 
#' may be assigned later. Most of these distinctions are largely intuitive, but care must be taken 
#' in choosing the order of constructed attributes. 
#' 
#' This function provides a synthetic population with the following characteristics as well as each 
#' synthetic individual's probability of inclusion. The included characteristics are: age, gender, 
#' marital status, educational attainment, employment status, nativity, poverty status, geographic 
#' mobility in the prior year, individual income, and race. **Note** that these are INDIVIDUAL attributes; 
#' they are not at the HOUSEHOLD level. Additional attributes which interest the user may be added 
#' in a similar manner via \code{\link{synthetic_new_attribute}}.
#' 
#' @param macro_data A macro dataset list: the result of \code{\link{pull_synth_data}}.
#' @param parallel Logical, defaults to \code{TRUE}. Do you wish to run the operation in
#' parallel?
#' @param leave_cores How many cores do you wish to leave open to other processing?
#' @return A \code{list} of the input macro datasets produced by
#' \code{\link{pull_synth_data}} and a \code{list} of synthetic micro datasets for each geographical 
#' subset within the specified macro geography.
#' @seealso \code{\link{pull_synth_data}}, \code{\link[acs]{acs.fetch}}, \code{\link[acs]{geo.make}}
#' @references Birkin, Mark, and M. Clarke. "SYNTHESIS-a synthetic spatial information system 
#' for urban and regional analysis: methods and examples." Environment and planning A 
#' 20.12 (1988): 1645-1671.
#' @export
#' 
#' @examples \dontrun{
#' # make geography
#' la_geo <- acs::geo.make(state= "CA", county= "Los Angeles", tract= "*")
#' # pull data elements for creating synthetic data
#' la_dat <- pull_synth_data(2014, 5, la_geo)
#' # derive synthetic data
#' la_synthetic <- derive_synth_datasets(la_dat, leave_cores= 0)
#' }

derive_synth_datasets <- function(macro_data,
                                  parallel= TRUE, leave_cores= 2) {
  
  # 01. Do some preliminaries
  #--------------------------------------------------------
  if (!is.macroACS(macro_data)) stop("Must input an appropriate macro_data object.")
  if (parallel == TRUE) {
    if (leave_cores < 0 | leave_cores > parallel::detectCores() | leave_cores %% 1 != 0) {
      stop("leave_cores must be an integer between 0 (not recommended) 
           and ", parallel::detectCores())
    }
  }
  # 02. Pull needed macro data and separate sub geographies
  #--------------------------------------------------------
  macro_data2 <- disaggregate_md(macro_data$estimates)
  n <- nrow(macro_data$estimates[[1]])
  
  # 03. Create synthetic data for each geography
  #--------------------------------------------------------
  if (parallel) {
    # create cluster
    nnodes <- min(n, parallel::detectCores() - leave_cores)
    if (grepl("Windows", utils::sessionInfo()$running)) {cl <- parallel::makeCluster(nnodes, type= "PSOCK")}
    else {cl <- parallel::makeCluster(nnodes, type= "FORK")}
    
    # parallel load balanced option:
    synth_data <- parallel::parLapplyLB(cl, macro_data2, synthesize)
    parallel::stopCluster(cl)
  } else {
    synth_data <- lapply(macro_data2, synthesize)
  }
  
  # 04. Return
  #--------------------------------------------------------
  syn <- vector(mode= "list", length= n)
  names(syn) <- rownames(macro_data$estimates[[1]])
  
  for (i in 1:length(synth_data)) {
    syn[[i]][[1]] <- macro_data2[[i]]
    syn[[i]][[2]] <- synth_data[[i]]
  }
  syn <- lapply(syn, function(l) {
    names(l) <- c("macro_constraints", "synthetic_micro")
    class(l) <- c("list", "macro_micro")
    return(l)
  })
  class(syn) <- c("synthACS", "list")
  return(synthetic_data= syn)
}




#--------------------------------------------------------
# helper function for synthesizing a micro dataset from a list of macro data vectors.
synthesize <- function(l) {
  temp <- synth_data_ag(    age_gender_vec= l$age_by_sex)
  temp <- synth_data_mar(   temp, mar_status_vec= l$marital_status)
  temp <- synth_data_edu(   temp, edu_vec= l$edu)
  temp <- synth_data_emp(   temp, emp_status_vec= l$emp_status)
  temp <- synth_data_nativ( temp, nativity_vec= l$nativity)
  temp <- synth_data_pov(   temp, pov_ge_vec= l$pov_status1, total_pop= l$age_by_sex[1])
  temp <- synth_data_geomob(temp, geomob_vec= l$geo_mob_edu)
  temp <- synth_data_inc(   temp, inc_nat_vec= l$by_inc_12mo)
  temp <- synth_data_race(  temp, race_vec= l$pop_by_race)
  return(temp)
}


# wrapper function to disaggregate_mdCPP. 
# for breaking apart macro data lists per sub geography
disaggregate_md <- function(macro_data) {
  # get list and column names
  c_names <- sapply(macro_data, colnames)
  l_names <- names(macro_data)
  # disagregate in C++
  macro_data <- lapply(macro_data, as.matrix) # c++ file requres matrices, not DFs
  md2 <- disaggregate_mdCPP(macro_data)
  # re-apply names and exit
  md2 <- lapply(md2, function(l, m_names) {names(l) <- m_names; return(l)}, m_names= l_names)
  return(lapply(md2, function(l) {mapply(function(l, nms) {
    names(l) <- nms
    return(l)
  }, l=l, nms= c_names)}))
}


