library(methods)

#new wrapper to go at end

#' Internal class
#'
#' in internal class
#'
#' @keywords internal
m_SignificantFeaturesSearch <- setRefClass("m_SignificantFeaturesSearch",
    fields = c('.inst', '.alpha', '.lmax', '.result_available', '.file_loaded'),
    methods = list(

        initialize = function(set_defaults=TRUE)
        {
            .self$.inst <-.self$.create_instance()
            .self$.alpha <- NULL
            .self$.lmax <- NULL
            .self$.result_available <- FALSE
            .self$.file_loaded <- FALSE
            if (set_defaults) {
                .self$set_alpha(0.05)
                .self$set_lmax(0)
            }
        },

        .create_instance = function()
        {
            NULL
        },

        .delete_instance = function(inst)
        {
            NULL
        },

        finalize = function()
        {
            .self$.delete_instance(.self$.inst)
        },

        .check_if_read_is_allowed = function()
        {
        },

        .mark_read_done = function() {
            .self$.result_available <- FALSE
            .self$.file_loaded <- TRUE
        },

        read_eth_files = function(data_path, labels_path, cov_path=NULL, encoding="dominant")
        {
            .self$.check_if_read_is_allowed()
            .self$.do_read_eth_files(data_path, labels_path, cov_path, encoding)
            .self$.mark_read_done()
        },

        # TODO: split for _int and _iset due to multi-inheritance
        .do_read_eth_files = function(data_path, labels_path, cov_path=NULL, encoding="dominant") {
            lib_read_eth_files(.self$.inst, data_path, labels_path, encoding)
        },

        read_plink_files = function(base_path, cov_path=NULL, encoding="dominant")
        {
            .self$.check_if_read_is_allowed()
            .self$.do_read_plink_files(base_path, cov_path, encoding)
            .self$.mark_read_done()
        },

        # TODO: split for _int and _iset due to multi-inheritance
        .do_read_plink_files = function(base_path, cov_path=NULL, encoding="dominant") {
            lib_read_plink_files(.self$.inst, base_path, encoding)
        },

        .check_if_write_is_allowed = function() {
            .self$.check_if_files_are_loaded()
        },

        write_eth_files = function(x_file, y_file, ...)
        {
            .self$.check_if_write_is_allowed()
            .self$.write_eth_files(x_file, y_file, ...)
        },

        .check_if_alpha_value_is_allowed = function(x)
        {
            if (!(x>=0 && x<=1)) {
                stop("you need to set alpha to a value between 0 and 1")
            }
        },

        set_alpha = function(alpha)
        {
            .self$.check_if_alpha_value_is_allowed(alpha)
            .self$.alpha <- alpha
            .self$.result_available <- FALSE
        },

        get_alpha = function() {
            return(.self$.alpha)
        },

        .check_if_lmax_value_is_allowed = function(x)
        {
            if (!(x%%1==0 && x>=0)) {
                stop("you need to set lmax to a non-negative integer value")
            }
        },

        set_lmax = function(lmax)
        {
            .self$.check_if_lmax_value_is_allowed(lmax)
            .self$.lmax <- lmax
            .self$.result_available <- FALSE
        },

        get_lmax = function() {
            return(.self$.lmax)
        },

        execute = function()
        {
            .self$.check_if_execute_is_allowed()
            .self$.result_available <- FALSE
            .self$.execute()
            .self$.result_available <- TRUE
        },

        .check_if_result_available = function()
        {
            if (!.self$.result_available) {
                stop("you need to call the execute method first")
            }
        },

        get_result = function()
        {
            .self$.check_if_result_available()
            result <- .self$.get_result()
            result["target.fwer"] <-.self$.alpha
            return(result)
        },

        write_summary = function(file_path)
        {
            .self$.check_if_result_available()
            .self$.write_summary(file_path)
        },

        .check_if_files_are_loaded = function()
        {
            if (!.self$.file_loaded) {
                stop("you need to call the read_eth_files / read_plink_files method first")
            }
        },

        .check_if_alpha_is_set = function()
        {
            if (is.null(.self$.alpha)) {
                stop("you need to call the set_alpha method first")
            }
        },

        .check_if_lmax_is_set = function()
        {
            if (is.null(.self$.lmax)) {
                stop("you need to call the set_lmax method first")
            }
        },

        .check_if_execute_is_allowed = function()
        {
            .self$.check_if_alpha_is_set()
            .self$.check_if_lmax_is_set()
            .check_if_files_are_loaded()
        },

        #NOTE NEW: print method
        show = function(){
            mytype <- "unknown"
            myclasstype = toString(class(.self)[[1]])
            if (myclasstype=="SignificantIntervalSearchExact")
                mytype = "FAIS"
            if (myclasstype=="SignificantIntervalSearchFastCmh")
                mytype = "FastCMH"
            if (myclasstype=="SignificantItemsetSearchFacs")
                mytype = "FACS"

            #extract alpha and lmxax
            myalpha <- toString(.self$.alpha)
            mylmax <- toString(.self$.lmax)

            #output message to be returned
            message1 <- paste0(mytype, " object with:")
            message2 <- paste0(" * alpha = ", myalpha)
            message3 <- paste0(" * lmax = ", mylmax)
            cat(message1, "\n")
            cat(message2, "\n")
            cat(message3, "\n")
        },

        write_profile = function(file_path)
        {
            .self$.check_if_result_available()
            lib_profiler_write_to_file(.self$.inst, file_path)
        }

    )
)

#' Internal class
#'
#' An internal class
#'
#' @keywords internal
m_SignificantIntervalSearch <- setRefClass("m_SignificantIntervalSearch",
    contains = c("m_SignificantFeaturesSearch"),
    methods = list(

        # .do_read_eth_files = function(data_path, labels_path, ...)
        # {
        #     lib_read_eth_files_int(.self$.inst, data_path, labels_path)
        # },

        # .do_read_plink_files = function(base_path, ...)
        # {
        #     lib_read_plink_files_int(.self$.inst, base_path)
        # },

        .execute = function() {
            lib_execute_int(.self$.inst, .self$.alpha, .self$.lmax)
        },

        .write_eth_files = function(x_file, y_file, ...)
        {
            lib_write_eth_files_int(.self$.inst, x_file, y_file)
        },

        write_filtered_intervals = function(file_path)
        {
            .self$.check_if_result_available()
            lib_filter_intervals_write_to_file(.self$.inst, file_path)
        },

        write_pvals_testable_intervals = function(file_path)
        {
            .self$.check_if_result_available()
            lib_pvals_testable_ints_write_to_file(.self$.inst, file_path)
        },

        write_pvals_significant_intervals = function(file_path)
        {
            .self$.check_if_result_available()
            lib_pvals_significant_ints_write_to_file(.self$.inst, file_path)
        },

        get_significant_intervals = function()
        {
            .self$.check_if_result_available()
            return(lib_get_significant_intervals(.self$.inst))
        },

        get_filtered_intervals = function()
        {
            .self$.check_if_result_available()
            return(lib_get_filtered_intervals(.self$.inst))
        },

        .get_result = function()
        {
            return(lib_get_result_fais(.self$.inst))
        }
    )

)

#' Internal class
#'
#' An internal class.
#'
#' @keywords internal
m_SignificantIntervalSearchFais <- setRefClass("m_SignificantIntervalSearchFais",
    contains = c("m_SignificantIntervalSearch"),
    methods = list(

        .write_summary = function(file_path)
        {
            lib_summary_write_to_file_fais(.self$.inst, file_path)
        }

    )

)

#' Internal class for search for significant regions
#'
#' Please use the \code{CASMAP} constructor.
#'
#' @keywords internal
SignificantIntervalSearchExact <- setRefClass("SignificantIntervalSearchExact",
    contains = "m_SignificantIntervalSearchFais",
    methods = list(
        .create_instance = lib_new_search_e,
        .delete_instance = lib_delete_search_e
    )
)

#' Approximate fast significant interval search
#'
#' Class for approximate significant intervals search with Tarone correction for
#' bounding intermediate FWERs.
#'
#' @keywords internal
SignificantIntervalSearchChi <- setRefClass("SignificantIntervalSearchChi",
    contains = c("m_SignificantIntervalSearchFais"),
    methods = list(
        .create_instance = lib_new_search_chi,
        .delete_instance = lib_delete_search_chi
    )
)


#' Internal class
#'
#' @keywords internal
m_SignificantFeaturesSearchWithCovariates <- setRefClass("m_SignificantFeaturesSearchWithCovariates",
    contains = c("m_SignificantFeaturesSearch"),
    fields = c(".cov_loaded"),
    methods = list(

        initialize = function(...) {
            callSuper(...)
            .self$.cov_loaded <- FALSE
        },

        .do_read_eth_files = function(data_path, labels_path, cov_path=NULL, encoding="dominant") {
            if (is.null(cov_path)) {
                callSuper(data_path, labels_path, encoding)
                .self$.cov_loaded <- FALSE
            } else {
                .self$.read_eth_files_with_cov(data_path, labels_path, cov_path, encoding)
                .self$.cov_loaded <- TRUE
            }
        },

        .do_read_plink_files = function(base_path, cov_path=NULL, encoding="dominant") {
            if (is.null(cov_path)) {
                callSuper(base_path, encoding)
                .self$.cov_loaded <- FALSE
            } else {
                .self$.read_plink_files_with_cov(base_path, cov_path, encoding)
                .self$.cov_loaded <- TRUE
            }
        },

        update_covariates_file = function(cov_path)
        {
            .self$.check_if_files_are_loaded()
            .self$.update_covariates_file(cov_path)
            .self$.result_available <- FALSE
            .self$.cov_loaded <- TRUE
        },

        .write_eth_files = function(x_file, y_file, cov_path=NULL, ...)
        {
            if (is.null(cov_path)) {
                lib_write_eth_files(.self$.inst, x_file, y_file)
            } else {
                .self$.check_if_covariates_are_loaded()
                .self$.write_eth_files_with_cov(x_file, y_file, cov_path)
            }
        },

        .check_if_covariates_are_loaded = function() {
            if (!.self$.cov_loaded) {
                #warning("assuming one covariate for all observations; to change covariates call the update_covariates_file method first")
            }
        },

        .check_if_execute_is_allowed = function()
        {
            callSuper()
            .self$.check_if_covariates_are_loaded()
        },

        .get_result = function()
        {
            return(lib_get_result_int(.self$.inst))
        }

    )
)


#' Fast significant interval search with categorical covariates
#'
#' Internal class, please use \code{CASMAP} constructor. 
#'
#' @keywords internal
SignificantIntervalSearchFastCmh <- setRefClass("SignificantIntervalSearchFastCmh",
    # Beware: order matters for calling overloaded covariates methods
    contains = c("m_SignificantFeaturesSearchWithCovariates", "m_SignificantIntervalSearch"),
    methods = list(
        .create_instance = lib_new_search_fastcmh,
        .delete_instance = lib_delete_search_fastcmh,

        .read_eth_files_with_cov = function(x_file, y_file, cov_path, encoding)
        {
            lib_read_eth_files_with_cov_fastcmh(.self$.inst, x_file, y_file, cov_path, encoding)
        },

        .read_plink_files_with_cov = function(base_path, cov_path, encoding)
        {
            lib_read_plink_files_with_cov_fastcmh(.self$.inst, base_path, cov_path, encoding)
        },

        .write_eth_files_with_cov = function(x_file, y_file, cov_path, ...)
        {
            lib_write_eth_files_with_cov_fastcmh(.self$.inst, x_file, y_file, cov_path)
        },

        .update_covariates_file = function(covariates_path)
        {
            lib_read_covariates_file_fastcmh(.self$.inst, covariates_path)
        },

        .write_summary = function(file_path)
        {
            lib_summary_write_to_file_fastcmh(.self$.inst, file_path)
        }

    )
)



#' Internal class
#'
#' @keywords internal
m_SignificantItemsetSearch <- setRefClass("m_SignificantItemsetSearch",
    contains = c("m_SignificantFeaturesSearch"),
    methods = list(

        # .do_read_eth_files = function(data_path, labels_path)
        # {
        #     lib_read_eth_files_iset(.self$.inst, data_path, labels_path)
        # },

        # .do_read_plink_files = function(base_path)
        # {
        #     lib_read_plink_files_iset(.self$.inst, base_path)
        # },

        .execute = function() {
            lib_execute_iset(.self$.inst, .self$.alpha, .self$.lmax)
        },

        .write_eth_files = function(x_file, y_file, ...)
        {
            lib_write_eth_files_iset(.self$.inst, x_file, y_file)
        },

        write_pvals_testable_itemsets = function(file_path)
        {
            .self$.check_if_result_available()
            lib_pvals_testable_isets_write_to_file(.self$.inst, file_path)
        },

        write_pvals_significant_itemsets = function(file_path)
        {
            .self$.check_if_result_available()
            lib_pvals_significant_isets_write_to_file(.self$.inst, file_path)
        },

        get_significant_itemsets = function()
        {
            .self$.check_if_result_available()
            return(lib_get_significant_itemsets(.self$.inst))
        },

        .get_result = function()
        {
            return(lib_get_result_iset(.self$.inst))
        }

    )

)


#' Significant itemsets search with categorical covariates
#'
#' Internal class, please use \code{CASMAP} constructor. 
#'
#' @keywords internal
SignificantItemsetSearchFacs <- setRefClass("SignificantItemsetSearchFacs",
    # Beware: order matters for calling overloaded covariates methods
    contains = c("m_SignificantFeaturesSearchWithCovariates", "m_SignificantItemsetSearch"),
    methods = list(
        .create_instance = lib_new_search_facs,
        .delete_instance = lib_delete_search_facs,

        .read_eth_files_with_cov = function(x_file, y_file, cov_path, encoding)
        {
            lib_read_eth_files_with_cov_facs(.self$.inst, x_file, y_file, cov_path, encoding)
        },

        .read_plink_files_with_cov = function(base_path, cov_path, encoding)
        {
            lib_read_plink_files_with_cov_facs(.self$.inst, base_path, cov_path, encoding)
        },

        .write_eth_files_with_cov = function(x_file, y_file, cov_path)
        {
            lib_write_eth_files_with_cov_facs(.self$.inst, x_file, y_file, cov_path)
        },

        .update_covariates_file = function(cov_path)
        {
            lib_read_covariates_file_facs(.self$.inst, cov_path)
        },

        .get_result = function()
        {
            return(lib_get_result_facs(.self$.inst))
        },

        .write_summary = function(file_path)
        {
            lib_summary_write_to_file_facs(.self$.inst, file_path)
        }

    )
)


#' A method to check value is numeric and in open interval
#' 
#' Checks if a value is numeric and strictly between two other values.
#' 
#' @param x Value to be checked. Needs to be numeric.
#' 
#' @param lower Lower bound. Default value is \code{0}.
#' 
#' @param upper Upper bound. Default value is \code{1}.
#' 
#' @return If numeric, and  strictly greater than \code{lower} and 
#'         strictly smaller than \code{upper}, then return \code{TRUE}. 
#'         Else return \code{FALSE}.
#' @keywords internal
isInOpenInterval <- function(x, lower=0, upper=1){
    inInterval <- TRUE
    if (is.finite(x)){
        if ( (x <= lower) || (x >= upper) )
            inInterval <- FALSE
    } else {
        inInterval <- FALSE
    }
    return (inInterval)
}


#' Check if a variable is boolean or not
#'
#' Checks if a variable is boolean, if not throws an error, otherwise
#' returns boolean.
#'
#' @param var The variable to be checked (if boolean).
#'
#' @param name The name of the variable to appear in any error message.
#'
#' @return If not boolean (or \code{NA}), throws error.
#'         If \code{NA}, return \code{FALSE}. Otherwise return 
#'         boolean value of \code{var}.
#' @keywords internal
checkIsBoolean <- function(var, name){
    if (is.logical(var)) {
        # if NA, return FALSE
        if (is.na(var))
            return(FALSE)

        # otherwise, just return its value (it must be TRUE/FALSE
        # from above check
        return(var)
    } else {
        message <- paste0("Error: ", name, " is not a boolean.")
        stop(message)
    }
}



#' Constructor for CASMAP class object.
#'
#' @field mode Either \code{'regionGWAS'} or \code{'higherOrderEpistasis'}.
#'
#' @field alpha A numeric value setting the Family-wise Error Rate (FWER).
#'              Must be strictly between \code{0} and \code{1}. Default
#'              value is \code{0.05}.
#'
#' @field max_comb_size A numeric specifying the maximum length of 
#'                      combinations. For example, if set to \code{4}, 
#'                      then only combinations of size between \code{1}
#'                      and \code{4} (inclusive) will be considered.
#'                      To consider combinations of arbitrary (maximal)
#'                      length, use value \code{0}, which is the default 
#'                      value.
#'
#' @details 
#' Constructor for CASMAP class object, which needs the \code{mode}
#' parameter to be set by the user. Please see the examples.
#'
#' 
#' @section Base method, for both modes:
#' \describe{
#'   \item{\code{readFiles}}{Read the data, label and possibly covariates 
#'                          files. Parameters are \code{genotype_file}, 
#'                          for the data, \code{phenotype_file} for the
#'                          labels and (optional) \code{covariates_file}
#'                          for the covariates. The option 
#'                          \code{plink_file_root} is not supported
#'                          in the current version, but will be supported
#'                          in future versions.}
#'
#'   \item{\code{setMode}}{Can set/change the mode, but note that any
#'                        data files will need to read in again using
#'                        the \code{readFiles} command.}
#'
#'   \item{\code{setTargetFWER}}{Can set/change the Family-wise 
#'                              Error Rate (FWER). Takes a numeric
#'                              parameter \code{alpha}, strictly between
#'                              \code{0} and \code{1}.}
#'   
#'   \item{\code{execute}}{Once the data files have been read, can execute the
#'                        algorithm. Please note that, depending on the size
#'                        of the data files, this could take a long time.}
#'   
#'   \item{\code{getSummary}}{Returns a data frame with a summary of the 
#'                           results from the execution, but not any
#'                           significant regions/itemsets. See 
#'                           \code{getSignificantRegions}, 
#'                           \code{getSignificantInteractions}, and
#'                           \code{getSignificantClusterRepresentatives}. }
#'   
#'   \item{\code{writeSummary}}{Directly write the information
#'                              from \code{getSummary} to file.}
#'   
#' }
#' 
#' @section \code{regionGWAS} Methods:
#' \describe{
#'   \item{\code{getSignificantRegions}}{Returns a data frame with the 
#'                           the significant regions. Only valid when
#'                           \code{mode='regionGWAS'}.}
#'   
#'   \item{\code{getSignificantClusterRepresentatives}}{Returns a data 
#'                           frame with the 
#'                           the representatives of the significant 
#'                           clusters. This will be a subset of the regions
#'                           returned from \code{getSignificantRegions}. 
#'                           Only valid when \code{mode='regionGWAS'}.}
#'   
#'   \item{\code{writeSignificantRegions}}{Writes the data from 
#'                           \code{getSignificantRegions} to file, which
#'                           must be specified in the parameter 
#'                           \code{path}.
#'                           Only valid when \code{mode='regionGWAS'}.}
#'   
#'   \item{\code{writeSignificantClusterRepresentatives}}{Writes the data 
#'                           from 
#'                           \code{getSignificantClusterRepresentatives} to 
#'                           file, which must be specified in the parameter 
#'                           \code{path}.
#'                           Only valid when \code{mode='regionGWAS'}.}
#'   
#'    }
#'   
#' @section \code{higherOrderEpistasis} Methods:
#' \describe{
#'   \item{\code{getSignificantInteractions}}{Returns the frame 
#'                           from \code{getSignificantInteractions} to 
#'                           file, which must be specified in the parameter 
#'                           \code{path}. Only valid 
#'                           when \code{mode='higherOrderEpistasis'}.}
#'   
#'   \item{\code{writeSignificantInteractions}}{Writes a data frame with  
#'                           the significant interactions. Only valid 
#'                           when \code{mode='higherOrderEpistasis'}.}
#'   
#'    }
#'   
#' @section References:
#' A. Terada, M. Okada-Hatakeyama, K. Tsuda and J. Sese
#'  \emph{Statistical significance of combinatorial regulations},
#'  Proceedings of the National Academy of Sciences (2013) 110 
#'  (32): 12996-13001
#'
#' F. Llinares-Lopez, D. G. Grimm, D. Bodenham,
#' U. Gieraths, M. Sugiyama, B. Rowan and K. Borgwardt,
#' \emph{Genome-wide detection of intervals of genetic heterogeneity 
#'       associated with complex traits},
#' ISMB 2015, Bioinformatics (2015) 31 (12): i240-i249
#'
#' L. Papaxanthos, F. Llinares-Lopez, D. Bodenham,
#' K .Borgwardt,
#' \emph{Finding significant combinations of features in the 
#'  presence of categorical covariates}, Advances
#'  in Neural Information Processing Systems 29 (NIPS 2016), 2271-2279.
#'
#' F. Llinares-Lopez, L. Papaxanthos, D. Bodenham,
#' D. Roqueiro and K .Borgwardt,
#' \emph{Genome-wide genetic heterogeneity discovery 
#'  with categorical covariates}.
#' Bioinformatics 2017, 33 (12): 1820-1828.
#'
#' @export
#' @examples
#'
#' ## An example using the "regionGWAS" mode
#' fastcmh <- CASMAP(mode="regionGWAS")      # initialise object
#'
#' datafile <- getExampleDataFilename()      # file name of example data
#' labelsfile <- getExampleLabelsFilename()  # file name of example labels
#' covfile <- getExampleCovariatesFilename() # file name of example covariates 
#'
#' # read the data, labels and covariate files
#' fastcmh$readFiles(genotype_file=getExampleDataFilename(),
#'                   phenotype_file=getExampleLabelsFilename(), 
#'                   covariate_file=getExampleCovariatesFilename() )
#'
#' # execute the algorithm (this may take some time)
#' fastcmh$execute()
#'
#' #get the summary results
#' summary_results <- fastcmh$getSummary()
#'
#' #get the significant regions
#' sig_regions <- fastcmh$getSignificantRegions()
#'
#' #get the clustered representatives for the significant regions
#' sig_cluster_rep <- fastcmh$getSignificantClusterRepresentatives()
#'
#'
#' ## Another example of regionGWAS
#' fais <- CASMAP(mode="regionGWAS")      # initialise object
#'
#' # read the data and labels, but no covariates
#' fastcmh$readFiles(genotype_file=getExampleDataFilename(),
#'                   phenotype_file=getExampleLabelsFilename())
#'
#'
#' ## Another example, doing higher order epistasis search
#' facs <- CASMAP(mode="higherOrderEpistasis")      # initialise object
#'
CASMAP <- setRefClass("CASMAP",
    fields = c('.mode', '.alpha', '.max_comb_size', '.core', '.use_covariates'),
    methods = list(
        initialize = function(mode, alpha=0.05, max_comb_size=0) {

            if (missing(mode)){
                message <- modeErrorMessage()
                stop(message)
            }
            #if gets past those checks, then setMode
            setMode(mode)
            setTargetFWER(alpha)
            .self$setMaxCombinationSize(max_comb_size)

            .self$.core <- NULL
            .self$.use_covariates <- NULL
        },

        getMode = function() {
            return(.self$.mode)
        },

        getTargetFWER = function() {
            return(.self$.alpha)
        },

        getMaxCombinationSize = function() {
            return(.self$.max_comb_size)
        },

        isInitialized = function() {
            return(!is.null(.self$.core))
        },

        .checkInitialized = function() {
            if (is.null(.self$.core)){
                stop("Object not initialized or hyperparameters changed since last execution. Please call method 'readFiles' prior to execute.")
            }
        },

        setMode = function(mode){
            if (!is.character(mode)){
                message <- modeErrorMessage()
                stop(message)
            } 
            if (modeNeedsMoreChars(mode)){
                message <- modeLengthErrorMessage()
                stop(message)
            }

            #it is a character string, so check:
            if (!isRegionGWASString(mode) && !isHigherOrderEpistasisString(mode)){
                message <- modeErrorMessage()
                stop(message)
            }
            #if reaches this stage, it is correctly set:
            if (isRegionGWASString(mode)){
                mode <- getRegionGWASString()
            }
            #if reaches this stage, it is correctly set:
            if (isHigherOrderEpistasisString(mode)){
                mode <- getHigherOrderEpistasisString()
            }
            .self$.mode <- mode
            .self$.core <- NULL
            .self$.use_covariates <- NULL
        },

        setTargetFWER = function(alpha=0.05) {
            #check alpha
            if (!isInOpenInterval(alpha)){
                stop("Target FWER 'alpha' needs to be a value strictly between 0 and 1.")
            }
            .self$.alpha <- alpha
            .self$.core <- NULL
            .self$.use_covariates <- NULL
        },

        setMaxCombinationSize = function(max_comb_size=0) {
            if (is.finite(max_comb_size)){
                max_comb_size <- floor(max_comb_size)
                if (.self$.mode == 'higherOrderEpistasis' & max_comb_size > 0){
                    print("The current implementation of higher-order epistasis analyses does not support a limited maximum number of interacting variants. The analysis will be carried out for an unlimited order.")
                    max_comb_size <- 0
                }
                if (max_comb_size < 0){
                    max_comb_size <- 0
                }
            } else {
                stop("Maximum combination size 'max_comb_size' needs to be either 0 (unlimited) or a positive integer.")
            }
            .self$.max_comb_size <- max_comb_size
            .self$.core <- NULL
            .self$.use_covariates <- NULL
        },

        .createCore = function() {
            if (!is.null(.self$.use_covariates)){
                # Instantiate object of the appropriate depending on options
                if (.self$.mode == 'regionGWAS' & !.self$.use_covariates){
                    .self$.core <- SignificantIntervalSearchChi()
                } else if (.self$.mode == 'regionGWAS' & .self$.use_covariates){
                    .self$.core <- SignificantIntervalSearchFastCmh()
                } else if (.self$.mode == 'higherOrderEpistasis'){
                    .self$.core <- SignificantItemsetSearchFacs()
                }

                # Set parameters of the object
                .self$.core$set_alpha(.self$.alpha)
                .self$.core$set_lmax(.self$.max_comb_size)

            } else {
                .self$.core <- NULL
                .self$.use_covariates <- NULL
            }
        },

        readFiles = function(genotype_file=NULL, phenotype_file=NULL, covariate_file=NULL, plink_file_root=NULL, encoding="dominant") {
            # Check whether user decided to use tab-separated text files (binary_format) or PLINK formatted files (plink_format)
            if (!missing(plink_file_root)){
                stop("plink format not currently supported. Please use binary format and specify 'genotype_file' and 'phenotype_file'.")
            }

            # At least one of the two must be two, otherwise raise an error
            binary_format <- !is.null(genotype_file) & !is.null(phenotype_file)
            if (missing(genotype_file) || missing(phenotype_file)){
                stop("'genotype_file' and 'phenotype_file' must both be specified as arguments.")
            }

            binary_format <- TRUE
            plink_format <- FALSE
            # Check that encoding type is correct
            if (!is.element(encoding, c('dominant', 'recessive'))){
                stop("Currently implemented encodings: 'dominant' and 'recessive' >")
            }

            # If an additional covariates file was specified, set the object into "CMH mode"
            .self$.use_covariates <- !missing(covariate_file)

            # Create appropriate "core" object
            .self$.createCore()

            # Give preference to plink_format over binary_format if, by any reason, a user decides to mess around and
            # specify both
            #plink_forat not currently supported
            #just to make sure it is set to FALSE
            plink_format <- FALSE
            if (plink_format){
                if(.self$.use_covariates){
                    .self$.core$read_plink_files(plink_file_root, 
                                                 covariate_file, 
                                                 encoding)
                } else {
                    .self$.core$read_plink_files(plink_file_root, encoding)
                }
            } else if(binary_format){
                if(.self$.use_covariates){
                    .self$.core$read_eth_files(data_path=genotype_file, 
                                               labels_path=phenotype_file, 
                                               cov_path=covariate_file, 
                                               encoding=encoding)
                } else {
                    .self$.core$read_eth_files(data_path=genotype_file, 
                                               labels_path=phenotype_file, 
                                               encoding=encoding)
                }
            } else{
                stop("'genotype_file' and 'phenotype_file' must both be specified as arguments.")
            }
        },

        execute = function(){
            .self$.checkInitialized()
            .self$.core$execute()
        },

        writeSummary = function(path){
            .self$.checkInitialized()
            .self$.core$write_summary(path)
        },

        writeProfile = function(path){
            .self$.checkInitialized()
            .self$.core$write_profile(path)
        },

        writeSignificantRegions = function(path){
            if (.self$.mode != 'regionGWAS'){
                stop("Method 'writeSignificantRegions' only available for region-based GWAS analyses.")
            }
            .self$.checkInitialized()
            .self$.core$write_pvals_significant_intervals(path)
        },

        writeSignificantClusterRepresentatives = function(path){
            if (.self$.mode != 'regionGWAS'){
                stop("Method 'writeSignificantClusterRepresentatives' only available for region-based GWAS analyses.")
            }
            .self$.checkInitialized()
            .self$.core$write_filtered_intervals(path)
        },

        writeSignificantInteractions = function(path){
            if (.self$.mode != 'higherOrderEpistasis'){
                stop("Method 'writeSignificantInteractions' only available for higher-order epistasis analyses.")
            }
            .self$.checkInitialized()
            .self$.core$write_pvals_significant_itemsets(path)
        },

        getSummary = function(){
            .self$.checkInitialized()
            return(.self$.core$get_result())
        },

        getSignificantRegions = function(){
            if (.self$.mode != 'regionGWAS'){
                stop("Method 'getSignificantRegions' only available for region-based GWAS analyses.")
            }
            .self$.checkInitialized()
            return(.self$.core$get_significant_intervals())
        },

        getSignificantClusterRepresentatives = function(){
            if (.self$.mode != 'regionGWAS'){
                stop("Method 'getSignificantClusterRepresentatives' only available for region-based GWAS analyses.")
            }
            .self$.checkInitialized()
            return(.self$.core$get_filtered_intervals())
        },

        getSignificantInteractions = function(){
            if (.self$.mode != 'higherOrderEpistasis'){
                stop("Method 'getSignificantInteractions' only available for higher-order epistasis analyses.")
            }
            .self$.checkInitialized()
            return(.self$.core$get_significant_itemsets())
        },

        show = function(){
            cat("CASMAP object with:", "\n")
            cat(paste(" * Mode =", .self$.mode), "\n")
            cat(paste(" * Target FWER =", .self$.alpha), "\n")
            cat(paste(" * Maximum combination size =", .self$.max_comb_size), "\n")
            if (!is.null(.self$.core)){
                cat(" * Input files read", "\n")
                cat(paste(" * Covariate =", .self$.use_covariates), "\n")
            } else{
                cat(" * No input files read", "\n")
            }
        }
    )
)

