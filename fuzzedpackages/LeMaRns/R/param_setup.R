#' An S4 class representing the inputs of the LeMans model
#'
#' @slot nsc A numeric value representing the number of length classes in the model.
#' @slot nfish A numeric value representing the number of fish species in the model.
#' @slot phi_min A numeric value representing the time step of the model.
#' @slot l_bound A numeric vector of length \code{nsc} representing the lower bounds of the length classes.
#' @slot u_bound A numeric vector of length \code{nsc} representing the upper bounds of the length classes.
#' @slot mid A numeric vector of length \code{nfish} representing the mid-point of the length classes.
#' @slot species_names A numeric or character vector of length \code{nfish} that denotes the names of the species in the model.
#' @slot Linf A numeric vector of length \code{nfish} representing the mid-point of the length classes.
#' @slot W_a A numeric vector representing the parameter \code{a} in the length-weight conversion.
#' @slot W_b A numeric vector representing the parameter \code{b} in the length-weight conversion.
#' @slot k k A numeric vector of length \code{nfish} representing the von Bertalanffy growth parameter \code{(1/yr)} for each species.
#' @slot Lmat A numeric vector of length \code{nsc} representing the length at which 50\% of the individuals are mature.
#' @slot mature A matrix with dimensions \code{nsc} and \code{nfish} with elements in the range 0-1 representing the proportion of mature individuals of each species in each length class.
#' @slot sc_Linf A numeric vector of length \code{nsc} representing the length class at which each species reaches its asymptotic length.
#' @slot wgt A matrix with dimensions \code{nsc} and \code{nfish} representing the weight of each species in each length class.
#' @slot phi A matrix with dimensions \code{nsc} and \code{nfish} representing the proportion of individuals that leave each length class.
#' @slot ration A matrix with dimensions \code{nsc} and \code{nfish} representing the amount of food required for fish of a given species and length class to grow according to the von Bertalanffy growth curve in a time step.
#' @slot other A numeric value representing the amount of other food (g) available from prey that is not explicitly represented in the model.
#' @slot M1 A matrix of dimensions \code{nsc} and \code{nfish} representing the natural mortality of each species for each length class.
#' @slot suit_M2 A list object of length \code{nfish}. Each element in the list is an array of dimensions \code{nsc}, \code{nsc} and \code{nfish} containing a value between 0 and 1 representing prey preference and prey suitability for each species and length class.
#' @slot Qs An array of dimensions \code{nsc}, \code{nfish} and \code{length(gear)} representing the catchability of each species by each of the fishing gears.
#' @slot stored_rec_funs A list object of length \code{rec_fun} where each element includes the stock recruitment function for a given species. If an invalid recruitment function is selected, \code{NULL} is returned and a warning message is shown.
#' @slot eps A numeric value specifying a numerical offset. The default value is \code{1e-5}.
#' @slot recruit_params A list object of length \code{nfish} specifying the parameters for the recruitment function.
#' @details Sets the class of a LeMansParam object, which is then used to run the LeMans model.
#' @export
setClass("LeMans_param", slots=
           c(nsc="numeric",
             nfish="numeric",
             phi_min="numeric",
             l_bound="numeric",
             u_bound="numeric",
             mid="numeric",
             species_names="character",
             Linf="numeric",
             W_a="numeric",
             W_b="numeric",
             k="numeric",
             Lmat="numeric",
             mature="matrix",
             sc_Linf="numeric",
             wgt="matrix",
             phi="matrix",
             ration="matrix",
             other="numeric",
             M1="matrix",
             suit_M2="list",
             Qs="array",
             stored_rec_funs="list",
             eps="numeric",
             recruit_params="list"))

#' A constructor for the \code{LeMansParam} class
#'
#' @description A constructor for the \code{\link{LeMansParam}} class.
#' @param df A data frame with \code{nfish} rows and the following columns: \code{Linf}, \code{W_a}, \code{W_b}, \code{k} and \code{Lmat}. See below for definitions of each of these parameters.
#' @param gdf A data frame with \code{nfish} rows and the following columns: \code{curve}, \code{catch_species}, \code{max_catchability}, \code{k} and \code{gear_name}. See below for definitions of each of these parameters.
#' @param nsc A numeric value representing the number of length classes in the model.
#' @param nfish A numeric value representing the number of fish species in the model.
#' @param pred_mu A numeric value representing the preferred predator-prey mass ratio.
#' @param pred_sigma A numeric value representing the width of the weight preference function.
#' @param other A numeric value representing the amount of other food (g) available from prey that is not explicitly represented in the model. The default is \code{1e12}.
#' @param bounds An optional argument specifying the bounds of the length classes.
#' @param calc_phi_min A logical statement indicating whether \code{phi_min} should be calculated within the function. The default is \code{FALSE}.
#' @param phi_min A fixed numeric value of \code{phi_min}, which represents the time step of the model. This parameter is only required if \code{calc_phi_min=FALSE}. The default is \code{0.1}.
#' @param vary_growth A logical statement indicating whether growth efficiency should vary for each species (\code{vary_growth=TRUE}) or be fixed at the value given by \code{fixed_growth} (\code{vary_growth=FALSE}). The default is \code{TRUE}.
#' @param growth_eff If \code{vary_growth==TRUE}, \code{growth_eff} is a numeric representing the growth efficiencies of a fish of length 0. If \code{vary_growth==FALSE}, \code{growth_eff} is a numeric value of length \code{1} representing a fixed growth efficiency for all fish. The default is 0.5.
#' @param growth_eff_decay A numeric value specifying the rate at which growth efficiency decreases as length increases to \code{Linf}. The default is 0.11.
#' @param eps A numeric value specifying a numerical offset. The default value is \code{1e-5}.
#' @param force_mature A logical statement indicating whether to force maturity for all fish in the largest length class. The default is \code{TRUE}.
#' @param species_names A numeric or character vector of length \code{nfish} that denotes the names of the species in the model.
#' @param Linf A numeric vector of length \code{nfish} representing the mid-point of the length classes.
#' @param W_a A numeric vector representing the parameter \code{a} in the length-weight conversion.
#' @param W_b A numeric vector representing the parameter \code{b} in the length-weight conversion.
#' @param k k A numeric vector of length \code{nfish} representing the von Bertalanffy growth parameter \code{(1/yr)} for each species.
#' @param Lmat A numeric vector of length \code{nsc} representing the length at which 50\% of the individuals are mature.
#' @param kappa A numeric vector of length \code{nfish} representing the rate of change from immature to mature fish.
#' @param tau A matrix with dimensions \code{nfish} and \code{nsc}. Row indices represent predators and column indices represent prey. A value of 1 at location \code{i}, \code{j} indicates prey \code{j} is eaten by predator \code{i}.
#' @param rec_fun A character vector representing the stock recruitment function to be applied to each species. The default value is \code{"hockey-stick"} but \code{rec_fun} can take a value of \code{"Ricker"}, \code{"Beverton-Holt"}, \code{"constant"}, or \code{"linear"} for each species.
#' @param recruit_params A list object of length \code{nfish} specifying the parameters for the recruitment function.
#' @param natmort_opt A character vector of length \code{1} describing the mortality function to be used for the species. The default value is \code{"std_RNM"} but can take a value of \code{"constant"} or \code{"linear"}. See \code{calc_M1} for more information.
#' @param Nmort A numeric vector of length \code{nfish} representing the maximum background mortality of each species. The default is \code{0.8} for all species.
#' @param prop A numeric vector of length \code{nfish} representing the proportion of length classes that have a non-zero background mortality. This is required only when the \code{natmort_opt} mortality function is used. The default is \code{3/4}.
#' @param curve A character vector of almost any length describing the type of curve to be used to determine the catchability of each species by fishing gear. By default, \code{curve} is a character vector of length \code{nfish} that takes the value \code{"logistic"} in each element, but it can also take a value of \code{"log-gaussian"} or \code{"knife-edge"}. If a custom curve is required for a particular species and/or fishing gear, the curve must be specified in \code{custom}.
#' @param catch_species A numeric value or character string describing the species to apply the catchability curve to.
#' @param max_catchability A numeric vector of length \code{curve} describing the maximum catchability for each catchability curve.
#' @param gear_name A character vector of the same length as \code{curve} and \code{species} describing the fishing gear that each element of \code{curve} and \code{species} relates to. By default, \code{gear_name} is a character vector of length \code{curve} that takes the value \code{paste("gear_", 1:length(curve), sep = "")}.
#' @param custom An array with dimensions \code{nsc}, \code{nfish} and the number of custom catchability curves that are required. \code{custom} represents the catchability of each species by the gears specified using custom catchability curves. By default, \code{custom} is set to \code{NULL}.
#' @param ... Additional arguments for calculating catchability curves. See \code{\link{calc_Q}} for more details.
#' @return An object of class \code{LeMansParam} for use in the LeMans model.
#' @details Converts objects of class data frame or vector to class LeMansParams for use in the LeMans model. \code{Linf}, \code{W_a}, \code{W_b}, \code{k} and \code{Lmat} are required as either a data frame or as vectors.
#' @seealso \code{\link{run_LeMans}}, \code{\linkS4class{LeMans_param}}
#' @examples
#' # To run the model with all inputs specified explicitly:
#' # Set up species-specific parameters
#' Linf <- NS_par$Linf # the von-Bertalanffy asymptotic length of each species (cm).
#' W_a <- NS_par$W_a # length-weight conversion parameter.
#' W_b <- NS_par$W_b # length-weight conversion parameter.
#' k <- NS_par$k # the von-Bertalnaffy growth parameter.
#' Lmat <- NS_par$Lmat # the length at which 50% of individuals are mature (cm).
#'
#' NS_params <- LeMansParam(species_names=NS_par$species_names, Linf=Linf, k=k, W_a=W_a, W_b=W_b,
#' Lmat=Lmat, tau=NS_tau, recruit_params=list(a=NS_par$a, b=NS_par$b), eta=rep(0.25, 21), L50=Lmat,
#' other=NS_other)
#'
#' ###############################################
#' # Alternatively:
#' NS_params <- LeMansParam(NS_par, tau=NS_tau, eta=rep(0.25, 21), L50=NS_par$Lmat, other=NS_other)
#' @export
setGeneric("LeMansParam",function(df, gdf, ...)
  standardGeneric('LeMansParam'))

# Takes df and separates it into constituant parts, before passing it on to the next function
#' @rdname LeMansParam
setMethod("LeMansParam", signature(df="ANY", gdf="ANY"),
          function(df, gdf, nfish=nrow(df), nsc=32,
                   # Global options
                   pred_mu=-2.25, pred_sigma=0.5, other=1e12,
                   # Options that will be rarely used
                   bounds=NULL, calc_phi_min=FALSE, phi_min=0.1, vary_growth = TRUE, growth_eff=0.5, growth_eff_decay=0.11, eps=1e-5, force_mature=TRUE,
                   # Species params
                   species_names=paste("species", 1:nfish, sep="_"),
                   kappa=rep(10, nfish),
                   tau=matrix(1, nfish, nfish),
                   rec_fun=rep("hockey-stick", nfish), recruit_params=list(a=18.835-4.133*df$Linf, b=rep(1e3/nfish, nfish)),
                   natmort_opt=rep("std_RNM", nfish), Nmort=rep(0.8, nfish), prop=rep(3/4, nfish),
                   # Fishing
                   curve=rep("logistic", nfish), catch_species=((0:(length(curve)-1))%%nfish)+1, max_catchability=rep(1, length(curve)), gear_name=paste("gear_", 1:length(curve), sep=""),
                   custom=NULL, ...) {
            if (is.data.frame(df)==F) {
              stop("df must be a data frame")
            }

            tmp <- setdiff(c("Linf", "k", "W_a", "W_b", "Lmat"), names(df))
            if (length(tmp)>0) {
              stop(paste("The following columns are essential but do not currently exist. Please check the column names and try again\n:", paste(tmp, collapse=", "), sep=""))
            }
            Linf <- df$Linf
            k <- df$k
            W_a <- df$W_a
            W_b <- df$W_b
            Lmat <- df$Lmat

            if (any(names(df)=="species_names")) {
              species_names <- df$species_names
            }
            if (any(names(df)=="kappa")) {
              kappa <- df$kappa
            }
            if (any(names(df)=="rec_fun")) {
              rec_fun <- df$rec_fun
            }
            if (any(names(df)=="natmort_opt")) {
              natmort_opt <- df$natmort_opt
            }
            if (any(names(df)=="rec_fun")) {
              rec_fun <- df$rec_fun
            }
            if (any(names(df)=="Nmort")) {
              Nmort <- df$Nmort
            }
            if (any(names(df)=="prop")) {
              prop <- df$prop
            }

            # Warning block
            tmp <- setdiff(names(df), c("Linf", "k", "W_a", "W_b", "Lmat", "species_names", "kappa", "rec_fun", "natmort_opt", "Nmort", "prop"))

            if (length(tmp)>0) {
              recruit_params[tmp] <- df[tmp]
              warning(paste("The following columns of df do not match any of the species arguments and were therefore added to recruit_params:\n", paste(tmp, collapse = ", "), sep=""))
            }

            if (missing("gdf")) {
              LeMansParam(nfish=nfish, nsc=nsc,
                          # Global options
                          pred_mu=pred_mu, pred_sigma=pred_sigma, other=other,
                          # Options that will be rarely used
                          bounds=bounds, calc_phi_min=calc_phi_min, phi_min=phi_min, vary_growth=vary_growth, growth_eff=growth_eff, growth_eff_decay=growth_eff_decay, eps=eps, force_mature=force_mature,
                          # Species params
                          species_names=species_names,
                          Linf=Linf,
                          k=k,
                          W_a=W_a, W_b=W_b,
                          Lmat=Lmat, kappa=kappa,
                          tau=tau,
                          rec_fun=rec_fun, recruit_params=recruit_params,
                          natmort_opt=natmort_opt, Nmort=Nmort, prop=prop,
                          # Fishing
                          curve=curve, catch_species=catch_species, max_catchability=max_catchability, gear_name=gear_name, custom=custom, ...=...)
            } else {
              LeMansParam(gdf=gdf, nfish=nfish, nsc=nsc,
                          #### global options
                          pred_mu=pred_mu, pred_sigma=pred_sigma, other=other,
                          #### options that will be rarely used
                          bounds=bounds,calc_phi_min=calc_phi_min, phi_min=phi_min, vary_growth=vary_growth, growth_eff=growth_eff, growth_eff_decay=growth_eff_decay, eps=eps, force_mature=force_mature,
                          #### possible data frame
                          species_names=species_names,
                          Linf=Linf,
                          k=k,
                          W_a=W_a, W_b=W_b,
                          Lmat=Lmat, kappa=kappa,
                          ### possible species
                          tau=tau,
                          rec_fun=rec_fun, recruit_params=recruit_params,
                          natmort_opt=natmort_opt, Nmort=Nmort, prop=prop,
                          ## fishing
                          curve=curve, catch_species=catch_species, max_catchability=max_catchability, gear_name=gear_name, custom=custom, ...=...)
            }
          })



#############################this ones takes the gdf data frame and does the same as above if  it exits
#' @rdname LeMansParam
setMethod("LeMansParam", signature(df="missing", gdf="ANY"),
          function(gdf, nfish=length(Linf), nsc=32,
                   #### global options
                   pred_mu=-2.25, pred_sigma=0.5,other=1e12,
                   #### options that will be rarely used
                   bounds=NULL, calc_phi_min=FALSE, phi_min=0.1, vary_growth = FALSE, growth_eff=0.5, growth_eff_decay=0.11, eps=1e-5, force_mature=TRUE,
                   Linf,
                   k,
                   W_a, W_b,
                   Lmat,
                   #### possible data frame
                   species_names=paste("species", 1:nfish, sep="_"),
                   kappa=rep(10,nfish),
                   ### possible species
                   tau=matrix(1, nfish, nfish),
                   rec_fun=rep("hockey-stick", nfish), recruit_params=list(a=18.835-4.133*Linf, b=rep(1e3/nfish, nfish)),
                   natmort_opt=rep("std_RNM", nfish), Nmort=rep(0.8, nfish), prop=rep(3/4, nfish),
                   ## fishing
                   curve=rep("logistic", nfish), catch_species=((0:(length(curve)-1))%%nfish)+1, max_catchability=rep(1, length(curve)), gear_name=paste("gear_", 1:length(curve), sep=""),
                   custom=NULL, ...) {
            ####check that gdf is data frame
            if (is.data.frame(gdf)==F) {
              stop("gdf should be a data frame")
            }
            tmp <- setdiff(c("curve", "catch_species", "max_catchability", "gear_name"), names(gdf))
            if (length(tmp)>0) {
              warning(paste("The following columns are not in the dataframe, default values will be used unless a seperate vector has been supplied\n:", paste(tmp, collapse=", "), sep=""))
            }
            if ("curve" %in% colnames(gdf)) {
              curve <- gdf$curve
            }
            if (is.numeric(gdf$catch_species) | is.integer(gdf$catch_species)) {
              catch_species <- gdf$catch_species
            } else {
              if (is.factor(gdf$catch_species)) {
                gdf$catch_species <- as.character(gdf$catch_species)
              }
              if (is.character(gdf$catch_species)) {
                catch_species_tmp <- sapply(gdf$catch_species, function(x, y) {which(y==x)}, y=species_names)
                if (is.integer(catch_species_tmp) & length(catch_species_tmp)==nrow(gdf)) {
                  catch_species <- catch_species_tmp
                } else {
                  stop("Incorrect species in catch_species")
                }
              } else {
                stop("catch_species is not a factor, numeric, integer or a character")
            }}

            if ("max_catchability" %in% colnames(gdf)) {
              max_catchability <- gdf$max_catchability
            }
            if ("gear_name" %in% colnames(gdf)) {
              gear_name <- gdf$gear_name
            }

            ###### setupMSnJB function (feed in the other bits)
            LeMansParam(nfish=nfish, nsc=nsc,
                  #### global options
                  pred_mu=pred_mu, pred_sigma=pred_sigma, other=other,
                  #### options that will be rarely used
                  bounds=bounds, calc_phi_min=calc_phi_min, phi_min=phi_min, vary_growth=vary_growth, growth_eff=growth_eff, growth_eff_decay=growth_eff_decay, eps=eps, force_mature=force_mature,
                  #### possible data frame
                  species_names=species_names,
                  Linf=Linf,
                  k=k,
                  W_a=W_a, W_b=W_b,
                  Lmat=Lmat, kappa=kappa,
                  ### possible species
                  tau=tau,
                  rec_fun=rec_fun, recruit_params=recruit_params,
                  natmort_opt=natmort_opt, Nmort=Nmort, prop=prop,
                  ## fishing
                  curve=curve, catch_species=catch_species, max_catchability=max_catchability, gear_name=gear_name, custom=custom, ...=...) ## ... is catchability parameters
            })

#################################### the final bit
#' @rdname LeMansParam
setMethod('LeMansParam', signature(df="missing", gdf="missing"),
          function(df, gdf, nfish=length(Linf), nsc=32,
                   #### global options
                   pred_mu=-2.25, pred_sigma=0.5, other=1e12,
                   #### options that will be rarely used
                   bounds=NULL, calc_phi_min=TRUE, phi_min=0.1, vary_growth=FALSE, growth_eff=0.5, growth_eff_decay=0.11, eps=1e-5, force_mature=TRUE,
                   #### possible data frame
                   species_names=paste("species", 1:nfish, sep="_"),
                   Linf, k, W_a, W_b, Lmat, kappa=rep(10, nfish),
                   ### possible species
                   tau = matrix(1, nfish, nfish),
                   rec_fun = rep("hockey-stick", nfish), recruit_params=list(a= 18.835-4.133*Linf, b=rep(1e3/nfish, nfish)),
                   natmort_opt=rep("std_RNM", nfish), Nmort=rep(0.8, nfish), prop=rep(3/4, nfish),
                   ## fishing
                   curve=rep("logistic", nfish), catch_species=((0:(length(curve)-1))%%nfish)+1, max_catchability=rep(1, length(curve)), gear_name=paste("gear_", 1:length(curve), sep=""), custom=NULL, ... ## ... is catchability parameters
                   ){

            # if (pred_mu <= 0) {
            #   stop("predator prey mass ratio cannot be 0 or negative")
            # }
            if (pred_sigma <= 0) {
              stop("pred_sigma must take a positive value")
            }
            #if (class(tau) != "matrix") {
            if (is(tau,"matrix")==FALSE){
             stop("tau must be a matrix of dimensions nfish and nfish")
            }
            #if (length(catch_species) < 0 | length(catch_species) > nfish) {
            #  stop("The length of catch_species is not equal to nfish")
            #}
            if (phi_min <= 0) {
              stop("phi_min must take a positive value")
            }
            if (growth_eff <= 0) {
              stop("growth_eff must take a positive value")
            }
            if (any(Nmort <= 0)) {
              warning("Nmort must take a positive value")
            }
            if (any(max_catchability < 0)) {
              warning("max_catchability cannot take a negative value")
            }
            if (length(species_names) < nfish | length(species_names) > nfish) {
              stop("The length of species_names is not equal to nfish")
            }
            if (is.character(species_names)==FALSE) {
              stop("species_names is not a character vector")
            }
            if (length(Linf) < nfish | length(Linf) > nfish) {
              stop("The length of Linf is not equal to nfish")
            }
            if (any(Linf <= 0)) {
              stop("Linf must take a positive value")
            }
            if (is.numeric(Linf)==FALSE) {
              stop("Linf is not a numeric vector")
            }
            if (length(k) < nfish | length(k) > nfish) {
              stop("The length of k is not equal to nfish")
            }
            if (is.numeric(k)==FALSE) {
              stop("k is not a numeric vector")
            }
            if (any(k < 0)) {
              stop("k cannot take a negative value")
            }
            if (length(W_a) < nfish | length(W_a) > nfish) {
              stop("The length of W_a is not equal to nfish")
            }
            if (is.numeric(W_a)== FALSE) {
              stop("W_a is not a numeric value")
            }
            if (any(W_a <= 0)) {
              stop("W_a must take a positive value")
            }
            if (length(W_b) < nfish | length(W_b) > nfish) {
              stop("The length of W_b is not equal to nfish")
            }
            if (is.numeric(W_b)== FALSE) {
              stop("W_b is not a numeric value")
            }
            if (any(W_b <= 0)) {
              warning("W_b must take a positive value")
            }
            if (length(Lmat) < nfish | length(Lmat) > nfish) {
              stop("The length of Lmat is not equal to nfish")
            }
            if (is.numeric(Lmat)== FALSE) {
              stop("Lmat is not a numeric vector")
            }
            if (any(Lmat <= 0)) {
              stop("Lmat must take a positive value")
            }
            if (any(Lmat >= Linf)) {
              warning("Lmat cannot be greater than Linf")
            }
            if (any(kappa <= 0)) {
              stop("kappa must take a positive value")
            }
            if (length(kappa) < nfish | length(kappa) > nfish ) {
              stop("The length of kappa is not equal to nfish")
            }
            if (is.numeric(kappa)== FALSE) {
              stop("kappa is not a numeric vector")
            }

            ### calculating
            if (is.null(bounds)) {
              maxsize <- max(Linf)*1.01 ## biggest size is 1% bigger than
              l_bound <- seq(0, maxsize, maxsize/nsc); l_bound <- l_bound[-length(l_bound)]
              u_bound <- seq(maxsize/nsc, maxsize, maxsize/nsc)
              mid <- l_bound+(u_bound-l_bound)/2
            } else {
              l_bound <- bounds[-length(bounds)]
              u_bound <- bounds[-1]
              mid <- l_bound+(u_bound-l_bound)/2
              nsc <- length(u_bound)
            }

            tmp <- calc_phi(k, Linf, nsc, nfish, u_bound, l_bound, calc_phi_min=calc_phi_min, phi_min=phi_min)
            phi <- tmp$phi
            phi_min <- tmp$phi_min
            tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min, vary_growth=vary_growth, growth_eff=growth_eff, growth_eff_decay=growth_eff_decay)
            ration <- tmp$ration
            sc_Linf <- tmp$sc_Linf
            wgt <- tmp$wgt
            g_eff <- tmp$g_eff

            mature <- calc_mature(Lmat, nfish, mid, kappa, sc_Linf, eps=eps, force_mature=force_mature)

            stored_rec_funs <- get_rec_fun(rec_fun)
            recruit_params <- do.call("Map", c(c, recruit_params))

            M1 <- calc_M1(nsc, sc_Linf, phi_min, natmort_opt=natmort_opt, Nmort=Nmort, prop=prop)

            prefs <- calc_prefs(pred_mu, pred_sigma, wgt, sc_Linf)

            suit_M2 <- calc_suit_vect(nsc, nfish, sc_Linf, prefs, tau)

            Qs <- calc_Q(curve=curve, species=catch_species, max_catchability=max_catchability, gear_name=gear_name, custom=custom, nsc, nfish, mid, l_bound, u_bound, species_names=species_names, ...)

            new("LeMans_param",
                nsc = nsc,
                nfish = nfish,
                phi_min = phi_min,
                l_bound = l_bound,
                u_bound = u_bound,
                mid = mid,
                species_names = species_names,
                Linf=Linf,
                W_a=W_a,
                W_b=W_b,
                k=k,
                Lmat=Lmat,
                mature=mature,
                sc_Linf=sc_Linf,
                wgt=wgt,
                phi=phi,
                ration=ration,
                other=other,
                M1=M1,
                suit_M2=suit_M2,
                Qs=Qs,
                stored_rec_funs=stored_rec_funs,
                eps=eps,
                recruit_params=recruit_params)
          })
