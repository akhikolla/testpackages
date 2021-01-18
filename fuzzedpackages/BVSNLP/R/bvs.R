#' High dimensional Bayesian variable selection using nonlocal priors
#'
#' @description This function performs Bayesian variable selection for high
#' dimensional design matrix using iMOM prior for non-zero coefficients. It
#' also performs adaptive hyperparameter selection for iMOM prior. Cleaning
#' the data in a preprocessing step and before any data analysis is done by
#' user preference. This function is for binary and survival time response
#' datasets. In the former, MCMC is used to search in the model space while for
#' the latter a stochastic search does that job. This function has the option
#' to do all the mentioned tasks in a parallel fashion, exploiting hundreds of
#' CPUs. It is highly recommended to use a cluster for this purpose. This
#' function also supports fixing covariates in variable selection process,
#' thus making them included in the final selected model with probability 1.
#' Categorical variable are also supported by this function as input covariates
#' to the selection process. They need to be well defined factor variables as
#' part of the input data frame. For the output, this function reports
#' necessary measurements that is common in Bayesian variable selection
#' algorithms. They include Highest Posterior Probability model, median
#' probability model and posterior inclusion probability for each of the
#' covariates in the design matrix.
#'
#' @param X The \code{n} times \code{p} input data frame containing the
#' covariates in the design matrix. The columns should represent genes and rows
#' represent the observed samples. The column names are used as gene names so
#' they should not be left as \code{NULL}. Moreover, the minimum number of
#' columns allowed is 3. The input data frame can also contain categorical
#' covariates that are appropriately defined as factor variables in R.
#' @param resp For logistic regression models it is the binary response
#' vector which could be either numeric or factor variable in R. For the Cox
#' proportional hazard models this is a two column matrix where the first
#' column contains survival time vector and the second column is the censoring
#' status for each observation.
#' @param prep A boolean variable determining if the preprocessing step should
#' be performed on the design matrix or not. That step contains removing
#' columns that have \code{NA}'s or all their elements are equal to 0, along
#' with standardizing non-binary columns. This step is recommended and thus the
#' default value is \code{TRUE}.
#' @param logT A boolean variable determining if log transform should be done
#' on continuous columns before scaling them in the preprocessing step.
#' Note that those columns should not contain any zeros or negative values.
#' @param fixed_cols A vector of indices showing the columns in the input
#' data frame that are not subject to the the selection procedure. These
#' columns are always in the final selected model. Note that if any of these
#' columns contain \code{NA}, they will be removed. Moreover, if a categorical
#' variable with \code{k} levels is chosen to be fixed, all \code{k-1} dummy
#' variables associated with it will be selected in the final model. 
#' @param eff_size This is the expected effect size in the model for a
#' standardized design matrix, which is basically the coefficient value that is
#' expected to occur the most based on some prior knowledge.
#' @param family Determines the type of data analysis. \code{logistic} is for
#' binary outcome data where logistic regression modeling is used, whereas
#' \code{survival} is for survival outcome data using Cox proportional
#' hazard model.
#' @param hselect A boolean variable indicating whether the automatic procedure
#' for hyperparameter selection should be run or not. The default value is
#' \code{TRUE}. Note that in this setting, \code{r} is always chosen to be 1.
#' @param nlptype Determines the type of nonlocal prior that is used in the
#' analyses. It can be "piMOM" for product inverse moment prior, or "pMOM" for
#' product moment prior. The default is set to piMOM prior.
#' @param r The paramter \code{r} of the iMOM prior, when no automatic
#' procedure for hyperparameter selection is done. As a result, this is
#' relevant only when \code{hselect = FALSE}, otherwise it is ignored.
#' @param tau The paramter \code{tau} of the iMOM prior, when no automatic
#' procedure for hyperparameter selection is done. As a result, this is
#' relevant only when \code{hselect = FALSE}, otherwise it is ignored.
#' @param niter Number of iterations. For binary response data, this
#' determines the number of MCMC iterations per CPU. For survival response data
#' this is the number of iterations per temperature schedule in the stochastic
#' search algorithm.
#' @param mod_prior Type of prior used for the model space. \code{unif} is
#' for a uniform binomial and \code{beta} is for a beta binomial prior. In the
#' former case, both hyper parameters in the beta prior are equal to \code{1},
#' but in the latter case those two hyper parameters are chosen as explained in
#' the reference papers. The default choice for this variable is the uniform
#' prior.
#' @param inseed The input seed for making the parallel processing
#' reproducible. This parameter is ignored in logistic regression models when
#' \code{cplng = FALSE}. The default value is \code{NULL} which means that each
#' time the search for model space is started from different starting points.
#' In case it is set to a number, it initializes the RNG for the first task and
#' subsequent tasks to get separate substreams.
#' @param cplng This parameter is only used in logistic regression models, and
#' indicating if coupling algorithm for MCMC output should be performed or not.
#' @param ncpu This is the number of cpus used in parallel processing. For
#' logistic regression models this is the number of parallel coupled chains
#' run at the same time. For survival outcome data this is the number of cpus
#' doing stochastic search at the same time to increase th enumber of visited
#' models.
#' @param parallel.MPI A boolean variable determining if MPI is used for
#' parallel processing or not. Note that in order to use this feature, your
#' system should support MPI and \code{Rmpi} and \code{doMPI} packages should
#' already be installed. The default is set to \code{FALSE} but in case your
#' system have the requirements, it is recommended to set this parameter to
#' \code{TRUE} as it is more efficient and results in faster run-time.
#' @return It returns a list containing different objects that depend on the
#' family of the model and the coupling flag for logistic regression models.
#' The following describes the objects in the output list based on different
#'combinations of those two input arguments.\cr \cr
#' \strong{1) } \code{family = logistic && cplng = FALSE}
#' \item{num_vis_models}{Number of unique models visited throughout the search
#' of the model space.}
#' \item{max_prob}{Maximum unnormalized probability among all visited models}
#' \item{HPM}{The indices of the model with highest posterior
#' probability among all visited models, with respect to the columns in
#' the output \code{des_mat}. This is not necessarily the same as the input
#' design matrix due to some changes to categorical variables. The names of
#' the selected columns can be checked using \code{gene_names}.
#' The corresponding design matrix is also one of the outputs that can be
#' checked in \code{des_mat}. If the output is \code{character[0]} it means
#' none of the variables of the design matrix is selected in the HPM and
#' HPM contains only the intercept.}
#' \item{beta_hat}{The coefficient vector for the selected model. The first
#' component is always for the intercept.}
#' \item{MPM}{The indices of median probability model. According to the paper
#' Barbieri et. al., this is defined to be the model consisting of those
#' variables whose posterior inclusion probability is at least 1/2. The order
#' of columns is similar to that is explained for \code{HPM}.}
#' \item{max_prob_vec}{A \code{1000} by \code{1} vector of unnormalized
#' probabilities of the first 1000 models with highest posterior probability
#' among all visited models. If the total number of visited models is less than
#' 1000, then the length of this vector would be equal to \code{num_vis_models}
#' . Note that the intercept is always used in calculating the probabilities
#' in this vector.}
#' \item{max_models}{A list containing models corresponding to
#' \code{max_prob_vec} vector. Each entry of this list contains the indices of
#' covariates for the model with posterior probability reported in the
#' corresponding entry in \code{max_prob_vec}. The intercept column is not
#' shown in this list as it is present in all of the models.}
#' \item{inc_probs}{A vector of length \code{p+1} containing the posterior
#' inclusion probability for each covariate in the design matrix. The order of
#' columns is with respect to processed design matrix, \code{des_mat}.}
#' \item{nlptype}{The type of nonlocal prior used in the analyses.}
#' \item{des_mat}{The design matrix used in the analysis where fixed columns
#' are moved to the beginning of the matrix and if \code{prep=TRUE}, the
#' columns containing \code{NA} are all removed. The reported indices in
#' selected models are all with respect to the columns of this matrix.}
#' \item{gene_names}{Names of the genes extracted from the design matrix.}
#' \item{r}{The hyperparameter for iMOM prior density function, calculated
#' using the proposed algorithm for the given dataset.}
#' \item{tau}{The hyperparameter for iMOM prior density function, calculated
#' using the proposed algorithm for the given dataset.}
#' \strong{2) } \code{family = logistic && cplng = TRUE}
#' \item{cpl_percent}{Shows what percentage of pairs of chains are coupled.}
#' \item{margin_probs}{A \code{k} by \code{1} vector of marginal probabilities
#' where element \code{i} shows the maximum marginal probability of the
#' data under the maximum model for the \eqn{i^{th}} pair of chains. \code{k}
#' is the number of paired chains which is the same as number of CPUs.}
#' \item{chains}{A \code{k} by \code{p} binary matrix, where each row is the
#' model for the \eqn{i^{th}} pair of chains. Note that the index of nonzero
#' elements are not necessarily in the same order as the input design matrix,
#' \code{X}, depending on existence of fixed columns in selection procedure.
#' As a result, always match the indices to the columns of the design matrix
#' that is reported as an output in \code{des_mat}.}
#' \item{nlptype}{The type of nonlocal prior used in the analyses.}
#' \item{cpl_flags}{A \code{k} by \code{1} binary vector, showing which pairs
#' are coupled, (=\code{1}) and which are not, (= \code{0}).}
#' \item{beta_hat}{A \code{k} by \code{(p+1)} matrix where each row is the
#' estimated coefficient for each modelin the rows of \code{Chains} variable.}
#' \item{uniq_models}{A list showing unique models with the indices of the
#' included covariates at each model.}
#' \item{freq}{Frequency of each of the unique models. It is used to find
#' the highest frquency model.}
#' \item{probs}{Unnormalized probability of each of the unique models.}
#' \item{des_mat}{The design matrix used in the analysis where fixed columns
#' are moved to the beginning of the matrix and if \code{prep=TRUE}, the
#' columns containing \code{NA} are all removed. The reported indices in
#' selected models are all with respect to the columns of this matrix.}
#' \item{gene_names}{Names of the genes extracted from the design matrix.}
#' \item{r}{The hyperparameter for iMOM prior density function, calculated
#' using the proposed algorithm for the given dataset.}
#' \item{tau}{The hyperparameter for iMOM prior density function, calculated
#' using the proposed algorithm for the given dataset.}
#' \strong{3) } \code{family = survival}
#' \item{num_vis_models}{Number of visited models during the whole process.}
#' \item{max_prob}{The unnormalized probability of the maximum model among
#' all visited models.}
#' \item{HPM}{The indices of the model with highest posterior
#' probability among all visited models, with respect to the columns in
#' \code{des_mat}. As a result, always look at the names of the selected
#' columns using \code{gene_names}. The corresponding design matrix is one of
#' the outputs that can be checked in \code{des_mat}.}
#' \item{MPM}{The indices of median probability model. According to the paper
#' Barbieri et. al., this is defined to be the model consisting of those
#' variables whose posterior inclusion probability is at least 1/2. The order
#' of columns is similar to that is explained for \code{HPM}.}
#' \item{beta_hat}{The coefficient vector for the selected model reported in
#' \code{HPM}.}
#' \item{max_prob_vec}{A \code{1000} by \code{1} vector of unnormalized
#' probabilities of the first 1000 models with highest posterior probability
#' among all visited models. If the total number of visited models is less than
#' 1000, then the length of this vector would be equal to \code{num_vis_models}
#' .}
#' \item{max_models}{A list containing models corresponding to
#' \code{max_prob_vec} vector. Each entry of this list contains the indices of
#' covariates for the model with posterior probability reported in the
#' corresponding entry in \code{max_prob_vec}.}
#' \item{inc_probs}{A \code{p} by \code{1} vector containing the posterior
#' inclusion probability for each covariate in the design matrix. The order of
#' columns is with respect to processed design matrix, \code{des_mat}.}
#' \item{nlptype}{The type of nonlocal prior used in the analyses.}
#' \item{des_mat}{The design matrix used in the analysis where fixed columns
#' are moved to the beginning of the matrix and if \code{prep=TRUE}, the
#' columns containing \code{NA} are all removed. The reported indices in
#' selected models are all with respect to the columns of this matrix.}
#' \item{start_models}{A \code{k} by \code{3} matrix showing the starting model
#' for each worker CPU. Obviously \code{k} is equal to the number of CPUs.}
#' \item{gene_names}{Names of the genes extracted from the design matrix.}
#' \item{r}{The hyperparameter for iMOM prior density function, calculated
#' using the proposed algorithm for the given dataset.}
#' \item{tau}{The hyperparameter for iMOM prior density function, calculated
#' using the proposed algorithm for the given dataset.}
#' @author Amir Nikooienejad
#' @references Nikooienejad, A., Wang, W., and Johnson, V. E. (2016). Bayesian
#' variable selection for binary outcomes in high dimensional genomic studies
#' using nonlocal priors. Bioinformatics, 32(9), 1338-1345.\cr\cr
#' Nikooienejad, A., Wang, W., & Johnson, V. E. (2020). Bayesian variable
#' selection for survival data using inverse moment priors. Annals of Applied
#' Statistics, 14(2), 809-828. \cr\cr
#' Johnson, V. E. (1998). A coupling-regeneration scheme for
#' diagnosing convergence in Markov chain Monte Carlo algorithms. Journal of
#' the American Statistical Association, 93(441), 238-248.\cr\cr
#' Shin, M., Bhattacharya, A., and Johnson, V. E. (2017). Scalable Bayesian
#' variable selection using nonlocal prior densities in ultrahigh dimensional
#' settings. Statistica Sinica.\cr\cr
#' Johnson, V. E., and Rossell, D. (2010). On the use of non-local prior
#' densities in Bayesian hypothesis tests. Journal of the Royal Statistical
#' Society: Series B (Statistical Methodology), 72(2), 143-170.\cr\cr
#' Barbieri, M. M., and Berger, J. O. (2004). Optimal predictive model
#' selection. The annals of statistics, 32(3), 870-897.
#' @seealso \code{\link{ModProb}}, \code{\link{CoefEst}}
#' @examples
#' ### Simulating Logistic Regression Data
#' n <- 200
#' p <- 40
#' set.seed(123)
#' Sigma <- diag(p)
#' full <- matrix(c(rep(0.5, p*p)), ncol=p)
#' Sigma <- full + 0.5*Sigma
#' cholS <- chol(Sigma)
#' Beta <- c(-1.9,1.3,2.2)
#' X <- matrix(rnorm(n*p), ncol=p)
#' X <- X%*%cholS
#' beta <- numeric(p)
#' beta[c(1:length(Beta))] <- Beta
#' XB <- X%*%beta
#' probs <- as.vector(exp(XB)/(1+exp(XB)))
#' y <- rbinom(n,1,probs)
#' colnames(X) <- paste("gene_",c(1:p),sep="")
#' X <- as.data.frame(X)
#'
#' ### Running 'bvs' function without coupling and with hyperparamter selection
#' ### procedure
#' bout <- bvs(X, y, family = "logistic", nlptype = "piMOM",
#'             mod_prior = "beta", niter = 50)
#'             
#' ### Highest Posterior Model
#' bout$HPM
#'
#'### Estimated Coefficients:
#' bout$beta_hat
#' 
#' ### Number of Visited Models:
#' bout$num_vis_models
bvs <- function(X, resp, prep = TRUE, logT = FALSE, fixed_cols = NULL,
                eff_size = 0.5, family = c("logistic", "survival"),
                hselect = TRUE, nlptype = "piMOM", r = 1, tau = 0.25,
                niter = 30, mod_prior = c("unif", "beta"), inseed = NULL,
                cplng = FALSE, ncpu = 4, parallel.MPI=FALSE){
  
  if(!class(X)=="data.frame") stop("input X should be a data frame!") 
  ol <- matprep(X, fixed_cols, prep, logT)
  X <- ol$fulmat
  gnames <- ol$gnames
  nf <- ol$nf
  
  if(family=="logistic"){
    y <- as.numeric(as.character(resp))
    dx <- dim(X)
    n <- dx[1]
    p <- dx[2]
    
    X <- cbind(rep(1, n), X)
    gnames <- c("Intercept", gnames)
    # ======================= Hyperparameters ========================
    
    cons <- 0;
    prp <- p / n
    ar <- 2 ^ n
    if (prp > 4 && ar < Inf){
      ac <- 0
      cons <- 0
      while (ar > ac) {
        cons <- cons + 1
        ac <- choose(p, cons)
      }
    }else{
      cons <- ceiling(log(p))
    }
    
    cons <- min(cons,ceiling(log(p)))
    if (mod_prior == "beta"){
      a <- cons; b <- p - a;
    }
    if (mod_prior == "unif"){
      a <- 1; b <- 1;
    }
    if (hselect){
      hyper <- HyperSelect(X, y, eff_size, nlptype, 20000, mod_prior,family)
      tau <- hyper$tau
      r <- 1
    }
    
    initProb <- cons / p
    exmat <- cbind(y, X)
    if(nlptype=="piMOM") nlptype_int <- 0
    if(nlptype=="pMOM") nlptype_int <- 1
    # =========================== Main ===============================
    
    if (!cplng){
      schain <- p
      while (schain > cons || schain == 0) {
        chain1 <- rbinom(p-nf, 1, initProb)
        schain <- sum(chain1)
      }
      chain1 <- as.numeric(c(rep(1,nf+1), chain1)) # one for the intercept.
      chain2 <- chain1
      Lregout <- logreg_bvs(exmat, chain1, nf, tau, r, nlptype_int, a, b, cons, niter,
                            cplng, chain2)
      
      #============ reading outputs ======================
      Hash_Key <- Lregout$hash_key; all_probs <- Lregout$hash_prob; 
      VisCovs <- Lregout$vis_covs;
      inds <- which(all_probs!=0);
      Hash_Key <- Hash_Key[inds]; all_probs <- all_probs[inds]; VisCovs <- VisCovs[inds];
      
      nvm <- length(unique(Hash_Key))
      uinds <- which(!duplicated(Hash_Key))
      all_probs <- all_probs[uinds]
      list_vis_covs <- VisCovs[uinds]
      
      outnum <- min(nvm,1000)
      sout <- sort(all_probs,decreasing = T,index.return=T)
      MaxMargs <- sout$x[1:outnum]
      minds <- sout$ix[1:outnum]
      max_marg <- MaxMargs[1]; indmax <- minds[1]
      sel_model <- list_vis_covs[[indmax]]
      
      gnames2 <- gnames[sel_model + 1];
      beta_hat <- Lregout$beta_hat; names(beta_hat) <- gnames2;
      gnames <- gnames[-1]
      sel_model <- sel_model[-1]
      
      MaxModels <- list(NULL)
      for (i in 1:outnum){
        MaxModels[[i]] <- list_vis_covs[[minds[i]]][-1]
      }
      inc_probs <- inc_prob_calc(all_probs,list_vis_covs,p+1)
      inc_probs <- inc_probs[-1]
      median_model <- which(inc_probs >= 0.5)
      
      #========================================#
      
      return(list(max_prob = max_marg, HPM = sel_model,
                  beta_hat = beta_hat, MPM = median_model, inc_probs= inc_probs,
                  max_prob_vec = MaxMargs, max_models = MaxModels,
                  num_vis_models = nvm, nlptype = nlptype, des_mat = X,
                  gene_names = gnames, r = r, tau = tau))
      
    } else {
      comb <- function(x, ...) {
        lapply(seq_along(x),
               function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
      }
      
      if(parallel.MPI){
        if (!requireNamespace("doMPI", quietly = TRUE)) {
          stop("Package doMPI needed for this function to work. Please install it.",
               call. = FALSE)
        } else {
          cl <- doMPI::startMPIcluster(count = ncpu)
          doMPI::registerDoMPI(cl)
          parout <- foreach(j = 1:ncpu, .combine = "comb", .multicombine = TRUE,
                            .init = list(list(), list(), list(), list()),
                            .packages = 'BVSNLP',
                            .options.mpi = list(seed = inseed)) %dopar% {
                              schain <- p
                              while (schain > cons || schain == 0) {
                                chain1 <- rbinom(p-nf, 1, initProb)
                                schain <- sum(chain1)
                              }
                              chain1 <- as.numeric(c(rep(1,nf+1), chain1))
                              schain <- p
                              while (schain > cons || schain == 0) {
                                chain2 <- rbinom(p-nf, 1, initProb)
                                schain <- sum(chain2)
                              }
                              chain2 <- as.numeric(c(rep(1,nf+1), chain2))
                              
                              Lregout <- logreg_bvs(exmat, chain1, nf, tau, r, nlptype_int, a, b,
                                                    cons, niter, cplng, chain2)
                              
                              maxChain <- as.logical(Lregout$max_chain)
                              maxMarg <- Lregout$max_prob
                              cflag <- Lregout$cplng_flag
                              bhat <- numeric(p + 1)
                              bhat[maxChain] <- Lregout$beta_hat
                              list(maxChain, maxMarg, cflag, bhat)
                            }
          doMPI::closeCluster(cl)
        }
        
      } else {
        cl <- makeCluster(ncpu)
        registerDoParallel(cl)
        opts <- list(preschedule=TRUE)
        if (!is.null(inseed)) {clusterSetRNGStream(cl, inseed)}
        ParOut <- foreach(j = 1:ncpu, .combine = "comb", .multicombine = TRUE,
                          .init = list(list(), list(), list(), list()),
                          .packages = 'BVSNLP',
                          .options.snow = opts ) %dopar% {
                            schain <- p
                            while (schain > cons || schain == 0) {
                              chain1 <- rbinom(p-nf, 1, initProb)
                              schain <- sum(chain1)
                            }
                            chain1 <- as.numeric(c(rep(1,nf+1), chain1))
                            schain <- p
                            while (schain > cons || schain == 0) {
                              chain2 <- rbinom(p-nf, 1, initProb)
                              schain <- sum(chain2)
                            }
                            chain2 <- as.numeric(c(rep(1,nf+1), chain2))
                            
                            Lregout <- logreg_bvs(exmat, chain1, nf, tau, r, nlptype_int, a, b,
                                                  cons, niter, cplng, chain2)
                            
                            maxChain <- as.logical(Lregout$max_chain)
                            maxMarg <- Lregout$max_prob
                            cflag <- Lregout$cplng_flag
                            bhat <- numeric(p + 1)
                            bhat[maxChain] <- Lregout$beta_hat
                            list(maxChain, maxMarg, cflag, bhat)
                          }
        stopCluster(cl)
      }
      
      MaxChain <- matrix(unlist(ParOut[[1]]), ncol = (p + 1), byrow = T)
      MaxMarg <- unlist(ParOut[[2]])
      cpl_flag <- unlist(ParOut[[3]])
      bhat <- matrix(unlist(ParOut[[4]]), ncol = (p + 1), byrow = T)
      
      cpl_percent <- sum(cpl_flag) / ncpu
      Final_Marg <- MaxMarg
      Final_Chains <- MaxChain
      
      D <- as.data.frame(cbind(Final_Chains, Final_Marg))
      Counts <- rep(1, length(Final_Marg))
      A <- aggregate(Counts, by = as.list(D), FUN = sum)
      Freq <- A[, p + 3]
      Probs <- A[, p + 2]
      UniqModels <- apply(A[, 1:(p + 1)], 1, function(x) which(x > 0))
      
      return(list(cpl_percent = cpl_percent, margin_probs = Final_Marg,
                  chains = Final_Chains, cpl_flags = cpl_flag, beta_hat = bhat,
                  freq = Freq, probs = Probs, uniq_models = UniqModels, nlptype = nlptype,
                  gene_names = gnames, r = r, tau = tau))
    }
  }
  
  ### ================================================================= 
  if(family=="survival"){
    
    TS <- resp
    time <- TS[, 1]
    status <- TS[, 2]
    
    sfidx <- nf+1
    
    dx <- dim(X)
    n <- dx[1]
    p <- dx[2]
    exmat <- cbind(time, status, X)
    if(nlptype=="piMOM") nlptype_int <- 0
    if(nlptype=="pMOM") nlptype_int <- 1
    # ======================= Hyperparameters ========================
    cons <- 1+nf
    
    if (mod_prior == "beta"){
      a <- cons; b <- p - a;
    }
    if (mod_prior == "unif"){
      a <- 1; b <- 1;
    }
    if (hselect){
      hyper <- HyperSelect(X, TS, eff_size, nlptype, 5000, mod_prior, family)
      tau <- hyper$tau
      r <- 1
    }
    
    ntimes <- 10
    d <- 2 * ceiling(log(p))
    temps <- seq(3, 1, length.out = ntimes)
    
    L <- ntimes
    J <- niter
    
    # =========================== Main ===============================
    
    comb <- function(x, ...) {
      lapply(seq_along(x),
             function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
    }
    
    if(parallel.MPI){
      if (!requireNamespace("doMPI", quietly = TRUE)) {
        stop("Package doMPI needed for this function to work. Please install it.",
             call. = FALSE)
      } else {
        cl <- doMPI::startMPIcluster(count = ncpu)
        doMPI::registerDoMPI(cl)
        parout <- foreach(j = 1:ncpu, .combine = "comb", .multicombine = TRUE,
                          .init = list(list(), list(), list(), list(), list(), list()),
                          .packages = 'BVSNLP',
                          .options.mpi = list(seed = inseed)) %dopar% {
                            cur_model <- sample(sfidx:p, 3)
                            if (nf > 0) cur_model <- c(1:nf,cur_model)
                            coxout <- cox_bvs(exmat, cur_model, nf, tau, r, nlptype_int, a, b,
                                              d, L, J, temps)
                            
                            maxmod <- coxout$max_model
                            maxprob <- coxout$max_prob
                            hashkey <- coxout$hash_key
                            allprobs <- coxout$all_probs
                            viscovs <- coxout$vis_covs_list
                            list(maxmod, maxprob, hashkey, allprobs, cur_model, viscovs)#,vismodels)
                          }
        doMPI::closeCluster(cl)
      }
    } else {
      cl <- makeCluster(ncpu)
      registerDoParallel(cl)
      opts <- list(preschedule=TRUE)
      if (!is.null(inseed)) {clusterSetRNGStream(cl, inseed)}
      parout <- foreach(j = 1:ncpu, .combine = "comb", .multicombine = TRUE,
                        .init = list(list(), list(), list(), list(), list(), list()),
                        .packages = 'BVSNLP',
                        .options.snow = opts ) %dopar% {
                          cur_model <- sample(sfidx:p, 3);
                          if (nf > 0) cur_model <- c(1:nf,cur_model)
                          coxout <- cox_bvs(exmat, cur_model, nf, tau, r, nlptype_int, a, b,
                                            d, L, J, temps)
                          
                          maxmod <- coxout$max_model
                          maxprob <- coxout$max_prob
                          hashkey <- coxout$hash_key
                          allprobs <- coxout$all_probs
                          viscovs <- coxout$vis_covs_list
                          list(maxmod, maxprob, hashkey, allprobs, cur_model, viscovs)#,vismodels)
                        }
      stopCluster(cl)
    }
    
    Hash_Key <- unlist(parout[[3]])
    All_Probs <- unlist(parout[[4]])
    CurModel <- matrix(unlist(parout[[5]]), ncol = 3, byrow = T)
    VisCovs <- NULL
    for (i in 1:ncpu){
      VisCovs <- c(VisCovs,parout[[6]][[i]])
    }
    
    num_vis_models <- length(unique(Hash_Key))
    uinds <- which(!duplicated(Hash_Key))
    all_probs <- All_Probs[uinds]
    list_vis_covs <- VisCovs[uinds]
    
    outnum <- min(num_vis_models,1000);
    sout <- sort(all_probs,decreasing = T,index.return=T)
    MaxMargs <- sout$x[1:outnum]
    minds <- sout$ix[1:outnum]
    max_marg <- MaxMargs[1]; indmax <- minds[1]
    sel_model <- list_vis_covs[[indmax]] + 1
    
    MaxModels <- list(NULL)
    for (i in 1:outnum){
      MaxModels[[i]] <- list_vis_covs[[minds[i]]] + 1
    }
    
    inc_probs <- inc_prob_calc(all_probs,list_vis_covs,p)
    median_model <- which(inc_probs >= 0.5)
    beta_hat <- CoefEst(X,TS,sel_model,nlptype,tau,r,"survival")
    
    return(list(num_vis_models = num_vis_models,
                max_prob = max_marg, HPM = sel_model, MPM = median_model,
                beta_hat = beta_hat, max_prob_vec = MaxMargs, max_models = MaxModels,
                inc_probs = inc_probs, nlptype = nlptype, des_mat = X, start_models = CurModel,
                r = r, tau = tau, gene_names = gnames))
  }
}

