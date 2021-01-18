#' Ordinal forests
#'
#' Constructs prediction rules using the ordinal forest (OF) method presented in Hornung (2020). \cr
#' The following tasks can be performed using OF: 1) Predicting the values of an ordinal target variable for new observations based on covariate values (see \code{\link{predict.ordfor}});
#' 2) Ranking the importances of the covariates with respect to predicting the values of the ordinal target variable. \cr
#' The default values for the hyperparameters \code{nsets}, \code{ntreeperdiv}, \code{ntreefinal}, \code{npermtrial}, and \code{nbest}
#' were found to be in a reasonable range in Hornung (2020) and it should not be necessary to alter these values in most situations. \cr
#' For details on OFs see the 'Details' section below. \cr
#' NOTE: Starting with package version 2.4, it is also possible to obtain class probability 
#' predictions in addition to the class point predictions and variable importance values
#' based on the class probabilities through using the (negative) ranked probability score (Epstein, 1969)
#' as performance function (\code{perffunction="probability"}, new default). Using the ranked probability score in the variable importance can be expected to deliver more stable variable rankings, because the ranked probability score accounts for the ordinal scale of the dependent variable.  In situations in which there is no need for predicting class probabilities, but simply
#' class predictions are sufficient, other performance functions may be more suitable. See the subsection "Performance functions" in the "Details" section below for further details.
#'
#' @param depvar character. Name of the dependent variable in \code{data}.
#' @param data data.frame. Data frame containing the covariates and a factor-valued ordinal target variable. The order of the levels of the latter
#' has to correspond to the order of the ordinal classes of the target variable.
#' @param nsets integer. Number of score sets tried prior to the approximation
#' of the optimal score set.
#' @param ntreeperdiv integer. Number of trees in the smaller regression forests constructed 
#' for each of the \code{nsets} different score sets tried.
#' @param ntreefinal integer. Number of trees in the larger regression forest 
#' constructed using the optimized score set (i.e., the OF).
#' @param perffunction character. Performance function. The default is \code{"probability"}. See 'Details', subsection 'Performance functions' below and \code{\link{perff}}.
#' @param classimp character. Class to priorize if \code{perffunction="oneclass"}.
#' @param classweights numeric. Needed if \code{perffunction="custom"}: vector of length equal to the number of classes. Class weights - the
#' higher the weight w_j assigned to class j is chosen, the higher the accuracy of the OF with respect to discerning observations in class j from observations
#' not in class j will tend to be.
#' @param nbest integer. Number of best score sets used to calculate the optimized score set.
#' @param naive boolean. If set to \code{TRUE}, a naive ordinal forest is constructed, that is, the score set used for the 
#' classes of the target variable is not optimized, but instead the following (naive) scores
#' are used: 1,2,3,... Note that it is strongly recommended to set \code{naive=FALSE} (default). The only advantage of choosing \code{naive=TRUE} is
#'  that the computational burden is reduced. However, the precision of the predictions of a prediction rule obtained using 
#' naive ordinal forest can be considerably worse than that of a corresponding prediction rule obtained using ordinal forest.
#' @param num.threads integer. Number of threads. Default is number of CPUs available (passed to the modified \code{ranger} code).
#' @param npermtrial integer. Number of permutations of the class width ordering to try for
#' the 2th to the \code{nsets}th score set tried prior to the calculation of
#' the optimized score set.
#' @param permperdefault boolean. If set to \code{TRUE}, \code{npermtrial} different permutations will per default 
#' be tried for the 2th to the \code{nsets}th score set used during the optimization - also for J! < \code{nsets}. Default is \code{FALSE}.
#' @param mtry integer. Number of variables to possibly split at in each node. Default is the (rounded down) square root of the number variables. 
#' @param min.node.size integer. Minimal node size. Default is 5, except if \code{perffunction="probability"}, in which case the default is 10.
#' @param replace boolean. Sample with replacement. Default is \code{TRUE}.
#' @param sample.fraction numeric. Fraction of observations to sample. Default is 1 for sampling with replacement and 0.632 for sampling without replacement.
#' @param always.split.variables character. Character vector with variable names to be always selected in addition to the \code{mtry} variables tried for splitting.
#' @param keep.inbag boolean. Save how often observations are in-bag in each tree. Default is \code{FALSE}.
#'
#' @details
#'
#' \subsection{Introduction}{
#' The ordinal forest (OF) method allows ordinal regression with high-dimensional and low-dimensional 
#' data. After having constructed an OF prediction rule using a training dataset, it can be used to predict the values of the ordinal target variable for new observations.
#' Moreover, by means of the (permutation-based) variable importance measure of OF, it is also 
#' possible to rank the covariates with respect to their importances in the prediction of the values
#' of the ordinal target variable. \cr
#' OF is presented in Hornung (2020). See the latter publication for details on the method. In the
#' following, a brief, practice-orientated introduction to OF is provided.
#' }
#'
#' \subsection{Methods}{
#' The concept of OF is based on the following assumption: There exists a (possibly latent) refined continuous 
#' variable y* underlying the observed ordinal target variable y (y in \{1,...,J\}, J number of classes), where y* determines the
#' values of y. The functional relationship between y* and y takes the form of a monotonically increasing
#' step function. Depending on which of J intervals ]c_1,\verb{ }c_2], \verb{ }]c_2,\verb{ }c_3], \verb{ } ..., \verb{ } ]c_J,\verb{ }c_\{J+1\}[
#' contains the value of y*, the ordinal target variable y takes a different value.
#'
#' In situations in which the values of the continuous target variable y* are known, they can be used in regression techniques for continuous response variables.
#' The OF method is, however, concerned with settings in which only the values of the classes of the ordinal target variable are given.
#' The main idea of OF is to optimize score values s_1,...,s_J to be used in place of the class values 1,...,J of the ordinal target variable
#' in standard regression forests by maximizing the out-of-bag (OOB) prediction performance measured by a performance function g (see section 
#' "Performance functions").
#' 
#' The approximation of the optimal score set consists of two steps:\cr
#' 1) Construct a large number of regression forests (b in 1,...,\code{nsets}) featuring limited numbers 
#' of trees, where each of these uses as the values of the target variable a randomly generated score set 
#' s_\{b,1\},...,s_\{b,J\}. For each forest constructed, calculate the value of the performance function g 
#' using the OOB estimated predictions of the values of the ordinal target variable and the corresponding 
#' true values.\cr
#' 2) Calculate the approximated optimal score set s_1,...,s_J as a summary over the \code{nbest} best score sets generated in 1),
#' that is, those \code{nbest} score sets that were associated with the highest values of the performance function g.
#'
#' After calculating the optimized score set, a larger regression forest is constructed using this optimized score set
#' s_1,...,s_J for the class values 1,...,J of the target variable. This regression forest is the OF prediction rule.
#'
#' Except in the case of using the (negative) ranked probabilty score as performance function, prediction is performed by majority voting of the predictions of the individual trees in the OF.
#' If the (negative) ranked probabilty score is used as performance function, both class predictions and predicted class probabilities are provided: The class probabilities
#' are obtained by averaging over the class probabilities predicted by the individual trees and the class predictions are obtained as the classes with maximum class probabilites.
#'
#' OF features a permutation variable importance measure that, if \code{perffunction} is set to \code{"probability"}, uses the ranked probability score as error measure
#' and the misclassification error else.
#' }
#'
#' \subsection{Hyperparameters}{
#' There are several hyperparameters, which do, however, not have to be optimized by the user in general, because the default values
#' used for these hyperparameters were seen to be in a reasonable range and the results seem to be quite robust with respect to the 
#' choices of the hyperparameter values.
#' 
#' These hyperparameters are described in the following:
#' \itemize{
#'   \item \code{nsets} \verb{   } Default value: 1000. The default value of the number of considered score sets in the approximation of the optimal score set 
#' is quite large. A large number of considered score sets is necessary to attain a high chance that some of the score sets are close enough to the optimal score set,
#' that is, the score set that leads to the optimal OOB prediction performance with respect to the considered performance function (provided with the argument \code{perffunction}).
#' \item \code{ntreeperdiv} \verb{   } Default value: 100. A very small number of trees considered per tried
#' score set might lead to a too strong variability in the assessments of the performances achieved for the individual score sets. For ultra-high dimensional
#' covariate data it might be necessary to choose a higher value for \code{ntreeperdiv} than the default value 100.
#' \item \code{ntreefinal} \verb{   } Default value: 5000. The number of trees \code{ntreefinal}
#' plays the same role as in conventional regression forests.
#' \item \code{npermtrial} \verb{   } Default value: 500. As stated above it is necessary to consider a large number of tried score sets \code{nsets} in the
#' optimization in order to increase the chance that the best of the considered score sets are close to the optimal score set.
#' To further increase this chance, it is in addition necessary that the collection of score sets tried is heterogeneous enough across 
#' the iterations. OF uses a particular algorithm for sampling the score sets tried that leads to a strongly heterogeneous collection of sets.
#' This algorithm features the hyperparameter \code{npermtrial}, where it has been seen in Hornung (2020) that the results are quite 
#' robust with respect to the choice of the value of this parameter.
#' \item \code{nbest} \verb{   } Default value: 10. In the case of a relatively small value of \code{nsets}, it is important that the 
#' number \code{nbest} of best score sets used to calculate the optimized score set is not strongly misspecified. A too large value of \code{nbest} 
#' leads to including suboptimal score sets into the calculation of the optimized score set that are too distinct from the optimal score set.
#' Conversely, a too small value of \code{nbest} leads to a high variance of the optimized score set. The combination \code{nsets=1000} and \code{nbest=10}
#' should lead to a good trade-off between the heterogeneity of the considered score sets and the variance in the estimation.
#' In Hornung (2020) this combination delivered good results and it was seen that using a very large value of \code{nbest} can lead to worse results.
#' }
#' }
#'
#' \subsection{Performance functions}{
#' As noted above, the different score sets tried during the estimation of the optimal score set are assessed with respect to their OOB prediction performance.
#' The choice of the specific performance function used in these assessments determines the specific kind of performance the ordinal forest should feature:
#' \itemize{
#' \item \code{perffunction="probability"} \verb{   } This choice should be made if it is of interest to predict class probabilties for the observations.
#' The ranked probability score is calculated between the predicted probabilities for the J classes and the observed class values. Because smaller values
#' of the ranked probability score correspond to a better prediction, the negative ranked probability score is considered as performance functions.
#' \item \code{perffunction="equal"} \verb{   } This choice should be made if it is of interest to classify observations from each class with the same accuracy independent of the class sizes. 
#' Youden's J statistic is calculated with respect to each class ("observation/prediction in class j" vs. "observation/prediction NOT in class j" (j=1,...,J))
#' and the simple average of the J results taken.
#' \item \code{perffunction="proportional"} \verb{   } This choice should be made if the main goal is to classify
#' correctly as many observations as possible. The latter is associated with a preference for larger classes at the 
#' expense of a lower classification accuracy with respect to smaller classes.
#' Youden's J statistic is calculated with respect to each class and subsequently a weighted average of these values is taken - with weights 
#' proportional to the number of observations representing the respective classes in the training data.
#' \item \code{perffunction="oneclass"} \verb{   } This choice should be made if it is merely relevant that observations 
#' in class \code{categ} can be distinguished as reliably as possible from observations not in class \code{categ}.
#' Class \code{categ} must be passed to \code{ordfor} via the argument \code{categ}.
#' Youden's J statistic is calculated with respect to class \code{categ}.
#' \item \code{perffunction="custom"} \verb{   } This choice should be made if there is a particular ranking of the classes with respect to their importance. 
#' Youden's J statistic is calculated with respect to each class. Subsequently, a weighted average
#' with user-specified weights (provided via the argument \code{classweights}) is taken. In this way, classes with 
#' higher weights are prioritized by the OF algorithm over classes with smaller weights.
#' }
#' }
#'
#' @return
#' \code{ordfor} returns an object of class \code{ordfor}.
#' An object of class "\code{ordfor}" is a list containing the following components: 
#' \item{forestfinal}{ object of class \code{"ranger"}. Regression forest constructed using the optimized score set (i.e., the OF). Required by \code{\link{predict.ordfor}}.  }
#' \item{bordersbest}{ vector of length J+1. Average over the \code{nbest} best partitions of [0,1]. Required by \code{\link{predict.ordfor}}. }
#' \item{forests}{ list of length \code{nsets}. The regression forests constructed for the \code{nsets} different score sets tried prior to the approximation of the optimal score set. }
#' \item{perffunctionvalues}{ vector of length \code{nsets}. Performance function values for all score sets tried prior to the approximation of the optimal score set. }
#' \item{bordersb}{ matrix of dimension \code{nsets} x (J+1). All \code{nsets} partitions of [0,1] considered. }
#' \item{classes}{ character vector of length J. Classes of the target variable. }
#' \item{nsets}{ integer. Number of score sets tried prior to the approximation of the optimal score set. }
#' \item{ntreeperdiv}{ integer. Number of trees per score set considered. }
#' \item{ntreefinal}{ integer. Number of trees of the OF prediction rule. }
#' \item{perffunction}{ character. Performance function used.  }
#' \item{classimp}{ character. If \code{perffunction="oneclass"}: class to priorize, NA else. }
#' \item{nbest}{ integer. Number of best score sets used to approximate the optimal score set. }
#' \item{classfreq}{ table. Class frequencies. }
#' \item{varimp}{ vector of length p. Permutation variable importance for each covariate. If \code{perffunction="probability"}, the ranked probability score is used as error measure in the variable importance. For all other choices of the performance function, the misclassification error is used. }
#'
#' @references
#' \itemize{
#'   \item Hornung R. (2020) Ordinal Forests. Journal of Classification 37, 4–17. <\doi{10.1007/s00357-018-9302-x}>.
#'   \item Epstein E.S. (1969) A scoring system for probability forecasts of ranked categories, Journal of Applied Meteorology. 8(6), 985-987.
#'   }
#'
#' @examples
#' data(hearth)
#'
#' set.seed(123)
#' hearthsubset <- hearth[sort(sample(1:nrow(hearth), size=floor(nrow(hearth)*(1/2)))),]
#' ordforres <- ordfor(depvar="Class", data=hearthsubset, nsets=50, nbest=5, ntreeperdiv=100, 
#'   ntreefinal=1000)
#' # NOTE: nsets=50 is not enough, because the prediction performance of the resulting 
#' # ordinal forest will be suboptimal!! In practice, nsets=1000 (default value) or a 
#' # larger number should be used.
#'
#' ordforres
#'
#' sort(ordforres$varimp, decreasing=TRUE)
#'
#' @export
ordfor <- 
  function(depvar, data, nsets=1000, ntreeperdiv=100, ntreefinal=5000, perffunction = c("probability", "equal", "proportional", "oneclass", "custom"), classimp, classweights, nbest=10, naive=FALSE, num.threads = NULL, npermtrial=500, permperdefault = FALSE, mtry = NULL, min.node.size = NULL, replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), always.split.variables = NULL, keep.inbag = FALSE) {
    
    if (is.null(num.threads)) {
      num.threads = 0
    } else if (!is.numeric(num.threads) | num.threads < 0) {
      stop("Error: Invalid value for num.threads")
    }
    
    perffunction <- perffunction[1]
    
    # Extract the target variable:
    y <- eval(parse(text=paste("data$", depvar, sep="")))
    
    if(!is.factor(y))
      stop("Error: dependent variable must be factor valued.")
    
    
    # Number of classes of target variable:
    J <- length(levels(y))
    
    
    # Make a version of the dataset that does not contain the
    # factor-valued target variable (in the iterations numeric
    # target variables are used and the factor-valued target
    # variable is remove here in order to avoid that it is
    # used as a covariate in rangerordfor):
    datait <- data
    eval(parse(text=paste("datait$", depvar, " <- NULL", sep="")))
    datait$ymetric <- NA
    
    ynum <- as.numeric(y)       
    
	
	if(is.null(min.node.size)) {
	
	if(perffunction=="probability")
	min.node.size <- 10
	else
	min.node.size <- 5
	
	} else {
	
		if(perffunction=="probability") {
	if(min.node.size > 10)
	  min.node.size <- 10
    warning("'min.node.size' must not be smaller than 10, if perffunction = 'probability'. -> 'min.node.size.' set to 10.")
    }
	
	}
	
    
    # The partition optimization is only performed if naive=FALSE:	
    
    if(!naive) {
      
      if(missing(classimp) & (perffunction=="oneclass"))
        stop("Class to priorize has to be provided if perffunction = 'oneclass'.")
      
      if(missing(classimp) & (perffunction!="oneclass"))
        classimp <- levels(y)[1]
      
      if(missing(classweights)) {
        classweights <- rep(1/J, J)
        if(perffunction=="custom")
          warning("perffunction = 'custom', but argument 'classweights' was missing. 'classweights' set to 1/J,...,1/J.")
      }
      
      if(missing(classimp) & (perffunction!="oneclass"))
        classimp <- levels(y)[1]
      
      # Assign performance function:
      
      if(perffunction=="proportional")
        perff <- perff_proportional
      
      if(perffunction=="equal")
        perff <- perff_equal
      
      if(perffunction=="oneclass")
        perff <- perff_oneclass
      
      if(perffunction=="custom")
        perff <- perff_custom
      
      
      # Optimization of the partition:
      
      # List which will contain the regression forests in the form
      # of ranger objects:
      forests <- list()
      # Initiate vector of performance function values:
      perffunctionvalues <- 0
      
      # Initiate matrix that will contain the metric values - the scores - used in
      # place of the the ordinal classes in the regression trees:
      cb <- matrix(nrow=nsets, ncol=J)
      
      # Simulate borders:  
      bordersb <- simulateborders(J=J, nsets=nsets, npermtrial=npermtrial, permperdefault=permperdefault)
      
      # Loop constructing a regression forest for each considered partition:
      for(b in 1:nsets) {
        
        # Get borders for the b-th forest:
        borders <- bordersb[b,]
        
        # Get values of the target variable:
        cb[b,] <- (borders[-1] + borders[-length(borders)])/2
        
        # Generate target variable for b-th forest:
        datait$ymetric <- qnorm(cb[b,][ynum])
        
        
        # Construct b-th regression tree:
        if(perffunction!="probability") {
          forests[[b]] <- rangerordfor(dependent.variable.name = "ymetric", data = datait, 
                                       num.trees = ntreeperdiv, num.threads=num.threads, 
                                       mtry=mtry, min.node.size=min.node.size, replace=replace, 
                                       sample.fraction=sample.fraction, 
                                       always.split.variables=always.split.variables,
                                       write.forest=FALSE)
          # Obtain OOB predictions:
          allpred <- forests[[b]]$predictions
          
          # Get OOB predictions of the classes:
          ynumtestpred <- sapply(allpred, function(x) max(which(x >= qnorm(bordersb[b,]))))
          
          # Remove erratic predictions:
          keepbool <- !is.infinite(ynumtestpred)
          ynumtestpred <- ynumtestpred[keepbool]
          ynumwithoutinf <- ynum[keepbool]
          
        }				   
        else {
          
          forests[[b]] <- rangerordfor(dependent.variable.name = "ymetric", data = datait, 
                                       num.trees = ntreeperdiv, num.threads=num.threads, 
                                       mtry=mtry, min.node.size=min.node.size, replace=replace, 
                                       sample.fraction=sample.fraction, 
                                       always.split.variables=always.split.variables,
                                       write.forest=FALSE, borders=qnorm(bordersb[b,][-c(1,length(bordersb[b,]))]), userps=TRUE)

									   ynumtestpred <- lapply(forests[[b]]$predictionsrps, function(x) x[[1]])
									   
          # Remove erratic predictions:
          keepbool <- sapply(ynumtestpred, function(x) !(NaN %in% x))
          ynumtestpred <- ynumtestpred[keepbool]
          ynumtestpred <- do.call("rbind", ynumtestpred)
          ynumwithoutinf <- ynum[keepbool]
          
        }
        
        
        # Calculate value of the performance function:
        if(perffunction!="probability")
          perffunctionvalues[b] <- perff(ynumwithoutinf, ynumtestpred, classimp, classweights)
        else
          perffunctionvalues[b] <- -verification::rps(obs=ynumwithoutinf, pred=ynumtestpred)$rps
        
      }
      
      # Take the average over the nbest best borders:
      bordersbest <- c(0, colMeans(bordersb[order(perffunctionvalues, decreasing=TRUE)[1:nbest],,drop=FALSE][,c(-1,-ncol(bordersb)),drop=FALSE]), 1)
      
      # Scores of the target variable:
      cbest <- (bordersbest[-1] + bordersbest[-length(bordersbest)])/2
      
      # Generate target variable for the large forest with
      # the optimized target variable:
      datait$ymetric <- qnorm(cbest[ynum])
      
    }  else {
      
      # For naive=FALSE the scores of the target variable are simply 1,2,3,...:
      cbest <- 1:J
      
      # The information on the optimization of the partitions is 'NA':
      bordersbest <- NA
      forests <- NA
      perffunctionvalues <- NA
      bordersb <- NA
      nsets <- NA
      ntreeperdiv <- NA
      perffunction <- NA
      nbest <- NA
      
      # Metric target variable:
      datait$ymetric <- cbest[ynum]
      
    }
    
    
    
    # Construct ordinal forest:
    if(!is.na(bordersbest[1])) {
      if(perffunction!="probability") {
        forestfinal <- rangerordfor(dependent.variable.name = "ymetric", data = datait, 
                                    num.trees = ntreefinal, importance="permutation", 
                                    num.threads=num.threads, borders=qnorm(bordersbest[-c(1,length(bordersbest))]), 
                                    mtry=mtry, min.node.size=min.node.size, replace=replace, sample.fraction=sample.fraction, always.split.variables=always.split.variables, 
                                    keep.inbag=keep.inbag)
      } else {
        forestfinal <- rangerordfor(dependent.variable.name = "ymetric", data = datait, 
                                    num.trees = ntreefinal, importance="permutation", 
                                    num.threads=num.threads, borders=qnorm(bordersbest[-c(1,length(bordersbest))]), 
                                    mtry=mtry, min.node.size=min.node.size, replace=replace, sample.fraction=sample.fraction, always.split.variables=always.split.variables, 
                                    keep.inbag=keep.inbag, userps=TRUE)
      }
    } else
      forestfinal <- rangerordfor(dependent.variable.name = "ymetric", data = datait, 
                                  num.trees = ntreefinal, importance="permutation", num.threads=num.threads, borders=(2:J) - 0.5,
                                  mtry=mtry, min.node.size=min.node.size, replace=replace, sample.fraction=sample.fraction, always.split.variables=always.split.variables, 
                                  keep.inbag=keep.inbag)
    
    # Ordinal classes of the target variable:
    classes <- levels(y)
    
    # Class frequencies:
    classfreq <- table(y)
    
    # Output informations:
    res <- list(forestfinal=forestfinal, bordersbest=bordersbest, forests=forests,
                perffunctionvalues=perffunctionvalues, bordersb=bordersb, 
                classes=classes, nsets=nsets, ntreeperdiv=ntreeperdiv, ntreefinal=ntreefinal, 
                perffunction = perffunction, classimp=ifelse(!is.na(perffunction) & perffunction=="oneclass", classimp, NA), nbest=nbest, classfreq=classfreq, varimp=forestfinal$variable.importance)
    class(res) <- "ordfor"
    
    # Output results:
    return(res)
    
  }
