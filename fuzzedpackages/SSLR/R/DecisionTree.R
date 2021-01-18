#' An S4 class to represent a class with more types values: null, numeric or character
#' @export
setClassUnion("nullOrNumericOrCharacter", c("NULL", "numeric","character"))

setClassUnion("nullOrNode", c("NULL"))

setClassUnion("nullOrTable", c("NULL","table"))

#' @title Class Node for Decision Tree
#' @description Class Node for Decision Tree
#' Slots: gini, num_samples, num_samples_per_class, predicted_class_value, feature_index
#' threshold, left, right, probabilities
setClass("Node",
         slots = c(
           gini = "numeric",
           num_samples = "numeric",
           num_samples_per_class = "nullOrTable",
           predicted_class_value = "nullOrNumericOrCharacter",
           feature_index = "character",
           threshold = "nullOrNumericOrCharacter",
           left = "nullOrNode",
           right = "nullOrNode",
           probabilities = "numeric"
         ),
         prototype = list(
           feature_index = NA_character_,
           num_samples = NA_integer_,
           threshold = NULL,
           gini = NA_integer_,
           left = NULL,
           right = NULL,
           probabilities = NA_integer_,
           predicted_class_value = NULL
         )
)

## add class "Node" to the union "nullOrNode"
setIs("Node", "nullOrNode")


#' @title Class DecisionTreeClassifier
#' @description Class DecisionTreeClassifier
#' Slots: max_depth, n_classes_, n_features_,  tree_, classes,  min_samples_split,
#' min_samples_leaf
#' @export
setClass("DecisionTreeClassifier",
         slots = c(
           max_depth = "numeric",
           n_classes_ = "numeric",
           n_features_ = "numeric",
           tree_ = "nullOrNode",
           classes = "character",
           min_samples_split = "numeric",
           min_samples_leaf = "numeric"
         ),
         prototype = list(
           max_depth = NA_integer_,
           n_classes_ = NA_integer_,
           n_features_ = NA_integer_,
           min_samples_split = NA_integer_,
           min_samples_leaf = NA_integer_,
           tree_ = NULL
         )
)

#' @title Gini or Variance by column
#' @description function used to calculate the gini coefficient or
#' variance according to the type of the column. This function is called
#' for the creation of the decision tree
#' @param X column to calculate variance or gini
gini_or_variance <- function(X){

  #If it is numeric, calculate var
  if(is.numeric(X)){
    return(var(X, na.rm = TRUE))
  }

  #If it is character, calculate gini
  else{
    return(calculate_gini(X))
  }

}

#' An S4 method to fit decision tree.
#' @param object DecisionTree object
#' @param ... This parameter is included for compatibility reasons.
#' @export
setGeneric("fit_decision_tree", function(object,...)
  standardGeneric("fit_decision_tree") )

#' @title Fit decision tree
#' @description method in class DecisionTreeClassifier used to build a Decision Tree
#' @param object A DecisionTreeClassifier object
#' @param X A object that can be coerced as data.frame. Training instances
#' @param y A vector with the labels of the training instances. In this vector
#' the unlabeled instances are specified with the value \code{NA}.
#' @param w weight parameter ranging from 0 to 1
#' @param min_samples_split the minimum number of observations to do split
#' @param min_samples_leaf the minimum number of any terminal leaf node
#' @export
setMethod(f="fit_decision_tree",
          signature="DecisionTreeClassifier",
          definition=function(object,X,y,min_samples_split = 20,
                              min_samples_leaf = ceiling(min_samples_split/3),
                              w = 0.5)
          {

            #W between 0 and 1
            if(w < 0 | w > 1){
              stop("W CAN BE BETWEEN [0,1]")
            }

            #min_samples_split >= 1
            if(min_samples_split <= 0){
              stop("min_samples_split should be greater than 1")
            }

            #min_samples_leaf >= 1
            if(min_samples_leaf <= 0){
              stop("min_samples_leaf should be greater than 1")
            }

            #min_samples_split >= min_samples_leaf
            if(min_samples_split < min_samples_leaf){
              stop("min_samples_split should be greater than min_samples_leaf")
            }


            #Calculate gini in entire labeled train
            gini_entire_labeled <- calculate_gini(y)

            #Calculate gini or variance according to column
            variances <- as.vector(sapply(X,gini_or_variance))

            #Parms list for recursive
            parms <- list(gini_entire_labeled = gini_entire_labeled , w = w, variances = variances,
                          min_samples_split = min_samples_split,
                          min_samples_leaf = min_samples_leaf)

            #Get num classes and num features
            object@n_classes_ <- as.numeric(length(unique(y)))
            object@n_features_ <- as.numeric(dim(X)[2])
            object@min_samples_leaf <- min_samples_leaf
            object@min_samples_split <- min_samples_split

            #Get names of classes
            if(is.factor(y))
              object@classes <- levels(y)

            #Call grow tree to build the tree
            object@tree_ <- grow_tree(object,X, y, parms)

            #If not root node
            if(is.null(object@tree_)){

              print(paste("NROW LABELED DATA IS ", length(y[!is.na(y)])))
              print(paste("min_samples_split is",min_samples_split))
              print(paste(" min_samples_leaf is", min_samples_leaf))

              stop("Root node is NULL. You can try with another min_samples_split OR min_samples_leaf")

            }

            #Return Decision Tree object
            return(object)
          }
)

#' @title Check value in leaf
#' @description Function to check value in leaf from numeric until character
#' @param value is the value in leaf node
#' @param threshold in leaf node
#' @return TRUE if <= in numeric or %in% in factor
check_value <- function(value,threshold){

  #Value returned
  result <- FALSE

  if(is.na(value))
    return(FALSE)

  if(is.numeric(threshold)){

    if(value <= threshold){
      result <- TRUE
    }

  }

  else if(is.character(threshold)){
    if(value %in% threshold){
      result <- TRUE
    }
  }
  result
}

#' An S4 method to predict inputs.
#' @param object DecisionTree object
#' @param ... This parameter is included for compatibility reasons.
setGeneric("predict_inputs", function(object,...)
  standardGeneric("predict_inputs") )

#' @title Predict inputs Decision Tree
#' @description Function to predict one input in Decision Tree
#' @param object DecisionTree object
#' @param inputs inputs to be predicted
#' @param type type prediction, class or prob
#' @export
setMethod(f="predict_inputs",
          signature="DecisionTreeClassifier",
          definition=function(object,inputs,type = "class")
          {

            `%notin%` <- Negate(`%in%`)

            if(type %notin% c("prob","class","numeric"))
              stop("Type should be prob, class or numeric")

            #Get first node in tree
            node = object@tree_

            #If node dont have left node OR right node
            while(!is.null(node@left) | !is.null(node@right)){

              #Break if conditions
              if(is.na(node@feature_index) | (is.null(node@left) && is.null(node@right)))
                break

              #Node is left if check value is True
              if(!is.null(node@left) & check_value(inputs[1,node@feature_index, drop = FALSE],node@threshold)){
                node = node@left
              }

              #else node is right
              else if(!is.null(node@right)){
                node = node@right
              }

              else
                break

            }

            #if type is class, return predicted class
            if(type == "class" | type == "numeric")
              return(node@predicted_class_value)

            #if type is prob, return probabilities
            else if(type == "prob")
              return(node@probabilities)
          }
)

#' An S4 method to best split
#' @param object DecisionTree object
#' @param ... This parameter is included for compatibility reasons.
setGeneric("best_split", function(object,...)
  standardGeneric("best_split") )

#' @title Best Split function
#' @description Function to get best split in Decision Tree.
#' Find the best split for node. "Beast" means that the mean
#' of impurity is the least possible.
#' To find the best division. Let's iterate through all the features.
#' All threshold / feature pairs will be computed in the numerical
#' features. In the features that are not numerical,
#' We get the best group of possible values
#' will be obtained based on an algorithm with the function
#' get_levels_categoric
#' @param object DecisionTree object
#' @param X is data
#' @param y is class values
#' @param parms parms in function
#' @return A list with: best_idx name of the feature with the best split or Null if it not be found
#' best_thr: threshold found in the best split, or Null if it not be found
setMethod(f="best_split",
          signature="DecisionTreeClassifier",
          definition=function(object,X,y, parms)
          {
            #Get the length from training instances
            m = length(y)

            #If m is 1 return list of nulls
            if(m <= 1){
              return(list(best_idx = NULL,
                          best_thr = NULL))
            }

            #Get actual gini
            best_gini <- calculate_gini(y)

            #Default values
            best_idx <- NULL
            best_thr <- NULL

            #Get num features
            n_features <- ncol(X)

            #sum_variances_total <- calculate_sum_variances(X,parms$variances)

            #Loop all features
            for(idx in 1:n_features){

              #IF FEATURE IS NUMERIC
              if(is.numeric(X[,idx])){

                #Get pairs values by classes
                myDf <- data.frame(X[,idx],y)

                #Delete duplicated
                myDf <- myDf[!duplicated(t(apply(myDf, 1, sort))),]

                #Order by value
                myDf <- myDf[order(myDf[,1]),]

                #Loop all pairs values
                for(j in 2:nrow(myDf)){

                  #Get threshold
                  threshold <- (myDf[j,1] + myDf[j -1,1]) / 2

                  #Get index of left and right node
                  index.left <- which(X[,idx] <= threshold)
                  index.right <- which(X[,idx] > threshold)

                  #If length of index with labeled data is greater than 0
                  if(length(y[!is.na(y[index.left])]) > 0 && length(y[!is.na(y[index.right])]) > 0){

                    #Get gini from left and right
                    gini_left <- calculate_gini(y[index.left])
                    gini_right <- calculate_gini(y[-index.left])

                    #Get length of index left
                    num_elements_left <- length(index.left)

                    #Calculate gini from labeled data in split: gini left + gini right
                    gini_labeled <- (num_elements_left/m) * gini_left + ((m - num_elements_left)/m) * gini_right

                    #Impurity in labeled is gini_labeld / gini in the entire labeled train (in parms)
                    gini_labeled <- gini_labeled/parms$gini_entire_labeled

                    #We need to calculate gini in feature
                    #In this case we need to calculate the variance in left and right

                    #Left
                    gini_variance_left <- (num_elements_left/m) * var(X[index.left,idx], na.rm = TRUE)

                    if(is.na(gini_variance_left))
                      gini_variance_left = 0

                    #Right
                    gini_variance_right <- ((m - num_elements_left)/m) * var(X[-index.left,idx], na.rm = TRUE)

                    if(is.na(gini_variance_right))
                      gini_variance_right = 0

                    #gini_variance_total <- sum((gini_variance_left + gini_variance_right)/parms$variances)

                    #Impurity in feature is:
                    impurity_feature <- (gini_variance_left + gini_variance_right)/parms$variances[idx]

                    #Gini in semi supervised is w * gini_labeled + (1 - w)/n_features * impurity_feature
                    gini <- parms$w * gini_labeled + ((1-parms$w)/n_features) * impurity_feature

                    #If gini < best_gini
                    if(gini < best_gini){
                      #Best gini is actual gini
                      best_gini <- gini
                      #Best idx is the actuyal colname
                      best_idx <- colnames(X)[idx]
                      #Best thr is actual threshold
                      best_thr <- threshold

                    }
                  }


                }

              }

              #IF FEATURE IS CHARACTER OR FACTOR
              else if(is.character(X[,idx]) | is.factor(X[,idx])){

                #Wee need best threshold in get_levels_categoric function
                threshold <- get_levels_categoric(X[,idx],y)

                #Get index of left and right node
                index.left <- which(X[,idx] %in% threshold)

                `%nin%` = Negate(`%in%`)

                #calculate index right
                index.right <- which(X[,idx] %nin% threshold)

                #If length of index with labeled data is greater than 0
                if(length(y[!is.na(y[index.left])]) > 0 && length(y[!is.na(y[index.right])]) > 0){

                  #Get gini from left and right
                  gini_left <- calculate_gini(y[index.left])
                  gini_right <- calculate_gini(y[-index.left])

                  #Get length of index left
                  num_elements_left <- length(index.left)

                  #Calculate gini from labeled data in split: gini left + gini right
                  gini_labeled <- (num_elements_left/m) * gini_left + ((m - num_elements_left)/m) * gini_right


                  #We need to calculate gini in feature
                  #In this case we need to calculate the gini in left and right

                  #Gini left
                  gini_variance_left <- (num_elements_left/m) * calculate_gini(X[index.left,idx])

                  if(is.na(gini_variance_left))
                    gini_variance_left = 0

                  #Gini right
                  gini_variance_right <- ((m - num_elements_left)/m) * calculate_gini(X[-index.left,idx])

                  if(is.na(gini_variance_right))
                    gini_variance_right = 0

                  #gini_variance_total <- sum((gini_variance_left + gini_variance_right)/parms$variances)

                  #Impurity in feature is:
                  impurity_feature <- (gini_variance_left + gini_variance_right)/parms$variances[idx]


                  #Gini in semi supervised is w * gini_labeled + (1 - w)/n_features * impurity_feature
                  gini <- parms$w * gini_labeled + ((1-parms$w)/n_features) * impurity_feature


                  #If gini < best_gini
                  if(gini < best_gini){

                    #Best gini is actual gini
                    best_gini <- gini
                    #Best idx is the actuyal colname
                    best_idx <- colnames(X)[idx]
                    #Best thr is actual threshold
                    best_thr <- threshold

                  }

                }

              }

            }

            #Return list
            return(list(best_idx = best_idx,
                        best_thr = best_thr))


          }
)

#' An S4 method to grow tree.
#' @param object DecisionTree object
#' @param ... This parameter is included for compatibility reasons.
setGeneric("grow_tree", function(object,...)
  standardGeneric("grow_tree") )

#' @title Function grow tree
#' @description Function to grow tree in Decision Tree
#' @param object DecisionTree instance
#' @param X data values
#' @param y classes
#' @param parms parameters for grow tree
#' @param depth depth in tree
#' @importFrom methods new
setMethod(f="grow_tree",
          signature="DecisionTreeClassifier",
          definition=function(object,X,y,parms, depth = 0)
          {
            #Check if exists elements labeled
            #&
            #nrow(X) < min_samples_leaf
            #It not can be leaf node
            if(length(y[!is.na(y)]) > 0 & length(y) >= object@min_samples_leaf){


              #Get num samples by class (table)

              num_samples_per_class <- NULL

              if(!is.numeric(y))
              num_samples_per_class <- table(y)

              #Get predicted class
              predicted_class_value <- names(which.max(table(y)))

              #Get probabilities

              probabilities <- NaN

              if(!is.numeric(y))
              probabilities <- as.numeric(table(y[!is.na(y)])/length(y[!is.na(y)]))

              #Create actual node
              node = new("Node",gini = calculate_gini(y), num_samples = length(y),
                         num_samples_per_class = num_samples_per_class,
                         predicted_class_value=predicted_class_value, probabilities = probabilities)

              #Check if depth <= max_depth
              if(as.numeric(depth) < as.numeric(object@max_depth) & calculate_gini(y)!=0
                 & length(y[!is.na(y)]) >= object@min_samples_split){

                #Get best split
                results_node <- best_split(object,X, y, parms = parms)

                #If not null, we nned to call recursively
                if(!is.null(results_node$best_idx)){
                  node@feature_index <- results_node$best_idx
                  node@threshold <- results_node$best_thr

                  #We need index left in threshold
                  indices_left <- NULL

                  #If feature is numeric
                  if(is.numeric(X[,node@feature_index])){
                    indices_left <- which(X[,node@feature_index,drop=FALSE] <= node@threshold)
                  }

                  #If feature is factor or character
                  else if(is.factor(X[,node@feature_index]) | is.character(X[,node@feature_index])){
                    indices_left <- which(X[,node@feature_index,drop=FALSE] %in% node@threshold)
                  }

                  if(ncol(X) > 0){

                    #Get train left and right nodes
                    X_left <- X[indices_left,]
                    y_left <- y[indices_left]

                    X_right <- X[-indices_left,]
                    y_right <- y[-indices_left]

                    #Call recursively to get left and right nodes
                    node@left <- grow_tree(object,X_left,y_left,parms, depth + 1)
                    node@right <- grow_tree(object,X_right,y_right, parms, depth + 1)
                  }

                }

              }

              #Return actual node
              return(node)

            }

            #If not exists labeled data
            else{
              return(NULL)
            }

          }

)


calculate_sum_variances <- function(X,variances){

  results <- as.vector(apply(X, 2 ,var)) / variances

  sum(results)

}

#' @title Function calculate gini
#' @description Function to calculate gini index.
#' Formula is: 1 - n:num_classes sum probabilitie_class ^ 2
#' @param column_factor class values
calculate_gini <- function(column_factor){

  #Get value from labeled
  column_factor <- column_factor[!is.na(column_factor)]

  #Get size of labeled data
  m <- length(column_factor)

  gini <- NA

  #If classification
  if(is.factor(column_factor) | is.character(column_factor)){

    #Get frequencies per class
    num_parent <- as.vector(table(column_factor))
    #Calculate gini
    #1 - n:num_classes Σ probabilitie_class ^ 2
    gini <- 1- sum( (num_parent / m) ** 2)
  }

  #Is numeric, calculate SSE
  else{
    media <- mean(column_factor,na.rm = TRUE)
    gini <- sum((as.numeric(column_factor) - media)**2,na.rm = TRUE)
  }


  #If it is na, return 1 (maximum)
  if(is.nan(gini))
    gini <- 1

  return(gini)
}



#' @title Function to compute Gini index
#' @description Function to compute Gini index
#' From: https://freakonometrics.hypotheses.org/20736
#' @param y values
#' @param classe classes
gini_prob=function(y,classe){
  n_values <- length(unique(y))

  T=table(y,classe)
  nx=apply(T,2,sum)
  n=sum(T)
  pxy=T/matrix(rep(nx,each=n_values),nrow=n_values)
  omega=matrix(rep(nx,each=n_values),nrow=n_values)/n
  g=-sum(omega*pxy*(1-pxy))
  return(g)
}

#' @title Function to get gtoup from gini index
#' @description Function to get group from gini index.
#' Used in categorical variable
#' From: https://freakonometrics.hypotheses.org/20736
#' @param column is the column
#' @param Y values
get_levels_categoric <- function(column,Y){

  Y <- as.numeric(Y)

  if(is.factor(Y) | is.character(Y))
  {
    df = data.frame(column,Y)
    probs  = prop.table(table(df));
    cond_prob = apply(as.matrix(probs),1,median)
    v_gini=rep(NA,length(cond_prob))
    group = names(cond_prob[order(cond_prob)])

    for(v in 1:length(cond_prob)){
      CLASSE=column %in% group[1:v]
      v_gini[v]=gini_prob(y=Y,classe=CLASSE)
    }
  }

  else{
    index.elements <- !is.na(Y)
    Y <- Y[index.elements]
    column <- column[index.elements]

    #FIXME: MEAN BY MEDIAN
    cond_prob=aggregate(Y,by=list(column),mean,na.rm = TRUE)
    group=cond_prob[order(cond_prob$x),1]

    v_gini=rep(NA,nrow(cond_prob))

    for(v in 1:nrow(cond_prob)){
      CLASSE=column %in% group[1:v]
      v_gini[v]=gini_prob(y=Y,classe=CLASSE)
    }

  }

  as.character(sort(group[1:which.max(v_gini)]))

}

#' @title  General Interface Decision Tree model
#' @description Decision Tree is a simple and effective semi-supervised
#' learning method.
#' Based on the article "Semi-supervised classification trees".
#' It also offers many parameters to modify the behavior of this method.
#' It is the same as the traditional Decision Tree
#' algorithm, but the difference is how the gini coefficient is calculated (classification).
#' In regression we use SSE metric (different from the original investigation)
#' It can be used in classification or regression. If Y is numeric is for regression, classification in another case
#' @details In this model we can make predictions with prob type
#' @param max_depth A number from 1 to Inf.
#' Is the maximum number of depth in Decision Tree
#' Default is 30
#' @param w weight parameter ranging from 0 to 1. Default is 0.5
#' @param min_samples_split the minimum number of observations to do split. Default is 20
#' @param min_samples_leaf the minimum number of any terminal leaf node. Default is ceiling(min_samples_split/3)
#' @references
#' Jurica Levati, Michelangelo Ceci, Dragi Kocev, Saso Dzeroski.\cr
#' \emph{Semi-supervised classification trees.}\cr
#' Published online: 25 March 2017
#' © Springer Science Business Media New York 2017
#' @example demo/DecisionTree.R
#' @importFrom methods new
#' @export
SSLRDecisionTree <- function(max_depth = 30, w = 0.5,
                             min_samples_split = 20,
                             min_samples_leaf = ceiling(min_samples_split/3)){

  train_function <- function(x,y){

    tree <- new("DecisionTreeClassifier",max_depth = max_depth)
    x <- as.data.frame(x)
    tree <- fit_decision_tree(tree,x,y , w = w,min_samples_split = min_samples_split,
                              min_samples_leaf = min_samples_leaf)

    result <- list(
      model = tree
    )

    if(is.factor(y))
    result$classes = levels(y)

    if(is.factor(y) | is.character(y))
      result$mode = "classification"

    else
      result$mode = "regression"

    if(result$mode == "classification")
    result$pred.params = c("class","raw","prob")

    else
      result$pred.params = c("numeric","raw")


    class(result) <- "SSLRDecisionTree_fitted"

    return(result)
  }

  args <- list(
    max_depth = max_depth,
    w = w,
    min_samples_split = min_samples_split,
    min_samples_leaf = min_samples_leaf
  )

  new_model_sslr(train_function, "DecisionTreeClassifier" ,args)

}

setGeneric("predict", function(object,...)
  standardGeneric("predict") )

#' @title Function to predict inputs in Decision Tree
#' @description Function to predict inputs in Decision Tree
#' @param object The Decision Tree object
#' @param inputs data to be predicted
#' @param type Is param to define the type of predict.
#' It can be "class", to get class labels
#' Or "prob" to get probabilites for class in each input.
#' Default is "class"
#' @export
setMethod(f="predict",
          signature="DecisionTreeClassifier",
          definition=function(object,inputs, type = "class")
          {
            #If type class, we need to create a vector
            if(type == "class" | type == "numeric"){
              results <- c()
            }

            #If type class, we need to create a data frame
            else if(type == "prob"){
              results <- data.frame()
            }

            #We need to iterate all inputs
            for(i in 1:nrow(inputs)){

              #If type class, we need to push predicted class in input
              if(type == "class" | type == "numeric"){
                results <- c(results, predict_inputs(object, inputs[i,, drop = FALSE], type))
              }

              #If type prob, we need to push probabilities per class
              else if(type == "prob"){
                results <- rbind(results,predict_inputs(object, inputs[i,, drop = FALSE], type))
              }
            }

            #If type class, we return tibble with values of vector
            if(type == "class"){
              results <- factor(results)
              levels(results) <- object@classes

              results
            }


            if(type == "numeric"){
              results <- as.numeric(results)

              results
            }

            #If type class, we return tibble with colnames are
            #the names of classes
            else if(type == "prob"){
              colnames(results) <- object@classes
            }

            results
          }
)


#' @title Function to create DecisionTree
#' @param max_depth max depth in tree
#' @importFrom methods new
#' @export
newDecisionTree <- function(max_depth){
  new("DecisionTreeClassifier",max_depth=max_depth)
}



#' @title Predictions of the SSLRDecisionTree_fitted method
#' @description Predicts the label of instances SSLRDecisionTree_fitted model.
#' @param object model SSLRDecisionTree_fitted.
#' @param x A object that can be coerced as matrix.
#' Depending on how was the model built, \code{x} is interpreted as a matrix
#' with the distances between the unseen instances and the selected training instances,
#' or a matrix of instances.
#' @param ... This parameter is included for compatibility reasons.
#' @param type of predict in principal model
#' @method predict SSLRDecisionTree_fitted
#' @return Vector with the labels assigned.
#' @importFrom stats predict
#' @importFrom magrittr %>%
predict.SSLRDecisionTree_fitted <- function(object, x, type = "class", ...) {

  if(type == "raw" & object$mode == "classification")
    type <- "class"

  else if(type == "raw" & object$mode == "regression")
    type <- "numeric"


  #With the format of chosen model
  result <- object$model %>% predict(x, type = type)

  return(result)
}

