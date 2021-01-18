\name{Custom Methods}
\alias{Custom Methods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Using Custom Depth Functions and Classifiers
}
\description{
To use a custom depth function or classifier one has to implement three functions: parameters validator, learning and calculating functions.

}
%\usage{
%ddalpha.train(depth = "MyDepth", separator = "MyClassifier", ...)
%}
%- maybe also 'usage' for other objects documented here.
%\arguments{
%  \item{x}{
%%     ~~Describe \code{x} here~~
%}
%}
\details{

\bold{To define a depth function:}
\describe{
  \item{.<NAME>_validate}{
  validates parameters passed to \code{\link{ddalpha.train}} and passes them to the \code{ddalpha} object.
  
  \tabular{ll}{IN:\cr
    \code{ddalpha} \tab {the ddalpha object, containing the data and settings (mandatory)}\cr
    \code{<custom parameters>} \tab  {parameters that are passed to the user-defined method}\cr
    \code{...} \tab {other parameters (mandatory)}
    %}  \tabular{ll}{
    \cr
    OUT:\cr 
    \code{list()} \tab {list of output parameters, after the validation is finished, these parameters are stored in the \code{ddalpha} object}
    }
  }

  \item{.<NAME>_learn}{
  trains the depth
  
  \tabular{ll}{IN:\cr
    \code{ddalpha} \tab {the ddalpha object, containing the data and settings}
    %}  \tabular{ll}{
    \cr
  MODIFIES:\cr
  \code{ddalpha} \tab store the calculated statistics in the \code{ddalpha} object \cr
  depths \tab calculate the depths of each pattern, e.g. \cr \tab \code{for (i in 1:ddalpha$numPatterns) ddalpha$patterns[[i]]$depths = .<NAME>_depths(ddalpha, ddalpha$patterns[[i]]$points)}
  %}  \tabular{ll}{
  \cr
  OUT:\cr
    \code{ddalpha} \tab {the updated \code{ddalpha} object}
    }
  
  }
  \item{.<NAME>_depths}{
  calculates the depths
  
    \tabular{ll}{IN:\cr
    \code{ddalpha} \tab {the ddalpha object, containing the data and settings} \cr
    \code{objects} \tab {the objects for which the depths are calculated} \cr
    \cr
    OUT:\cr
      \code{depths} \tab {the calculated depths for each object (rows), with respect to each class (cols)}
      }
  }
}
Usage: \code{ddalpha.train(data, depth = "<NAME>", <custom parameters>, ...)}

\verb{
#### Custom depths ####

.MyDepth_validate  <- function(ddalpha, mydepth.parameter = "value", ...){ 
  print("MyDepth validating")
  # validate the parameters
  if (!is.valid(mydepth.parameter)){
    warning("Argument \"mydepth.parameter\" not specified correctly. Default value is used")
    mydepth.parameter = "value"
    # or stop("Argument \"mydepth.parameter\" not specified correctly.")
  }
  
  # the values from the return list will be stored in the ddalpha object
  return (list(mydepthpar = mydepth.parameter))
}

.MyDepth_learn <- function(ddalpha){
  print("MyDepth learning")
  
  #1. Calculate any statistics based on data that .MyDepth_depths needs
  #   and store them to the ddalpha object: 
  ddalpha$mydepth.statistic = "some value"
  
  #2. Calculate depths for each pattern
  for (i in 1:ddalpha$numPatterns){
    ddalpha$patterns[[i]]$depths = .MyDepth_depths(ddalpha, ddalpha$patterns[[i]]$points)
  }
  
  return(ddalpha)
}

.MyDepth_depths <- function(ddalpha, objects){
  print("MyDepth calculating")
  depths <- NULL
  
  # The depth parameters are accessible in the ddalpha object:
  mydepth.parameter = ddalpha$mydepth.parameter
  mydepth.statistic = ddalpha$mydepth.statistic
  
  #calculate the depths of the objects w.r.t. each pattern
  for (i in 1:ddalpha$numPatterns){
    depth_wrt_i = #calculate depths of the objects, as vector
    
    depths <- cbind(depths, depth_wrt_i)
  }
    
  return (depths)
}

ddalpha.train(data, depth = "MyDepth", ...)
}

\bold{To define a classifier:}
\describe{
  \item{.<NAME>_validate}{
  validates parameters passed to \code{\link{ddalpha.train}} and passes them to the \code{ddalpha} object
    \tabular{ll}{IN:\cr
    \code{ddalpha} \tab {the ddalpha object, containing the data and settings (mandatory)}\cr
    \code{<custom parameters>} \tab  {parameters that are passed to the user-defined method}\cr
    \code{...} \tab {other parameters (mandatory)}
    \cr
    OUT:\cr 
    \code{list()} \tab {list of output parameters, after the validation is finished, these parameters are stored in the \code{ddalpha} object. 
    
    In case of a multiclass classifier the validator must return \code{methodSeparatorBinary = F} and/or pass \code{aggregation.method = "none"} to \code{ddalpha.train}}
    }
  }
  \item{.<NAME>_learn}{
  trains the classifier. Is different for binnary and mylticlass classifiers.
    \tabular{ll}{IN:\cr
    \code{ddalpha} \tab {the ddalpha object, containing the data and settings}\cr
    \code{index1} \tab {(only for binary) index of the first class}\cr
    \code{index2} \tab {(only for binary) index of the second class}\cr
    \code{depths1} \tab {(only for binary) depths of the first class w.r.t. all classes}\cr
    \code{depths2} \tab {(only for binary) depths of the second class w.r.t. all classes}\cr\cr
    \tab depths w.r.t. only given classes are received by \code{depths1[,c(index1, index2)]}\cr\cr
    \tab for the multiclass classifiers the depths are accessible via \code{ddalpha$patterns[[i]]$depths}
    \cr
    OUT:\cr 
    \code{classifier} \tab { the trained \code{classifier} object}
    }
  }
  \item{.<NAME>_classify}{
  classifies the objects
   \tabular{ll}{IN:\cr
    \code{ddalpha} \tab {the ddalpha object, containing the data and global settings} \cr
    \code{classifier} \tab {the previously trained classifier}\cr
    \code{objects} \tab {the objects (depths) that are classified} \cr
    \cr
    OUT:\cr
      \code{result} \tab { a vector with classification results}\cr
        \tab {(binary) the objects from class \code{"classifier$index1"} get positive values} \cr
        \tab {(multiclass) the objects get the numbers of patterns in \code{ddalpha}}
    }
      
  }
}
Usage: \code{ddalpha.train(data, separator = "<NAME>", ...)}

\verb{
#### Custom classifiers ####


.MyClassifier_validate  <- function(ddalpha, my.parameter = "value", ...){
  print("MyClassifier validating")
  # validate the parameters
  ...
  
  # always include methodSeparatorBinary. 
  # TRUE for the binary classifier, FALSE otherwise
  return(list(methodSeparatorBinary = T, 
              my.parameter = my.parameter
              ))
}

# a binary classifier 
# the package takes care of the voting procedures. Just train it as if there are only two classes
.MyClassifier_learn <- function(ddalpha, index1, index2, depths1, depths2){
  print("MyClassifier (binary) learning")
  
  # The parameters are accessible in the ddalpha object:
  my.parameter = ddalpha$my.parameter
  
  #depths w.r.t. only given classes are received by
  depths1[,c(index1, index2)]
  depths2[,c(index1, index2)]
  
  # train the classifier
  classifier <-  ...
  
  #return the needed values in a list, e.g.
  return(list(
    coefficients  = classifier$coefficients,
    ... ))
}

# a multiclass classifier 
.MyClassifier_learn <- function(ddalpha){
  print("MyClassifier (multiclass) learning")
  
  # access the data through the ddalpha object, e.g.
  for (i in 1:ddalpha$numPatterns){
    depth  <- ddalpha$patterns[[i]]$depths
    number <- ddalpha$patterns[[i]]$cardinality
    ...
  }
  
  # train the classifier
  classifier <-  ...
  
  # return the classifier
  return(classifier)
}

# the interface of the classify function is equal for binary and multiclass classifiers
.MyClassifier_classify <- function(ddalpha, classifier, depths){
  print("MyClassifier classifying")
  
  # The global parameters are accessible in the ddalpha object:
  my.parameter = ddalpha$my.parameter
  # The parameters generated by .MyClassifier_learn are accessible in the classifier object:
  classifier$coefficients
  
  # here are the depths w.r.t. the first class
  depths[,classifier$index1]
  # here are the depths w.r.t. the second class
  depths[,classifier$index2]
  
  # fill results in a vector, so that:
  # (binary) the objects from class "classifier$index1" get positive values
  # (multiclass) the objects get the numbers of patterns in ddalpha
  result <- ...
  
  return(result)
}

ddalpha.train(data, separator = "MyClassifier", ...)
}
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ddalpha.train}}
}
\examples{

\dontrun{

#### example:  Euclidean depth ####

#.Euclidean_validate is basically not needed

.Euclidean_learn <- function(ddalpha){
  print("Euclidean depth learning")
  
  #1. Calculate any statistics based on data that .MyDepth_depths needs
  #   and store them to the ddalpha object: 
  for (i in 1:ddalpha$numPatterns){
    ddalpha$patterns[[i]]$center <- colMeans(ddalpha$patterns[[i]]$points)
  }
  
  #2. Calculate depths for each pattern
  for (i in 1:ddalpha$numPatterns){
    ddalpha$patterns[[i]]$depths = .Euclidian_depths(ddalpha, ddalpha$patterns[[i]]$points)
  }
  
  return(ddalpha)
}

.Euclidean_depths <- function(ddalpha, objects){
  print("Euclidean depth calculating")
  depths <- NULL
  
  #calculate the depths of the objects w.r.t. each pattern
  for (i in 1:ddalpha$numPatterns){
    # The depth parameters are accessible in the ddalpha object:
    center = ddalpha$patterns[[i]]$center
    
    depth_wrt_i <- 1/(1 + colSums((t(objects) - center)^2))
    depths <- cbind(depths, depth_wrt_i)
  }
  
  return (depths)
}

#### example:  binary decision tree ####

library(rpart)

.tree_validate  <- function(ddalpha, ...){
  print("tree validating")
  return(list(methodSeparatorBinary = T))
}

# a binary classifier 
# the package takes care of the voting procedures. Just train it as if there are only two classes
.tree_learn <- function(ddalpha, index1, index2, depths1, depths2){
  print("tree learning")
  
  # prepare the data
  data = as.data.frame(cbind( (rbind(depths1, depths2)), 
                              c(rep(1, times = nrow(depths1)), rep(-1, times = nrow(depths1)))))
  names(data) <- paste0("V",seq_len(ncol(data)))
  names(data)[ncol(data)] <- "O"
  # train the classifier
  classifier <- rpart(O~., data = data)
  
  #return the needed values in a list, e.g.
  return(classifier)
}

# the interface of the classify function is equal for binary and multiclass classifiers
.tree_classify <- function(ddalpha, classifier, depths){
  print("tree classifying")
  
  # fill results in a vector, so that the objects from class "classifier$index1" get positive values
  data = as.data.frame(depths)
  names(data) <- paste0("V",seq_len(ncol(data)))
  result <- predict(classifier, as.data.frame(depths), type = "vector")
  
  return(result)
}




#### checking ####

library(ddalpha)
data = getdata("hemophilia")

ddalpha = ddalpha.train(data = data, depth = "Euclidean", separator = "tree")
c = ddalpha.classify(ddalpha, data[,1:2])
errors = sum(unlist(c) != data[,3])/nrow(data)
print(paste("Error rate: ",errors))

# building the depth contours using the custom depth
depth.contours.ddalpha(ddalpha, asp = T, levels = seq(0.5, 1, length.out = 10))

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ depth }
\keyword{ classif }
\keyword{ custom }
