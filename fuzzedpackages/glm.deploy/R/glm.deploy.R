# glm.deploy:  Generates the source code of the scoring algorithm (predict) of a GLM object to C or JAVA to deploy/operationalize it outside R.
#
# Copyright (c) 2018, Oscar J. Castro-Lopez, Ines F. Vega-Lopez
#
# This file is part of the glm.deploy package for R.
#
# The glm.deploy: package is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The glm.deploy: package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details (http://www.gnu.org/licenses/).
######################################################################################

glm.deploy <- function(model, filename = NULL, language, path = NULL) {
  if (!inherits(model, "glm"))
    stop("ERROR: Not a glm object")

  allcoefs = dummy.coef(model)
  Fields = c()
  Coefficients = c()
  Datatypes = attributes(model$terms)$dataClasses
  Types = c()
  Arguments = attributes(model$terms)$dataClasses
  Factors = c()
  Factorsname = c()
  Intercept_label = paste0("", names(attributes(model$terms)$dataClasses[1]))
  i = 1
  j = 1
  hasfactor = FALSE
  for (x in 1:length(allcoefs)) {
    if (length(allcoefs[[x]]) == 1) {
      Fields[i] = names(allcoefs[x])
      Factors[i] = NA
      Factorsname[i] = NA
      Coefficients[i] = allcoefs[[x]]
      Types[i] = unname(Datatypes[grep(names(allcoefs[x]), names(Datatypes))])
      i = i + 1
      j = j + 1
    } else{
      ##If data is factor
      hasfactor = TRUE
      j = j + 1
      for (y in 1:length(allcoefs[[x]])) {
        Fields[i] = paste0(names(allcoefs[x]), names(allcoefs[[x]][y]))
        Factors[i] = names(allcoefs[[x]][y])
        Factorsname[i] = names(allcoefs[x])
        Coefficients[i] = allcoefs[[x]][y]
        Types[i] = unname(Datatypes[grep(names(allcoefs[x]), names(Datatypes))])
        i = i + 1
      }
    }
  }
  #Substitute invalid characters for variable names, and transform variable names to lower case
  names(Arguments) = tolower(gsub("\\.", "_", names(Arguments)))
  Fields = tolower(gsub("\\.", "_", Fields))
  Intercept_label = tolower(gsub("\\.", "_", Intercept_label))
  #UTF8 ENCODING
  Fields = enc2utf8(Fields)
  Intercept_label = enc2utf8(Intercept_label)
  Arguments = enc2utf8(Arguments)

  tbl = data.frame(
    Field = Fields,
    Coefficient = Coefficients,
    Type = Types,
    Factor = Factors,
    Factorname = Factorsname
  )
  if (is.null(filename)) {
    filename = paste0('glm_',Intercept_label)
  }

  glmdeploy_cpp(
    tbl,
    Arguments,
    model$family$family,
    model$family$link,
    Intercept_label,
    hasfactor,
    filename,
    language,
    path
  )
}

#' @name glm2c
#' @title C source code generator for rapid deployment of glm predictive models
##' @description The \code{glm2c()} function is used to generate source code in C language
##' implementing a given glm predictive model. It implements the following two functions;
##'  the glm_xxx_response() and glm_xxx_link(), where xxx stands for the name of the target variable
##'  of the glm object.\cr \cr
##' After the invocation of the \code{glm2c()} function two files are generated:\cr
##' \itemize{
#'   \item A .c file with the two scoring functions.\cr
#'   \item An .h file with the prototypes of the two functions of the .c file.\cr
#'  }
##'
##' @param model A fitted object of class "glm".
##' @param filename OPTIONAL The name of the output file(s), the default filenames are "glm_xxx.c" and "glm_xxx.h", where xxx is the target variable's name.
##' @param path The directory path where files are going to be saved.
##' @note All numeric variables used as input to the glm object are treated as doubles, whereas factor variables are treated as strings.
##' @seealso \code{\link{glm2java}}
##' @author Oscar Castro-Lopez, Ines Vega-Lopez
##' @examples
##'
#'  # Example with the iris dataset with a Logical target and numeric variables,
#'  # using the binomial family and the logit link function
#'  data(iris)
#'  iristest = iris
#'  iristest$Virginica = ifelse(iristest$Species == 'virginica', TRUE,FALSE)
#'  iristest$Species = NULL
#'
#'  # Load Package
#'  library(glm.deploy)
#'  # For repeatable results
#'  set.seed(123)
#'  # Generate the fitted glm object
#'  m = glm(Virginica ~ ., family = binomial(logit), data=iristest)
#'  # Call the glm2c() function with default filename
#'  glm2c(m,,tempdir())
#'
#'  # Call the glm2c() function with custom filename
#'  glm2c(m,'my_glm_virginica', tempdir())
#'
#'  # The glm2c() function generates the files: "glm_virginica.c" and
#'  # "glm_virginica.h"
#'
#'\dontrun{
#'---------------Contents of the "glm_virgninica.c" file---------------
#'
#' #include <stdlib.h>
#' #include <stdio.h>
#' #include <string.h>
#' #include <math.h>
#'
#' double glm_virginica_link(double sepal_length,
#'                           double sepal_width,
#'                           double petal_length,
#'                           double petal_width){
#'   double new_sepal_length = -2.46522019518341 * sepal_length;
#'   double new_sepal_width = -6.68088701405762 * sepal_width;
#'   double new_petal_length = 9.4293851538836 * petal_length;
#'   double new_petal_width = 18.2861368877881 * petal_width;
#'
#'   return -42.6378038127854+new_sepal_length+
#'                            new_sepal_width+
#'                            new_petal_length+
#'                            new_petal_width;
#' }
#' double glm_virginica_response(double sepal_length,
#'                               double sepal_width,
#'                               double petal_length,
#'                               double petal_width){
#'   return 1/(1+exp(-glm_virginica_link(sepal_length,
#'                                       sepal_width,
#'                                       petal_length,
#'                                       petal_width)));
#' }
#'----End of Contents of the "glm_virgninica.c" file------------------
#'--------------------------------------------------------------------
#'
#'-----Contents of the "glm_virgninica.h" file------------------------
#' double glm_virginica_link(double sepal_length,
#'                           double sepal_width,
#'                           double petal_length,
#'                           double petal_width);
#' double glm_virginica_response(double sepal_length,
#'                               double sepal_width,
#'                               double petal_length,
#'                               double petal_width);
#'-----End of Contents of the "glm_virgninica.h" file-----------------
#'--------------------------------------------------------------------
#'
##' Usage of the functions in another programs;
##' 1) We need to add an include line #include "virginica_glm.h" to all
##' source files that use library definitions.
##' 2) Link the .c file with the library object file.
##'     gcc -c glm_virginica.c
##' 3) The following is an example file "test.c" to call the functions
##' and print the result:
#'
#'-------------------"test.c"---------------
#' #include <stdio.h>
#' #include "glm_virgnica.h" //Added to call the scoring functions.
#'
#' int main(int argc, char *argv[]){
#'   printf("%f\n",glm_virginica_link(5.7,2.5,5.0,2.0));
#'   printf("%f\n",glm_virginica_response(5.7,2.5,5.0,2.0));
#'   return 0;
#' }
#'---------------End of "test.c"---------------
#'---------------------------------------------
#'
##' 4) Compile the "test.c" file and link it to the glm_virginica shared
##' library, we also need to add the "-lm" option to link it to the
##' math.h library:
##' gcc test.c -o test glm_virginica.o -lm
##'
##' 5) Finally Run the test.o program in linux:
##' ./test
##' }
glm2c <- function(model, filename = NULL, path = NULL) {
  if(is.null(path))
    stop("ERROR: A directory path must be provided")

  glm.deploy(model, filename, 0, path)
}

#' @name glm2java
#' @title Java source code generator for rapid deployment of glm predictive models
##' @description The \code{glm2java()} function is used to generate source code in Java language
##' implementing a given glm predictive model. It implements the following two methods;
##'  the glm_xxx_response() and glm_xxx_link(), where xxx stands for the name of the target variable
##'  of the glm object.\cr \cr
##' After invocation of the \code{glm2java()}, a .java file is generated containing the two predict methods which are declared as public static inside a java class called "glm_xxx_class".
##' @param model A fitted object of class "glm".
##' @param filename OPTIONAL The name of the output file, the default file name is  "glm_xxx_class.java", where xxx is the target variable's name.
##' @param path The directory path where files are going to be saved.
##' @note All numeric variables used as input to the glm object are treated as doubles, whereas factors variables are treated as strings.
##' @seealso \code{\link{glm2java}}
##' @author Oscar Castro-Lopez, Ines Vega-Lopez
##' @examples
#'  # Example with the iris dataset with a Logical target and numeric
#'  # variables, using the binomial family and the logit link function
#'  data(iris)
#'  iristest = iris
#'  iristest$Virginica = ifelse(iristest$Species == 'virginica', TRUE,FALSE)
#'  iristest$Species = NULL
#'
#'  # Load Package
#'  library(glm.deploy)
#'  # For repeatable results
#'  set.seed(123)
#'  # Generate the fitted glm object
#'  m = glm(Virginica ~ ., family = binomial(logit), data=iristest)
#'  # Call the glm2java() function with default filename
#'  glm2java(m,, tempdir())
#'  # Call the glm2java() function with custom filename
#'  glm2java(m,'my_glm_virginica', tempdir())
#'
#'  # The glm2java() function generates the file "glm_virginica_class.java".
#'
#'\dontrun{
#'----------Contents of the "glm_virgninica_class.java" file-------
#'   package test;
#'   public class glm_virginica_class{
#'
#'   public static double glm_virginica_link(double sepal_length,
#'                                           double sepal_width,
#'                                           double petal_length,
#'                                           double petal_width){
#'       double new_sepal_length = -2.46522019518341 * sepal_length;
#'       double new_sepal_width = -6.68088701405762 * sepal_width;
#'       double new_petal_length = 9.4293851538836 * petal_length;
#'       double new_petal_width = 18.2861368877881 * petal_width;
#'
#'       return -42.6378038127854+new_sepal_length+
#'                                new_sepal_width+
#'                                new_petal_length+
#'                                new_petal_width;
#'     }
#'     public static double glm_virginica_response(double sepal_length,
#'                                                 double sepal_width,
#'                                                 double petal_length,
#'                                                 double petal_width){
#'       return 1/(1+Math.exp(-glm_virginica_link(sepal_length,
#'                                                sepal_width,
#'                                                petal_length,
#'                                                petal_width)));
#'     }
#'
#'   }
#'---------------End of "glm_virgninica_class.java"---------------
#'----------------------------------------------------------------
#' To use these methods in another class just add
#' the "import glm_virginica_class.*;"
#' }
glm2java <- function(model, filename = NULL, path = NULL) {
  if(is.null(path))
    stop("ERROR: A directory path must be provided")
  glm.deploy(model, filename, 1, path)
}
