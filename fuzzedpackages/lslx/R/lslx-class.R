#' R6 class for semi-confirmatory structural equation modeling via penalized likelihood
#'
#' @docType class
#' @useDynLib lslx
#' @import stats
#' @import ggplot2
#' @importFrom Rcpp sourceCpp
#' @importFrom R6 R6Class
#' @keywords NULL
#' @format NULL
#' @usage NULL
#' @return Object of \code{lslx} R6 class for fitting semi-confirmatory structural equation modeling (SEM) with penalized likelihood (PL).
#'
#'
#'
#' @section Usage:
#' \code{lslx} is an \code{R6ClassGenerator} for constructing an \code{lslx} object that has methods for fitting semi-confirmatory SEM.
#' In a simpliest case, the use of \code{lslx} involves three major steps
#' \enumerate{
#' \item {Initialize a new \code{lslx} object by specifying a model and importing a data set.
#'
#' \code{r6_lslx <- lslx$new(model, data)}
#' }
#'
#' \item {
#' Fit the specified model to the imported data with given fitting control.
#'
#' \code{r6_lslx$fit(penalty_method, lambda_grid, delta_grid)}
#' }
#'
#' \item{
#' Summarize the fitting results with specified selector.
#'
#' \code{r6_lslx$summarize(selector)}
#' }
#' }
#' 
#' To cite \pkg{lslx} in publications use:
#' 
#' Po-Hsien Huang (in press). lslx: Semi-Confirmatory Structural Equation Modeling via Penalized Likelihood. Journal of Statistical Software.
#' 
#'
#'
#' @section Overview:
#' \pkg{lslx} is a package for fitting semi-confirmatory structural equation modeling (SEM) via penalized likelihood (PL) developed by Huang, Chen, and Weng (2017).
#' In this semi-confirmatory method, an SEM model is distinguished into two parts: a confirmatory part and an exploratory part.
#' The confirmatory part includes all of the freely estimated parameters and fixed parameters that are allowed for theory testing.
#' The exploratory part is composed by a set of penalized parameters describing relationships that cannot be clearly determined by available substantive theory.
#' By implementing a sparsity-inducing penalty and choosing an optimal penalty level, the relationships in the exploratory part can be efficiently identified by the sparsity pattern of these penalized parameters.
#' After Version 0.6.7, \pkg{lslx} also supports penalized least squares for SEM with ordianl data under delta parameterization. 
#' The technical details of \pkg{lslx} can be found in its JSS paper (Huang, 2020) <doi:10.18637/jss.v093.i07> or Vignette for Package lslx (https://cran.r-project.org/web/packages/lslx/vignettes/vignette-lslx.pdf).
#'
#' The main function \code{lslx} generates an object of \code{lslx} R6 class.
#' R6 class is established via package \pkg{R6} (Chang, 2017) that facilitates encapsulation object-oriented programming in \pkg{R} system.
#' Hence, the \code{lslx} object is self-contained.
#' On the one hand, \code{lslx} object stores model, data, and fitting results.
#' On the other hand, it has many built-in methods to respecify model, fit the model to data, and test goodness of fit and coefficients.
#' The initialization of a new \code{lslx} object requires importing a model and a data set to be analyzed.
#' After an \code{lslx} object is initialized, build-in methods can be used to modify the object, find the estimates, and summarize fitting result.
#' Details of object initialization is described in the section of \emph{Initialize Method}.
#'
#' In the current semi-confirmatory approach, the model specification is quite similar to the traditional practice of SEM except that some parameters can be set as penalized.
#' Model specification in \code{lslx} mainly relies on the argument model when creating a new \code{lslx} object.
#' After a \code{lslx} object is initialized, the initialized model can be still modified by set-related methods.
#' These set-related methods may hugely change the initialized model by just one simple command.
#' This two-step approach allows users specifying their own models flexibly and efficiently.
#' Details of the model specification can be found in the sections of \emph{Model Syntax} and \emph{Set-Related Methods}.
#'
#' Given a penalty level, \pkg{lslx} finds a PL estimate by minimizing a penalized maximum likelihood (ML) loss or a least squares loss functions (including OLS, DWLS, and WLS).
#' The penalty function can be set as lasso (Tibshirani, 1996), ridge (Hoerl & Kennard, 1970), elastic net (Zou & Hastie, 2005), or mcp (minimax concave penalty; Zhang, 2010).
#' \pkg{lslx} solves the optimization problem based on an improved \pkg{glmnet} method (Friedman, Hastie, & Tibshirani, 2010) made by Yuan, Ho, and Lin (2012).
#' The underlying optimizer is written by using \pkg{Rcpp} (Eddelbuettel & Francois, 2011) and \pkg{RcppEigen} (Bates & Eddelbuettel, 2013).
#' Our experiences show that the algorithm can efficiently find a local minimum provided that (1) the starting value is reasonable, and (2) the saturated covariance matrix is not nearly singular.
#' Details of optimization algorithm and how to implement the algorithm can be found in the sections of \emph{Optimization Algorithm} and \emph{Fit-Related Methods}.
#'
#' When conducting SEM, missing data are easily encountered.
#' \pkg{lslx} can handle missing data problem by listwise deletion and two-step methods.
#' Details of the methods for missing data can be found in the section of \emph{Missing Data}.
#'
#' After fitting the specified model to data under all of the considered penalty levels, an optimal penalty level should be chosen.
#' A naive method for penalty level selection is using information criteria.
#' Huang, Chen, and Weng (2017) have shown the asymptotic properties of Akaike information criterion (AIC) and Bayesian information criterion (BIC) in selecting the penalty level.
#' In \pkg{lslx}, information criteria other an AIC and BIC can be also used.
#' However, the empirical performances of these included criteria should be further studied.
#' Details of choosing an optimal penalty level can be found in the section of \emph{Penalty Level Selection}.
#'
#' Given a penalty level, it is important to evaluate the goodness-of-fit of selected model and coefficients.
#' In \pkg{lslx}, it is possible to make statistical inferences for goodness-of-fit and coefficients.
#' However, the inference methods assume that no model selection is conducted, which is not true in the case of using PL.
#' After version 0.6.4, several post-selection mehtods are available (Huang, 2019b).
#' Details of statistical inference can be found in the sections of \emph{Model Fit Evaluation} and \emph{Coefficient Evaluation}.
#' Implementations of these methods can be found in the sections of \emph{Summarize Method} and \emph{Test-Related Methods}.
#'
#' Besides making statistical inference, \code{lslx} has methods for plotting the fitting results, include visualizing quality of optimization and the values of information criteria, fit indices, and coefficient estimates.
#' Details of methods for plotting can be found in the section of \emph{Plot-Related Methods}.
#'
#' An object of \code{lslx} R6 class is composed by three R6 class objects: \code{lslxModel}, \code{lslxData}, and \code{lslxFitting}.
#' \code{lslxModel} contains the specified model and \code{lslxData} object stores the imported data to be analyzed.
#' When fitting the model to data, a reduced model and data will be sent to \code{lslxFitting}.
#' After the underlying optimizer finishes its job, the fitting results will be also stored in \code{lslxFitting}.
#' Since the three members are set as private, they can be only assessed by defined member functions.
#' Other than the three members, quantities that are crucial for SEM can be also extracted, such as model-implied moments, information matrix, and etc..
#' Details of methods for obtaining private members and SEM-related quantities can be found in the sections of \emph{Get-Related Methods} and \emph{Extract-Related Methods}.
#'
#'
#'
#' @section Model Syntax:
#' With \pkg{lslx} the relationships among observed variables and latent factors are mainly specified via equation-like syntax.
#' The creation of syntax in \pkg{lslx} is highly motivated by \pkg{lavaan} (Rosseel, 2012), a successful package for fitting SEM.
#' However, \pkg{lslx} utilizes slightly more complex, but still intuitive operators to describe relations among variables.
#'
#'
#'
#'
#' \bold{Example 1: Multiple Regression Model}
#'
#' Consider the first example of model that specifies a multiple regression model
#'
#'  \code{y <= x1 + x2}
#'
#' In this example, an dependent variable \code{y} is predicted by \code{x1} and \code{x2}.
#' These three variables are all observed variables and should appear in the given data set.
#' The operator \code{<=} means that the regression coefficients from the right-hand side (RHS) variables to the left-hand side (LHS) variables should be freely estimated.
#'
#' Although it is a very simple example, at least four important things behind this example should be addressed:
#' \enumerate{
#' \item{For any endogenous variable (i.e., an variable that is influenced by any other variable), its residual term is not required to be specified.
#' In this example, \code{lslx} recognizes \code{y} is an endogenous variable and hence the variance of corresponding residual will be set as freely estimated parameter.
#' It is also possible to explicitly specify the variance of residual of \code{y} by \code{y <=> y}.
#' Here, the operator \code{<=>} indicates the covariance of the RHS and LHS variables should be freely estimated.}
#' \item{If all of the specified equations do not contain the intercept variable \code{1}, then the intercept of each endogenous and observed variable will be freely estimated.
#' Since the intercept variable \code{1} doesn't appear in this example, the intercept for \code{y} will be set as a freely estimated parameter.
#' We can explicitly set the intercept term by \code{y <= 1}.
#' However, under the situation that many equations are specified, once the intercept variable \code{1} appears in some equation, intercept terms in other equations should be explicitly specified.
#' Otherwise, the intercepts of endogenous and observed variables in other equations will be fixed at zero.
#' }
#' \item{For any set of exogeneous variables (i.e., variables that are not influenced by any other variable in the specified system),
#' not only their variances will be freely estimated, but also their pairwise covariances will be set as freely estimated parameters.
#' In this example, \code{x1} and \code{x2} are both exogeneous variables.
#' Hence, their variance and pairwise covariances will be automatically set as freely estimated parameters.
#' These covariances can be explicitly stated by simply \code{x1 + x2 <=> x1 + x2}.
#' The syntax parser in \code{lslx} will consider variance/covariance of each combination of LHS and RHS variables.}
#' \item{The intercepts (or means) of exogeneous and observed variables are always set as freely estimated parameters.
#' In this example, the intercepts of \code{x1} and \code{x2} will be freely estimated.
#' It can be stated explicitly by \code{x1 + x2 <= 1}.
#' Also, the \code{lslx} parser will know that the intercept variable \code{1} has effect on all of \code{x1} and \code{x2}.}
#' }
#' The previous regression example can be equivalently represented by \code{x1 + x2 => y}.
#' In \pkg{lslx} all of the directed operators can be reversed under the stage of model specification.
#' Users can choose the directions of operators according to their own preference.
#'
#' The unique feature of \pkg{lslx} is that all of the parameters can be set as penalized.
#' To penalize all of the regression coefficients, the equation can be modified as
#'
#'  \code{y <~ x1 + x2}
#'
#' Here, the operator \code{<~} means that the regression coefficients from the RHS variables to the LHS variables should be estimated with penalization.
#' If only the coefficient of \code{x1} should be penalized, we can use prefix to partly modify the equation
#'
#'  \code{y <= pen() * x1 + x2}
#'
#' or equivalently
#'
#'  \code{y <~ x1 + free() * x2}
#'
#' Both \code{pen()} and \code{free()} are prefix to modify the \pkg{lslx} operators.
#' \code{pen()} makes the corresponding parameter to be penalized and \code{free()} makes it to be freely estimated.
#' Inside the parentheses, starting values can be specified.
#' Any prefix must present before some variable name and divided by asterisk \code{*}.
#' Note that prefix can appear in the either RHS or LHS of operators and its function can be 'distributed' to the variables in the other side.
#' For example, \code{free() * y <~ x1 + x2} will be interpreted as that all of the coefficients should be freely estimated.
#' However, any prefix cannot simultaneously appear on both sides of operators, which may result in an ambiguity specification.
#'
#'
#'
#'
#' \bold{Example 2: Factor Analysis Model}
#'
#' Now, we consider another example of equation specification.
#'
#'
#'  \code{y1 + y2 + y3 <=: f1}
#'
#'  \code{y4 + y5 + y6 <=: f2}
#'
#'  \code{y7 + y8 + y9 <=: f3}
#'
#' This example is a factor analysis model with nine observed variables and three latent factors.
#' In \pkg{lslx}, defining a latent factor can be through the operator \code{<=:} which means that the RHS factor is defined by LHS observed variables.
#' The observed variables must be presented in the given data set.
#' Of course, \code{f1} can be equivalently defined by \code{f1 :=> y1 + y2 + y3}.
#'
#' As addressed in the first example, the \code{lslx} parser will automatically set many parameters that are not directly presented in these equations.
#' \enumerate{
#' \item {
#' In this example, all of the observed variables are directed by some latent factor and hence they are endogenous.
#' The variances of their residuals will be set as freely estimated parameters.
#' Also, their intercepts will be freely estimated since no intercept variable \code{1} presents in the specified equations.
#' }
#' \item {
#' The three latent factors \code{f1}, \code{f2}, and \code{f3} are exogenous variables.
#' Their pairwise covariances will be also set as freely estimated parameters.
#' However, because they are latent but not observed, their intercepts will be fixed at zero.
#' If user hope to estimate the latent factor means, they should add an additional equation \code{f1 + f2 + f3 <= 1}.
#' After adding this equation, on the one hand, the latent intercepts will be set as free as indicated by that equation.
#' On the other hand, since intercept variable \code{1} now presents in the specified equations, the intercepts for the endogenous and observed variables will be then fixed at zero.
#' }
#' }
#' So far, the specification for the factor analysis model is not complete since the scales of factors are not yet determined.
#' In SEM, there are two common ways for scale setting.
#' The first way is to fix some loading per factor.
#' For example, we may respecify the model via
#'
#'  \code{fix(1) * y1 + y2 + y3 <=: f1}
#'
#'  \code{fix(1) * y4 + y5 + y6 <=: f2}
#'
#'  \code{fix(1) * y7 + y8 + y9 <=: f3}
#'
#' The prefix \code{fix(1)} will fix the corresponding loadings to be one.
#' Simply using \code{1 * y1} will be also interpreted as fixing the loading of \code{y1} at one to mimic \pkg{lavaan}.
#' The second way for scale setting is fixing the variance of latent factors, which can be achieved by specifying additional equations
#'
#'  \code{fix(1) * f1 <=> f1}
#'
#'  \code{fix(1) * f2 <=> f2}
#'
#'  \code{fix(1) * f3 <=> f3}
#'
#' Note that in the current version of \pkg{lslx}, scale setting will be not made automatically.
#' Users must accomplish it manually.
#'
#' When conducting factor analysis, we may face the problem that each variable may not be influenced by only one latent factor.
#' The semi-confirmatory factor analysis, which penalizes some part of loading matrix, can be applied in this situation.
#' One possible model specification for the semi-confirmatory approach is
#'
#'  \code{y1 + y2 + y3 <=: f1}
#'
#'  \code{y4 + y5 + y6 <=: f2}
#'
#'  \code{y7 + y8 + y9 <=: f3}
#'
#'  \code{y4 + y5 + y6 + y7 + y8 + y9 <~: f1}
#'
#'  \code{y1 + y2 + y3 + y7 + y8 + y9 <~: f2}
#'
#'  \code{y4 + y5 + y6 + y7 + y8 + y9 <~: f3}
#'
#'  \code{fix(1) * f1 <=> f1}
#'
#'  \code{fix(1) * f2 <=> f2}
#'
#'  \code{fix(1) * f3 <=> f3}
#'
#' In this specification, loadings in the non-independent cluster will be also estimated but with penalization.
#' 
#' After version 0.6.3, \pkg{lslx} supports basic \pkg{lavaan} operators. 
#' The previous model can be equivalently specified as
#'
#'  \code{f1 =~ y1 + y2 + y3}
#'
#'  \code{f2 =~ y4 + y5 + y6}
#'
#'  \code{f3 =~ y7 + y8 + y9}
#'
#'  \code{pen() * f1 =~ y4 + y5 + y6 + y7 + y8 + y9}
#'
#'  \code{pen() * f2 =~ y1 + y2 + y3 + y7 + y8 + y9}
#'
#'  \code{pen() * f3 =~ y4 + y5 + y6 + y7 + y8 + y9}
#'
#'  \code{f1 ~~ 1 * f1}
#'
#'  \code{f2 ~~ 1 * f2}
#'
#'  \code{f3 ~~ 1 * f3}
#'
#'
#' \bold{Example 3: Path Models with both Observed Variables and Latent Factors}
#'
#' In the third example, we consider a path model with both observed variables and latent factors
#'
#'  \code{fix(1) * y1 + y2 + y3 <=: f1}
#'
#'  \code{fix(1) * y4 + y5 + y6 <=: f2}
#'
#'  \code{fix(1) * y7 + y8 + y9 <=: f3}
#'
#'  \code{f3 <= f1 + f2}
#'
#'  \code{f1 + f2 + f3 <~ x1 + x2}
#'
#'  \code{f1 <~> f2}
#'
#' The first three equations specify the measurement model for \code{y1} - \code{y9} and \code{f1} - \code{f3}.
#' The forth equation describes the relations among latent factor.
#' The fifth equation sets all the coefficients from \code{x1} - \code{x2} to \code{f1} - \code{f3} to be penalized.
#' The final equation states that the covariance of residuals of \code{f1} and \code{f2} is estimated with penalization, which is achieved by the operator \code{<~>}.
#'
#' Like Example 1 and 2, many parameters in the current example are automatically set by \pkg{lslx}.
#' \enumerate{
#' \item {
#' Because \code{y1} - \code{y9} and \code{f1} - \code{f3} are all endogenous, the variances of their residuals will be treated as freely estimated parameters.
#' Also, due to the non-presence of intercept variable \code{1}, the intercept of \code{y1} - \code{y9} will be set as free parameters
#' and the intercept of \code{f1} - \code{f3} will be set as zero.
#' }
#' \item {
#' The variance, intercepts, and pairwise covariances of exogenous and observed variables \code{x1} - \code{x2} will be all estimated freely.
#' }
#' }
#'
#' In this example, we can see that model specification in \pkg{lslx} is quite flexible.
#' Like usual SEM, users can specify their models according some substantive theory.
#' If no theory is available to guide the relationships in some part of the model,
#' the semi-confirmatory approach can set this part as exploratory by setting the corresponding parameters as penalized.
#'
#'
#'
#'
#' \bold{Example 4: Multi-Group Factor Analysis Model}
#'
#' In the fourth example, we consider a multi-group factor analysis model.
#'
#'  \code{fix(1) * y1 + y2 + y3 <=: f1}
#'
#'  \code{fix(1) * y4 + y5 + y6 <=: f2}
#'
#'  \code{fix(1) * y7 + y8 + y9 <=: f3}
#'
#' The syntax specifies a factor analysis model with nine observed variables and three latent factors.
#' Loadings for \code{y2}, \code{y3}, \code{y5}, \code{y6}, \code{y8}, and \code{y9} are freely estimated in both groups.
#' Loadings for \code{y1}, \code{y4}, and \code{y7} are set as fixed for scale setting in both groups.
#' You may observe that the syntax for multi-group analysis is the same as that for single group analysis.
#' That is true because in \pkg{lslx} a multi-group analysis is mainly identified by specifying a group variable.
#' If the imported data can be divided into several samples based on some group variable (argument \code{group_variable} in \code{new} method, please see the section of \emph{Initialize Method}) for group labeling, \pkg{lslx} will automatically conduct multi-group analysis
#' (see example of \emph{Semi-Confirmatory Multi-Group Factor Analysis} in the Section of \emph{Examples}).
#'
#' Sometimes, we may hope to specify different model structures for the two groups.
#' It can be achieved by using vector version of prefix, which is also motivated by the syntax in \pkg{lavaan}.
#' For example, if we hope to restrict the loading for \code{y2} to be 1 in the first group but set it as freely estimate parameter in the second group.
#' Then we may use
#'
#' \code{fix(1) * y1 + c(fix(1), free()) * y2 + y3 <=: f1}
#'
#' Note that the order of groups is important here.
#' Since \pkg{lslx} treats the group variable as \code{factor}, its order is determined by the sorted name of groups.
#' For example, if three groups \code{c}, \code{a}, and \code{b} are considered, then the first group is \code{a}, the second is \code{b}, and the third is \code{c}.
#'
#'
#' In the current version of \pkg{lslx}, coefficient constraints cannot be imposed.
#' It seems that testing coefficient invariance across groups is impossible in \pkg{lslx}.
#' However, the present package parameterizes group coefficients in different way compared to other SEM software (Huang, 2018).
#' Under \pkg{lslx}, each group coefficient is decomposed into a sum of a reference component and an increment component.
#' If the reference component is assumed to be zero, the increment component represents the group coefficient, which is equivalent to the usual parameterization in other software solutions.
#' On the other hand, if some group is set as reference (argument \code{reference_group} in \code{new} method, please see the section of \emph{Initialize Method}),
#' then the reference component now represents the group coefficient of the reference group and other increment components represent the differences from the reference group.
#' The coefficient invariance across groups can be evaluated by examining the value or sparsity of the corresponding increment component.
#'
#'@section Optimization Algorithm:
#' Let \eqn{\theta} denote the vector of model parameter.
#' \pkg{lslx} tries to find a PL estimate for \eqn{\theta} by minimizing the objective function
#' \eqn{objective(\theta, \lambda) = loss(\theta) + regularizer(\theta, \lambda)}
#' where \eqn{loss} is the ML loss function, \eqn{regularizer} is a regularizer, possibly lasso (Tibshirani, 1996) or mcp (Zhang, 2010), and \eqn{\lambda} is a regularization parameter.
#' The optimization algorithm for minimizing the PL criterion is based on an improved \pkg{glmnet} method (Friedman, Hastie, & Tibshirani, 2010) made by Yuan, Ho, and Lin (2012).
#' The algorithm can be understood as a quasi-Newton method with inner loop and outer loop.
#' The inner loop of the algorithm derives a quasi-Newton direction by minimizing a quadratic approximated objective function via coordinate descent.
#' To save the computation time, the Hessian matrix for the quadratic term is approximated by the identity matrix, the Broyden-Fletcher-Goldfarb-Shanno (BFGS) method, or the expected Hessian (Fisher scoring).
#' Although the computational cost of BFGS approximation is much smaller than calculating expected hessian,
#' our experience shows that the two methods perform similarly in terms of computation time because more outer iterations are required for BFGS.
#' The inner loop stops if the change of the derived direction is quite small.
#' The outer loop of the algorithm updates the value of parameter estimate via the derived quasi-Newton direction and Armijo's rule.
#' The outer loop stops if the maximal absolute element of subgradient of objective function is smaller than the specified tolerance.
#' The minimizer is the so-called PL estimates.
#' Note that the PL estimates is a function of penalty level, i.e., PL estimates can vary under different penalty levels.
#' An optimal penalty level can be chosen by using model selection criterion.
#'
#' In \pkg{lslx}, PL estimates under each penalty level and convexity level specified by user will be calculated.
#' The convexity levels will be sorted from large to small based on the suggestion of Mazumder (2011).
#' The previous obtained PL estimate will be used as warm start for further minimization problem.
#' Since the solution path is continuous, the warm start can speed up the convergence of minimization (see Friedman, Hastie, & Tibshirani, 2010).
#'
#'
#' @section Missing Data:
#' When conducting SEM, it is easy to encounter the problem of missing data.
#' In \pkg{lslx}, missing data can be handled by the listwise deletion method and the two-stage method (Yuan & Bentler, 2000).
#' The listwise deletion method only uses fully complete observations for further SEM analysis.
#' If the missing mechanism is missing completely at random (MCAR; Rubin, 1976), the listwise deletion method can yield a consistent estimator.
#' The two-stage method first calculates the saturated moments by minimizing the likelihoods based on all of the available observations and then use the obtained saturated moment estimates for further SEM analysis.
#' Under the assumption of missing at random (MAR; Rubin, 1976), it has been shown that the two-stage method can yield a consistent estimate.
#' In addition, the standard errors of coefficients can be also consistently estimated if a correct asymptotic covariance of saturated moments is used.
#' Because the two-stage approach is generally valid and efficient compared to the listwise deletion method,
#' \pkg{lslx} set the two-stage method as default for handling the missing data problem.
#' The current version also supports the use of auxiliary variables (see Savalei & Bentler, 2008).
#' If the two-stage method is implemented, the standard error formula will be corrected for the presence of missing data (see Yuan & Lu, 2008 for technical details).
#'
#' So far, \pkg{lslx} doesn't include the full-information maximum likelihood (FIML) method for missing values.
#' One reason is that PL can be computationally intensive if many penalty levels are considered.
#' The additional E-step in each iteration of FIML makes the problem worse.
#' Another reason is that the two-step method has been shown to outperform FIML in simulation settings (Savalei & Falk, 2014).
#' Therefore, we tend to believe that the implementation of FIML in PL may not bring further advantages over the two-step method.
#'
#'
#'@section Penalty Level Selection:
#' Penalty level selection in \pkg{lslx} is based on optimizing the value of some information criterion.
#' Many information criteria are available for this task.
#' In the current version, available information criteria are
#'
#' \describe{
#'   \item{\code{aic}}{Akaike Information Criterion (Akaike, 1974)
#'   \deqn{AIC(\theta)=loss(\theta) - (2 / N) * df(\theta) } }
#'   \item{\code{aic3}}{Akaike Information Criterion with Penalty Being 3 (Sclove, 1987)
#'   \deqn{AIC3(\theta)=loss(\theta) - (3 / N) * df(\theta) }}
#'   \item{\code{caic}}{Consistent Akaike Information Criterion (Bozdogan, 1987)
#'   \deqn{CAIC(\theta)=loss(\theta) - ((log(N) + 1) / N) * df(\theta) }}
#'   \item{\code{bic}}{Bayesian Information Criterion (Schwarz, 1978)
#'   \deqn{BIC(\theta)=loss(\theta) - (log(N) / N) * df(\theta) } }
#'   \item{\code{abic}}{Adjusted Bayesian Information Criterion (Sclove, 1987)
#'   \deqn{ABIC(\theta)=loss(\theta) - (log((N + 2) / 24 ) / N) * df(\theta) }}
#'   \item{\code{hbic}}{Haughton Bayesian Information Criterion (Haughton, 1997)
#'   \deqn{HBIC(\theta)=loss(\theta) - (log(N / \pi ) / N) * df(\theta) }}
#' }
#'
#' where
#' \itemize{
#' \item{
#' \eqn{N}: total number of sample size;
#' }
#' \item{
#' \eqn{G}: total number of group;
#' }
#' \item{
#' \eqn{loss(\theta)}: the loss value under estimate \eqn{\theta};
#' }
#' \item{
#' \eqn{df(\theta)}: the degree of freedom defined as (1) \eqn{G * P * (P + 3) / 2 - e(\theta)} with \eqn{e(\theta)} being the number of non-zero elements in \eqn{\theta} for Lasso and MCP; 
#' or (2) the expectation of likelihood ratio statistics with ridge for ridge and elastic net.
#' }
#' }
#' Note the formula for calculating the information criteria in \pkg{lslx} are different to other software solutions.
#' The loss function value is used to replace the likelihood function value and hence the penalty term is also divided by sample size \eqn{N}.
#' For each information criterion, a robust version is calculated if raw data is available.
#' Their corresponding names are \code{raic}, \code{raic3}, \code{rcaic}, \code{rbic}, \code{rabic}, and \code{rhbic} with "r" standing for "robust".
#' These robust criteria use the Satorra-Bentler scaling factor for correcting degree of freedom.
#' For the case of normal data and correctly specified model, the two versions will be the same asymptotically.
#'
#'
#' Huang, Chen, and Weng (2017) have study the asymptotic behaviors of \code{aic} and \code{bic} under penalized estimation.
#' They show that under suitable conditions, \code{aic} can select a model with minimum expected loss and \code{bic} can choose the most parsimonious one from models that attain the minimum expected loss.
#' By the order of penalty term, we may expect: (1) the large sample behaviors of \code{aic3} and \code{tic} will be similar to \code{aic};
#' and (2) the asymptotic behaviors of \code{caic}, \code{abic}, and \code{hbic} will be similar to \code{bic}.
#' However, their small-sample performances require further studies.
#'
#'
#'@section Model Fit Evaluation:
#' Given a chosen penalty level, we may evaluate the overall model fit by using fit indices.
#' In the current version, available fit indices for model evaluation are
#' \describe{
#' \item{\code{rmsea}}{Root Mean Square Error of Approximation (Steiger, 1998; Steiger & Lind, 1980)
#' \deqn{ RMSEA(\theta)=\sqrt(G * max(loss(\theta) / df(\theta) - 1 / N, 0)) }}
#' \item{\code{cfi}}{Comparative Fit Index (Bentler, 1990)
#' \deqn{ CFI(\theta)=(max(loss_0 - df_0 / N, 0) - max(loss(\theta) - df(\theta) / N, 0)) / max(loss_0 - df_0 / N, 0)}}
#' \item{\code{nnfi}}{Non-Normed Fit Index (Tucker & Lewis, 1973)
#' \deqn{ NNFI(\theta)=(loss_0 / df_0 - loss(\theta) / df(\theta)) / (loss_0 / df_0 - 1 /N) }}
#' \item{\code{srmr}}{Standardized Root Mean of Residual (Bentler, 1995)
#' \deqn{ SRMR(\theta)=\sqrt(\sum_g w_g \sum_i \sum_{j \leq i} ((\sigma_{gij} - s_{gij})^2 / (\sigma_{gii} * \sigma_{gjj})) / (G * P * (P + 1) / 2)}
#' \deqn{+\sum_g w_g \sum_i ((\mu_{gi} - m_{gi}) ^ 2 / \sigma_{gii}) / (G * P)) }}
#' }
#' where
#' \itemize{
#' \item{
#' \eqn{N}: total number of sample size;
#' }
#' \item{
#' \eqn{G}: total number of groups;
#' }
#' \item{
#' \eqn{P}: number of observed variables;
#' }
#' \item{
#' \eqn{w_g}: sample weight of group \eqn{g};
#' }
#' \item{
#' \eqn{loss(\theta)}: the loss value under estimate \eqn{\theta};
#' }
#' \item{
#' \eqn{\#(\theta)}: the number of non-zero elements in \eqn{\theta};
#' }
#' \item{
#' \eqn{df(\theta)}: the degree of freedom defined by \eqn{G * P * (P + 3) / 2 -\#(\theta)};
#' }
#' \item{
#' \eqn{loss_0}: the loss value under baseline model;
#' }
#' \item{
#' \eqn{df_0}: the degree of freedom under baseline model;
#' }
#' \item{
#' \eqn{\sigma_{gij}}: the \eqn{(i,j)} element of model implied covariance at group \eqn{g};
#' }
#' \item{
#' \eqn{s_{gij}}: the \eqn{(i,j)} element of sample covariance at group \eqn{g};
#' }
#' \item{
#' \eqn{\mu_{gi}}: the \eqn{i} element of model implied mean at group \eqn{g};
#' }
#' \item{
#' \eqn{m_{gi}}: the \eqn{i} element of sample mean at group \eqn{g};
#' }
#' }
#' In \pkg{lslx}, the baseline model is the model that assumes a diagonal covariance matrix and a saturated mean.
#' Hence, the baseline model may not be appropriate if users hope to evaluate the goodness-of-fit of mean structure.
#'
#'
#' It is also possible to test overall model fit by formal statistical test.
#' In the current version, statistical tests for likelihood ratio (LR) and root mean square error of approximation (RMSEA) can be implemented.
#' If raw data is available, \pkg{lslx} calculates mean-adjusted versions of LR statistic (Satorra & Bentler, 1994) and RMSEA intervals (Brosseau-Liard, Savalei & Li, 2012; Li & Bentler, 2006).
#' It should be noted that the classical tests may not be valid after penalty level selection because the task of penalty level selection may destroy the sampling distribution of test statistics (see Pötscher, 1991 for discussion).
#' Valid post model selection inference methods require further development.
#'
#'
#' @section Coefficient Evaluation:
#' Given a chosen penalty level, we may evaluate the significance of coefficients (or parameters).
#' In the current version, standard errors based on the expected/observed Fisher information matrix and the sandwich formula are available (see Yuan & Hayashi, 2006 for discussion).
#' Because the sandwich formula is generally valid compared to the approaches based on Fisher information, \pkg{lslx} uses the sandwich formula as default whenever raw data is available.
#' Note that sandwich covariance matrix in \pkg{lslx} is calculated based on Equation (14) in Yuan and Hayashi (2006) but not Equation (2.12a) in Browne (1984) to accommodate the potential model misspecification.
#' Again, the significance tests may not be valid after penalty level selection.
#' After version 0.6.4, several post-selection mehtods are available (Huang, in press), inlcuding PoSI method with Scheffe constant (Berk, Brown, Buja, Zhang, & Zhao, 2013) and Polyhedral method (Lee, Sun, Sun, & Taylor, 2016).
#'
#'
#' @section Initialize Method:
#'
#' \preformatted{$new(model, data, numeric_variable, ordered_variable, 
#'   group_variable, reference_group, weight_variable, auxiliary_variable, 
#'   sample_cov, sample_mean, sample_size, sample_moment_acov, verbose = TRUE)}
#' \describe{
#' \item{\bold{Arguments}}{
#'
#' }
#' \item{\code{model}}{A \code{character} with length one to represent the model specification.}
#' \item{\code{data}}{A \code{data.frame} of raw data.
#' It must contains variables specified in \code{model} (and possibly the variables specified by \code{group_variable} and \code{weight_variable}).}
#' \item{\code{numeric_variable}}{A \code{character} to specify which response variables should be transfromed into \code{numeric}.}
#' \item{\code{ordered_variable}}{A \code{character} to specify which response variables should be transfromed into \code{ordered}.}
#' \item{\code{weight_variable}}{A \code{character} with length one to specify what variable is used for sampling weight.}
#' \item{\code{auxiliary_variable}}{A \code{character} to specify what variable(s) is used as auxiliary variable(s) for estimating saturated moments when missing data presents and two-step method is implemented.
#' Auxiliary variable(s) must be numeric. If any categorical auxiliary is considered, please transform it into dummy variables before initialization.}
#' \item{\code{group_variable}}{A \code{character} with length one to specify what variable is used for labeling group.}
#' \item{\code{reference_group}}{A \code{character} with length one to specify which group is set as reference. }
#' \item{\code{sample_cov}}{A numeric \code{matrix} (single group case) or a \code{list} of numeric \code{matrix} (multi-group case) to represent sample covariance matrixs. It must have row and column names that match the variable names specified in \code{model}.}
#' \item{\code{sample_mean}}{A \code{numeric} (single group case) or a \code{list} of \code{numeric} (multi-group case) to represent sample mean vectors.}
#' \item{\code{sample_size}}{A \code{numeric} (single group case) with length one or a \code{list} of \code{numeric} (multi-group case) to represent the sample sizes.}
#' \item{\code{sample_moment_acov}}{A numeric \code{matrix} (single group case) or a \code{list} of numeric \code{matrix} (multi-group case) to represent asymptotic covariance for moments.}
#' \item{\code{verbose}}{A \code{logical} to specify whether messages made by \code{lslx} should be printed.}
#' }
#' \bold{Details}
#'
#' \code{$new()} initializes a new object of \code{lslx} R6 class for fitting semi-confirmatory structural equation modeling (SEM).
#' In most cases, a new \code{lslx} object is initialized by supplying \code{model} and \code{data}.
#' For details of syntax for model specification, see the section of \emph{Model Syntax}.
#' By default, types of response variables (\code{numeric} versus \code{ordered}) will be infered according to the imported \code{data.frame}.
#' If users hope to transform variables types in \code{lslx}, \code{numeric_variable} and \code{ordered_variable} can be used. 
#' When multi-group analysis is desired, argument \code{group_variable} should be given to specify what variable is used for labeling group.
#' Argument \code{reference_group} can be used to set reference group.
#' Note that if some group is set as reference, the coefficients in other groups will represent increments from the reference.
#' When the missingness of data depends on some other variables, \code{auxiliary_variable} can be used to specify auxiliary variables for estimate saturated moments.
#' For details of missing data, see the section of \emph{Missing Data}.
#' If raw data is not available, \code{lslx} also supports initialization via sample moments.
#' In that case, \code{sample_cov} and \code{sample_size} are required.
#' If \code{sample_mean} is missing under moment initialization, it is assumed to be zero.
#'
#'
#' @section Set-Related Methods:
#' \preformatted{$free_coefficient(name, start, verbose = TRUE)
#' $penalize_coefficient(name, start, verbose = TRUE)
#' $fix_coefficient(name, start, verbose = TRUE)
#'
#' $free_directed(left, right, group, verbose = TRUE)
#' $penalize_directed(left, right, group, verbose = TRUE)
#' $fix_directed(left, right, group, verbose = TRUE)
#'
#' $free_undirected(both, group, verbose = TRUE)
#' $penalize_undirected(both, group, verbose = TRUE)
#' $fix_undirected(both, group, verbose = TRUE)
#'
#' $free_block(block, group, type, verbose = TRUE)
#' $penalize_block(block, group, type, verbose = TRUE)
#' $fix_block(block, group, type, verbose = TRUE)
#'
#' $free_heterogeneity(block, group, verbose = TRUE)
#' $penalize_heterogeneity(block, group, verbose = TRUE)
#' $fix_heterogeneity(block, group, verbose = TRUE)}
#' \describe{
#' \item{\bold{Arguments}}{
#'
#' }
#' \item{\code{name}}{A \code{character} to indicate which coefficients should be reset.}
#' \item{\code{start}}{A \code{numeric} to specify starting values.
#' The length of \code{start} should be one or match the length of \code{name} to avoid ambiguity.
#' If \code{start} is missing, the starting value will be set as
#' (1) \code{NA} for free or penalized coefficient; and (2) \code{0} for fixed coefficient.}
#' \item{\code{left}}{A \code{character} to indicate variable names in the left-hand side of operator \code{"<-"}.}
#' \item{\code{right}}{A \code{character} to indicate variable names in the right-hand side of operator \code{"<-"}.}
#' \item{\code{both}}{A \code{character} to indicate variable names in both side of operator \code{"<->"}.}
#' \item{\code{group}}{A \code{character} to indicate group names that the specified relations belong to.}
#' \item{\code{block}}{A \code{character} with length one to indicate a block such that the corresponding target coefficient will be reset.
#' Its value must be \code{"f<-1"}, \code{"y<-1"}, \code{"f<-f"}, \code{"f<-y"}, \code{"y<-f"}, \code{"y<-y"}, \code{"f<->f"}, \code{"f<->y"}, \code{"y<->f"}, or \code{"y<->y"}.}
#' \item{\code{type}}{A \code{character} to indicate which type of parameter should be changed in the given block.
#' Its value must be \code{"free"}, \code{"fixed"}, or \code{"pen"}. If \code{type} is not specified, the types of all parameters in the given \code{block} will be modified.}
#' \item{\code{verbose}}{A \code{logical} to specify whether messages made by \code{lslx} should be printed.}
#' }
#' \bold{Details}
#'
#' Set-related methods include several member functions that can be used to modify the initialized model specification.
#' Like most encapsulation objects, set-related function is used to modify the inner members of object.
#' So far, the set-related methods are all established to modify the model.
#' The data are protected without any modification.
#'
#' \code{$free_coefficient()} / \code{$penalize_coefficient()} / \code{$fix_coefficient()} sets the coefficient named \code{name} as FREE / PENALIZED / FIXED with starting value \code{start}.
#' In the case of single group analysis, argument \code{name} can be replaced by relations, i.e., the group name can be omitted.
#'
#' \code{$free_directed()} / \code{$penalize_directed()} / \code{$fix_directed()} sets all the regression coefficients from variables in \code{right} to variables in \code{left} at groups in \code{group} as FREE / PENALIZED / FIXED.
#'
#' \code{$free_undirected()} / \code{$penalize_undirected()} / \code{$fix_undirected()} sets all the covariances among variables in \code{both} at groups in \code{group} as FREE / PENALIZED / FIXED.
#' Note that this method all always not modify the variance of variables specified in \code{both}.
#'
#' \code{$free_block()} / \code{$penalize_block()} / \code{$fix_block()} sets all the parameters belonging to \code{block} at \code{group} with \code{type} as FREE / PENALIZED / FIXED.
#'
#' \code{$free_heterogeneity()} / \code{$penalize_heterogeneity()} / \code{$fix_heterogeneity()} sets every target coefficient as FREE / PENALIZED / FIXED.
#' A target coefficient should satisfy that (1) it belongs to \code{block} at \code{group}, and
#' (2) it has either free or penalized reference component.
#' The method is only available in the case of multi-group analysis with specified reference group.
#'
#' Inside \code{lslx} object, every coefficient (or parameter) has its own name and belongs to some block.
#' \itemize{
#' \item{The coefficient name is constructed by a relation and a group name.
#' A relation is defined by combining a left variable name, an operator, and a right variable name.
#' For example, \code{"y1<-f1"} is a relation to represent the coefficient from \code{"f1"} to \code{"y1"}.
#' Note that in relation we only use operators \code{"<-"} and \code{"<->"}.
#' \code{"->"} is illegal and no distinction between \code{"="} and \code{"~"} are made.
#' In multi-group analysis cases, a coefficient name requires explicitly specifying group name.
#' For example, \code{"y1<-f1/G1"} is the parameter name for the path coefficient from \code{"f1"} to \code{"y1"} in group \code{"G1"}.}
#' \item{
#' Block is defined by the types of left variable and right variable, and the operator.
#' In \pkg{lslx}, \code{"y"} is used to indicate observed response, \code{"f"} is used for latent factor, and \code{"1"} is for intercept.
#' Hence, \code{"y<-f"} is the block that contains all the coefficients from latent factors to observed responses.
#' There are 10 possible distinct blocks: \code{"f<-1"}, \code{"y<-1"}, \code{"f<-f"}, \code{"f<-y"}, \code{"y<-f"}, \code{"y<-y"}, \code{"f<->f"}, \code{"f<->y"}, \code{"y<->f"}, and \code{"y<->y"}
#' }
#' }
#' Arguments in set-related methods may rely on the naming rule of coefficent name and block to modify model specification.
#'
#' @section Fit-Related Methods:
#'
#' \preformatted{$fit(penalty_method = "lasso", lambda_grid = "default", delta_grid = "default",
#'   step_grid = "default", loss = "default", algorithm = "default", 
#'   missing_method = "default", start_method = "default", 
#'   lambda_direction = "default", lambda_length = 50L, delta_length = 3L,
#'   threshold_value = 0.3, iter_out_max = 100L, iter_in_max = 50L, 
#'   iter_other_max = 500L, iter_armijo_max = 100L, tol_out = 1e-3, 
#'   tol_in = 1e-3, tol_other = 1e-7, step_size = 0.5, momentum = 0, 
#'   armijo = 1e-5, ridge_cov = 0, ridge_hessian = 1e-4, warm_start = TRUE, 
#'   positive_variance = TRUE, minimum_variance = 1e-4, enforce_cd = FALSE, 
#'   random_update = TRUE, weight_matrix = NULL, verbose = TRUE)
#' $fit_none(...)
#' $fit_lasso(lambda_grid = "default", ...)
#' $fit_ridge(lambda_grid = "default", ...)
#' $fit_elastic_net(lambda_grid = "default", delta_grid = "default", ...)
#' $fit_mcp(lambda_grid = "default", delta_grid = "default", ...)
#' $fit_forward(step_grid = "default", ...)
#' $fit_backward(step_grid = "default", ...)}
#' 
#'\describe{
#'\item{\bold{Arguments}}{
#'
#'}
#'\item{\code{penalty_method}}{A \code{character} to specify the penalty method.
#' There are two class penalty methods can be specified.
#' For penalized estimation with a regularizer, \code{"lasso"}, \code{"ridge"}, \code{"elastic_net"}, and \code{"mcp"} can be used.
#' For penalized estimation with stepwise search, \code{"forward"} and \code{"backward"} can be implemented. 
#' If no penalty is considered, we can set \code{penalty_method} as \code{"none"}.}
#'\item{\code{lambda_grid}}{A non-negative \code{numeric} to specify penalty levels for the regularizer.
#'   If it is set as \code{"default"}, its value will be generated automatically based on the variable scales.}
#' \item{\code{delta_grid}}{A non-negative \code{numeric} to specify the combination weight for \code{"elastic_net"} or the convexity level for \code{"mcp"}.
#'   If it is set as \code{"default"}, its value will be generated automatically.}
#' \item{\code{step_grid}}{A non-negative \code{numeric} to specify a grid for stepwise search with \code{"forward"} and \code{"backward"}.}
#'\item{\code{loss}}{A \code{character} to determine the loss function.
#'   The current version supports \code{"ml"} (maximum likelihood), \code{"uls"} (unweighted least squares), 
#'   \code{"dwls"} (diagonal weighted least squres), and \code{"wls"} (weighted least squares).
#'   The maximum likelihood is only available for all continuous response variables.
#'   If the argument is set as \code{"default"}, then (1) \code{"ml"} will be implemented for all continuous response variables; 
#'   (2) \code{"dwls"} will be implemented for ordinal or mixed types response variables.}
#'\item{\code{algorithm}}{A \code{character} to determine the method of optimization.
#'   The current version supports \code{"gd"} (gradient descent), \code{"bfgs"} (Broyden-Fletcher-Goldfarb-Shanno), \code{"fisher"} (Fisher scoring), and \code{"dynamic"} (an adaptive algorithm).
#'   If the argument is set as \code{"default"}, then \code{"dynamic"} will be implemented.}
#'\item{\code{missing_method}}{A \code{character} to determine the method for handling missing data (or \code{NA}).
#'   The current version supports \code{"two_stage"} and \code{"listwise_deletion"}.
#'   If the argument is set as \code{"default"} and a raw data set is available, the \code{"two_stage"} will be implemented.
#'   If the argument is set as \code{"default"} and only moment data is available, the \code{"listwise_deletion"} will be used
#'   (actually, in this case no missing presences).}
#'\item{\code{start_method}}{A \code{character} to determine the method for calculating unspecified starting values.
#'   The current version supports \code{"mh"} (McDonald & Hartmann, 1992) and \code{"heuristic"}.
#'   If the argument is set as \code{"default"}, the \code{"mh"} will be implemented.}
#'\item{\code{lambda_direction}}{A \code{character} to determine the "direction" of \code{lambda_grid}.
#'   \code{"decrease"} sorts \code{lambda_grid} from large to small. 
#'   On the contrary, \code{"increase"} sorts \code{lambda_grid} from small to large. 
#'   If the argument is set as \code{"default"}, \code{"increase"} will be used when the smallest element of \code{lambda_grid} is larger than zero;
#'   otherwise, \code{"decrease"} is assumed.}
#'\item{\code{lambda_length}}{A \code{numeric} to specify the length of automatically generated \code{lambda_grid} under \code{lambda_grid = "default"}.}
#'\item{\code{delta_length}}{A \code{numeric} to specify the length of automatically generated \code{delta_grid} under \code{delta_grid = "default"}.}
#'\item{\code{threshold_value}}{A \code{numeric} to specify the "largest threshold value" for \code{lambda_grid} initialization.}
#'\item{\code{iter_out_max}}{A positive \code{integer} to specify the maximal iterations for outer loop of the modified \code{glmnet} algorithm.}
#'\item{\code{iter_in_max}}{A positive \code{integer} to specify the maximal iterations for inner loop of the modified \code{glmnet} algorith.}
#'\item{\code{iter_other_max}}{A positive \code{integer} to specify the maximal iterations for other loop.}
#'\item{\code{iter_armijo_max}}{A positive \code{integer} to specify the maximal iterations for searching step-size via Armijo rule.}
#'\item{\code{tol_out}}{A small positive \code{numeric} to specify the tolerance (convergence criterion) for outer loop of the modified \code{glmnet} algorithm.}
#'\item{\code{tol_in}}{A small positive \code{numeric} to specify the tolerance (convergence criterion) for inner loop of the modified \code{glmnet} algorithm.}
#'\item{\code{tol_other}}{A small positive \code{numeric} to specify the tolerance (convergence criterion) for other loop.}
#'\item{\code{step_size}}{A positive \code{numeric} smaller than one to specify the step-size.}
#'\item{\code{momentum}}{A \code{numeric} between 0 and 1 for momentum parameter.}
#'\item{\code{armijo}}{A small positive \code{numeric} for the constant in Armijo rule.}
#'\item{\code{ridge_cov}}{A small positive \code{numeric} for the ridge of sample covariance matrix.}
#'\item{\code{ridge_hessian}}{A small positive \code{numeric} for the ridge of approximated hessian in optimization.}
#'\item{\code{ridge_weight}}{A small positive \code{numeric} for the ridge of weight matrix in weighted least squares.}
#'\item{\code{warm_start}}{A \code{logical} to specify whether the warm start approach should be used.}
#'\item{\code{positive_variance}}{A \code{logical} to specify whether the variance estimate should be constrained to be larger than \code{minimum_variance}.}
#'\item{\code{minimum_variance}}{A \code{numeric} to specify the minimum value of variance if \code{positive_variance = TRUE}.}
#'\item{\code{enforce_cd}}{A \code{logic} to specify whether coordinate descent should be used when no penalty function is used. Its default value is TRUE.}
#'\item{\code{random_update}}{A \code{logic} to specify whether coordinate descent should be conducted in a random order. Its default value is TRUE.}
#'\item{\code{weight_matrix}}{A \code{list} with length equaling to the number of groups to specify a user-defined weight matrix for least squares loss.}
#'\item{\code{verbose}}{A \code{logical} to specify whether messages made by \code{lslx} should be printed.}
#'\item{\code{...}}{Other passing arguments for calling \code{$fit()}.}
#'}
#' \bold{Details}
#'
#' Fit-related methods are used for fitting model to data.
#' The success of these methods may depend on the specified fitting control.
#' For details of optimization algorithm, see the section of \emph{Optimization Algorithm}.
#'
#' \code{$fit()} fits the specified model to data by minimizing a penalized loss function.
#' It is the most comprehensive fit method and hence many arguments can be specified.
#'
#' \code{$fit_none()} fits the specified model to data by minimizing a loss function without penalty.
#' It is a user convinient wrapper of \code{$fit()} with \code{penalty_method = "none"}.
#'
#' \code{$fit_lasso()} fits the specified model to data by minimizing a loss function with lasso penalty (Tibshirani, 1996).
#' It is a user convinient wrapper of \code{$fit()} with \code{penalty_method = "lasso"}.
#'
#' \code{$fit_ridge()} fits the specified model to data by minimizing a loss function with ridge penalty (Hoerl & Kennard, 1970).
#' It is a user convinient wrapper of \code{$fit()} with \code{penalty_method = "ridge"}.
#' 
#' \code{$fit_elastic_net()} fits the specified model to data by minimizing a loss function with elastic net penalty (Zou & Hastie, 2005).
#' It is a user convinient wrapper of \code{$fit()} with \code{penalty_method = "elastic_net"}.
#' 
#' \code{$fit_mcp()} method fits the specified model to data by minimizing a loss function with mcp (Zhang, 2010).
#' It is a user convinient wrapper of \code{$fit()} with \code{penalty_method = "mcp"}.
#' 
#' \code{$fit_forward()} method fits the specified model to data by minimizing a loss function with forward searching.
#' It is a user convinient wrapper of \code{$fit()} with \code{penalty_method = "forward"}.
#'
#' \code{$fit_backward()} method fits the specified model to data by minimizing a loss function with backward searching.
#' It is a user convinient wrapper of \code{$fit()} with \code{penalty_method = "backward"}.
#' 
#'
#' @section Summarize Method:
#' \preformatted{$summarize(selector, lambda, delta, step, standard_error = "default", 
#'   debias = "default", inference = "default", alpha_level = .05, 
#'   include_faulty = FALSE, style = "default", mode = "default", digit = 3, 
#'   interval = TRUE, output)}
#'\describe{
#'\item{\bold{Arguments}}{
#'
#'}
#'\item{\code{selector}}{A \code{character} to specify a selector for determining an optimal penalty level.
#'   Its value can be any one in \code{"aic"}, \code{"aic3"}, \code{"caic"}, \code{"bic"}, \code{"abic"}, \code{"hbic"},
#'   or their robust counterparts \code{"raic"}, \code{"raic3"}, \code{"rcaic"}, \code{"rbic"}, \code{"rabic"}, \code{"rhbic"} if raw data is available.}
#'\item{\code{lambda}}{A \code{numeric} to specific a chosen optimal penalty level. 
#'   If the specified \code{lambda} is not in \code{lambda_grid}, a nearest legitimate value will be used. }
#'\item{\code{delta}}{A \code{numeric} to specific a chosen optimal convexity level.
#'   If the specified \code{delta} is not in \code{delta_grid}, a nearest legitimate value will be used.}
#'\item{\code{step}}{A \code{numeric} to specific a chosen step for stepwise searching.
#'   If the specified \code{step} is not in \code{step_grid}, a nearest legitimate value will be used.}
#'\item{\code{standard_error}}{A \code{character} to specify the standard error to be used for hypothesis testing.
#'   The argument can be either \code{"sandwich"}, \code{"expected_information"}, and \code{"observed_information"}.
#'   If it is specified as \code{"default"}, it will be set as
#'   (1) \code{"sandwich"} when raw data is available; (2) \code{"observed_information"} when only moment data is available.}
#'\item{\code{debias}}{A \code{character} to specify a debias method for obtaining a debiased estimator.
#'   Its value can be either \code{"none"} or \code{"one_step"}. 
#'   If it is specified as \code{"default"}, \code{"none"} will be used unless \code{post = "polyhedral"} is used.}
#'\item{\code{inference}}{A \code{character} to specify the method for post selection inference.
#'   The current version supports \code{"naive"}, \code{"polyhedral"}, and \code{"scheffe"}.
#'   If it is specified as \code{"default"}, \code{"naive"} will be used.}
#'\item{\code{alpha_level}}{A \code{numeric} to specify the alpha level for constructing 1 - alpha confidence intervals.}
#'\item{\code{include_faulty}}{A \code{logical} to specify whether non-convergence or non-convexity results should be removed for penalty level selection.
#'   Non-convergence result determined by examining the maximal elements of absolute objective gradient and the number of iterations.
#'   non-convexity result is determined by checking the minimum of univariate approximate hessian.}
#'\item{\code{style}}{A \code{character} to specify whether the style of summary.
#'   Its value must be either \code{"default"}, \code{"manual"}, \code{"minimal"}, or \code{"maximal"}.}
#'\item{\code{mode}}{A \code{character} to specify the mode of summary. 
#'   Its value must be \code{"print"} or \code{"return"}. 
#'   \code{"print"} will print the summary result and \code{"return"} will return a \code{list} of the summary result.
#'   If it is specified as \code{"default"}, it will be set as \code{"print"}.}
#'\item{\code{digit}}{An \code{integer} to specify the number of digits to be displayed.}
#'\item{\code{interval}}{A \code{logical} to specify whether the confidence interval should be printed.}
#'}
#' \bold{Details}
#'
#' \code{$summarize()} prints a summary for the fitting result under a selected peanlty/convexity level.
#' It requires users to specify which selector should be used or a combination of \code{gamma} and \code{lambda}..
#' By default, the summary includes model information, numerical conditions, fit indices, coefficient estimates, and related statistical inference.
#' It can be modified by the \code{style} argument. 
#' For details of evaluation and inference methods, see the sections of \emph{Model Fit Evaluation} and \emph{Coefficient Evaluation}.
#'
#' @section Test-Related Methods:
#'
#' \preformatted{$test_lr(selector, lambda, delta, step, include_faulty = FALSE)
#' $test_rmsea(selector, lambda, delta, step, 
#'   alpha_level = .05, include_faulty = FALSE)
#' $test_coefficient(selector, lambda, delta, step, 
#'   standard_error = "default", debias = "default", inference = "default", 
#'   alpha_level = .05, include_faulty = FALSE)}
#'\describe{
#'\item{\bold{Arguments}}{
#'
#'}
#'\item{\code{selector}}{A \code{character} to specify a selector for determining an optimal penalty level.
#'   Its value can be any one in \code{"aic"}, \code{"aic3"}, \code{"caic"}, \code{"bic"}, \code{"abic"}, \code{"hbic"},
#'   or their robust counterparts \code{"raic"}, \code{"raic3"}, \code{"rcaic"}, \code{"rbic"}, \code{"rabic"}, \code{"rhbic"} if raw data is available.}
#'\item{\code{delta}}{A \code{numeric} to specific a chosen optimal weight for elastic net or convexity level for mcp.
#'   If the specified \code{delta} is not in \code{delta_grid}, a nearest legitimate value will be used.}
#'\item{\code{step}}{A \code{numeric} to specific a chosen step for stepwise searching.
#'   If the specified \code{step} is not in \code{step_grid}, a nearest legitimate value will be used.}
#'\item{\code{standard_error}}{A \code{character} to specify the standard error to be used for hypothesis testing.
#'   The argument can be either \code{"sandwich"}, \code{"expected_information"}, and \code{"observed_information"}.
#'   If it is specified as \code{"default"}, it will be set as
#'   (1) \code{"sandwich"} when raw data is available; (2) \code{"observed_information"} when only moment data is available.}
#'\item{\code{debias}}{A \code{character} to specify a debias method for obtaining a debiased estimator.
#'   Its value can be either \code{"none"} or \code{"one_step"}. 
#'   If it is specified as \code{"default"}, \code{"none"} will be used unless \code{post = "polyhedral"} is used.}
#'\item{\code{inference}}{A \code{character} to specify the method for post selection inference.
#'   The current version supports \code{"naive"}, \code{"polyhedral"}, and \code{"scheffe"}.
#'   If it is specified as \code{"default"}, \code{"naive"} will be used.}
#'\item{\code{alpha_level}}{A \code{numeric} to specify the alpha level for constructing 1 - alpha confidence intervals.}
#'\item{\code{include_faulty}}{A \code{logical} to specify whether non-convergence or non-convexity results should be removed for penalty level selection.
#'   Non-convergence result determined by examining the maximal elements of absolute objective gradient and the number of iteration.
#'   non-convexity result is determined by checking the minimum of univariate approximate hessian.}
#'}
#' \bold{Details}
#'
#' Test-related methods are used to obtain the result of specific statistical test.
#' So far, only tests for likelihood ratio (LR), root mean square error of approximation (RMSEA), and coefficients are available.
#'
#' \code{$test_lr()} returns a \code{data.frame} of result for likelihood ratio test.
#' If raw data is available, it also calculates a mean-adjusted statistic.
#' For details of significance test method for LR, see the section of \emph{Model Fit Evaluation}.
#'
#' \code{$test_rmsea()} returns a \code{data.frame} of result for rmsea confidence intervals.
#' If raw data is available, it also calculates a mean-adjusted confidence interval (Brosseau-Liard, Savalei & Li, 2012; Li & Bentler, 2006).
#' For details of confidence interval construction for RMSEA, see the section of \emph{Model Fit Evaluation}.
#'
#' \code{$test_coefficient()} returns a \code{data.frame} of result for coefficient significance and confidence interval.
#' For details of standard error formula for coefficients, see the section of \emph{Coefficient Evaluation}.
#'
#' @section Plot-Related Methods:
#' \preformatted{$plot_numerical_condition(condition, x_scale = "default", 
#'   mode = "default")
#' $plot_information_criterion(criterion, x_scale = "default", 
#'   mode = "default")
#' $plot_fit_index(index, x_scale = "default", mode = "default")
#' $plot_coefficient(block, left, right, both, x_scale = "default", 
#'   mode = "default")}
#' \describe{
#' \item{\bold{Arguments}}{
#'
#' }
#' \item{\code{condition}}{A \code{character} to specify which numerical conditions should be plotted.
#' Its value must be \code{"objective_value"}, \code{"objective_gradient_abs_max"}, \code{"objective_hessian_convexity"},
#' \code{"n_iter_out"}, \code{"loss_value"}, \code{"n_nonzero_coefficient"}, or their combination.}
#'  \item{\code{criterion}}{A \code{character} to specify which information criteria should be plotted.
#'  Its value must be \code{"aic"}, \code{"aic3"}, \code{"caic"}, \code{"bic"}, \code{"abic"}, \code{"hbic"},
#'  \code{"raic"}, \code{"raic3"}, \code{"rcaic"}, \code{"rbic"}, \code{"rabic"}, \code{"rhbic"}, or their combination.}
#'  \item{\code{index}}{A \code{character} to specify which fit indices should be plotted.
#'  Its value must be \code{"rmsea"}, \code{"cfi"}, \code{"nnfi"}, \code{"srmr"}, or their combination. }
#'  \item{\code{block}}{A \code{character} with length one to indicate a block such that the corresponding target coefficient will be reset.
#'  Its value must be one of \code{"f<-1"}, \code{"y<-1"}, \code{"f<-f"}, \code{"f<-y"}, \code{"y<-f"}, \code{"y<-y"}, \code{"f<->f"}, \code{"f<->y"}, \code{"y<->f"}, or \code{"y<->y"}.}
#'  \item{\code{left}}{A \code{character} to specify the variables in the left-hand side of operator in \code{block}.}
#'  \item{\code{right}}{A \code{character} to specify the variables in the right-hand side of operator in \code{block}.}
#'  \item{\code{both}}{A \code{character} to specify the variables in both sides of operator in \code{block}.}
#'  \item{\code{lambda_scale}}{A \code{character} to specify the scale of lambda (x-axis) for \code{coord_trans()} in \pkg{ggplot2}.}
#'  \item{\code{mode}}{A \code{character} to specify the mode of plot. 
#'   Its value must be \code{"plot"} or \code{"return"}. 
#'   \code{"plot"} will plot the result and \code{"return"} will return a \code{data.frame} for plot.
#'   If it is specified as \code{"default"}, it will be set as \code{"plot"}.}
#' }
#' \bold{Details}
#'
#' Plot-related methods are used for visualizing the fitting results.
#'
#' \code{$plot_numerical_condition()} plots the values of selected numerical conditions.
#' It can be used to assess the quality of optimization.
#' By default, \code{"n_iter_out"}, \code{"objective_gradient_abs_max"}, and \code{"objective_hessian_convexity"} across the given penalty levels are plotted.
#'
#' \code{$plot_information_criterion()} shows how the values of information criteria vary with penalty levels.
#' By default, \code{"aic"}, \code{"aic3"}, \code{"caic"}, \code{"bic"}, \code{"abic"}, and \code{"hbic"} are plotted.
#'
#' \code{$plot_fit_index()} shows how the values of fit indices vary with penalty levels.
#' By default, \code{"rmsea"}, \code{"cfi"}, \code{"nnfi"}, and \code{"srmr"} are plotted.
#'
#' \code{$plot_coefficient()} visualizes the solution paths of coefficients belonging to the intersection of \code{block}, \code{left}, \code{right}, and \code{both} arguments.
#' By default, all of the coefficients are plotted.
#'
#' @section Get-Related Methods:
#' \preformatted{$get_model()
#' $get_data()
#' $get_fitting()
#'
#' }
#' \bold{Details}
#'
#' Get-related methods are defined to obtain a deep copy of members inside \code{lslx}.
#' Note that all of the data members of \code{lslx} are set as private to protect the inner data structure.
#' They cannot be assessed directly via \code{$}.
#'
#' \code{$get_model()} returns a deep copy of \code{model} member in the current \code{lslx} object.
#'
#' \code{$get_data()} returns a deep copy of \code{data} member in the current \code{lslx} object.
#'
#' \code{$get_fitting()} returns a deep copy of \code{fitting} member in the current \code{lslx} object.
#'
#' @section Extract-Related Methods:
#'
#' \preformatted{$extract_specification()
#' $extract_saturated_cov()
#' $extract_saturated_mean()
#' $extract_saturated_moment_acov()
#'
#' $extract_penalty_level(selector, lambda, delta, step, 
#'   include_faulty = FALSE)
#' 
#' $extract_numerical_condition(selector, lambda, delta, step, 
#'   include_faulty = FALSE)
#' $extract_information_criterion(selector, lambda, delta, step, 
#'   include_faulty = FALSE)
#' $extract_fit_index(selector, lambda, delta, step, include_faulty = FALSE)
#' $extract_cv_error(selector, lambda, delta, step, include_faulty = FALSE)
#' 
#' $extract_coefficient(selector, lambda, delta, step, type = "default", 
#'   include_faulty = FALSE)
#' $extract_debiased_coefficient(selector, lambda, delta, step, type = "default", 
#'   include_faulty = FALSE)
#'
#' $extract_implied_cov(selector, lambda, delta, step, 
#'   include_faulty = FALSE)
#' $extract_implied_mean(selector, lambda, delta, step, 
#'   include_faulty = FALSE)
#' $extract_residual_cov(selector, lambda, delta, step, 
#'   include_faulty = FALSE)
#' $extract_residual_mean(selector, lambda, delta, step, 
#'   include_faulty = FALSE)
#'
#' $extract_coefficient_matrix(selector, lambda, delta, step, block, 
#'   include_faulty = FALSE)
#' 
#' $extract_moment_jacobian(selector, lambda, delta, step, 
#'   type = "default", include_faulty = FALSE)
#'
#' $extract_expected_information(selector, lambda, delta, step, 
#'   type = "default", include_faulty = FALSE)
#' $extract_observed_information(selector, lambda, delta, step, 
#'   type = "default", include_faulty = FALSE)
#' 
#' $extract_bfgs_hessian(selector, lambda, delta, step, type = "default", 
#'   include_faulty = FALSE)
#' $extract_score_acov(selector, lambda, delta, step, type = "default", 
#'   include_faulty = FALSE)
#' $extract_coefficient_acov(selector, lambda, delta, step, 
#'   standard_error = "default", type = "default", include_faulty = FALSE)
#'
#' $extract_loss_gradient(selector, lambda, delta, step, type = "default", 
#'   include_faulty = FALSE)
#' $extract_regularizer_gradient(selector, lambda, delta, step, 
#'   type = "default", include_faulty = FALSE)
#' $extract_objective_gradient(selector, lambda, delta, step, 
#'   type = "default", include_faulty = FALSE)}
#' \describe{
#' \item{\bold{Arguments}}{
#'
#' }
#'\item{\code{selector}}{A \code{character} to specify a selector for determining an optimal penalty level.
#'   Its value can be any one in \code{"aic"}, \code{"aic3"}, \code{"caic"}, \code{"bic"}, \code{"abic"}, \code{"hbic"},
#'   or their robust counterparts \code{"raic"}, \code{"raic3"}, \code{"rcaic"}, \code{"rbic"}, \code{"rabic"}, \code{"rhbic"} if raw data is available.}
#'\item{\code{lambda}}{A \code{numeric} to specific a chosen optimal penalty level. 
#'   If the specified \code{lambda} is not in \code{lambda_grid}, a nearest legitimate value will be used. }
#'\item{\code{delta}}{A \code{numeric} to specific a chosen optimal weight for elastic net or convexity level for mcp.
#'   If the specified \code{delta} is not in \code{delta_grid}, a nearest legitimate value will be used.}
#'\item{\code{step}}{A \code{numeric} to specific a chosen step for stepwise searching.
#'   If the specified \code{step} is not in \code{step_grid}, a nearest legitimate value will be used.}
#'\item{\code{standard_error}}{A \code{character} to specify the standard error to be used for hypothesis testing.
#'   The argument can be either \code{"sandwich"}, \code{"expected_information"}, and \code{"observed_information"}.
#'   If it is specified as \code{"default"}, it will be set as
#'   (1) \code{"sandwich"} when raw data is available; (2) \code{"observed_information"} when only moment data is available.}
#' \item{\code{block}}{A \code{character} with length one to indicate a block such that the corresponding target coefficient will be reset.
#'  Its value must be \code{"f<-1"}, \code{"y<-1"}, \code{"f<-f"}, \code{"f<-y"}, \code{"y<-f"}, \code{"y<-y"}, \code{"f<->f"}, \code{"f<->y"}, \code{"y<->f"}, or \code{"y<->y"}.}
#'  \item{\code{type}}{A \code{character} to specify the type of parameters that will be used to compute the extracted quantity. 
#'  The argument can be either \code{"all"}, \code{"fixed"}, \code{"free"}, \code{"pen"}, \code{"effective"} (include \code{"free"} + \code{"selected"}), and \code{"selected"} (non-zero element of \code{"pen"}).
#'  If it is specified as \code{"default"}, it will be set as \code{all}.}
#'\item{\code{include_faulty}}{A \code{logical} to specify whether non-convergence or non-convexity results should be removed for penalty level selection.
#'   Non-convergence result determined by examining the maximal elements of absolute objective gradient and the number of iterations.
#'   non-convexity result is determined by checking the minimum of univariate approximate hessian.}
#'}
#' \bold{Details}
#'
#' Many extract-related methods are defined to obtain quantities that can be used for further SEM applications or model diagnosis.
#' Some of these quantities only depend on data (e.g., saturated sample covariance matrix), but some of them relies on a penalty level (e.g., gradient of objective function).
#' An optimal penalty level can be determined by specifying a \code{selector} or a combination of \code{gamma} and \code{lambda}.
#' When implementing fit-related methods, \code{lslx} only records necessary results for saving memory.
#' Therefore, the extract-related methods not only extract objects but may possibly re-compute some of them.
#' If the extracted quantity is a function of model coefficients, note that fixed coefficient is still considered as a valid variable.
#' For example, the sub-gradient of objective function is calculated by considering both estimated coefficients and fixed coefficients.
#' Hence, it is not appropriate to use all the elements of the sub-gradient to evaluate the optimality of solution because the values of fixed coefficients are not optimized.
#' Fortunately, the \code{type} argument can be used to choose desired parameters by their types.
#'
#' \code{$extract_specification()} returns a \code{data.frame} of model specification.
#'
#' \code{$extract_saturated_cov()} returns a \code{list} of saturated sample covariance matrix(s).
#'
#' \code{$extract_saturated_mean()} returns a \code{list} of saturated sample mean vector(s).
#'
#' \code{$extract_saturated_moment_acov()} returns a \code{list} of asymptotic covariance matrix(s) of saturated moments.
#' Note that if raw data is not available, asymptotic covariance matrix is calculated by assuming normality for data.
#'
#' \code{$extract_penalty_level()} returns a \code{character} of the index name of the optimal penalty level.
#'
#' \code{$extract_numerical_condition()} returns a \code{numeric} of the numerical conditions.
#'
#' \code{$extract_information_criterion()} returns a \code{numeric} of the values of information criteria.
#'
#' \code{$extract_fit_indice()} returns a \code{numeric} of the values of fit indices.
#'
#' \code{$extract_coefficient()} returns a \code{numeric} of estimates of the coefficients.
#'
#' \code{$extract_implied_cov()} returns a \code{list} of model-implied covariance matrix(s).
#'
#' \code{$extract_implied_mean()} returns a \code{list} of model-implied mean vector(s).
#'
#' \code{$extract_residual_cov()} returns a \code{list} of residual matrix(s) of covariance.
#'
#' \code{$extract_residual_mean()} returns a \code{list} of residual vector(s) of mean.
#'
#' \code{$extract_coefficient_matrix()} returns a \code{list} of coefficient matrix(s) specified by \code{block}.
#'
#' \code{$extract_moment_jacobian()} returns a \code{matrix} of Jacobian of moment structure.
#'
#' \code{$extract_expected_information()} returns a \code{matrix} of the expected Fisher information matrix.
#'
#' \code{$extract_observed_information()} returns a \code{matrix} of the observed Fisher information matrix.
#' Note that the observed information matrix is calculated via numerical differentiation for the gradient of loss.
#'
#' \code{$extract_bfgs_hessian()} returns a \code{matrix} of the BFGS Hessian matrix.
#'
#' \code{$extract_score_acov()} returns a \code{matrix} of the asymptotic covariance of scores.
#'
#' \code{$extract_coefficient_acov()} returns a \code{matrix} of the asymptotic covariance of coefficients.
#' For details of standard error formula, see the section of \emph{Coefficient Evaluation}.
#'
#' \code{$extract_loss_gradient()} returns a \code{matrix} of the gradient of loss function.
#'
#' \code{$extract_regularizer_gradient()} returns a \code{matrix} of the sub-gradient of regularizer.
#'
#' \code{$extract_objective_gradient()} returns a \code{matrix} of the sub-gradient of objective function.
#'
#'
#' @references
#'
#' Akaike, H. (1974). A new look at the statistical model identification. IEEE Transactionson Automatic Control, 19(6), 716–723.
#'
#' Bates, D., & Eddelbuettel, D. (2013). Fast and Elegant Numerical Linear Algebra Using the {RcppEigen} Package. Journal of Statistical Software, 52(5), 1–24.
#'
#' Bentler, P. M. (1995). EQS structural equations program manual. Encino, CA: Multivariate Software.
#'
#' Bentler, P. (1990). Comparative fit indices in structural models. Psychological Bulletin, 107(2), 238–246.
#'
#' Berk, R., Brown, L., Buja, A., Zhang, K., & Zhao, L. (2013). Valid postselection inference. The Annals of Statistics, 41(2), 802–837. 
#'
#' Bozdogan, H. (1987). Model selection and Akaike’s Information Criterion (AIC): The general theory and its analytical extensions. Psychometrika, 52(3), 345–370.
#'
#' Browne, M. W. (1984). Asymptotic distribution-free methods for the analysis of covariance structures. British Journal of Mathematical and Statistical Psychology, 37(1), 62–83.
#'
#' Chang, W. (2017). R6: Classes with Reference Semantics.
#'
#' Eddelbuettel, D., & François, R. (2011). {Rcpp}: Seamless {R} and {C++} Integration. Journal of Statistical Software, 40(8), 1–18.
#'
#' Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1–22.
#'
#' Haughton, D. M. A., Oud, J. H. L., & Jansen, R. A. R. G. (1997). Information and other criteria in structural equation model selection. Communications in Statistics - Simulation and Computation, 26(4), 1477–1516.
#'
#' Hoerl, A. E., & Kennard, R. W. (1970). Ridge Regression: Biased Estimation for Nonorthogonal Problems. Technometrics, 12(1), 55–67.
#' 
#' Huang, P. H. (2018). A Penalized Likelihood Method for Multi-Group Structural Equation Modeling. British Journal of Mathematical and Statistical Psychology, 71(3),  499-522.
#' 
#' Huang, P. H. (2020). lslx: Semi-Confirmatory Structural Equation Modeling via Penalized Likelihood. Journal of Statistical Software. 93(7), 1-37.
#' 
#' Huang, P. H. (in press). Post-selection inference in Structural Equation Modeling. Multivariate Behavioral Research.
#'
#' Huang, P. H., Chen, H., & Weng, L. J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354.
#'
#' Brosseau-Liard, P. E., Savalei, V., & Li, L. (2012). An Investigation of the Sample Performance of Two Nonnormality Corrections for RMSEA. Multivariate Behavioral Research, 47(6), 904-930.
#' 
#' Lee, J. D., Sun, D. L., Sun, Y., & Taylor, J. E. (2016). Exact postselection inference, with application to the lasso. The Annals of Statistics, 44(3), 907–927.
#' 
#' Li, L., & Bentler, P. M. (2006). Robust statistical tests for evaluating the hypothesis of close fit of misspecified mean and covariance structural models.
#' UCLA Statistics Preprint #506. Los Angeles: University of California.
#'
#' Mazumder, R., Friedman, J. H., & Hastie, T. (2011). SparseNet: Coordinate Descent With Nonconvex Penalties. Journal of the American Statistical Association, 106(495), 1125–1138.
#'
#' McDonald, R. P., & Hartmann, W. M. (1992). A procedure for obtaining initial values of parameters in the RAM model. Multivariate Behavioral Research, 27(1), 57–76.
#'
#' Pötscher, B. M. (1991). Effects of model selection on inference. Econometric Theory, 7(2), 163-185.
#'
#' Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1–36.
#'
#' Rubin, D. B. (1976). Inference and Missing Data. Biometrika, 63(3), 581–592.
#'
#' Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and standard errors in covariance structure analysis.
#' In A. von Eye & C. C. Clogg (Eds.), Latent variable analysis: Applications to developmental research (pp. 399–419). Thousand Oaks, CA: Sage.
#'
#' Savalei, V. & Falk, C. F. (2014). Robust two-stage approach outperforms robust full information maximum likelihood with incomplete nonnormal data. Structural Equation Modeling: A Multidisciplinary Journal, 21(2), 280-302.
#'
#' Savalei, V. & Bentler, P. M. (2009). A Two-Stage Approach to Missing Data: Theory and Application to Auxiliary Variables, Structural Equation Modeling: A Multidisciplinary Journal, 16(3), 477-497.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. The Annals of Statistics, 6(2), 461–464.
#'
#' Sclove, S. L. (1987). Application of model-selection criteria to some problems in multivariate analysis. Psychometrika, 52(3), 333–343.
#'
#' Steiger, J. H. (1998). A Note on Multiple Sample Extensions of the RMSEA Fit Index. Structural Equation Modeling-a Multidisciplinary Journal, 5(4), 411–419.
#'
#' Steiger, J. H., & Lind, J. C. (1980). Statistically-based tests for the number of common factors. In Paper presented at the annual meeting of the Psychometric Society.
#'
#' Tibshirani, R. (1996). Regression Selection and Shrinkage via the Lasso. Journal of the Royal Statistical Society B, 58(1), 267–288.
#'
#' Tucker, L. R., & Lewis, C. (1973). A reliability coefficient for maximum likelihood factor analysis. Psychometrika, 38(1), 1–10.
#'
#' Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for mean and covariance structure analysis with nonnormal missing data. Sociological Methodology, 30(1), 165–200.
#'
#' Yuan, K. H., & Hayashi, K. (2006). Standard errors in covariance structure models: Asymptotics versus bootstrap. British Journal of Mathematical and Statistical Psychology, 59(2), 397–417.
#'
#' Yuan, K.-H., & Lu, L. (2008). SEM with missing data and unknown population distributions using two-stage ML: Theory and its application. Multivariate Behavioral Research, 43(4), 621–652.
#'
#' Yuan, G. X., Ho, C. H., & Lin, C. J. (2012). An Improved GLMNET for L1-regularized Logistic Regression. Journal of Machine Learning Research, 13(1), 1999–2030.
#'
#' Zhang, C. H. (2010). Nearly unbiased variable selection under minimax concave penalty. Annals of Statistics, 38(2), 894–942.
#'
#' Zou, H., & Hastie, T. (2005). Regularization and Variable Selection via the Elastic Net. Journal of the Royal Statistical Society B, 67(2), 301–320.
#'
#' @examples
#' \donttest{
#' ## EXAMPLE: Regression Analysis with Lasso Penalty ##
#' # run `vignette("regression-analysis")` to see the vignette
#' # generate data for regression analysis
#' set.seed(9487)
#' x <- matrix(rnorm(2000), 200, 10)
#' colnames(x) <- paste0("x", 1:10)
#' y <- matrix(rnorm(200), 200, 1)
#' data_reg <- data.frame(y, x)
#'
#' # specify regression model with penalized coefficients
#' model_reg <- "y <= x1 + x2 + x3 + x4
#'               y <~ x5 + x6 + x7 + x8 + x9 + x10"
#'
#' # initialize lslx object via specified model and raw data
#' lslx_reg <- lslx$new(model = model_reg, 
#'                      data = data_reg)
#'
#' # fit specified model to data with lasso under specified penalty levels
#' lslx_reg$fit(penalty_method = "lasso",
#'              lambda_grid = seq(.00, .30, .02))
#'
#' # summarize fitting result under penalty level selected by 'aic'
#' lslx_reg$summarize(selector = "aic")
#'
#'
#' ## EXAMPLE: Semi-Confirmatory Factor Analysis ##
#' # run `vignette("factor-analysis")` to see the vignette
#' # specify semi-confirmatory factor analysis model
#' model_fa <- "visual  :=> x1 + x2 + x3
#'              textual :=> x4 + x5 + x6
#'              speed   :=> x7 + x8 + x9
#'              visual  :~> x4 + x5 + x6 + x7 + x8 + x9
#'              textual :~> x1 + x2 + x3 + x7 + x8 + x9
#'              speed   :~> x1 + x2 + x3 + x4 + x5 + x6
#'              visual  <=> 1 * visual
#'              textual <=> 1 * textual
#'              speed   <=> 1 * speed"
#'           
#' # initialize lslx object via specified model and raw data
#' lslx_fa <- lslx$new(model = model_fa,
#'                     data = lavaan::HolzingerSwineford1939)
#'                     
#' # fit with mcp under specified penalty levels and convexity levels
#' lslx_fa$fit(penalty_method = "mcp", 
#'             lambda_grid = seq(.02, .60, .02), 
#'             delta_grid = c(1.5, 3.0, Inf))
#' 
#' # summarize fitting result under penalty level selected by 'bic'
#' lslx_fa$summarize(selector = "bic")
#'
#'
#' ## EXAMPLE: Semi-Confirmatory Structural Equation Modeling ##
#' # run `vignette("structural-equation-modeling")` to see the vignette
#' # specify structural equation modeling model
#' model_sem <- "fix(1) * x1 + x2 + x3   <=: ind60
#'               fix(1) * y1 + y2 + y3 + y4 <=: dem60
#'               fix(1) * y5 + y6 + y7 + y8 <=: dem65
#'               dem60 <= ind60
#'               dem65 <= ind60 + dem60"
#'
#' # initialize lslx object via specified model and sample moments
#' lslx_sem <- lslx$new(model = model_sem,
#'                      sample_cov = cov(lavaan::PoliticalDemocracy),
#'                      sample_size = nrow(lavaan::PoliticalDemocracy))
#'
#' # set some covariances of errors as penalized 
#' lslx_sem$penalize_coefficient(name = c("y1<->y5",
#'                                        "y2<->y4",
#'                                        "y2<->y6",
#'                                        "y3<->y7",
#'                                        "y4<->y8",
#'                                        "y6<->y8"))
#'
#' # fit with lasso under default penalty levels
#' lslx_sem$fit_lasso(lambda_length = 25)
#' 
#' # summarize fitting result under penalty level selected by 'abic'
#' lslx_sem$summarize(selector = "abic")
#'
#'
#' ## EXAMPLE: Factor Analysis with Missing Data ##
#' # run `vignette("missing-data-analysis")` to see the vignette
#' # create missing values for x5 and x9 by the code in package semTools
#' data_miss <- lavaan::HolzingerSwineford1939
#' data_miss$x5 <- ifelse(data_miss$x1 <= quantile(data_miss$x1, .3), 
#'                        NA, data_miss$x5)
#' data_miss$age <- data_miss$ageyr + data_miss$agemo/12
#' data_miss$x9 <- ifelse(data_miss$age <= quantile(data_miss$age, .3), 
#'                        NA, data_miss$x9)
#' 
#' # specify confirmatory factor analysis model
#' model_miss <- "visual  :=> x1 + x2 + x3 
#'                textual :=> x4 + x5 + x6
#'                speed   :=> x7 + x8 + x9
#'                visual  <=> 1 * visual
#'                textual <=> 1 * textual
#'                speed   <=> 1 * speed"
#'
#' # "ageyr" and "agemo" are set as auxiliary variables
#' lslx_miss <- lslx$new(model = model_miss,
#'                       data = data_miss,
#'                       auxiliary_variable = c("ageyr", "agemo"))
#'                       
#' # penalize all covariances among residuals
#' lslx_miss$penalize_block(block = "y<->y", 
#'                          type = "fixed",
#'                          verbose = FALSE)
#' 
#' # fit with lasso under default penalty levels
#' lslx_miss$fit_lasso(lambda_length = 25)
#' 
#' # summarize fitting result under penalty level selected by 'raic'
#' lslx_miss$summarize(selector = "raic")
#' 
#' 
#' ## EXAMPLE: Multi-Group Factor Analysis ##
#' # run `vignette("multi-group-analysis")` to see the vignette
#' # specify multi-group factor analysis model
#' model_mgfa <- "visual  :=> 1 * x1 + x2 + x3
#'                textual :=> 1 * x4 + x5 + x6
#'                speed   :=> 1 * x7 + x8 + x9"
#' 
#' # "school" is set as group variable and "Pasteur" is specified as reference
#' lslx_mgfa <- lslx$new(model = model_mgfa,
#'                       data = lavaan::HolzingerSwineford1939,
#'                       group_variable = "school",
#'                       reference_group = "Pasteur")
#'
#' # penalize increment components of loadings and intercepts in 'Grant-White'
#' lslx_mgfa$penalize_heterogeneity(block = c("y<-1", "y<-f"), 
#'                                  group = "Grant-White")
#'
#' # free increment components of means of latent factors in 'Grant-White'
#' lslx_mgfa$free_block(block = "f<-1", 
#'                      group = "Grant-White")
#'                      
#' # fit with mcp under default penalty levels and specified convexity levels
#' lslx_mgfa$fit_mcp(lambda_length = 25)
#' 
#' # summarize fitting result under penalty level selected by 'hbic'
#' lslx_mgfa$summarize(selector = "hbic")
#' }
#' 
#' 
#' @export
lslx <-
  R6::R6Class(
    classname = "lslx",
    inherit = prelslx,
    private = list(
      model = "lslxModel",
      data = "lslxData",
      fitting = "lslxFitting"
    )
  )
