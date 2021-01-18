
######################################################################################################*
######################################################################################################*
#' @title Minimax and Standardized Maximin D-Optimal Designs
#' @description
#'
#'  Finds minimax and standardized maximin D-optimal designs for linear and nonlinear models.
#'  It should be used when the user assumes the unknown parameters belong to a parameter region
#'  \eqn{\Theta}, which is called ``region of uncertainty'',  and the
#'    purpose is to protect the experiment from the worst case scenario
#'   over \eqn{\Theta}.
#'
#'
#'@param formula A linear or nonlinear model \code{\link[stats]{formula}}.
#' A symbolic description of the model consists of predictors and the unknown model parameters.
#' Will be coerced to a \code{\link[stats]{formula}} if necessary.
#'@param predvars A vector of characters. Denotes the predictors in the \code{\link[stats]{formula}}.
#'@param parvars A vector of characters. Denotes the unknown parameters in the \code{\link[stats]{formula}}.
#' @param family A description of the response distribution and the link function to be used in the model.
#'  This can be a family function, a call to a family function or a character string naming the family.
#'   Every family function has a link argument allowing to specify the link function to be applied on the response variable.
#'    If not specified, default links are used. For details see \code{\link[stats]{family}}.
#'     By default, a linear gaussian model \code{gaussian()} is applied.
#' @param lx Vector of lower bounds for the predictors. Should be in the same order as \code{predvars}.
#' @param ux Vector of upper bounds for the predictors. Should be in the same order as \code{predvars}.
#' @param lp Vector of lower bounds for the model parameters. Should be in the same order as \code{parvars} or \code{param} in the argument \code{fimfunc}.
#' @param up Vector of upper bounds for the model parameters. Should be in the same order as \code{parvars} or \code{param} in the argument \code{fimfunc}.
#' When a parameter is known (has a fixed value), its associated lower and upper bound values  in \code{lp} and \code{up}  must be set equal.
#' @param iter Maximum number of iterations.
#' @param k Number of design points. Must be at least equal to the number of model parameters to avoid singularity of the FIM.
#' @param n.grid Only required when the parameter space is
#' going to be discretized.
#' The total number of grid points from the parameter space is \code{n.grid^p}.
#'  When \code{n.grid > 0}, optimal design protects the experimenter against the worst case scenario only over the grid points, and not over the continuous parameter space.
#'   The resulting designs may not be globally optimal.
#'   In some literature, this type of designs has been used as a compromise to
#'   the minimax type designs to avoid continuous optimization problem over
#'    the parameter space and simplify the minimax design problems.
#'     Especially when the design criterion is convex with respect to the
#'   given parameter space at every given design from the design space,
#'    the obtained design may also be globally optimal (because the maximum of a convex function is attained on the bounds, and here, are included in the grid points).
#'   See 'Details' of \code{\link{minimax}}.
#' @param fimfunc A function. Returns the FIM as a \code{matrix}. Required when \code{formula} is missing. See 'Details' of \code{\link{minimax}}.
#' @param ICA.control ICA control parameters. For details, see \code{\link{ICA.control}}.
#' @param sens.minimax.control Control parameters to construct the answering set required for verify the general equivalence theorem and calculating the ELB. For details, see the function \code{\link{sens.minimax.control}}.
#' @param sens.control Control Parameters for Calculating the ELB. For details, see \code{\link{sens.control}}.
#' @param crt.minimax.control Control parameters to optimize the minimax or standardized maximin criterion at a given design over a \strong{continuous} parameter space (when \code{n.grid = 0}).
#'  For details, see the function \code{\link{crt.minimax.control}}.
#' @param standardized  Maximin standardized design? When \code{standardized = TRUE}, the argument \code{localdes} must be given.
#'  Defaults to \code{FALSE}. See 'Details' of \code{\link{minimax}}.
#' @param initial A matrix of the  initial design points and weights that will be inserted into the initial solutions (countries) of the algorithm.
#'  Every row is a design, i.e.  a concatenation of \code{x} and \code{w}. Will be coerced to a \code{matrix} if necessary.  See 'Details' of \code{\link{minimax}}.
#' @param localdes A function that takes the parameter values  as inputs and returns the design points and weights of the locally optimal design.
#'  Required when \code{standardized = "TRUE"}. See 'Details' of \code{\link{minimax}}.
#' @param npar Number of model parameters.  Used when \code{fimfunc} is given instead of \code{formula} to specify the number of model parameters.
#'   If not specified truly, the sensitivity (derivative) plot may be shifted below the y-axis. When \code{NULL} (default), it is set to \code{length(lp)}.
#' @param plot_3d Which package should be used to plot the sensitivity (derivative) function for two-dimensional design space. Defaults to \code{"lattice"}.
#' @param x A vector of candidate design (support) points.
#' When is not set to \code{NULL} (default),
#'  the algorithm only finds the optimal weights for the candidate points in  \code{x}.
#'    Should be set when the user has a finite number of candidate design points  and the purpose
#'    is to find the optimal weight for each of them (when zero, they will be excluded from the design).
#' For design points with more than one dimension, see 'Details' of \code{\link{sensminimax}}.
#' @param crtfunc (Optional) a function that specifies an arbitrary criterion. It must have especial arguments and output. See 'Details' of \code{\link{minimax}}.
#' @param sensfunc (Optional) a function that specifies the sensitivity function for \code{crtfunc}. See 'Details' of \code{\link{minimax}}.
#' @details
#'
#' Let \eqn{\Xi} be the space of all  approximate designs with
#'  \eqn{k} design points (support points) at \eqn{x_1, x_2, ...,  x_k}
#'  from  the design space \eqn{\chi} with
#'  corresponding weights  \eqn{w_1, . . . ,w_k}.
#'  Let \eqn{M(\xi, \theta)} be the Fisher information
#'   matrix (FIM) of a \eqn{k-}point design \eqn{\xi} and \eqn{\theta} be the vector of
#'    unknown parameters.
#'  A  minimax D-optimal design \eqn{\xi^*}{\xi*} minimizes over \eqn{\Xi}
#'   \deqn{\max_{\theta \in \Theta} -\log|M(\xi, \theta)|.}{
#'    max over \Theta -log|M(\xi, \theta)|.}
#'
#'  A standardized maximin D-optimal design \eqn{\xi^*}{\xi*} maximizes over \eqn{\Xi}
#'  \deqn{\inf_{\theta \in \Theta} \left[\left(\frac{|M(\xi, \theta)|}{|M(\xi_{\theta}, \theta)|}\right)^\frac{1}{p}\right],}{
#'   inf over \Theta {|M(\xi, \theta)| / |M(\xi_\theta, \theta)|}^p,}
#'   where \eqn{p} is the number of model parameters and \eqn{\xi_\theta} is the locally D-optimal design with respect to \eqn{\theta}.\cr
#'
#' A minimax criterion (cost function or objective function) is evaluated at each design (decision variables) by maximizing the criterion over the parameter space.
#' We call the optimization problem over the parameter space as \emph{inner optimization problem}.
#' Two different strategies may be
#'  applied to solve the inner problem at a given design (design points and weights):
#' \enumerate{
#' \item \strong{Continuous inner problem}: we optimize the criterion
#'  over a continuous parameter space using the function \code{\link[nloptr]{nloptr}}.
#' In this case, the tuning parameters can be regulated
#'  via the argument \code{\link{crt.minimax.control}}, when the most influential one
#'  is \strong{\code{maxeval}}.
#' \item  \strong{Discrete inner problem}: we map the parameter space to
#'  the grid points and optimize the criterion over a discrete parameter space.
#' In this case, the number of grid points can be regulated via \code{n.grid}.
#' This strategy is quite efficient (ans fast) when  the maxima most likely attain the vertices of the continuous parameter space at any given design.
#' The output design here protects the experiment from the worst scenario
#' over the grid points.
#' }
#'
#'
#' The \code{formula} is used to automatically create the Fisher information matrix (FIM)
#'  for a linear or nonlinear model provided that the distribution of the
#'   response variable belongs to the natural exponential family.
#' Function \code{minimax} also provides an
#'  option to assign a user-defined FIM directly via the argument  \code{fimfunc}.
#'  In this case, the
#' argument \code{fimfunc} takes a \code{function} that has three arguments as follows:
#'  \enumerate{
#'   \item \code{x} a vector of design points. For design points with more than one dimension,
#'    it is a concatenation of the design points, but dimension-wise.
#'    For example, let the model has three predictors   \eqn{(I, S, Z)}.
#'     Then,  a two-point design is of the format
#'    \eqn{\{\mbox{point1} = (I_1, S_1, Z_1), \mbox{point2} = (I_2, S_2, Z_2)\}}{{point1 = (I1, S1, Z1), point2 = (I2, S2, Z2)}}.
#'     and the argument \code{x} is equivalent to
#'     \code{x = c(I1, I2, S1, S2, Z1, Z2)}.
#'   \item \code{w} a vector that includes the design weights associated with \code{x}.
#'   \item \code{param} a vector of parameter values associated with \code{lp} and \code{up}.
#' }
#'  The output must be the Fisher information matrix with number of rows equal to \code{length(param)}. See 'Examples'.
#'
#'
#'
#'  Minimax optimal designs can have very different criterion values depending on the nominal set of parameter values.
#'  Accordingly, it is desirable to standardize the criterion and control for the potentially widely varying magnitude of the criterion (Dette, 1997).
#'  Evaluating a standardized maximin criterion requires knowing locally optimal designs.
#' We strongly advise setting \code{standardized = 'TRUE'} only when analytical solutions for the locally D-optimal designs is available.
#' In this case, for any initial estimate of the unknown parameters,
#' an analytical solution  for the locally optimal design, i.e, the support points \code{x}
#'  and the corresponding weights \code{w}, must be provided via the argument \code{localdes}.
#' Therefore, depending on how the model is specified, \code{localdes} is a \code{function} with the following arguments (input).
#' \itemize{
#' \item If \code{formula} is given (\code{!missing(formula)}):
#'     \itemize{
#'     \item The parameter names given by \code{parvars} in the same order.
#'
#'     }
#' \item If FIM is given via the argument \code{fimfunc} (\code{missing(formula)}):
#'    \itemize{
#'    \item \code{param}: A vector of the parameters equal to the argument \code{param} in \code{fimfunc}.
#'     }
#' }
#' This function must return a list with the components \code{x} and \code{w} (they match the same arguments in the function \code{fimfunc}).
#'   See 'Examples'.\cr
#'
#' The standardized D-criterion is equal to the  D-efficiency and it must be between 0 and 1.
#'  However, in practice, when running the algorithm,  it may be the case that
#'  the criterion takes a value larger than one.
#'   This may happen  because the user-function that is  given via \code{localdes} does not
#'   return the true (accurate) locally optimal designs for some
#'    requested initial estimates of the parameters  from \eqn{\Theta}.
#'  In this case, the function \code{minimax}
#' throw an error where the error message helps the user
#'  to debug her/his function.
#'
#'
#' Each row of \code{initial} is one design, i.e. a concatenation of values for design (support) points  and the associated design weights.
#' Let \code{x0} and \code{w0} be the vector of initial values with exactly the same length and order as \code{x} and \code{w} (the arguments of \code{fimfunc}).
#'  As an example, the first row  of the matrix \code{initial} is equal to \code{initial[1, ] = c(x0, w0)}.
#'   For models with more than one predictors, \code{x0} is a concatenation of the initial values for design points, but \strong{dimension-wise}.
#'   See the details of the argument \code{fimfunc}, above.
#'
#'  To verify the optimality of the output design by the general equivalence theorem,
#'  the user can either \code{plot} the results or set  \code{checkfreq} in \code{\link{ICA.control}}
#'  to \code{Inf}. In either way, the function \code{\link{sensminimax}} is called for verification.
#'   Note that  the function \code{\link{sensminimax}} always verifies the optimality of a design assuming a continues parameter space.
#' See 'Examples'.
#'
#'  \code{crtfunc} is a function that is used
#'  to specify a new criterion.
#'   Its arguments are:
#'   \itemize{
#'  \item design points \code{x} (as a \code{vector}).
#'  \item design weights \code{w} (as a \code{vector}).
#'  \item model parameters as follows.
#'      \itemize{
#'          \item If \code{formula} is specified:
#'                they should be the same parameter specified by \code{parvars}.
#'          \item If FIM is specified via the argument \code{fimfunc}:
#'              \code{param} that is a vector of the parameters in \code{fimfunc}.
#'               }
#'        \item \code{fimfunc} is a \code{function} that takes the other arguments of \code{crtfunc}
#'        and returns the computed Fisher information matrix as a \code{matrix}.
#'      }
#'    The \code{crtfunc} function must return the criterion value.
#'     \code{crtfunc}. It has one more argument than  \code{crtfunc},
#'      which is  \code{xi_x}. It denotes the design point of the degenerate design
#'      and  must be a vector with the same length as the number of  predictors.
#'     For more details, see 'Examples'.
#'
#'
#' @return
#'
#'  an object of class \code{minimax} that is a list including three sub-lists:
#' \describe{
#'   \item{\code{arg}}{A list of design and algorithm parameters.}
#'   \item{\code{evol}}{A list of length equal to the number of iterations that stores
#'    the information about the best design (design with least criterion value)
#'     of each iteration. \code{evol[[iter]]} contains:
#'     \tabular{lll}{
#'       \code{iter}                   \tab      \tab Iteration number.\cr
#'       \code{x}                      \tab      \tab Design points. \cr
#'       \code{w}                      \tab      \tab Design weights. \cr
#'       \code{min_cost}               \tab      \tab Value of the criterion for the best imperialist (design).  \cr
#'       \code{mean_cost}              \tab      \tab Mean of the criterion values of all the imperialists. \cr
#'       \code{sens}                   \tab      \tab An object of class \code{'sensminimax'}. See below. \cr
#'       \code{param}                  \tab      \tab Vector of parameters.\cr
#'     }
#'   }
#'
#'   \item{\code{empires}}{A list of all the  empires of the last iteration.}
#'   \item{\code{alg}}{A list with following information:
#'     \tabular{lll}{
#'       \code{nfeval}           \tab      \tab Number of function evaluations.  It does not count the function evaluations from checking the general equivalence theorem.\cr
#'       \code{nlocal}           \tab      \tab Number of successful local searches. \cr
#'       \code{nrevol}           \tab      \tab Number of successful revolutions. \cr
#'       \code{nimprove}         \tab      \tab Number of successful movements toward the imperialists in the assimilation step. \cr
#'       \code{convergence}      \tab      \tab Stopped by \code{'maxiter'} or \code{'equivalence'}?\cr
#'     }
#'   }
#'   \item{\code{method}}{A type of optimal designs used.}
#'   \item{\code{design}}{Design points and weights at the final iteration.}
#'   \item{\code{out}}{A data frame of design points, weights, value of the criterion for the best imperialist (min_cost), and Mean of the criterion values of all the imperialistsat each iteration (mean_cost).}
#' }
#'
#' The list \code{sens} contains information about the design verification by the general equivalence theorem. See \code{sensminimax} for more details.
#' It is given  every \code{ICA.control$checkfreq} iterations
#' and also the last iteration if   \code{ICA.control$checkfreq >= 0}. Otherwise, \code{NULL}.
#'
#'  \code{param} is a vector of parameters that is the global minimum of
#'   the minimax criterion or the global maximum of the standardized maximin criterion over the parameter space, given  the current \code{x}, \code{w}.
#'
#' @note
#' For larger parameter space or model with more number of unknown parameters,
#'  it is always important to increase the value of  \code{ncount} in \code{ICA.control}
#' and \code{optslist$maxeval} in \code{crt.minimax.control} to produce very accurate designs.
#'
#' Although standardized criteria have been preferred theoretically,
#' in practice, they should be applied only
#'   when  an analytical solution for
#'   the locally D-optimal designs is available for the model
#'   of interest.
#'    Otherwise, we encounter a three-level nested-optimization algorithm, which is very slow.
#'
#' @references
#' Atashpaz-Gargari, E, & Lucas, C (2007). Imperialist competitive algorithm: an algorithm for optimization inspired by imperialistic competition. In 2007 IEEE congress on evolutionary computation (pp. 4661-4667). IEEE.\cr
#' Dette, H. (1997). Designing experiments with respect to 'standardized' optimality criteria. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 59(1), 97-110. \cr
#' Masoudi E, Holling H, Wong WK (2017). Application of Imperialist Competitive Algorithm to Find Minimax and Standardized Maximin Optimal Designs. Computational Statistics and Data Analysis, 113, 330-345. <doi:10.1016/j.csda.2016.06.014>\cr
#' @example inst/examples/minimax_examples.R
#' @importFrom utils capture.output
#' @export
#' @seealso \code{\link{sensminimax}}
minimax <- function(formula, predvars, parvars, family = gaussian(),
                    lx, ux, lp, up, iter, k,
                    n.grid = 0,
                    fimfunc = NULL,
                    ICA.control = list(),
                    sens.control = list(),
                    sens.minimax.control = list(),
                    crt.minimax.control = list(),
                    standardized = FALSE,
                    initial = NULL,
                    localdes = NULL,
                    npar = length(lp),
                    plot_3d = c("lattice", "rgl"),
                    x = NULL,
                    crtfunc = NULL,
                    sensfunc = NULL){
  ### how to control ICA.control sens.minimax.control
  if (!is.numeric(n.grid) || n.grid < 0)
    stop("value of 'n.grid' must be >= 0")
  # if (!is.numeric(maxeval) || param_maxeval <= 0)
  #   stop("value of 'param_maxeval' must be > 0")
  if (!is.logical(standardized))
    stop("'standardized' must be logical")
  if (standardized)
    type <- "standardized" else
      type <- "minimax"

    # The number of the grid points used is length(lp)^n.grid

    #n.grid^length(lp) - (n.grid -2)^length(lp)

    crt.minimax.control <- do.call("crt.minimax.control", crt.minimax.control)
    # check if the length of x0 be equal to the length of lp and up!!
    #minimax.control$inner_maxeval <-  param_maxeval

    if (n.grid>0){
      crt.minimax.control$param_set <- make_grid(lp = lp, up = up, n.grid = n.grid)
      crt.minimax.control$inner_space <- "discrete"
    }else
      crt.minimax.control$inner_space <- "continuous"
    if (is.null(crtfunc))
      crt_type = "D" else
        crt_type = "user"
    if (!is.null(x))
      k <- length(x)/length(lx)


        out <- minimax_inner(formula = formula,
                         predvars = predvars, parvars = parvars, family = family,
                         lx = lx, ux = ux, lp = lp, up = up, iter = iter,
                         k = k,
                         fimfunc = fimfunc,
                         ICA.control = ICA.control,
                         sens.control = sens.control,
                         sens.minimax.control = sens.minimax.control,
                         crt.minimax.control = crt.minimax.control,
                         type = type,
                         initial = initial,
                         localdes = localdes,
                         npar = npar,
                         crt_type = crt_type,
                         multipars = list(),
                         plot_3d = plot_3d[1],
                         only_w_varlist = list(x = x),
                         user_crtfunc = crtfunc,
                         user_sensfunc = sensfunc)

        out$method = 'minimax' # 06202020@seongho
        if (!missing(formula)){
      out$call = formula # 06202020@seongho
    }else{
      out$call = NULL
    }

    dout = NULL
    fout = NULL
    fiter = length(out$evol)
    if(fiter>0){
      tout = c()
      for(i in 1:fiter){
        sout = out$evol[[i]]
        tout = rbind(tout,c(i,sout$x,sout$w,sout$min_cost,sout$mean_cost))
      }
      dout = as.data.frame(tout)
      #dimnames(dout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost")

      if(length(sout$x)==length(sout$w)){
        dimnames(dout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost")
      }else{
        tlen = length(sout$x)/length(sout$w)
        tdimx = c()
        for(i in 1:tlen){
          tdimx = c(tdimx,paste(paste("x",i,sep=""),1:k,sep=""))
        }
        dimnames(dout)[[2]] = c("iter",tdimx,paste("w",1:k,sep=""),"min_cost","mean_cost")
      }

      # extract sens outcomes
      sout = out$evol[[fiter]]
      tsens = capture.output(sout$sens)
      tpos.max = grep("Maximum",tsens)
      tpos.elb = grep("ELB",tsens)
      tpos.time = grep("Verification",tsens)
      tmax = NA
      telb = NA
      ttime = NA
      if(length(tpos.max)>0){
        tmax = as.numeric(strsplit(gsub("\\s+","",tsens[tpos.max]),"is")[[1]][2])
      }
      if(length(tpos.elb)>0){
        telb = as.numeric(strsplit(gsub("\\s+","",tsens[tpos.elb]),"is")[[1]][2])
      }
      if(length(tpos.time)>0){
        ttime = as.numeric(strsplit(gsub("seconds!","",gsub("\\s+","",tsens[tpos.time])),"required")[[1]][2])
      }
      t2out = c(fiter,sout$x,sout$w,sout$min_cost,sout$mean_cost,tmax,telb,ttime)
      fout = as.data.frame(t(t2out))
      #dimnames(fout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")

      if(length(sout$x)==length(sout$w)){
        dimnames(fout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
      }else{
        tlen = length(sout$x)/length(sout$w)
        tdimx2 = c()
        for(i in 1:tlen){
          tdimx2 = c(tdimx2,paste(paste("x",i,sep=""),1:k,sep=""))
        }
        dimnames(fout)[[2]] = c("iter",tdimx2,paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
      }

    }

    out$out = dout
    #out$out2 = out$evol
    out$design = fout

    return(out)
}

######################################################################################################*
######################################################################################################*
#' @title Verifying Optimality of The Minimax and Standardized maximin D-optimal Designs
#'
#' @description
#'
#' It plots the sensitivity (derivative) function of the minimax or
#'  standardized maximin D-optimal criterion
#' at a given approximate (continuous) design and also
#'  calculates its efficiency lower bound (ELB) with respect
#' to the optimality criterion.
#' For an approximate (continuous) design, when the design space is one or two-dimensional,
#'  the user can visually verify the optimality of the design by observing the
#' sensitivity plot. Furthermore, the proximity of the design to the optimal design
#'  can be measured by the  ELB without knowing the latter.
#'  See, for more details, Masoudi et al. (2017).
#'
#' @inheritParams minimax
#' @param x Vector of the design (support) points. See 'Details' of \code{\link{sensminimax}} for models with more than one predictors.
#' @param w Vector of the corresponding design weights for \code{x}.
#' @param calculate_criterion Calculate the optimality criterion? See 'Details' of \code{\link{sensminimax}}.
#' @param plot_3d Which package should be used to plot the sensitivity (derivative) function for models with two predictors.
#'   Either \code{"rgl"} or \code{"lattice"} (default).
#' @param plot_sens Plot the sensitivity (derivative) function? Defaults to \code{TRUE}.
#' @param crt.minimax.control Control parameters to calculate the value of the  minimax or standardized maximin optimality criterion  over the  continuous parameter space.
#'   Only applicable when \code{calculate_criterion = TRUE}.
#'   For more details, see \code{\link{crt.minimax.control}}.
#' @param silent Do not print anything? Defaults to \code{FALSE}.
#' @param crtfunc (Optional) a function that specifies an arbitrary criterion. It must have especial arguments and output. See 'Details' of \code{\link{minimax}}.
#' @param sensfunc (Optional) a function that specifies the sensitivity function for \code{crtfunc}. See 'Details' of \code{\link{minimax}}.
#' @details
#' Let the unknown parameters belong to \eqn{\Theta}.
#' A design \eqn{\xi^*}{\xi*} is minimax D-optimal among all designs on \eqn{\chi} if and only if there exists a probability measure \eqn{\mu^*}{\mu*} on
#'    \deqn{A(\xi^*) = \left\{\nu \in \Theta \mid -log|M(\xi^*, \nu)| = \max_{\theta \in \Theta} -log|M(\xi^*, \theta)| \right\},}{
#'      A(\xi*) = {\nu belongs to \Theta where -log|M(\xi*, \nu) = maxima of function -log|M(\xi*, \theta)| with respect to \theta over \Theta} ,}
#'        such that the following inequality holds for all \eqn{\boldsymbol{x} \in \chi}{x belong to \chi}
#'         \deqn{c(\boldsymbol{x}, \mu^*, \xi^*) = \int_{A(\xi^*)} tr M^{-1}(\xi^*, \nu)I(\boldsymbol{x}, \nu)\mu^* d(\nu)-p \leq 0,}{
#'          c(x, \mu*, \xi*) = integration over A(\xi*) with integrand  tr M^-1(\xi*, \nu)I(x, \nu)\mu* d(\nu)-p <= 0,}
#'           with equality at all support points of \eqn{\xi^*}{\xi*}.
#'            Here, \eqn{p} is the number of model parameters. \eqn{c(\boldsymbol{x}, \mu^*, \xi^*)}{c(x, \mu*, \xi*)} is called \strong{sensitivity} or \strong{derivative} function.
#' The set \eqn{A(\xi^*)}{A(\xi*)} is sometimes called \bold{answering set} of
#'  \eqn{\xi^*}{\xi*} and the measure \eqn{\mu^*}{\mu*} is a sub-gradient of the
#'    non-differentiable criterion evaluated at \eqn{M(\xi^*,\nu)}{M(\xi*,\nu)}.\cr
#' For the standardized maximin D-optimal designs, the answering set \eqn{N(\xi^*)}{N(\xi*)} is
#'    \deqn{N(\xi^*) = \left\{\boldsymbol{\nu} \in \Theta \mid \mbox{eff}_D(\xi^*, \boldsymbol{\nu}) = \min_{\boldsymbol{\theta} \in \Theta} \mbox{eff}_D(\xi^*, \boldsymbol{\theta}) \right\}.
#'      }{N(\xi*) = \{\nu belongs to \Theta  where eff_D(\xi*, \nu) = minima of function eff_D(\xi*, \theta)  with respect to \theta over \Theta\},} where
#'      \eqn{\mbox{eff}_D(\xi, \boldsymbol{\theta}) =  (\frac{|M(\xi, \boldsymbol{\theta})|}{|M(\xi_{\boldsymbol{\theta}}, \boldsymbol{\theta})|})^\frac{1}{p}}{
#'      eff_D(\xi, \theta) =  (|M(\xi, \theta)|/|M(\xi_\theta, \theta)|)^(1/p)} and \eqn{\xi_\theta} is the locally D-optimal design with respect to \eqn{\theta}.
#'      See 'Details' of \code{\link{sens.minimax.control}} on how we find the answering set.
#'
#' The argument  \code{x} is the vector of design points.
#'  For design points with more than one dimension (the models with more than one predictors),
#'    it is a concatenation of the design points, but \strong{dimension-wise}.
#'    For example, let the model has three predictors   \eqn{(I, S, Z)}.
#'     Then,  a two-point optimal design has the following points:
#'    \eqn{\{\mbox{point1} = (I_1, S_1, Z_1), \mbox{point2} = (I_2, S_2, Z_2)\}}{{point1 = (I1, S1, Z1), point2 = (I2, S2, Z2)}}.
#'     Then, the argument \code{x} is equal to
#'     \code{x = c(I1, I2, S1, S2, Z1, Z2)}.
#'
#' ELB is a measure of  proximity of a design to the optimal design without knowing the latter.
#' Given a design, let \eqn{\epsilon} be the global maximum
#'  of the sensitivity (derivative) function with respect \eqn{x} where \eqn{x \in \chi}{x belong to \chi}.
#' ELB is given by \deqn{ELB = p/(p + \epsilon),}
#' where \eqn{p} is the number of model parameters. Obviously,
#' calculating ELB requires finding \eqn{\epsilon} and
#' another optimization problem to be solved.
#' The tuning parameters of this optimization can be regulated via the argument \code{\link{sens.minimax.control}}.
#' See, for more details, Masoudi et al. (2017).
#'
#'   The criterion value for the minimax D-optimal design is (global maximum over \eqn{\Theta})
#'  \deqn{\max_{\theta \in \Theta} -\log|M(\xi, \theta)|;}{max -log|M(\xi, \theta)|;}
#'  for the standardized maximin D-optimal design is (global minimum over \eqn{\Theta})
#'  \deqn{\inf_{\theta \in \Theta} \left[\left(\frac{|M(\xi, \theta)|}{|M(\xi_{\theta}, \theta)|}\right)^\frac{1}{p}\right].}{
#'   inf {|M(\xi, \theta)| / |M(\xi_\theta, \theta)|}^p.}
#'
#'  This function confirms the optimality assuming only a continuous parameter space \eqn{\Theta}.
#'
#' @return
#'  an object of class \code{sensminimax} that is a list with the following elements:
#'  \describe{
#'  \item{\code{type}}{Argument \code{type} that is required for print methods.}
#'  \item{\code{optima}}{A \code{matrix} that stores all the local optima over the parameter space.
#'   The cost  (criterion) values are stored in a column named \code{Criterion_Value}.
#'  The last column (\code{Answering_Set})
#'   shows if the optimum belongs to the answering set (1) or not (0). See 'Details' of \code{\link{sens.minimax.control}}.
#'    Only applicable for minimax or standardized maximin designs.}
#'  \item{\code{mu}}{Probability measure on the answering set.
#'   Corresponds to the rows of \code{optima} for which the associated row in column \code{Answering_Set} is equal to 1.
#'    Only applicable for minimax or standardized maximin designs.}
#'  \item{\code{max_deriv}}{Global maximum of the sensitivity (derivative) function (\eqn{\epsilon} in 'Details').}
#'  \item{\code{ELB}}{D-efficiency lower bound. Can not be larger than 1. If negative, see 'Note' in \code{\link{sensminimax}} or  \code{\link{sens.minimax.control}}.}
#'  \item{\code{merge_tol}}{Merging tolerance to create the answering set from the set of all local optima. See 'Details' in \code{\link{sens.minimax.control}}.
#'   Only applicable for minimax or standardized maximin designs.}
#'  \item{\code{crtval}}{Criterion value. Compare it with the column \code{Crtiterion_Value} in \code{optima} for minimax and standardized maximin designs.}
#'  \item{\code{time}}{Used CPU time (rough approximation).}
#'  }
#' @note
#' Theoretically, ELB can not be larger than 1. But if so, it may have one of the following reasons:
#' \itemize{
#' \item \code{max_deriv} is not a GLOBAL maximum.  Please increase  the value of the parameter \code{maxeval} in \code{\link{sens.minimax.control}} to find the global maximum.
#' \item The sensitivity function is shifted below the y-axis because
#' the number of model parameters has not been specified correctly (less value given).
#' Please specify the correct number of model parameters via argument \code{npar}.
#' }
#'  Please increase the value of the parameter
#'   \code{n_seg} in \code{\link{sens.minimax.control}}
#'    for  models with larger number of parameters or large parameter space to find the true
#'    answering set for minimax and standardized maximin designs. See \code{\link{sens.minimax.control}} for more details.
#' @references
#' Masoudi E, Holling H, Wong W.K. (2017). Application of Imperialist Competitive Algorithm to Find Minimax and Standardized Maximin Optimal Designs. Computational Statistics and Data Analysis, 113, 330-345. \cr
#' @example inst/examples/sensminimax_examples.R
#' @export
sensminimax <- function (formula, predvars, parvars,
                         family = gaussian(),
                         x, w,
                         lx, ux,
                         lp, up,
                         fimfunc = NULL,
                         standardized = FALSE,
                         localdes = NULL,
                         sens.control = list(),
                         sens.minimax.control = list(),
                         calculate_criterion = TRUE,
                         crt.minimax.control = list(),
                         plot_3d = c("lattice", "rgl"),
                         plot_sens = TRUE,
                         npar = length(lp),
                         silent = FALSE,
                         crtfunc = NULL,
                         sensfunc = NULL
){

  if (!is.logical(standardized))
    stop("'standardized' must be logical")
  if (standardized)
    type <- "standardized" else
      type <- "minimax"

    if(calculate_criterion){
      crt.minimax.control <- do.call("crt.minimax.control", crt.minimax.control)
      crt.minimax.control$inner_space <- "continuous"
    }
    if (is.null(crtfunc))
      crt_type = "D" else
        crt_type = "user"
      # check if the length of x0 be equal to the length of lp and up!!
      #minimax.control$inner_maxeval <-  param_maxeval

      # if (n.grid>0){
      #   crt.minimax.control$param_set <- make_grid(lp = lp, up = up, n.grid = n.grid)
      #   crt.minimax.control$inner_space <- "discrete"
      # }else


      out <- sensminimax_inner(formula = formula, predvars = predvars, parvars = parvars,
                               family = family,
                               x = x, w = w,
                               lx = lx, ux = ux,
                               lp = lp, up = up,
                               fimfunc = fimfunc,
                               sens.minimax.control =  sens.minimax.control,
                               sens.control = sens.control,
                               type = type,
                               localdes = localdes,
                               plot_3d = plot_3d[1],
                               plot_sens = plot_sens,
                               varlist = list(),
                               calledfrom = "sensfuncs",
                               npar = npar,
                               crt.minimax.control = crt.minimax.control,
                               calculate_criterion = calculate_criterion,
                               crt_type = crt_type,
                               multipars = list(),
                               silent = silent,
                               user_crtfunc = crtfunc,
                               user_sensfunc = sensfunc)
      out$method = "minimax" # 06222020@seongho

      return(out)
}
######################################################################################################*
######################################################################################################*
#' @title Locally D-Optimal Designs
#' @description
#'  Finds locally D-optimal designs for linear and nonlinear models.
#'  It should be used when a vector of initial estimates is available for the unknown model parameters.
#'  Locally optimal designs may not be efficient when the initial estimates are  far away from the true values of the parameters.
#' @inheritParams minimax
#' @param inipars A vector of initial estimates for the unknown parameters.
#' It must match \code{parvars} or the argument \code{param} of the function \code{fimfunc}, when provided.
#' @param npar Number of model parameters.  Used when \code{fimfunc} is given instead of \code{formula} to specify the number of model parameters.
#'   If not given, the sensitivity plot may be shifted below the y-axis.
#'   When \code{NULL}, it is set to \code{length(inipars)}.
#' @importFrom utils capture.output
#' @export
#' @seealso \code{\link{senslocally}}
#' @details
#'  Let \eqn{M(\xi, \theta_0)} be the Fisher information
#'   matrix (FIM) of a \eqn{k-}point design \eqn{\xi} and
#'    \eqn{\theta_0}  be the vector of the initial  estimates for the unknown parameters.
#'  A locally D-optimal design \eqn{\xi^*}{\xi*} minimizes over \eqn{\Xi}
#'   \deqn{-\log|M(\xi, \theta_0)|.}{-log|M(\xi, \theta_0)|.}
#'
#' One can adjust the tuning parameters in \code{\link{ICA.control}} to set a stopping rule
#' based on the general equivalence theorem. See "Examples" below.
#' @inherit minimax return
#'
#'
#' @references
#' Masoudi E, Holling H, Wong W.K. (2017). Application of Imperialist Competitive Algorithm to Find Minimax and Standardized Maximin Optimal Designs. Computational Statistics and Data Analysis, 113, 330-345. \cr
#' @example inst/examples/locally_examples.R
locally <- function(formula, predvars, parvars, family = gaussian(),
                    lx, ux,  iter, k,
                    inipars,
                    fimfunc = NULL,
                    ICA.control = list(),
                    sens.control = list(),
                    initial = NULL,
                    npar = length(inipars),
                    plot_3d = c("lattice", "rgl"),
                    x = NULL,
                    crtfunc = NULL,
                    sensfunc = NULL){


  if (!missing(formula)){
    if (length(inipars) != length(parvars))
      stop("lengtb of 'inipars' is not equal to the length of 'parvars'")
  }
  if (is.null(x))
    if (k < length(inipars))
      stop("\"k\" must be larger than the number of parameters to avoid singularity")
  if (!is.null(x))
    k <- length(x)/length(lx) # Bug!!! how you can understand the number of

  # if (!is.null(x)){
  #   #dim(x)[2] is the dimension of the predictors
  #   #dim(x)[1] is the number of the weights
  #
  #   ## lower bound and upper bound of the design space when only
  #   ## w is requested. Because lx and ux will be replcaed by 0 and 1,
  #   # when the x is given (w is only the objective variables)
  #   # only applied for equivalence theorem
  #   #lx_when_only_w <- rep(lx, length(x)/k)
  #   #ux_when_only_w <- rep(ux, length(x)/k)
  #   #lx_when_only_w <- lx
  #   #ux_when_only_w <- ux
  #   ## here lx and ux are not the lowerbound and upperbond of x, but w
  #   #lx <- rep(0, length(x)/k)
  #   #ux <- rep(1, length(x)/k)
  # }
  if (is.null(crtfunc))
    crt_type = "D" else
      crt_type = "user"


    #cat("#### BEFORE OUT ####\n")

    out <- minimax_inner(formula = formula,
                         predvars = predvars, parvars = parvars, family = family,
                         lx = lx, ux = ux, lp = inipars, up = inipars, iter = iter, k = k,
                         fimfunc = fimfunc,
                         ICA.control = ICA.control,
                         sens.control = sens.control,
                         crt.minimax.control = list(inner_space = "locally"),
                         type = "locally",
                         initial = initial,
                         localdes = NULL,
                         npar = npar,
                         robpars = list(),
                         crt_type = crt_type,
                         multipars = list(),
                         plot_3d = plot_3d[1],
                         only_w_varlist = list(x = x),
                         user_crtfunc = crtfunc,
                         user_sensfunc = sensfunc)

    out$method = 'locally' # 06202020@seongho
    if (!missing(formula)){
      out$call = formula # 06202020@seongho
    }else{
      out$call = NULL
    }

    #cat("#### AFTER OUT ####\n")

    dout = NULL
    fout = NULL
    fiter = length(out$evol)
    if(fiter>0){

      #cat("############-1-1-1-1-1-1- AT THE END OF OUT ##############\n")

      tout = c()
      for(i in 1:fiter){
        sout = out$evol[[i]]
        tout = rbind(tout,c(i,sout$x,sout$w,sout$min_cost,sout$mean_cost))
      }
      dout = as.data.frame(tout)
      #dimnames(dout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost")

      #cat("############----------11111111111100000 AT THE END OF OUT ##############\n")


      if(length(sout$x)==length(sout$w)){
        dimnames(dout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost")
      }else{
        tlen = length(sout$x)/length(sout$w)
        tdimx = c()
        for(i in 1:tlen){
          tdimx = c(tdimx,paste(paste("x",i,sep=""),1:k,sep=""))
        }
        dimnames(dout)[[2]] = c("iter",tdimx,paste("w",1:k,sep=""),"min_cost","mean_cost")
      }

      #cat("############----------11111111111122222222222200000 AT THE END OF OUT ##############\n")
      #cat("############ fiter = ",fiter," ############# -----\n")
      #print(out$evol[[fiter]])

      # extract sens outcomes
      sout = out$evol[[fiter]]
      tsens = capture.output(sout$sens)
      tpos.max = grep("Maximum",tsens)
      tpos.elb = grep("ELB",tsens)
      tpos.time = grep("Verification",tsens)
      tmax = NA
      telb = NA
      ttime = NA
      if(length(tpos.max)>0){
        tmax = as.numeric(strsplit(gsub("\\s+","",tsens[tpos.max]),"is")[[1]][2])
      }
      if(length(tpos.elb)>0){
        telb = as.numeric(strsplit(gsub("\\s+","",tsens[tpos.elb]),"is")[[1]][2])
      }
      if(length(tpos.time)>0){
        ttime = as.numeric(strsplit(gsub("seconds!","",gsub("\\s+","",tsens[tpos.time])),"required")[[1]][2])
      }
      t2out = c(fiter,sout$x,sout$w,sout$min_cost,sout$mean_cost,tmax,telb,ttime)
      fout = as.data.frame(t(t2out))
      #dimnames(fout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")

      #cat("############----------1111111111112222222222223333333333333300000 AT THE END OF OUT ##############\n")


      if(length(sout$x)==length(sout$w)){
        dimnames(fout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
      }else{
        tlen = length(sout$x)/length(sout$w)
        tdimx2 = c()
        for(i in 1:tlen){
          tdimx2 = c(tdimx2,paste(paste("x",i,sep=""),1:k,sep=""))
        }
        dimnames(fout)[[2]] = c("iter",tdimx2,paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
      }

      #cat("############----------00000 AT THE END OF OUT ##############\n")

    }

    #cat("############00000 AT THE END OF OUT ##############\n")

    out$out = dout
    #out$out2 = out$evol
    out$design = fout

    #cat("############ AT THE END OF OUT ##############\n")

    return(out)
}

######################################################################################################*
######################################################################################################*
#' @title Verifying Optimality of The Locally D-optimal Designs
#' @description
#' It plots the sensitivity (derivative) function of the
#'  locally D-optimal criterion
#' at a given approximate (continuous) design and also
#'  calculates its efficiency lower bound (ELB) with respect
#' to the optimality criterion.
#' For an approximate (continuous) design, when the design space is one or two-dimensional,
#'  the user can visually verify the optimality of the design by observing the
#' sensitivity plot. Furthermore, the proximity of the design to the optimal design
#'  can be measured by the  ELB without knowing the latter.
#'  See, for more details, Masoudi et al. (2017).
#' @inheritParams sensminimax
#' @inheritParams locally
#' @param npar Number of model parameters.  Used when \code{fimfunc} is given instead of \code{formula} to specify the number of model parameters.
#'   If not given, the sensitivity plot may be shifted below the y-axis.
#'   When \code{NULL}, it is set to \code{length(inipars)}.
#' @example inst/examples/senslocally_examples.R
#' @inherit sensminimax return
#' @export
#' @details
#'
#' Let \eqn{\theta_0} denotes the vector of initial estimates for the unknown parameters.
#' A design \eqn{\xi^*}{\xi*} is locally D-optimal among all designs on \eqn{\chi} if and only if
#'  the following inequality holds for all \eqn{\boldsymbol{x} \in \chi}{x belong to \chi}
#'         \deqn{c(\boldsymbol{x}, \xi^*, \theta_0) =  tr M^{-1}(\xi^*, \theta_0)I(\boldsymbol{x}, \theta_0)-p \leq 0,}{
#'          c(x, \xi*, \theta_0) =  tr M^-1(\xi*, \theta0)I(x, \theta_0)-p <= 0,}
#'           with equality at all support points of \eqn{\xi^*}{\xi*}.
#'            Here, \eqn{p} is the number of model parameters.
#'             \eqn{c(\boldsymbol{x},\xi^*, \theta_0)}{c(x,\xi*, \theta_0)} is called \strong{sensitivity} or \strong{derivative} function.
#'
#' ELB is a measure of  proximity of a design to the optimal design without knowing the latter.
#' Given a design, let \eqn{\epsilon} be the global maximum
#'  of the sensitivity (derivative) function over \eqn{x \in \chi}{x belong to \chi}.
#' ELB is given by \deqn{ELB = p/(p + \epsilon),}
#' where \eqn{p} is the number of model parameters. Obviously,
#' calculating ELB requires finding \eqn{\epsilon} and
#' another optimization problem to be solved.
#' The tuning parameters of this optimization can be regulated via the argument  \code{\link{sens.minimax.control}}.
#' See, for more details, Masoudi et al. (2017).
#'
#' @note
#'
#' Theoretically, ELB can not be larger than 1. But if so, it may have one of the following reasons:
#' \itemize{
#' \item \code{max_deriv} is not a GLOBAL maximum.  Please increase  the value of the parameter \code{maxeval} in \code{\link{sens.minimax.control}} to find the global maximum.
#' \item The sensitivity function is shifted below the y-axis because
#' the number of model parameters has not been specified correctly (less value given).
#' Please specify the correct number of model parameters via the argument \code{npar}.
#' }
#'
#' @references
#' Masoudi E, Holling H, Wong W.K. (2017). Application of Imperialist Competitive Algorithm to Find Minimax and Standardized Maximin Optimal Designs. Computational Statistics and Data Analysis, 113, 330-345. \cr
senslocally <- function (formula, predvars, parvars,
                         family = gaussian(),
                         x, w,
                         lx, ux,
                         inipars,
                         fimfunc = NULL,
                         sens.control = list(),
                         calculate_criterion = TRUE,
                         plot_3d = c("lattice", "rgl"),
                         plot_sens = TRUE,
                         npar = length(inipars),
                         silent = FALSE,
                         crtfunc = NULL,
                         sensfunc = NULL){


  if (!missing(formula)){
    if (length(inipars) != length(parvars))
      stop("lengtb of 'inipars' is not equal to the length of 'parvars'")
  }
  if (is.null(npar))
    npar <- length(inipars)

  if (is.null(crtfunc))
    crt_type = "D" else
      crt_type = "user"


    out <- sensminimax_inner(formula = formula, predvars = predvars, parvars = parvars,
                             family = family,
                             x = x, w = w,
                             lx = lx, ux = ux,
                             lp = inipars, up = inipars,
                             fimfunc = fimfunc,
                             sens.control = sens.control,
                             type = "locally",
                             localdes = NULL,
                             plot_3d = plot_3d[1],
                             plot_sens = plot_sens,
                             varlist = list(),
                             calledfrom = "sensfuncs",
                             npar = npar,
                             crt.minimax.control = list(inner_space = "locally"),
                             calculate_criterion = calculate_criterion,
                             robpars = list(),
                             crt_type = crt_type,
                             multipars = list(),
                             silent = silent,
                             user_crtfunc = crtfunc,
                             user_sensfunc = sensfunc)
    out$method = "locally" # 06222020@seongho

    return(out)
}

######################################################################################################*
######################################################################################################*
#' @title Robust D-Optimal Designs
#'
#' @description
#' Finds Robust designs or optimal  in-average designs for linear and nonlinear models.
#'  It is useful when a set of different vectors of initial estimates
#'   along with a discrete probability measure
#'   are available for the unknown model parameters.
#'   It is a discrete version of \code{\link{bayes}}.
#' @inheritParams senslocally
#' @inheritParams locally
#' @param prob A vector of the probability measure \eqn{\pi} associated with each row of \code{parset}.
#' @param parset A matrix that provides the vector of initial estimates for the model parameters, i.e. support of \eqn{\pi}.
#'  Every row is one vector  (\code{nrow(parset) == length(prob)}). See 'Details'.
#' @param npar Number of model parameters.  Used when \code{fimfunc} is given instead of \code{formula} to specify the number of model parameters.
#'   If not given, the sensitivity plot may be shifted below the y-axis.
#'    When \code{NULL}, it is set to \code{dim(parset)[2]}.
#' @inherit locally return
#' @details
#' Let \eqn{\Theta} be a set of initial estimates for the unknown parameters.
#' A robust criterion is evaluated at the elements of \eqn{\Theta} weighted by a probability measure
#' \eqn{\pi} as follows:
#' \deqn{B(\xi, \pi) = \int_{\Theta}|M(\xi, \theta)|\pi(\theta) d\theta.}{
#' B(\xi, \Pi) = intergation over \Theta \Psi(\xi, \theta)\pi(\theta) d\theta.}
#' A robust design \eqn{\xi^*}{\xi*}   maximizes \eqn{B(\xi, \pi)} over the space of all designs.
#'
#'  When the model is given via \code{formula},
#'   columns of \code{parset} must match the parameters introduced
#'   in \code{parvars}.
#'   Otherwise, when the model is introduced via \code{fimfunc},
#'   columns of \code{parset} must match the argument \code{param} in \code{fimfunc}.
#'
#'  To verify the optimality of the output design by the general equivalence theorem,
#'  the user can either \code{plot} the results or set  \code{checkfreq} in \code{\link{ICA.control}}
#'  to \code{Inf}. In either way, the function \code{\link{sensrobust}} is called for verification.
#' One can also adjust the tuning parameters in \code{\link{ICA.control}} to set a stopping rule
#' based on the general equivalence theorem. See 'Examples' below.
#'
#' @importFrom utils capture.output
#' @export
#' @note
#' When a continuous prior distribution for the unknown model parameters is available,  use \code{\link{bayes}}.
#' When only one initial estimates of the unknown model parameters is available (\eqn{\Theta} has only one element),  use  \code{\link{locally}}.
#' @seealso \code{\link{bayes}} \code{\link{sensrobust}}
#' @example inst/examples/robust_examples.R
robust <- function(formula, predvars, parvars, family = gaussian(),
                   lx, ux,  iter, k,
                   prob,
                   parset,
                   fimfunc = NULL,
                   ICA.control = list(),
                   sens.control = list(),
                   initial = NULL,
                   npar = dim(parset)[2],
                   plot_3d = c("lattice", "rgl"),
                   x = NULL,
                   crtfunc = NULL,
                   sensfunc = NULL){


  if (length(prob) != dim(parset)[1])
    stop("length of \"prior\" is not equal to the number of rows of \"param\"")
  if (!missing(formula)){
    if (dim(parset)[2] != length(parvars))
      stop("number of columns of  'parset' is not equal to the length of 'parvars'")
  }
  if (k < dim(parset)[2])
    stop("\"k\" must be larger than the number of parameters to avoid singularity")
  if (is.null(crtfunc))
    crt_type = "D" else
      crt_type = "user"
    # if (is.null(npar))
    #   npar <- dim(parset)[2]
    ## you must provide npar
    if (!is.null(x))
      k <- length(x)/length(lx)


        out <- minimax_inner(formula = formula,
                         predvars = predvars, parvars = parvars,
                         family = family,
                         lx = lx, ux = ux, lp = NA, up = NA, iter = iter, k = k,
                         fimfunc = fimfunc,
                         ICA.control = ICA.control,
                         sens.control = sens.control,
                         crt.minimax.control = list(inner_space = "robust_set"),
                         type = "robust",
                         initial = initial,
                         localdes = NULL,
                         npar = npar,
                         robpars = list(prob = prob, parset = parset),
                         crt_type = crt_type,
                         multipars = list(),
                         plot_3d = plot_3d[1],
                         only_w_varlist = list(x = x),
                         user_crtfunc = crtfunc,
                         user_sensfunc = sensfunc)

        out$method = 'robust' # 06202020@seongho
        if (!missing(formula)){
      out$call = formula # 06202020@seongho
    }else{
      out$call = NULL
    }

    dout = NULL
    fout = NULL
    fiter = length(out$evol)
    if(fiter>0){
      tout = c()
      for(i in 1:fiter){
        sout = out$evol[[i]]
        tout = rbind(tout,c(i,sout$x,sout$w,sout$min_cost,sout$mean_cost))
      }
      dout = as.data.frame(tout)
      #dimnames(dout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost")

      if(length(sout$x)==length(sout$w)){
        dimnames(dout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost")
      }else{
        tlen = length(sout$x)/length(sout$w)
        tdimx = c()
        for(i in 1:tlen){
          tdimx = c(tdimx,paste(paste("x",i,sep=""),1:k,sep=""))
        }
        dimnames(dout)[[2]] = c("iter",tdimx,paste("w",1:k,sep=""),"min_cost","mean_cost")
      }

      # extract sens outcomes
      sout = out$evol[[fiter]]
      tsens = capture.output(sout$sens)
      tpos.max = grep("Maximum",tsens)
      tpos.elb = grep("ELB",tsens)
      tpos.time = grep("Verification",tsens)
      tmax = NA
      telb = NA
      ttime = NA
      if(length(tpos.max)>0){
        tmax = as.numeric(strsplit(gsub("\\s+","",tsens[tpos.max]),"is")[[1]][2])
      }
      if(length(tpos.elb)>0){
        telb = as.numeric(strsplit(gsub("\\s+","",tsens[tpos.elb]),"is")[[1]][2])
      }
      if(length(tpos.time)>0){
        ttime = as.numeric(strsplit(gsub("seconds!","",gsub("\\s+","",tsens[tpos.time])),"required")[[1]][2])
      }
      t2out = c(fiter,sout$x,sout$w,sout$min_cost,sout$mean_cost,tmax,telb,ttime)
      fout = as.data.frame(t(t2out))
      #dimnames(fout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")

      if(length(sout$x)==length(sout$w)){
        dimnames(fout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
      }else{
        tlen = length(sout$x)/length(sout$w)
        tdimx2 = c()
        for(i in 1:tlen){
          tdimx2 = c(tdimx2,paste(paste("x",i,sep=""),1:k,sep=""))
        }
        dimnames(fout)[[2]] = c("iter",tdimx2,paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
      }

    }

    out$out = dout
    #out$out2 = out$evol
    out$design = fout

    return(out)
}
######################################################################################################*
######################################################################################################*
#' @title Verifying Optimality of The Robust Designs
#'
#' @description
#' It plots the sensitivity (derivative) function of the
#'  robust criterion
#' at a given approximate (continuous) design and also
#'  calculates its efficiency lower bound (ELB) with respect
#' to the optimality criterion.
#' For an approximate (continuous) design, when the design space is one or two-dimensional,
#'  the user can visually verify the optimality of the design by observing the
#' sensitivity plot. Furthermore, the proximity of the design to the optimal design
#'  can be measured by the  ELB without knowing the latter.
#'  See, for more details, Masoudi et al. (2017).
#' @inheritParams robust
#' @inheritParams senslocally
#' @inherit senslocally return
#' @details
#' Let \eqn{\Theta}  be the set initial estimates for the model parameters and \eqn{\pi} be a probability measure having support in  \eqn{\Theta}.
#' A design \eqn{\xi^*}{\xi*} is robust with respect to  \eqn{\pi}
#' if the following inequality holds for all \eqn{\boldsymbol{x} \in \chi}{x belong to \chi}:
#'  \deqn{c(\boldsymbol{x}, \pi, \xi^*) = \int_{\pi} tr M^{-1}(\xi^*, \theta)I(\boldsymbol{x}, \theta)\pi(\theta) d(\theta)-p \leq 0,}{
#'          c(x, \pi, \xi*)  = integration over \pi with integrand tr M^-1(\xi*, \theta)I(x, \theta)\pi(\theta) d(\theta)-p <= 0}
#'           with equality at all support points of \eqn{\xi^*}{\xi*}.
#'            Here, \eqn{p} is the number of model parameters.
#'
#'  ELB is a measure of  proximity of a design to the optimal design without knowing the latter.
#' Given a design, let \eqn{\epsilon} be the global maximum
#'  of the sensitivity (derivative) function over \eqn{x \in \chi}{x belong to \chi}.
#' ELB is given by \deqn{ELB = p/(p + \epsilon),}
#' where \eqn{p} is the number of model parameters. Obviously,
#' calculating ELB requires finding \eqn{\epsilon} and
#' another optimization problem to be solved.
#' The tuning parameters of this optimization can be regulated via the argument  \code{\link{sens.minimax.control}}.
#'
#' @note
#' Theoretically, ELB can not be larger than 1. But if so, it may have one of the following reasons:
#' \itemize{
#' \item \code{max_deriv} is not a GLOBAL maximum.  Please increase  the value of the parameter \code{maxeval} in \code{\link{sens.minimax.control}} to find the global maximum.
#' \item The sensitivity function is shifted below the y-axis because
#' the number of model parameters has not been specified correctly (less value given).
#' Please specify the correct number of model parameters via the argument \code{npar}.
#' }
#'
#'
#' @seealso \code{\link{bayes}} \code{\link{sensbayes}} \code{\link{robust}}
#' @export
#' @example inst/examples/sensrobust_examples.R
sensrobust <- function (formula, predvars, parvars, family = gaussian(),
                        x, w,
                        lx, ux,
                        prob,
                        parset,
                        fimfunc = NULL,
                        sens.control = list(),
                        calculate_criterion = TRUE,
                        plot_3d = c("lattice", "rgl"),
                        plot_sens = TRUE,
                        npar = dim(parset)[2],
                        silent = FALSE,
                        crtfunc = NULL,
                        sensfunc = NULL){


  if (length(prob) != dim(parset)[1])
    stop("length of \"prior\" is not equal to the number of rows of \"param\"")
  if (!missing(formula)){
    if (dim(parset)[2] != length(parvars))
      stop("number of columns of  'parset' is not equal to the length of 'parvars'")
  }
  if (is.null(crtfunc))
    crt_type = "D" else
      crt_type = "user"
    # if (is.null(npar))
    #   npar <- dim(parset)[2]
    ## you must provide npar


    out <- sensminimax_inner(formula = formula, predvars = predvars, parvars = parvars,
                             family = family,
                             x = x, w = w,
                             lx = lx, ux = ux,
                             lp = NA, up = NA,
                             fimfunc = fimfunc,
                             sens.control = sens.control,
                             type = "robust",
                             localdes = NULL,
                             plot_3d = plot_3d[1],
                             plot_sens = plot_sens,
                             varlist = list(),
                             calledfrom = "sensfuncs",
                             npar = npar,
                             crt.minimax.control = list(inner_space = "robust_set"),
                             calculate_criterion = calculate_criterion,
                             robpars = list(prob = prob, parset = parset),
                             crt_type = crt_type,
                             multipars = list(),
                             silent = silent,
                             user_crtfunc = crtfunc,
                             user_sensfunc = sensfunc)
    out$method = "robust" # 06222020@seongho

    return(out)
}
######################################################################################################*
######################################################################################################*
#' @title
#' Locally Multiple Objective Optimal Designs for the 4-Parameter Hill Model
#'
#' @description
#'  The 4-parameter Hill model is of the form
#'  \deqn{f(D) = c + \frac{(d-c)(\frac{D}{a})^b}{1+(\frac{D}{a})^b} + \epsilon,}{
#'  f(D) = c + (d-c)(D/a)^b/(1 + (D/a)^b) + \epsilon,}
#' where \eqn{\epsilon \sim N(0, \sigma^2)}{\epsilon ~ N(0, \sigma^2)},
#'  \eqn{D} is the dose level and the predictor,
#' \eqn{a} is the ED50,
#'  \eqn{d} is the upper limit of response,
#'   \eqn{c} is the lower limit of response and
#'    \eqn{b} denotes the Hill constant that
#'  control the flexibility in the slope of the response curve.\cr
#'  Sometimes, the Hill model is re-parameterized and written as
#'  \deqn{f(x) = \frac{\theta_1}{1 + exp(\theta_2 x + \theta_3)} + \theta_4,}{
#'  f(x)= \theta_1/(1 + exp(\theta_2*x + \theta_3)) + \theta_4,}
#'   where \eqn{\theta_1 = d - c}, \eqn{\theta_2 = - b},
#'   \eqn{\theta_3 = b\log(a)}{\theta_3 = b*log(a)}, \eqn{\theta_4 = c}, \eqn{\theta_1 > 0},
#'   \eqn{\theta_2 \neq 0}{\theta_2 not equal to 0}, and \eqn{-\infty < ED50 < \infty},
#'   where \eqn{x = log(D) \in [-M, M]}{x = log(D) belongs to [-M, M]}
#'   for some sufficiently large value of \eqn{M}.
#'   The new form is sometimes called  4-parameter logistic model.\cr
#'  The function \code{multiple} finds locally multiple-objective optimal designs for estimating the model parameters, the ED50, and the MED, simultaneously.
#'    For more details, see Hyun and  Wong (2015).
#'
#' @param minDose Minimum dose \eqn{D}. For the 4-parameter logistic model, i.e. when \code{Hill_par = FALSE}, it is the minimum of \eqn{log(D)}.
#' @param maxDose  Maximum dose \eqn{D}. For the 4-parameter logistic model, i.e. when \code{Hill_par = FALSE}, it is the maximum of \eqn{log(D)}.
#' @inheritParams minimax
#' @param lambda A vector of relative importance of each of the three criteria,
#'  i.e. \eqn{\lambda = (\lambda_1, \lambda_2, \lambda_3)}.
#'   Here \eqn{0 < \lambda_i < 1} and  s \eqn{\sum \lambda_i = 1}.
# user select weights, where \eqn{\lambda_1}{\lambda1} is the weight for estimating parameters,
# \eqn{\lambda_2}{\lambda2} is the weight for estimating median effective dose level (ED50), and \eqn{\lambda_3}{\lambda3} is the weight for estimating minimum effective dose level (MED).
#' @param delta   Predetermined meaningful value of the minimum effective dose MED.
#' When \eqn{\delta < 0 }, then \eqn{\theta_2 > 0} or when \eqn{\delta > 0}, then \eqn{\theta_2 < 0}.
#' @param inipars A vector of initial estimates for the vector of parameters  \eqn{(a, b, c, d)}.
#'  For the 4-parameter logistic model, i.e. when \code{Hill_par = FALSE},
#'  it is  a vector of initial estimates for \eqn{(\theta_1, \theta_2,\theta_3, \theta_4)}.
#' @param Hill_par Hill model parameterization? Defaults to \code{TRUE}.
#' @param tol Tolerance for finding the general inverse of the Fisher information matrix. Defaults to \code{.Machine$double.xmin}.
#' @references
#' Hyun, S. W., and Wong, W. K. (2015). Multiple-Objective Optimal Designs for Studying the Dose Response Function and Interesting Dose Levels. The international journal of biostatistics, 11(2), 253-271.
#' @details
#'  When \eqn{\lambda_1 > 0}, then the number of support points \code{k}
#'   must at least be four to avoid singularity of the Fisher information matrix.
#'
#' One can adjust the tuning parameters in \code{\link{ICA.control}} to set a stopping rule
#' based on the general equivalence theorem. See 'Examples' below.
#' @note
#' This function is NOT appropriate for  finding  c-optimal designs for estimating 'MED' or 'ED50' (single objective optimal designs)
#'  and  the results may not be stable.
#'  The reason is that for the c-optimal criterion
#'  the generalized inverse of the Fisher information matrix
#'   is not stable and depends
#'  on the tolerance value (\code{tol}).
#' @importFrom utils capture.output
#' @export
#' @inherit locally return
#' @seealso \code{\link{sensmultiple}}
#' @example inst/examples/multiple_examples.R
multiple <- function(minDose, maxDose,
                     iter, k,
                     inipars,
                     Hill_par = TRUE,
                     delta,
                     lambda,
                     fimfunc = NULL,
                     ICA.control = list(),
                     sens.control = list(),
                     initial = NULL,
                     tol = sqrt(.Machine$double.xmin),
                     x = NULL){

  lx <- minDose
  ux <- maxDose
  if (length(lambda) != 3)
    stop("length of 'lambda' must be 3")
  if (sum(lambda) != 1)
    stop("sum of 'lambda' must be 1")
  if(any(lambda >1) || any(lambda < 0))
    stop("each element of 'lambda' must be between 0 an 1")
  if (!is.numeric(delta))
    stop("'delta' must be numeric")
  if (!is.logical(Hill_par))
    stop("'Hill_param' must be logical")
  if (any(lx > ux))
    stop("'ux' must be greater than lx")
  if (length(inipars) != 4)
    stop("length of 'inipars' must be 4")
  if (k < length(inipars))
    stop("\"k\" must be larger than the number of parameters to avoid singularity")

  if (!Hill_par){
    names(inipars) <- paste("theta", 1:4, sep = "")
    if (inipars["theta2"] < 0)
      if (!delta>0)
        stop("'delta > 0' when theta2 < 0'")
    if (inipars["theta2"] > 0)
      if (!delta<0)
        stop("'delta < 0' when theta2 > 0'")
    if (round(inipars["theta2"], 5) == 0)
      stop("'theta2 can not be zero")
    if (inipars["theta1"]<= 0)
      stop("theta1 must be positive")
  }else{
    names(inipars) <- letters[1:4]
    if (inipars["a"] <= 0 || inipars["d"] <= 0 || inipars["c"] <= 0)
      stop("a, c, and d must be positive")
    if (inipars["c"] > inipars["d"])
      if (lx <= 0)
        stop("'lx' must be positive")
    lx <- log(lx)
    ux <- log(ux)
    if (any(is.nan(c(lx, ux))))
      stop("'NaN produced for 'lx' or 'ux' when taking logarithm. Provide 'lx' and 'ux' accroding to the Hill model parameterization")
    inipars <- c(theta1 = inipars["d"]-inipars["c"], theta2 = -inipars["b"], theta3 = inipars["b"] * log(inipars["a"]), theta4 = inipars["c"])
    if (!is.null(x)) # we should convert the x to the scale of the 4PL model
      x <- log(x)
  }

  ## you must provide npar
  out <- minimax_inner(lx = lx, ux = ux, lp = inipars, up = inipars,
                       iter = iter, k = k,
                       fimfunc = FIM_logistic_4par,
                       ICA.control = ICA.control,
                       sens.control = sens.control,
                       crt.minimax.control = list(inner_space = "locally"),
                       type = "locally",
                       initial = initial,
                       localdes = NULL,
                       npar = 4,
                       robpars = list(),
                       crt_type = "multiple",
                       multipars = list(delta = delta, lambda = lambda, tol = tol),
                       plot_3d = "lattice",
                       only_w_varlist = list(x = x))
  #if (Hill_par)
  #  for (j in 1:length(out$evol))
  #    out$evol[[j]]$x <- exp(out$evol[[j]]$x)
  out$method = 'locally' # 06202020@seongho
  #if (!missing(formula)){
  #  out$call = formula # 06202020@seongho
  #}else{
    out$call = NULL
  #}

  dout = NULL
  fout = NULL
  fiter = length(out$evol)
  if(fiter>0){
    tout = c()
    for(i in 1:fiter){
      sout = out$evol[[i]]
      tout = rbind(tout,c(i,sout$x,sout$w,sout$min_cost,sout$mean_cost))
    }
    dout = as.data.frame(tout)
    #dimnames(dout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost")

    if(length(sout$x)==length(sout$w)){
      dimnames(dout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost")
    }else{
      tlen = length(sout$x)/length(sout$w)
      tdimx = c()
      for(i in 1:tlen){
        tdimx = c(tdimx,paste(paste("x",i,sep=""),1:k,sep=""))
      }
      dimnames(dout)[[2]] = c("iter",tdimx,paste("w",1:k,sep=""),"min_cost","mean_cost")
    }

    # extract sens outcomes
    sout = out$evol[[fiter]]
    tsens = capture.output(sout$sens)
    tpos.max = grep("Maximum",tsens)
    tpos.elb = grep("ELB",tsens)
    tpos.time = grep("Verification",tsens)
    tmax = NA
    telb = NA
    ttime = NA
    if(length(tpos.max)>0){
      tmax = as.numeric(strsplit(gsub("\\s+","",tsens[tpos.max]),"is")[[1]][2])
    }
    if(length(tpos.elb)>0){
      telb = as.numeric(strsplit(gsub("\\s+","",tsens[tpos.elb]),"is")[[1]][2])
    }
    if(length(tpos.time)>0){
      ttime = as.numeric(strsplit(gsub("seconds!","",gsub("\\s+","",tsens[tpos.time])),"required")[[1]][2])
    }
    t2out = c(fiter,sout$x,sout$w,sout$min_cost,sout$mean_cost,tmax,telb,ttime)
    fout = as.data.frame(t(t2out))
    #dimnames(fout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")

    if(length(sout$x)==length(sout$w)){
      dimnames(fout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
    }else{
      tlen = length(sout$x)/length(sout$w)
      tdimx2 = c()
      for(i in 1:tlen){
        tdimx2 = c(tdimx2,paste(paste("x",i,sep=""),1:k,sep=""))
      }
      dimnames(fout)[[2]] = c("iter",tdimx2,paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
    }

  }

  out$out = dout
  #out$out2 = out$evol
  out$design = fout

  return(out)
}

######################################################################################################*
######################################################################################################*
#' @title Verifying Optimality of The Multiple Objective Designs for The 4-Parameter Hill Model
#'
#' @description
#' This function uses general equivalence theorem to verify
#' the optimality of a multiple objective optimal design found for
#'  the 4-Parameter Hill model and  the 4-parameter logistic model.
#' For more details, See Hyun and Wong (2015).
#' @inheritParams multiple
#' @param dose A vector of design points. It is  either dose values or logarithm of dose values when \code{Hill_par = TRUE}.
#' @param w A vector of design weights.
#' @param silent Do not print anything? Defaults to \code{FALSE}.
#' @param calculate_criterion Calculate the criterion? Defaults to \code{TRUE}.
#' @param plot_sens Plot the sensitivity (derivative) function? Defaults to \code{TRUE}.
#' @export
#' @inherit senslocally return
#' @inherit multiple details
#' @details
#' ELB is a measure of  proximity of a design to the optimal design without knowing the latter.
#' Given a design, let \eqn{\epsilon} be the global maximum
#'  of the sensitivity (derivative) function over \eqn{x \in \chi}{x belong to \chi}.
#' ELB is given by \deqn{ELB = p/(p + \epsilon),}
#' where \eqn{p} is the number of model parameters. Obviously,
#' calculating ELB requires finding \eqn{\epsilon} and
#' another optimization problem to be solved.
#' The tuning parameters of this optimization can be regulated via the argument  \code{\link{sens.minimax.control}}.
#' See, for more details, Masoudi et al. (2017).
#'
#' @note
#' DO NOT use this function to verify  c-optimal designs for estimating 'MED' or 'ED50' (verifying single objective optimal designs) because the results may be unstable.
#'  The reason is that for the c-optimal criterion the generalized inverse of the Fisher information matrix is not stable and depends
#'  on the tolerance value (\code{tol}).
#'
#'  Theoretically, ELB can not be larger than 1. But if so, it may have one of the following reasons:
#' \itemize{
#' \item \code{max_deriv} is not a GLOBAL maximum.  Please increase  the value of the parameter \code{maxeval} in \code{\link{sens.minimax.control}} to find the global maximum.
#' \item The sensitivity function is shifted below the y-axis because
#' the number of model parameters has not been specified correctly (less value given).
#' Please specify the correct number of model parameters via argument \code{npar}.
#' }
#' @references
#' Hyun, S. W., and Wong, W. K. (2015). Multiple-Objective Optimal Designs for Studying the Dose Response Function and Interesting Dose Levels. The international journal of biostatistics, 11(2), 253-271.
#' @seealso \code{\link{multiple}}
#' @example inst/examples/sensmultiple_examples.R
sensmultiple <- function (dose, w,
                          minDose, maxDose,
                          inipars,
                          lambda,
                          delta,
                          Hill_par = TRUE,
                          sens.control = list(),
                          calculate_criterion = TRUE,
                          plot_sens = TRUE,
                          tol = sqrt(.Machine$double.xmin),
                          silent = FALSE){
  x <- dose
  lx <- minDose
  ux <- maxDose
  if (length(x) != length(w))
    stop("length of 'x' and 'w' is not equal")
  if (length(lambda) != 3)
    stop("length of 'lambda' must be 3")
  if (sum(lambda) != 1)
    stop("sum of 'lambda' must be 1")
  if(any(lambda >1) || any(lambda < 0))
    stop("each element of 'lambda' must be between 0 an 1")
  if (!is.numeric(delta))
    stop("'delta' must be numeric")
  if (!is.logical(Hill_par))
    stop("'Hill_param' must be logical")
  if (any(lx > ux))
    stop("'ux' must be greater than lx")
  if (length(inipars) != 4)
    stop("length of 'inipars' must be 4")
  if (!Hill_par){
    names(inipars) <- paste("theta", 1:4, sep = "")
    if (inipars["theta2"] < 0)
      if (!delta>0)
        stop("'delta > 0' when theta2 < 0'")
    if (inipars["theta2"] > 0)
      if (!delta<0)
        stop("'delta < 0' when theta2 > 0'")
    if (round(inipars["theta2"], 5) == 0)
      stop("'theta2 can not be zero")
    if (inipars["theta1"]<= 0)
      stop("theta1 must be positive")
  }else{
    names(inipars) <- letters[1:4]
    if (inipars["a"] <= 0 || inipars["d"] <= 0 || inipars["c"] <= 0)
      stop("a, c, and d must be positive")
    if (inipars["c"] > inipars["d"])
      if (lx <= 0)
        stop("'lx' must be positive")
    lx <- log(lx)
    ux <- log(ux)
    x <- log(x)
    if (any(is.nan(c(lx, ux))))
      stop("'NaN produced for 'lx' or 'ux' when taking logarithm. Provide 'lx' and 'ux' accroding to the Hill model parameterization")
    inipars <- c(theta1 = inipars["d"]-inipars["c"], theta2 = -inipars["b"], theta3 = inipars["b"] * log(inipars["a"]), theta4 = inipars["c"])
  }

  out <- sensminimax_inner(x = x, w = w,
                           lx = lx, ux = ux,
                           lp = inipars, up = inipars,
                           fimfunc = FIM_logistic_4par,
                           sens.control = sens.control,
                           type = "locally",
                           localdes = NULL,
                           plot_3d = "lattice", # not used
                           plot_sens = plot_sens,
                           varlist = list(),
                           calledfrom = "sensfuncs",
                           npar = 4,
                           crt.minimax.control = list(inner_space = "locally"),
                           calculate_criterion = calculate_criterion,
                           robpars = list(),
                           crt_type = "multiple",
                           multipars = list(delta = delta, lambda = lambda, tol = tol),
                           silent = silent)
  out$method = "locally" # 06222020@seongho

  return(out)
}

######################################################################################################*
######################################################################################################*
#' @title Locally DP-Optimal Designs
#' @inheritParams minimax
#' @inheritParams bayescomp
#' @param inipars Vector. Initial values for the unknown parameters. It will be passed to the information matrix and also probability function.
#' @param k Number of design points. When \code{alpha = 0}, then \code{k} can be less than the number of parameters.
#' @param npar Number of model parameters.  Used when \code{fimfunc} is given instead of \code{formula} to specify the number of model parameters.
#'   If not given, the sensitivity plot may be shifted below the y-axis. When \code{NULL}, it will be set here to \code{length(inipars)}.
#' @importFrom utils capture.output
#' @export
#' @inherit locally return
#' @description
#'  Finds compound locally DP-optimal designs that meet the dual goal of parameter estimation and
#'   increasing the probability of a particular outcome in a binary response model.
#'A compound locally DP-optimal design maximizes  the product of the efficiencies of a design \eqn{\xi} with respect to D- and average P-optimality, weighted by a pre-defined mixing constant
#' \eqn{0 \leq \alpha \leq 1}{0 <= \alpha <= 1}.
#' @details
#' Let \eqn{\Xi} be the space of all  approximate designs with
#'  \eqn{k} design points (support points) at \eqn{x_1, x_2, ...,  x_k}
#'   from  design space \eqn{\chi} with
#'  corresponding weights  \eqn{w_1,... ,w_k}.
#'  Let \eqn{M(\xi, \theta)} be the Fisher information
#'   matrix (FIM) of a \eqn{k-}point design \eqn{\xi},
#'    \eqn{\theta_0} is a user-given vector of initial estimates for the  unknown parameters \eqn{\theta} and
#'    \eqn{p(x_i, \theta)} is the ith probability of success
#' given by \eqn{x_i} in a binary response model.
#'   A compound locally DP-optimal design   maximizes over \eqn{\Xi}
#' \deqn{ \frac{\alpha}{q}\log|M(\xi, \theta_0)| + (1- \alpha)
#'\log \left( \sum_{i=1}^k w_ip(x_i, \theta_0) \right).}{
#'  \alpha/q log|M(\xi, \theta_0)| + (1- \alpha)
#'log ( \sum w_i p(x_i, \theta_0)).
#'}
#'
#' Use \code{\link{plot}} function to verify the general equivalence theorem for the output design or change \code{checkfreq} in \code{\link{ICA.control}}.
#'
#' One can adjust the tuning parameters in \code{\link{ICA.control}} to set a stopping rule
#' based on the general equivalence theorem. See "Examples" in \code{\link{locally}}.
#' @example inst/examples/locallycomp_examples.R
#' @references  McGree, J. M., Eccleston, J. A., and Duffull, S. B. (2008). Compound optimal design criteria for nonlinear models. Journal of Biopharmaceutical Statistics, 18(4), 646-661.
locallycomp <- function(formula, predvars, parvars, family = gaussian(),
                        lx, ux,
                        alpha,
                        prob,
                        iter, k,
                        inipars,
                        fimfunc = NULL,
                        ICA.control = list(),
                        sens.control = list(),
                        initial = NULL,
                        npar = length(inipars),
                        plot_3d = c("lattice", "rgl")){


  if (!missing(formula)){
    if (length(inipars) != length(parvars))
      stop("lengtb of 'inipars' is not equal to the length of 'parvars'")
  }
  if (alpha !=0)
    if (k < length(inipars))
      stop("\"k\" must be larger than the number of parameters to avoid singularity")

  if (is.formula(prob)){
    prob <- create_prob(prob = prob, predvars = predvars, parvars = parvars)
  }else{
    if (!is.function(prob))
      stop("'prob' must be either a function or a formula")
    if (!formalArgs(prob) %in% c("x", "param"))
      stop("arguments of 'prob' must be 'x' and 'param'")
  }

    out <- minimax_inner(formula = formula,
                       predvars = predvars, parvars = parvars, family = family,
                       lx = lx, ux = ux, lp = inipars, up = inipars, iter = iter, k = k,
                       fimfunc = fimfunc,
                       ICA.control = ICA.control,
                       sens.control = sens.control,
                       crt.minimax.control = list(inner_space = "locally"),
                       type = "locally",
                       initial = initial,
                       localdes = NULL,
                       npar = npar,
                       robpars = list(),
                       crt_type = "DPA",
                       multipars = list(),
                       plot_3d = plot_3d[1],
                       compound = list(prob = prob, alpha = alpha, npar = npar))


    out$method = 'locally' # 06202020@seongho

    if (!missing(formula)){
    out$call = formula # 06202020@seongho
  }else{
    out$call = NULL
  }

  plen = length(predvars)

  dout = NULL
  fout = NULL
  if(!is.null(out$evol)){
  #if(F){
    fiter = length(out$evol)
    if(fiter>0){
      tout = c()
      for(i in 1:fiter){
        sout = out$evol[[i]]
        tout = rbind(tout,c(i,sout$x,sout$w,sout$min_cost,sout$mean_cost))
      }
      dout = as.data.frame(tout)
      if(plen==1){
        dimnames(dout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost")
      }else{
        tdimx = c()
        for(i in 1:plen){
          tdimx = c(tdimx,paste(paste("x",i,sep=""),1:k,sep=""))
        }
        dimnames(dout)[[2]] = c("iter",tdimx,paste("w",1:k,sep=""),"min_cost","mean_cost")
      }
      # extract sens outcomes
      sout = out$evol[[fiter]]
      tsens = capture.output(sout$sens)
      tpos.max = grep("Maximum",tsens)
      tpos.elb = grep("ELB",tsens)
      tpos.time = grep("Verification",tsens)
      tmax = NA
      telb = NA
      ttime = NA
      if(length(tpos.max)>0){
        tmax = as.numeric(strsplit(gsub("\\s+","",tsens[tpos.max]),"is")[[1]][2])
      }
      if(length(tpos.elb)>0){
        telb = as.numeric(strsplit(gsub("\\s+","",tsens[tpos.elb]),"is")[[1]][2])
      }
      if(length(tpos.time)>0){
        ttime = as.numeric(strsplit(gsub("seconds!","",gsub("\\s+","",tsens[tpos.time])),"required")[[1]][2])
      }
      t2out = c(fiter,sout$x,sout$w,sout$min_cost,sout$mean_cost,tmax,telb,ttime)
      fout = as.data.frame(t(t2out))
      if(plen==1){
        dimnames(fout)[[2]] = c("iter",paste("x",1:k,sep=""),paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
      }else{
        tdimx = c()
        for(i in 1:plen){
          tdimx = c(tdimx,paste(paste("x",i,sep=""),1:k,sep=""))
        }
        dimnames(fout)[[2]] = c("iter",tdimx,paste("w",1:k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
      }
    }
  }
  out$out = dout
  #out$out2 = out$evol
  out$design = fout

  return(out)
}
######################################################################################################*
######################################################################################################*
#' @title Verifying Optimality of The Locally DP-optimal Designs
#' @inheritParams sensminimax
#' @inheritParams sensbayescomp
#' @param inipars Vector of initial estimates for the unknown parameters.
#' It must match \code{parvars} or argument \code{param} of the function provided in \code{fimfunc}.
#' @param npar Number of model parameters.  Used when \code{fimfunc} is given instead of \code{formula} to specify the number of model parameters.
#'   If not given, the sensitivity plot may be shifted below the y-axis.
#'   When \code{NULL}, it is set  to \code{length(inipars)}.
#'@description
#'  This function plot the sensitivity (derivative) function given an approximate (continuous) design and calculate the efficiency lower bound (ELB) for locally DP-optimal designs.
#' Let \eqn{\boldsymbol{x}}{x} belongs to \eqn{\chi} that denotes the design space.
#' Based on the general equivalence theorem, generally, a design \eqn{\xi^*}{\xi*} is optimal if and only if the value of its sensitivity (derivative) function
#' be non-positive for all \eqn{\boldsymbol{x}}{x} in \eqn{\chi} and it only reaches zero
#' when \eqn{\boldsymbol{x}}{x} belong to the support of \eqn{\xi^*}{\xi*} (be equal to one of the design point).
#' Therefore, the user can look at the sensitivity plot and the ELB to decide whether the
#' design is optimal or close enough to the true optimal design.
#'
#' @export
#' @inherit senslocally return
#' @example inst/examples/senslocallycomp_examples.R
#' @references  McGree, J. M., Eccleston, J. A., and Duffull, S. B. (2008). Compound optimal design criteria for nonlinear models. Journal of Biopharmaceutical Statistics, 18(4), 646-661.
senslocallycomp <- function (formula, predvars, parvars,
                             alpha,
                             prob,
                             family = gaussian(),
                             x, w,
                             lx, ux,
                             inipars,
                             fimfunc = NULL,
                             sens.control = list(),
                             calculate_criterion = TRUE,
                             plot_3d = c("lattice", "rgl"),
                             plot_sens = TRUE,
                             npar = length(inipars),
                             silent = FALSE){


  if (!missing(formula)){
    if (length(inipars) != length(parvars))
      stop("lengtb of 'inipars' is not equal to the length of 'parvars'")
  }

  if (is.formula(prob)){
    prob <- create_prob(prob = prob, predvars = predvars, parvars = parvars)
  }else{
    if (!is.function(prob))
      stop("'prob' must be either a function or a formula")
    if (!formalArgs(prob) %in% c("x", "param"))
      stop("arguments of 'prob' must be 'x' and 'param'")
  }


  out <- sensminimax_inner(formula = formula, predvars = predvars, parvars = parvars,
                           family = family,
                           x = x, w = w,
                           lx = lx, ux = ux,
                           lp = inipars, up = inipars,
                           fimfunc = fimfunc,
                           sens.control =  sens.control,
                           type = "locally",
                           localdes = NULL,
                           plot_3d = plot_3d[1],
                           plot_sens = plot_sens,
                           varlist = list(),
                           calledfrom = "sensfuncs",
                           npar = npar,
                           crt.minimax.control = list(inner_space = "locally"),
                           calculate_criterion = calculate_criterion,
                           robpars = list(),
                           crt_type = "DPA",
                           multipars = list(),
                           silent = silent,
                           compound = list(prob = prob, alpha = alpha, npar = npar))

  out$method = "locally" # 06222020@seongho

  return(out)
}

######################################################################################################*
######################################################################################################*
#  roxygen
#' Plotting \code{minimax} Objects
#' @description
#'  This function plots the evolution of the ICA algorithm (iteration vs the best (minimum) criterion value at each iteration) and also verifies the optimality of the last obtained design
#'  using  the general equivalence theorem. It plots the sensitivity function and calculates the ELB for the  best design generated at iteration number  \code{iter}.
#' @param x An object of class \code{minimax}.
#' @param iter Iteration number. if \code{NULL} (default), it will be set to the last iteration.
#' @param sensitivity Logical. If \code{TRUE} (default), the general equivalence theorem is used to check the optimality if the best design in iteration number \code{iter} and the sensitivity function will be plotted.
#' @param calculate_criterion  Logical. Re-calculate the criterion value (maybe with a set of new tuning parameters to be sure of the globality of the maximum over the parameter space given the design)? It only assumes a continuous parameter space for the minimax and standardized maximin designs.  Defaults to \code{FALSE}. See 'Details'.
#' @param sens.control Control Parameters for Calculating the ELB. For details, see the function \code{\link{sens.control}}.
#' @param sens.minimax.control Control parameters to verify general equivalence theorem. For details, see \code{\link{sens.minimax.control}}.
#' If \code{NULL} (default), it will be set to the  tuning parameters used to create object \code{x}.
#' @param crt.minimax.control Control parameters to optimize the minimax or standardized maximin criterion at a given design over a \strong{continuous} parameter space.
#'  For details, see \code{\link{crt.minimax.control}}.
#' @param sens.bayes.control Control parameters to verify general equivalence theorem for the Bayesian optimal designs. For details, see \code{\link{sens.bayes.control}}. If \code{NULL} (default), it will be set to the  tuning parameters used to create object \code{x}.
#' @param crt.bayes.control Control parameters to optimize the integration in the Bayesian criterion at a given design over a \strong{continuous} parameter space. For details, see \code{\link{crt.bayes.control}}. If \code{NULL} (default), it will be set to the  tuning parameters used to create object \code{x}.
#'  If \code{NULL} (default), it will be set to the  tuning parameters used to create object \code{x}.
#' @param silent Do not print anything? Defaults to \code{TRUE}.
#' @param plot_3d Which package should be used to plot the sensitivity function for two-dimensional design space. Defaults to \code{plot_3d = "lattice"}.
#' Only applicable when \code{sensitivity = TRUE}.
#' @param evolution Plot Evolution? Defaults to \code{FALSE}.
#' @param ... Argument with no further use.
#' @seealso \code{\link{minimax}}, \code{\link{locally}}, \code{\link{robust}}
#' @details
#'
#' In addition to verifying the general equivalence theorem,
#'  this function makes it possible to re-calculated the criterion value
#'   for the output designs using a new set of tuning parameters, especially,
#'    a large value for \code{maxeval} in the function \code{\link{crt.minimax.control}}.
#'  This is useful for  minimax and standardized maximin optimal designs to assess
#'   the robustness of the
#'  criterion value with respect to different values of \code{maxeval}.
#'  To put it simple, for these designs, the user can re-calculate the
#'  criterion value (finds the global maximum over the parameter space given an output design in a minimax problem)
#'   with larger values for  \code{maxeval} in \code{\link{crt.minimax.control}}
#'  to be sure that the function \code{nloptr} finds global optima of the inner
#'  optimization problem over the parameter space using the default value
#'  (or the user-given value) of \code{maxeval}.
#'  If increasing the value of \code{maxeval} returns different criterion values,
#'  then the results can not be trusted and the algorithm should be repeated with a higher value for \code{maxeval}.
#' @export
plot.minimax <- function(x, 
                         iter = NULL,
                         sensitivity = TRUE,
                         calculate_criterion = FALSE,
                         sens.minimax.control = list(),
                         crt.minimax.control = list(),
                         sens.bayes.control = list(),
                         crt.bayes.control = list(),
                         sens.control = list(),
                         silent = TRUE,
                         plot_3d = c("lattice", "rgl"),
                         evolution = FALSE,
                         ...){
  if (!evolution & !sensitivity){
    warning("Both 'sensitivity' and 'evolution' set to be FALSE.\nNo action is done in plot function!")
    return(invisible(NULL))
  }
  #if (any(class(x) != c("list", "minimax"))) # 062020@seongho
  if (any(class(x) != c("minimax")))
    stop("'x' must be of class 'minimax'")

  ## to not be confused with design points
  obj <- x
  arg <- obj$arg
  if(is.null(iter))
    totaliter <- length(obj$evol) else
      totaliter <- iter
  if (totaliter > length(x$evol))
    stop("'iter' is larger than the maximum number of iterations")

  #cat("#### HERE HERE ####\n")

  if(x$method!="bayes"){
    #sens.minimax.control = sens.minimax.bayes.control
    #crt.minimax.control = crt.minimax.bayes.control
    if (calculate_criterion || sensitivity){
      if (is.null(sens.minimax.control)){
        sens.minimax.control <- arg$sens.minimax.control}
      else {
        sens.minimax.control <- do.call("sens.minimax.control", sens.minimax.control)
      }
      if (is.null(sens.control)){
        sens.control <- arg$sens.control}
      else {
        sens.control <- do.call("sens.control", sens.control)
      }
      if (is.null(crt.minimax.control)){
        crt.minimax.control <- arg$crt.minimax.control}
      else {
        crt.minimax.control <- do.call("crt.minimax.control", crt.minimax.control)
        if (arg$type == "minimax")
          crt.minimax.control$inner_space <- "continuous"
      }

      optim_starting <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id,  npred){
        if (!arg$is.only.w)
          q1 <- c(x, w) else
            q1 <- w
          out <- optim2(fn = fn, lower = lower, upper = upper,
                        n_seg = sens.minimax.control$n_seg,
                        q = q1,
                        fixedpar = fixedpar, fixedpar_id = fixedpar_id,
                        npred= npred)
          minima <- out$minima
          counts <- out$counts
          return(list(minima =minima, counts = counts))
      }
      sens_varlist <-list(fixedpar = arg$fixedpar, fixedpar_id = arg$fixedpar_id,
                          npred = length(arg$lx),
                          crfunc_sens = arg$crfunc_sens,
                          lp_nofixed = arg$lp_nofixed,
                          up_nofixed = arg$up_nofixed,
                          plot_3d = plot_3d,
                          npar = arg$npar,
                          optim_starting = optim_starting,
                          fimfunc_sens = arg$FIM_sens,
                          Psi_x_minus_minimax = arg$Psi_funcs$Psi_x_minus_minimax,
                          Psi_x = arg$Psi_funcs$Psi_x,
                          Psi_xy = arg$Psi_funcs$Psi_xy, Psi_Mu = arg$Psi_funcs$Psi_Mu)

      sens_res <- sensminimax_inner(x = obj$evol[[totaliter]]$x, w = obj$evol[[totaliter]]$w,
                                    lx = arg$lx, ux = arg$ux,
                                    lp = arg$lp_nofixed, up = arg$up_nofixed,
                                    fimfunc = arg$FIM,
                                    sens.minimax.control = sens.minimax.control,
                                    sens.control = sens.control,
                                    type = arg$type,
                                    localdes = arg$localdes,
                                    plot_sens = TRUE,
                                    varlist = sens_varlist,
                                    calledfrom = "plot",
                                    npar = arg$npar,
                                    calculate_criterion = calculate_criterion,
                                    crt.minimax.control = crt.minimax.control,
                                    robpars = arg$robpars,
                                    plot_3d = plot_3d[1],
                                    silent = silent,
                                    calculate_sens = sensitivity)
    }

    #cat("#### HOWABOUT ####\n")

    if (evolution){
      ## extract evolution data from the object
      mean_cost <- sapply(1:totaliter, FUN = function(j)obj$evol[[j]]$mean_cost)
      min_cost <- sapply(1:totaliter, FUN = function(j)obj$evol[[j]]$min_cost)
      if (calculate_criterion)
        min_cost[totaliter] <- sens_res$crtval


      type <- obj$arg$type

      # plot setting
      legend_place <- "topright"
      legend_text <- c( "Best Imperialist", "Mean of Imperialists")
      line_col <- c("firebrick3", "blue4")
      if (type == "minimax")
        title1 <- "cost value"
      if (type == "standardized")
        title1 <- "minimum efficiency"
      if (type == "locally" || type == "robust")
        title1 <- "log determinant of inverse of FIM"
      if (type == "multiple_locally")
        title1 <- "criterion value"

      ylim = switch(type,
                    "minimax" = c(min(min_cost) - .07, max(mean_cost[1:totaliter]) + .2),
                    "standardized" = c(min(mean_cost[1:totaliter]) - .07, max( min_cost) + .2),
                    "locally" = c(min(min_cost) - .07, max(mean_cost[1:totaliter]) + .2),
                    "robust" = c(min(min_cost) - .07, max(mean_cost[1:totaliter]) + .2))

      PlotEffCost(from = 1,
                  to = (totaliter),
                  AllCost = min_cost, ##all criterion up to now (all cost function)
                  UserCost = NULL,
                  DesignType = type,
                  col = line_col[1],
                  xlab = "Iteration",
                  ylim = ylim,
                  lty = 1,
                  title1 = title1,
                  plot_main = TRUE)
      ICAMean_line <-  mean_cost[1:(totaliter)]
      lines(x = 1:(totaliter),
            y = ICAMean_line,
            col = line_col[2], type = "s", lty = 5)
      legend(x = legend_place,  legend = legend_text,lty = c(1,5, 3), col = line_col, xpd = TRUE, bty = "n")
    }

    #cat("####****#### HERE????? #####******#####\n")

  }else{ # bayes
    #sens.bayes.control = sens.minimax.bayes.control
    #crt.bayes.control = crt.minimax.bayes.control
    if (calculate_criterion || sensitivity){
      if (is.null(sens.bayes.control)){
        sens.bayes.control <- arg$sens.bayes.control }else {
          sens.bayes.control <- do.call("sens.bayes.control", sens.bayes.control)
        }
      if (is.null(sens.control)){
        sens.control <- arg$sens.control }else {
          sens.control <- do.call("sens.control", sens.control)
        }
      if (is.null(crt.bayes.control)){
        crt.bayes.control <- arg$crt.bayes.control
        Psi_x_bayes  <- arg$Psi_funcs$Psi_x_bayes
        Psi_xy_bayes  <- arg$Psi_funcs$Psi_xy_bayes
      } else {

        crt.bayes.control <- do.call("crt.bayes.control", crt.bayes.control)
        temp_psi <- create_Psi_bayes(type = arg$type, prior = arg$prior, FIM = arg$FIM, lp = arg$prior$lower,
                                     up = arg$prior$upper, npar = arg$npar,
                                     truncated_standard = arg$truncated_standard,
                                     const = arg$const, sens.bayes.control = sens.bayes.control,
                                     compound = arg$compound,
                                     method =  sens.bayes.control$method,
                                     user_sensfunc = arg$user_sensfunc)
        Psi_x_bayes  <- temp_psi$Psi_x_bayes
        Psi_xy_bayes  <- temp_psi$Psi_xy_bayes
      }

      vertices_outer <- make_vertices(lower = arg$lx, upper = arg$ux)
      sens_varlist <-list(npred = length(arg$lx),
                          # plot_3d = "lattice",
                          npar = arg$npar,
                          fimfunc_sens = arg$FIM_sens,
                          Psi_x_bayes  = Psi_x_bayes,
                          Psi_xy_bayes  = Psi_xy_bayes,
                          crfunc = arg$crfunc,
                          vertices_outer = vertices_outer)

      sens_res <- sensbayes_inner(x = obj$evol[[totaliter]]$x, w = obj$evol[[totaliter]]$w,
                                  lx = arg$lx, ux = arg$ux,
                                  fimfunc = arg$FIM,
                                  prior = arg$prior,
                                  sens.control = sens.control,
                                  sens.bayes.control = sens.bayes.control,
                                  crt.bayes.control = crt.bayes.control,
                                  type = arg$type,
                                  plot_3d = plot_3d[1],
                                  plot_sens = TRUE,
                                  const = arg$const,
                                  compound = arg$compound,### you dont need compund here
                                  varlist = sens_varlist,
                                  calledfrom =  "plot",
                                  npar = arg$npar,
                                  calculate_criterion = calculate_criterion,
                                  silent = silent,
                                  calculate_sens = sensitivity)
    }
    if (evolution){
      ## extract evolution data from the object
      mean_cost <- sapply(1:totaliter, FUN = function(j)obj$evol[[j]]$mean_cost)
      min_cost <- sapply(1:totaliter, FUN = function(j)obj$evol[[j]]$min_cost)
      if (calculate_criterion)
        min_cost[totaliter] <- sens_res$crtval


      type <- obj$arg$type

      # plot setting
      legend_place <- "topright"
      legend_text <- c( "Best Imperialist", "Mean of Imperialists")
      line_col <- c("firebrick3", "blue4")
      title1 <- "Bayesian criterion"
      ylim = switch(type, "bayesian_D" = c(min(min_cost) - .07, max(mean_cost[1:(totaliter)]) + .2))


      PlotEffCost(from = 1,
                  to = (totaliter),
                  AllCost = min_cost, ##all criterion up to now (all cost function)
                  UserCost = NULL,
                  DesignType = type,
                  col = line_col[1],
                  xlab = "Iteration",
                  ylim = ylim,
                  lty = 1,
                  title1 = title1,
                  plot_main = TRUE)
      ICAMean_line <-  mean_cost[1:(totaliter)]
      lines(x = 1:(totaliter),
            y = ICAMean_line,
            col = line_col[2], type = "s", lty = 5)
      legend(x = legend_place,  legend = legend_text,lty = c(1,5, 3), col = line_col, xpd = TRUE, bty = "n")
    }
  }

  #cat("#### ************** THERE **************#####\n")

  if(sensitivity || calculate_criterion)
    return(sens_res) else
      return(invisible(NULL))
}
######################################################################################################*
######################################################################################################*
#' Printing \code{minimax} Objects
#'
#' Print method for an object of class \code{minimax}.
#' 
#' @param x An object of class \code{minimax}.
#' @param ... Argument with no further use.
#' @export
#' @seealso \code{\link{minimax}}, \code{\link{locally}}, \code{\link{robust}}, \code{\link{bayes}}
print.minimax <- function(x, ...){

  if (any(class(x) != c("minimax")))
    stop("'x' must be of class 'minimax'")
  object <- x
  totaliter <- length(object$evol)
  type <- x$arg$type

  # 06202020@seongho
  flab = ""
  if(x$method=="locally")
    flab = "locally optimal designs"
  else if(x$method=="minimax")
    flab = "minimax optimal designs"
  else if(x$method=="robust")
    flab = "robust or optimum-on-average optimal designs"
  else if(x$method=="bayes")
    flab = "Bayesian optimal designs"

  cat("\nFinding ",flab,"\n")
  #cat("\nCall:\n",
  #    paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  if(!is.null(x$call)){
    cat("\nCall:\n",
        paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }else{
    cat("\nCall:\n","FIM defined directly via the argument 'fimfunc'","\n\n",sep = "")
  }

  #cat("Optimal designs' profile:\n")
  tpos = seq(1,totaliter,length=min(10,totaliter))
  print(object$out[tpos,])
  #cat("")

  cat("\nOptimal designs (k=",x$arg$k,"):\n",sep="")
  cat(" ",print_xw_char(x = object$evol[[totaliter]]$x, w =  object$evol[[totaliter]]$w,
                       npred = length(object$arg$lx), is.only.w = object$arg$is.only.w,
                       equal_weight = object$arg$ICA.control$equal_weight)
      ,sep=""
  )
  cat("\n\n ICA iteration:", totaliter)
  cat("\n Criterion value: ", object$evol[[totaliter]]$min_cost)

  cat("\n Total number of function evaluations:", object$alg$nfeval)
  cat("\n Total number of successful local search moves:", object$alg$nlocal)
  cat("\n Total number of successful revolution moves:", object$alg$nrevol)

  cat("\n Convergence:", object$alg$convergence)

  if (object$arg$ICA.control$only_improve)
    cat("\n Total number of successful assimilation moves:", object$alg$nimprove)
  if (type == "minimax")
    cat( "\n Vector of maximum parameter values: ", object$evol[[totaliter]]$param)
  if (type == "standardized")
    cat( "\n Vector of minimum parameter values: ", object$evol[[totaliter]]$param)
  cat("\n CPU time:", object$arg$time[1], " seconds!\n")

  if (!is.null(object$evol[[totaliter]]$sens))
    print(object$evol[[totaliter]]$sens)

  return(invisible(NULL))
}

######################################################################################################*
######################################################################################################*
#' Printing \code{sensminimax} Objects
#'
#' Print method for an object of class \code{sensminimax}.
#' @param x An object of class \code{sensminimax}.
#' @param ... Argument with no further use.
#' @export
#' @seealso \code{\link{sensminimax}}, \code{\link{senslocally}}, \code{\link{sensrobust}}

print.sensminimax <- function(x,...){
  #if (any(class(x) != c("list", "sensminimax"))) # 06202020@seongho
  if (any(class(x) != c("sensminimax")))
    stop("'x' must be of class 'sensminimax'")

  #cat("#############************** SENSE HERE ################\n")
  #print(x)
  #print(names(x))
  #cat("################################ SENSE END ################################\n")

  #if(x$method!="bayes"){
    #cat("\n***********************************************************************")
    if (!is.null(x$max_deriv))
      cat("\n Maximum of the sensitivity function is ", x$max_deriv, "\n Efficiency lower bound (ELB) is ", x$ELB)
    if (!is.null(x$crtval))
      cat("\n Criterion value is ", x$crtval)
    if (x$type != "locally" & x$type != "robust"){
      cat("\n Verification required",x$time, "seconds!", "\n Adjust the control parameters in 'sens.minimax.control' ('n_seg') or in 'sens.bayes.control' for higher speed.", "\n\n")
    }else{
      cat("\n Verification required",x$time, "seconds!", "\n\n")
    }
  #}else{
  #  #cat( #"\n***********************************************************************")
  #  cat("\n Maximum of the sensitivity function is ", x$max_deriv, "\n Efficiency lower bound (ELB) is ", x$ELB,
  #    "\n Verification required",x$time, "seconds!",
  #    "\n Adjust the control parameters in 'sens.bayes.control' for higher speed\n\n")
    #   if (x$type == "minimax")
    #     optimchar <- "\nSet of all maxima over parameter space\n"
    #   if (x$type == "standardized")
    #     optimchar <- "\nSet of all minima over parameter space\n"
    #   cat("\nAnswering set:\n", answer, optimchar, optima, "\n")
    # }
    #
    # if (!is.null(x$crtval))
    #   cat("Criterion value found by 'nloptr':", x$crtval)
    #cat("***********************************************************************")
  #}

  return(invisible(NULL))
}

######################################################################################################*
######################################################################################################*
#' Returns Control Parameters for Optimizing Minimax Criteria Over The Parameter Space
#'
#'
#' The function \code{crt.minimax.control} returns a list of \code{\link[nloptr]{nloptr}} control parameters for optimizing the minimax criterion over the parameter space.\cr
#' The key tuning parameter for our application is \strong{\code{maxeval}.}
#' Its value should be increased when either the dimension or the size of the parameter space becomes larger
#'  to avoid pre-mature convergence in the inner optimization problem over the parameter space.
#' If the CPU time matters, the user should find an appropriate speed-accuracy trade-off  for her/his own design problem.
#'
#' @param x0 Vector of the starting values for the optimization problem (must be from the parameter space).
#' @param optslist A list. It will be passed to the argument \code{opts} of the function \code{\link[nloptr]{nloptr}}. See 'Details'.
#' @param ... Further arguments will be passed to \code{\link{nl.opts}} from package \code{\link[nloptr]{nloptr}}.
#' @importFrom nloptr nl.opts
#' @importFrom nloptr nloptr.print.options
#' @return A list of control parameters for the function \code{\link[nloptr]{nloptr}}.
#' @details
#'  Argument \code{optslist} will be passed to the argument \code{opts} of the function \code{\link[nloptr]{nloptr}}:
#'  \describe{
#'   \item{\code{stopval}}{Stop minimization when an objective value <= \code{stopval} is found. Setting \code{stopval} to \code{-Inf} disables this stopping criterion (default).}
#'   \item{\code{algorithm}}{Defaults to \code{NLOPT_GN_DIRECT_L}. DIRECT-L is a deterministic-search algorithm based on systematic division of the search domain into smaller and smaller hyperrectangles.}
#'   \item{\code{xtol_rel}}{Stop when an optimization step (or an estimate of the optimum) changes every parameter by less than \code{xtol_rel} multiplied by the absolute value of the parameter. Criterion is disabled if \code{xtol_rel} is non-positive. Defaults to \code{1e-5}.}
#'   \item{\code{ftol_rel}}{Stop when an optimization step (or an estimate of the optimum) changes the objective function value by less than \code{ftol_rel} multiplied by the absolute value of the function value. Criterion is disabled if \code{ftol_rel} is non-positive. Defaults to \code{1e-8}.}
#'   \item{\code{maxeval}}{Stop when the number of function evaluations exceeds \code{maxeval}. Criterion is disabled if \code{maxeval} is non-positive. Defaults to \code{1000}. See below.}
#' }
#'  A detailed explanation of all the options is shown by \code{nloptr.print.options()} in package \code{\link[nloptr]{nloptr}}.
#'
#' @export
#' @examples
#' crt.minimax.control(optslist = list(maxeval = 2000))
crt.minimax.control <- function (x0 = NULL,
                                 optslist = list(stopval = -Inf,
                                                 algorithm = "NLOPT_GN_DIRECT_L",
                                                 xtol_rel = 1e-6,
                                                 ftol_rel = 0,
                                                 maxeval = 1000), ...){

  optstlist2 <- do.call(c, list(optslist, list(...)))
  outlist <- suppressWarnings(nloptr::nl.opts(optstlist2))
  outlist["algorithm"] <- optslist$algorithm

  ## outlist has the defaut values of the nl.opts, when any component is null.
  # we play with that here
  if (is.null(optslist$algorithm))
    outlist$algorithm <- "NLOPT_GN_DIRECT_L"
  if (is.null(optslist$stopval))
    outlist$stopval<- -Inf
  if (is.null(optslist$xtol_rel))
    outlist$xtol_rel <- 1e-6
  if (is.null(optslist$ftol_rel))
    outlist$ftol_rel <- 0
  if (is.null(optslist$maxeval))
    outlist$maxeval <- 1000

  return(list(x0 = x0, optslist = outlist))
}
######################################################################################################*
######################################################################################################*
#' @title Returns Control Parameters for Verifying General Equivalence Theorem For Minimax Optimal Designs
#'
#'
#' @description
#' This function returns a list of control parameters that are used to find
#' the ``answering set'' for minimax and
#' standardized maximin designs.
#'  The answering set is required to  obtain the sensitivity (derivative) function in order to verify the optimality of
#'  a given design.
#'
#'
#' @param n_seg For a given design, the number of starting points in the local search to find all the local maxima of the minimax criterion over the parameter space is equal to \code{(n_seg + 1)^p}. Defaults to \code{6}.
#' Please increase its value when the parameter space is large. It is also applicable for standardized maximin designs. See 'Details' of \link{sens.minimax.control}.
#' @param merge_tol Merging tolerance. It is used  to  specify the elements of the answering set
#' by choosing only the local maxima (found by the local search) that are nearer to
#' the global maximum.  See 'Details' of \link{sens.minimax.control}. Defaults to \code{0.005}.
#' We advise to not change its default value because it has been successfully tested on many optimal design problems.
#' @return A list of control parameters for verifying the general equivalence
#'  theorem for minimax and standardized maximin optimal designs.
#' @details
#'  Given a design, an ``answering set'' is a subset of all the local optima
#'   of the optimality criterion over the parameter space.
#' Answering set is used to obtain the sensitivity function
#' of a minimax or standardized maximin criterion.
#'    Therefore, an invalid  answering set may result in a false
#'   sensitivity plot and ELB.
#'  Unfortunately, there is no theoretical rule on how to choose the number of elements of
#'   the answering set; and they  have to be found by trial and error.
#'  Given a design, the answering set for a minimax criterion is obtained as follows:
#'  \itemize{
#'  \item{Step 1: }{Find all the local maxima of the optimality criterion (minimax)  over the parameter space.
#'   For this purpose,  the parameter space is divided into \code{(n_seg + 1)^p} segments,
#'    where p is the number of unknown model parameters.
#'    Then, each boundary point of the resulted segments (intervals) is assigned to the argument
#'    \code{par} of the function \code{optim} in order to start a local search
#'     using the \code{"L-BFGS-B"} method.}
#'  \item{Step 2: }{Pick the ones nearest to the global minimum subject to a merging tolerance
#'   \code{merge_tol} (default \code{0.005}).}
#' }
#' Obviously, the answering set is a subset of all the local maxima over the parameter space (or local minima in case of standardized maximin criteria)
#' Therefore, it is very important to be able to find all the local maxima to create the true answering set with no missing elements.
#'  Otherwise, even when the design is optimal, the sensitivity (derivative) plot may not reveal its optimality.
#'
#'    Note that the minimax criterion (or standardized maximin criterion)
#'    is a multimodel function especially near the optimal design and
#'    this makes the job of finding all the locall maxima (minima) over the
#'    parameter space very complicated.
#'
#'
#' @export
#' @examples
#' sens.minimax.control()
#' sens.minimax.control(n_seg = 4)
sens.minimax.control <- function(n_seg = 6, merge_tol = .005){

  # the default value is tha same as the ones in the list, so no modification is
  if (!is.numeric(n_seg) || n_seg <= 0)
    stop("The value of 'n_seg' must be > 0")
  if (!is.numeric(merge_tol) || merge_tol <= 0)
    stop("The value of 'merge_tol' must be > 0")

  return(list(n_seg = n_seg, merge_tol = merge_tol))
}
######################################################################################################*
######################################################################################################*
# roxygen
#' @title Updating an Object of Class \code{minimax}
#'
#' @description  Runs the ICA optimization algorithm on an object of class \code{minimax} for more number of iterations  and updates the results.
#'
#' @param object An object of class \code{minimax}.
#' @param iter Number of iterations.
#' @param ... An argument of no further use.
#' @importFrom nloptr directL
#' @importFrom utils capture.output
#' @seealso \code{\link{minimax}}
#' @export
update.minimax <- function(object, iter, ...){
  if (all(class(object) != c("minimax")))
    stop("''object' must be of class 'minimax'")
  if (missing(iter))
    stop("'iter' is missing")

  #if(object$method!="bayes"){
  if(any(object$arg$type==c("D", "DPA", "DPM", "multiple", "user"))){
    bayes.update(object=object,iter=iter,...)
  }else{
    minimax.update(object=object,iter=iter,...)
  }
}

#' @importFrom utils capture.output
minimax.update <- function(object, iter, ...){ # 06212020@seongho
  # ... is an argument of no use. Only to match the generic update
  #if (all(class(object) != c("list", "minimax"))) # 06202020@seongho

  # blocked by 06212020@seongho
  #if (all(class(object) != c("minimax")))
  #    stop("''object' must be of class 'minimax'")
  #if (missing(iter))
  #  stop("'iter' is missing")

  arg <- object$arg
  ICA.control <- object$arg$ICA.control
  crt.minimax.control <- object$arg$crt.minimax.control
  sens.minimax.control <-  object$arg$sens.minimax.control
  sens.control <-  object$arg$sens.control
  evol <- object$evol
  evol.flag <- 0 # 1: updating, 0: not updating # 06212020@seongho
  #if (!arg$is.only.w)
  npred <- length(arg$lx) #else
  #    npred <- NA #dim(arg$x)[2]
  ## all fo the types for optim_on_average will be set to be equal to "robust"
  ## but the arg$type remains unchanged to be used in equivalence function!!
  # if( grepl("on_average", arg$type))
  #   type = "robust" else
  #     type <- arg$type
  type <- arg$type

  # warning: no arg$type must be used further

  if (type == "robust")
    npar <- dim(arg$param)[2] else
      npar <- length(arg$lp)
  if (ICA.control$equal_weight)
    w_equal <- rep(1/arg$k, arg$k)

  ###############################################################*
  ## multi_locally is the same as locally in update!

  if (type == "multiple_locally"){
    type <- "locally"
    ## rewuired for setting the title of plots
    multi_type <- TRUE
  }else
    multi_type <- FALSE

  # if (type == "multiple_minimax")
  #   type <- "minimax"

  if (!(type %in% c("minimax", "standardized", "locally", "robust")))
    stop("Bug: the type must be 'minimax' or 'standardized' or 'locally' or 'ave' in 'iterate.minimax\nset  'multiple_locally' to 'locally'")
  # because they have the same configuration. But we need to know the multi becasue of the verifying and plot methods!
  ##################################################################*

  # if (type == "locally")
  #   param_locally <- arg$up
  # if (type == "robust")
  #   param_set <- arg$robpars$parset


  ############################################################*
  ###finding if there is any fixed parameters.
  # only if type != "locally"
  #if (type != "locally" & control$inner_space != "vertices" & control$inner_space != "discrete"){
  # if (type != "locally" && type != "robust"){
  #   # here we search if one of the parameters are fixed. then we pass it to the optimization function in the inner problem because otherwise it may casue an error.
  #   any_fixed <- sapply(1:length(lp), function(i) lp [i] == up[i]) # is a vector
  #   if (any(any_fixed)){
  #     is_fixed <- TRUE
  #     fixedpar_id <- which(any_fixed)
  #     fixedpar <- lp[fixedpar_id]
  #     lp <- lp[-fixedpar_id]
  #     up <- up[-fixedpar_id]
  #     ## warning: lp and up are channged here if 'any_fix == TRUE'
  #   }else{
  #     fixedpar <- NA
  #     fixedpar_id <- NA
  #     is_fixed <- FALSE
  #   }
  # }else{
  #   fixedpar <- NA
  #   fixedpar_id <- NA
  #   is_fixed <- FALSE
  # }
  if(crt.minimax.control$inner_space == "discrete"){
    if(!is.na(arg$fixedpar))
      discrete_set <- crt.minimax.control$param_set[, -arg$fixedpar_id, drop = FALSE] else
        discrete_set <- crt.minimax.control$param_set
  }else
    discrete_set <-NULL
  ########################################################*

  ########################################################*
  # plot setting
  #plot_cost <- control$plot_cost
  #plot_sens <- control$plot_sens
  legend_place <- "topright"
  legend_text <- c( "Best Imperialist", "Mean of Imperialists")
  line_col <- c("firebrick3", "blue4")
  if (type == "minimax")
    title1 <- "cost value"
  if (type == "standardized")
    title1 <- "minimum efficiency"
  if (type == "locally" || type == "robust")
    title1 <-  "log determinant of inverse of FIM"
  if (multi_type)
    title1 <- "criterion value"
  ##################################################################*
  ## In last iteration the check functions should be applied??
  check_last <- ifelse(ICA.control$checkfreq == 0, FALSE, TRUE)

  # ##################################################################*
  # ### re-defimimg crfunc to handle fixed parameters.
  # crfunc <- arg$crfunc
  # if (is_fixed){
  #   crfunc2 <- function(param, q, fixedpar = NA, fixedpar_id = NA, npred){
  #     # if (any(!is.na(fixedpar))){
  #     #   if (any(is.na(fixedpar_id)))
  #     #     stop("'fixedpar' index is missing.")
  #     param_new <- rep(NA, length(param) + length(fixedpar))
  #     param_new[fixedpar_id] <- fixedpar
  #     param_new[-fixedpar_id] <- param
  #     param <- param_new
  #     #}
  #     out <- crfunc(param = param, q = q, npred = npred)
  #     return(out)
  #   }
  # }else{
  #   crfunc2 <- function(param, q, fixedpar = NA, fixedpar_id = NA, npred){
  #     # no use for fixedpar  and fixedpar_id = NA
  #     out <- crfunc(param = param, q = q, npred = npred)
  #     return(out)
  #   }
  # }
  # #####################################################################*

  ####################################################################*
  ### Psi as a function of x and x, y for plotting. Psi_x defined as minus psi to find the minimum
  ## Psi_x is mult-dimensional, x can be of two dimension.
  #
  #   Psi_x_minus <- function(x1, mu,  FIM,  x, w,  answering){
  #     ## mu and answering are only to avoid having another if when we want to check the maximum of sensitivity function
  #     Out <- arg$Psi_x(x1 = x1, mu =  mu, FIM = FIM,  x = x, w = w, answering = answering)
  #     return(-Out)
  #   }
  if(length(arg$lx) == 1)
    Psi_x_plot <-  arg$Psi_x ## for PlotPsi_x
  # it is necessary to distniguish between Psi_x for plotting and finding ELB becasue in plotting for models with two
  # explanatory variables the function should be defined as a function of x, y (x, y here are the ploints to be plotted)
  if(length(arg$lx) == 2)
    Psi_x_plot <- arg$Psi_xy
  #when length(lx) == 1, then Psi_x_plot = Psi_x
  #########################################################################*

  #########################################################################*
  # required for finding the answering set for verification
  #if (length(lp) <= 2)
  optim_starting <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id,  npred){
    if (!arg$is.only.w)
      q1 <- c(x, w) else
        q1 <- w
      out <- optim2(fn = fn, lower = lower, upper = upper,
                    n_seg = sens.minimax.control$n_seg,
                    q = q1,
                    fixedpar = fixedpar, fixedpar_id = fixedpar_id,
                    npred= npred)
      minima <- out$minima
      counts <- out$counts
      return(list(minima =minima, counts = counts))
  }
  optim_func <- create_optim_func(type = type, lp_nofixed = arg$lp_nofixed, up_nofixed = arg$up_nofixed,
                                  crt.minimax.control = crt.minimax.control,
                                  discrete_set = discrete_set, robpars = arg$robpars,
                                  inipars = arg$inipars, is.only.w = arg$is.only.w)

  ################################################################################*
  ## x_id, w_id are the index of x and w in positions
  #cost_id is the index of
  ## in symmetric case the length of x_id can be one less than the w_id if the number of design points be odd!

  if (!arg$is.only.w){
    if (ICA.control$sym)
      x_id <- 1:floor(arg$k/2) else
        x_id <- 1:(arg$k * npred)

      if (!ICA.control$equal_weight)
        w_id <- (x_id[length(x_id)] + 1):length(arg$ld) else
          w_id <- NA
  }else{
    w_id <- 1:length(arg$ld)
    x_id <- NA
  }

  ###column index of cost in  matrix output of the inner problem
  if (type != "robust")
    CostColumnId <- length(arg$lp_nofixed) + 1 else
      CostColumnId <- dim(arg$robpars$parset)[2] + 1


  ## warning: not the lp withot fixed param

  ##########################################################################*
  ## whenever Calculate_Cost is used, the fixed_arg list should be passed to
  ## fixed argumnet for function Calculate_Cost
  fixed_arg = list(x_id = x_id,
                   w_id = w_id,
                   sym = ICA.control$sym ,
                   sym_point = ICA.control$sym_point,
                   CostColumnId = CostColumnId,
                   crfunc = arg$crfunc,
                   lp = arg$lp_nofixed, ## NULL for locally and optim_on_average
                   up = arg$up_nofixed, ## NULL for locally and optim_on_average
                   fixedpar = arg$fixedpar,
                   fixedpar_id = arg$fixedpar_id,
                   optim_func = optim_func,
                   npred = npred,
                   type = type,
                   equal_weight = ICA.control$equal_weight,
                   k = arg$k,
                   Calculate_Cost = Calculate_Cost_minimax,
                   is.only.w = arg$is.only.w)

  if (type == "robust")
    fixed_arg$parset <- arg$robpars$parset

  ## for sensitivity checking
  sens_varlist <-list(fixedpar = arg$fixedpar, fixedpar_id = arg$fixedpar_id,
                      npred = npred,
                      crfunc_sens = arg$crfunc_sens,
                      lp_nofixed = arg$lp_nofixed,
                      up_nofixed = arg$up_nofixed,
                      plot_3d = "lattice",
                      npar = arg$npar,
                      optim_starting = optim_starting,
                      fimfunc_sens = arg$FIM_sens,
                      Psi_x_minus_minimax = arg$Psi_funcs$Psi_x_minus_minimax, Psi_x = arg$Psi_funcs$Psi_x,
                      Psi_xy = arg$Psi_funcs$Psi_xy, Psi_Mu = arg$Psi_funcs$Psi_Mu,
                      is.only.w = arg$is.only.w)
  ############################################################################*

  ############################################################################*
  # Initialization when evol is NULL
  ############################################################################*
  if (is.null(evol)){
    evol.flag = 0
    ## set the old seed if call is from minimax
    if (!is.null(ICA.control$rseed))
      set.seed(ICA.control$rseed)

    msg <- NULL
    revol_rate <- ICA.control$revol_rate
    maxiter <- iter
    totaliter <- 0
    #evol <- list()
    min_cost <- c() ## cost of the best imperialists
    mean_cost <- c() ## mean cost of all imperialists
    check_counter <- 0 ## counter to count the check
    total_nlocal  <- 0 ## total number of successful local search
    if (!ICA.control$lsearch)
      total_nlocal <- NA
    total_nrevol <- 0 ## total number of successful revolution
    total_nimprove <- 0 ##total number of improvements due to assimilation
    prev_time <- 0
    ############################################## Initialization for ICA
    InitialCountries <- GenerateNewCountry(NumOfCountries = ICA.control$ncount,
                                           lower = arg$ld,
                                           upper = arg$ud,
                                           sym = ICA.control$sym,
                                           w_id = w_id,
                                           x_id = x_id,
                                           npred= npred,
                                           equal_weight = ICA.control$equal_weight,
                                           is.only.w = arg$is.only.w)

    if (!is.null(arg$initial))
      InitialCountries[1:dim(arg$initial)[1], ] <- arg$initial
    InitialCost <- vector("double", ICA.control$ncount)
    temp <- fixed_arg$Calculate_Cost(mat = InitialCountries, fixed_arg = fixed_arg)
    total_nfeval <-  temp$nfeval
    InitialCost <-  temp$cost
    inparam <- temp$inner_optima
    ## waring inparam for optim_on_average does not have any meaning!
    temp <- NA # safety
    ##Now we should sort the initial countries with respect to their initial cost
    SortInd <- order(InitialCost)
    InitialCost <- InitialCost[SortInd] # Sort the cost in assending order. The best countries will be in higher places
    InitialCountries <- InitialCountries[SortInd,, drop = FALSE] #  Sort the population with respect to their cost. The best country is in the first column
    inparam <- inparam[SortInd, , drop = FALSE]
    # creating empires
    Empires <- CreateInitialEmpires(sorted_Countries = InitialCountries,
                                    sorted_Cost = InitialCost,
                                    Zeta = ICA.control$zeta,
                                    NumOfInitialImperialists = ICA.control$nimp,
                                    NumOfAllColonies = (ICA.control$ncount - ICA.control$nimp),
                                    sorted_InnerParam = inparam)
    best_imp_id<- 1 ## the index of list in which best imperialists is in.
    ########################################################################*
  }
  ##########################################################################*

  ##########################################################################*
  # when we are updating the object for more number of iterations
  ##########################################################################*
  if (!is.null(evol)){
    evol.flag = 1
    ## reset the seed!
    if (exists(".Random.seed")){
      GlobalSeed <- get(".Random.seed", envir = .GlobalEnv)
      #if you call directly from update and not minimax!
      on.exit(assign(".Random.seed", GlobalSeed, envir = .GlobalEnv))
    }

    msg <- object$best$msg
    prev_iter <- length(evol) ##previous number of iterationst
    maxiter <- iter + prev_iter
    totaliter <- prev_iter
    mean_cost <- sapply(1:(totaliter), FUN = function(j) evol[[j]]$mean_cost)
    min_cost <- sapply(1:(totaliter), FUN = function(j) evol[[j]]$min_cost)
    Empires <- object$empires
    prev_time <- arg$time

    check_counter <- arg$updating$check_counter
    total_nfeval <- object$alg$nfeval
    total_nlocal <- object$alg$nlocal
    total_nrevol <- object$alg$nrevol
    total_nimprove <- object$alg$nimprove

    imp_cost <- round(sapply(object$empires, "[[", "ImperialistCost"), 12)
    best_imp_id<- which.min(imp_cost)
    revol_rate <-  arg$updating$revol_rate
    ##updating the random seed
    if (!is.null(ICA.control$rseed)){
      do.call("RNGkind",as.list(arg$updating$oldRNGkind))  ## must be first!
      assign(".Random.seed", arg$updating$oldseed , .GlobalEnv)
    }
  }
  ##########################################################################*
  space_size <- arg$ud - arg$ld
  continue = TRUE

  ###########################################################################*
  ### start of the while loop until continue == TRUE
  ###########################################################################*
  while (continue == TRUE){
    totaliter <- totaliter + 1
    check_counter <- check_counter + 1
    # if (totaliter == 1058)
    #   browser()
    revol_rate <- ICA.control$damp * revol_rate
    ## revolution rate is increased by damp ration in every iter

    #########################################################################*
    ################################################################# for loop over all empires[ii]
    for(ii in 1:length(Empires)){

      ########################################## local search is only for point!
      if (ICA.control$lsearch){


        LocalSearch_res <- LocalSearch(TheEmpire =  Empires[[ii]],
                                       lower = arg$ld,
                                       upper = arg$ud,
                                       l = ICA.control$l,
                                       fixed_arg = fixed_arg)
        Empires[[ii]] <- LocalSearch_res$TheEmpire
        total_nfeval <- total_nfeval + LocalSearch_res$nfeval
        total_nlocal <- total_nlocal + LocalSearch_res$n_success
      }
      ###################################################################*

      ############################################################## Assimilation
      temp5 <- AssimilateColonies2(TheEmpire = Empires[[ii]],
                                   AssimilationCoefficient = ICA.control$assim_coeff,
                                   VarMin = arg$ld,
                                   VarMax = arg$ud,
                                   ExceedStrategy = "perturbed",
                                   sym = ICA.control$sym,
                                   AsssimilationStrategy = ICA.control$assim_strategy,
                                   MoveOnlyWhenImprove = ICA.control$only_improve,
                                   fixed_arg = fixed_arg,
                                   w_id = w_id,
                                   equal_weight = ICA.control$equal_weight)
      ##Warning: in this function the colonies position are changed but the imperialist and the
      ##cost functions of colonies are not updated yet!
      ##they will be updated after revolution
      Empires[[ii]] <- temp5$TheEmpire
      total_nfeval <- total_nfeval + temp5$nfeval
      total_nimprove <-  total_nimprove + temp5$nimprove
      ##########################################################################*

      ############################################################### Revolution
      temp4 <- RevolveColonies(TheEmpire = Empires[[ii]],
                               RevolutionRate = revol_rate,
                               NumOfCountries = ICA.control$ncount,
                               lower = arg$ld,
                               upper = arg$ud,
                               sym = ICA.control$sym,
                               sym_point = ICA.control$sym_point,
                               fixed_arg = fixed_arg,
                               w_id = w_id,
                               equal_weight = ICA.control$equal_weight)
      Empires[[ii]] <- temp4$TheEmpire
      total_nrevol <- total_nrevol + temp4$nrevol
      total_nfeval <- total_nfeval + temp4$nfeval
      ############################################################*
      Empires[[ii]] <- PossesEmpire(TheEmpire = Empires[[ii]])

      ##after updating the empire the total cost should be updated
      ## Computation of Total Cost for Empires
      Empires[[ii]]$TotalCost <- Empires[[ii]]$ImperialistCost + ICA.control$zeta * mean(Empires[[ii]]$ColoniesCost)

    }
    ############################################################ end of the loop for empires [[ii]]
    ###############################################################################################*

    #################################################### Uniting Similiar Empires
    if (length(Empires)>1){

      Empires <- UniteSimilarEmpires(Empires = Empires,
                                     Zeta = ICA.control$zeta,
                                     UnitingThreshold = ICA.control$uniting_threshold,
                                     SearchSpaceSize = space_size)
    }
    ############################################################################*
    # zeta is necessary to update the total cost!
    Empires <- ImperialisticCompetition(Empires = Empires, Zeta = ICA.control$zeta)

    ############################################################## save the seed
    # we get the seed here because we dont know if in cheking it wil be chaned
    #we save the seed when we exit the algorithm
    oldseed <- get(".Random.seed", envir = .GlobalEnv)
    oldRNGkind <- RNGkind()
    ############################################################################*
    ############################################################################*
    # extracing the best emperor and its position
    imp_cost <- round(sapply(Empires, "[[", "ImperialistCost"), 12)
    min_cost[totaliter] <- switch(type, "minimax" = min(imp_cost), "standardized" = -min(imp_cost),
                                  "locally" =  min(imp_cost), "robust" = min(imp_cost))
    mean_cost[totaliter] <- switch(type, "minimax" = mean(imp_cost), "standardized" = -mean(imp_cost),
                                   "locally" = mean(imp_cost), "robust" = mean(imp_cost))

    best_imp_id <- which.min(imp_cost) ## which list contain the best imp
    if (!ICA.control$equal_weight)
      w <- Empires[[best_imp_id]]$ImperialistPosition[, w_id] else
        w <- w_equal
    if (!arg$is.only.w)
      x <- Empires[[best_imp_id]]$ImperialistPosition[, x_id] else
        x <- arg$only_w_varlist$x
    inparam <- Empires[[best_imp_id]]$ImperialistInnerParam
    if (length(arg$lp_nofixed)==1)
      inparam <- t(inparam)


    ##modifying the answering set if there is any fixed parameters.
    ## does not applicable for locally and optim_on_average
    if (any(!is.na(arg$fixedpar))){
      fix_inparam <- c(arg$fixedpar, inparam)
      NumOfParam <- 1:length(fix_inparam)
      inparam <- fix_inparam[order( c(arg$fixedpar_id, setdiff(NumOfParam, arg$fixedpar_id)))]
    }

    if (ICA.control$sym){
      x_w <- ICA_extract_x_w(x = x, w = w, sym_point = ICA.control$sym_point)
      x <- x_w$x
      w <- x_w$w
    }
    ##sort Point

    if (!arg$is.only.w){
      # do not sort if you have only weights because they correspond the points
      if (npred == 1){
        w <- w[order(x)]
        x <- sort(x)
      }
    }
    ############################################################################*

    ################################################################ print trace
    if (ICA.control$trace){
      # if (type != "locally" && type != "robust")
      #   cat("\nICA iter:", totaliter, "\nDesign Points:\n", x, "\nWeights: \n", w,
      #       "\nCriterion value: ", min_cost[totaliter], "\nparam: ",
      #       inparam,"\n") else
      #         cat("\nICA iter:", totaliter, "\nDesign Points:\n", x, "\nWeights: \n", w,
      #             "\nbest criterion value: ", min_cost[totaliter],"\n")

      if (!arg$is.only.w)
        cat("\nIteration:", totaliter, "\nDesign Points:\n", x, "\nWeights: \n", w,
            "\nCriterion value: ", min_cost[totaliter],
            "\nTotal number of function evaluations:", total_nfeval, "\nTotal number of successful local search moves:", total_nlocal,
            "\nTotal number of successful revolution moves:", total_nrevol, "\n") else
              cat("\nIteration:", totaliter, "\nWeights: \n", w,
                  "\nCriterion value: ", min_cost[totaliter],
                  "\nTotal number of function evaluations:", total_nfeval, "\nTotal number of successful local search moves:", total_nlocal,
                  "\nTotal number of successful revolution moves:", total_nrevol, "\n")

      if (ICA.control$only_improve)
        cat("Total number of successful assimilation moves:", total_nimprove, "\n")
      if (type == "minimax")
        cat( "Vector of maximum parameter values: ", inparam,"\n")
      if (type == "standardized")
        cat( "Vector of minimum parameter values: ", inparam,"\n")
    }
    ############################################################################*

    if ( min_cost[totaliter] == 1e-24)
      warning("Computational issue! maybe the design is singular!\n")

    ################################################################### continue
    if (totaliter ==  maxiter){
      continue <- FALSE
      convergence = "Maximum_Iteration"
    }
    if(length(Empires) == 1 && ICA.control$stop_rule == "one_empire"){
      continue <- FALSE
      convergence = "one_empire"
    }
    ## the continue also can be changed in check
    ############################################################################*

    ################################################################################# plot_cost
    if (ICA.control$plot_cost) {
      #ylim for efficiency depends on the criterion type because
      ylim = switch(type,
                    "minimax" = c(min(min_cost) - .07, max(mean_cost[1:(totaliter)]) + .2),
                    "standardized" = c(min(mean_cost[1:(totaliter)]) - .07, max( min_cost) + .2),
                    "locally" = c(min(min_cost) - .07, max(mean_cost[1:(totaliter)]) + .2))
      PlotEffCost(from = 1,
                  to = (totaliter),
                  AllCost = min_cost, ##all criterion up to now (all cost function)
                  UserCost = NULL,
                  DesignType = type,
                  col = line_col[1],
                  xlab = "Iteration",
                  ylim = ylim,
                  lty = 1,
                  title1 = title1,
                  plot_main = TRUE)
      ICAMean_line <-  mean_cost[1:(totaliter)]
      lines(x = 1:(totaliter),
            y = ICAMean_line,
            col = line_col[2], type = "s", lty = 5)
      legend(x = legend_place,  legend = legend_text,lty = c(1,5, 3), col = line_col, xpd = TRUE, bty = "n")
    }
    ############################################################################*

    ####################################################################################*
    #check the equivalence theorem and find ELB
    ####################################################################################*
    ## we check the quvalence theorem in the last iteration anyway. but we may not plot it.
    if (check_counter == ICA.control$checkfreq || (check_last && !continue)){
      check_counter <- 0
      if (arg$ICA.control$trace){
        #cat("\n*********************************************************************")
        if (!continue)
          cat("\nOptimization is done!\n")
        cat("Requesting design verification by the general equivalence theorem\n")
      }
      # if (type == "robust")
      #   type1 <- "locally" else
      #     type1 <- type
      ## Note: we pass the localdes here but we dont use it
      sens_res <- sensminimax_inner(x = x, w = w, lx = arg$lx, ux = arg$ux, lp = arg$lp_nofixed, up = arg$up_nofixed,
                                    fimfunc = arg$FIM,
                                    sens.minimax.control = sens.minimax.control,
                                    sens.control = sens.control,
                                    #nloptr.control.sens = nloptr.control.sens,
                                    type = type, localdes = NULL, plot_sens = ICA.control$plot_sens,
                                    varlist = sens_varlist, calledfrom = "iter",
                                    npar = arg$npar,
                                    calculate_criterion = FALSE,
                                    robpars = arg$robpars,
                                    plot_3d = arg$plot_3d,
                                    silent = !arg$ICA.control$trace)


      ##########################################################################*
      GE_confirmation <- ( sens_res$ELB >= ICA.control$stoptol)
      # print trace that is related to checking
      # if (ICA.control$trace)
      #   cat("maximum of sensitivity:", sens_res$max_deriv, "\nefficiency lower bound (ELB):", sens_res$ELB, "\n")
      if (GE_confirmation && ICA.control$stop_rule == "equivalence"){
        continue <- FALSE
        convergence <- "equivalence"
      }
    }else
      sens_res <- NULL
    # max_deriv <- answering <- answering_cost <-all_optima <- all_optima_cost  <- mu <- ELB <- NA
    # if (type == "locally" || type == "robust"){
    #   answering <- NA # now we dont need answering. We required it before for checking so we set it to NA
    #   mu <- 1
    # }
    ####################################################################### end of check
    ####################################################################### save
    # evol[[totaliter]] <- list(iter = totaliter,
    #                           x = x,
    #                           w = w,
    #                           min_cost = min_cost[totaliter],
    #                           mean_cost = mean_cost[totaliter],
    #                           all_optima = sens_res$all_optima,
    #                           all_optima_cost = sens_res$all_optima_cost,
    #                           answering = sens_res$answering,
    #                           answering_cost = sens_res$answering_cost,
    #                           mu = sens_res$mu,
    #                           max_deriv = sens_res$max_deriv,
    #                           ELB = sens_res$ELB)
    evol[[totaliter]] <- list(iter = totaliter, x = x, w = w, min_cost = min_cost[totaliter], mean_cost = mean_cost[totaliter], sens = sens_res)

    if (type != "locally" && type != "robust"){
      evol[[totaliter]]$param = inparam
    } else
      evol[[totaliter]]$param = NA
    ############################################################################*

    ################################################################ print trace
    # if (ICA.control$trace){
    #   cat("total local search:", total_nlocal, "\n")
    #   cat("total revolution:", total_nrevol, "\n")
    #   if (ICA.control$only_improve)
    #     cat("total improve:", total_nimprove, "\n")
    # }
    ############################################################################*
  }
  ##############################################################################*
  ### end of the while loop over continue == TRUE
  ##############################################################################*

  if (!ICA.control$only_improve)
    total_nimprove <- NA

  #if (ELB >= control$stoptol && control$stop_rule == "equivalence")
  #  convergence = "equivalence" else

  ##############################################################################*
  # check the appropriateness of the maxeval
  # if (type != "locally" & type != "robust" & control$inner_space == "continuous"){
  #   if (control$check_inner_maxeval){
  #     check_temp <- check_maxeval(fn = crfunc2, lower = lp, upper = up, maxeval = control$inner_maxeval,
  #                                 fixedpar = fixedpar, fixedpar_id = fixedpar_id, npred = npred, q = c(x, w))
  #     msg <- check_temp$msg
  #   }else
  #     msg <- NULL
  # }

  ##################*
  msg <- NULL
  ##############################################################################*

  ######################################################################## saving
  ## we add the following to arg becasue dont want to document it in Rd files
  # updating parameters
  object$arg$updating$check_counter <- check_counter
  object$arg$updating$oldseed <- oldseed
  object$arg$updating$oldRNGkind <- oldRNGkind
  object$arg$updating$revol_rate = revol_rate ## different from revolrate
  object$evol <- evol
  object$empires <- Empires

  if(evol.flag==1){
    dout = NULL
    fout = NULL
    fiter = length(object$evol)
    if(fiter>0){
      tout = c()
      for(i in 1:fiter){
        sout = object$evol[[i]]
        tout = rbind(tout,c(i,sout$x,sout$w,sout$min_cost,sout$mean_cost))
      }
      dout = as.data.frame(tout)
      if(length(sout$x)==length(sout$w)){
        dimnames(dout)[[2]] = c("iter",paste("x",1:object$arg$k,sep=""),paste("w",1:object$arg$k,sep=""),"min_cost","mean_cost")
      }else{
        tlen = length(sout$x)/length(sout$w)
        tdimx = c()
        for(i in 1:tlen){
          tdimx = c(tdimx,paste(paste("x",i,sep=""),1:object$arg$k,sep=""))
        }
        dimnames(dout)[[2]] = c("iter",tdimx,paste("w",1:object$arg$k,sep=""),"min_cost","mean_cost")
      }
      # extract sens outcomes
      sout = object$evol[[fiter]]
      tsens = capture.output(sout$sens)
      tpos.max = grep("Maximum",tsens)
      tpos.elb = grep("ELB",tsens)
      tpos.time = grep("Verification",tsens)
      tmax = NA
      telb = NA
      ttime = NA
      if(length(tpos.max)>0){
        tmax = as.numeric(strsplit(gsub("\\s+","",tsens[tpos.max]),"is")[[1]][2])
      }
      if(length(tpos.elb)>0){
        telb = as.numeric(strsplit(gsub("\\s+","",tsens[tpos.elb]),"is")[[1]][2])
      }
      if(length(tpos.time)>0){
        ttime = as.numeric(strsplit(gsub("seconds!","",gsub("\\s+","",tsens[tpos.time])),"required")[[1]][2])
      }
      t2out = c(fiter,sout$x,sout$w,sout$min_cost,sout$mean_cost,tmax,telb,ttime)
      fout = as.data.frame(t(t2out))
      #dimnames(fout)[[2]] = c("iter",paste("x",1:object$arg$k,sep=""),paste("w",1:object$arg$k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
      if(length(sout$x)==length(sout$w)){
        dimnames(fout)[[2]] = c("iter",paste("x",1:object$arg$k,sep=""),paste("w",1:object$arg$k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
      }else{
        tlen = length(sout$x)/length(sout$w)
        tdimx2 = c()
        for(i in 1:tlen){
          tdimx2 = c(tdimx2,paste(paste("x",i,sep=""),1:object$arg$k,sep=""))
        }
        dimnames(fout)[[2]] = c("iter",tdimx2,paste("w",1:object$arg$k,sep=""),"min_cost","mean_cost","max_sens","elb","time_sec")
      }
    }

    object$out = dout
    #out$out2 = out$evol
    object$design = fout
  }

  object$alg <- list(
    nfeval = total_nfeval,
    nlocal = total_nlocal,
    nrevol = total_nrevol,
    nimprove = total_nimprove,
    convergence = convergence)
  #msg = msg

  ## so object is 'res', 'arg' and 'evol'
  ## arg has a list named update as well
  ###### end of saving
  ##############################################################################*
  object$arg$time <- proc.time() - arg$time_start + prev_time
  return(object)
}
######################################################################################################*
######################################################################################################*
