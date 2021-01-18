#include <r2sundials.h>
//' @title Solving ODE System and Sensitivity Equations
//'
//' @description \code{r2cvodes} sets up necessary structures and calls \code{cvodes()} from SUNDIALS library to solve user defined ODE system \eqn{y' = f(t, y, p)}, \eqn{y(t0) = y0}, where \eqn{p} is a constant parameter vector. If requested, corresponding forward sensitivity equations \eqn{s'[i] = df/dy s[i] + df/dp[i]}, \eqn{s[i](t0) = dy0(p)/dp[i]} (here \eqn{s[i](t)=dy(t)/dp[i]}) can be solved simultaneously with the original ODE system. Root finding and proceeding can be defined as well.
//'
//' @param yv const numeric vector, initial values of state vector (\eqn{y0}). Its length defines the number of equations in the ODE system and is refered hereafter as 'Neq'.
//' @param times const numeric vector, time point values at which the solution is stored
//' @param frhs R function or XPtr pointer to Rcpp function, defines rhs \eqn{f(t, y, p)}, i.e. returns (or fills in place) first derivative vector of length \code{Neq} (cf. details for argument list)
//' @param param any R object (default \code{NULL}), parameter object (can be a vector, list, environment, ...) passed through to user defined functions namely to \code{frhs}. It can be useful to get the code more readable as named parameters are authorized. (Starting from this one, the following parameters are optional)
//' @param tstop const numeric scalar (default \code{as.numeric(c())}) if non empty, defines a time point beyond which the calculation is not made (e.g. out of definition domain). If a vector is given, only the first value is used.
//' @param abstol const numeric scalar (default \code{1.e-8}), absolute tolerance in ODE solving (used for values which are close to 0)
//' @param reltol const double (default \code{1.e-8}), relative tolerance in ODE solving (used for values which are not so close to 0)
//' @param integrator integer scalar (default \code{as.integer(c())}), defines which integration scheme should be used by cvodes: implicit (\code{CV_BDF}) or explicit (\code{CV_ADAMS}). Default empty vector is equivalent to \code{CV_BDF}. Constants \code{CV_BDF} and \code{CV_ADAMS} (as other mentioned constants) are defined by the package and are available for user.
//' @param maxord const integer scalar (default \code{0}), defines maximal order of time scheme. Default 0 values is equivalent to 5 for \code{integrator=CV_BDF} or to 12 for \code{integrator=CV_ADAMS}
//' @param maxsteps const integer scalar (default \code{0}), maximum of internal steps before reaching next time point from 'times'. Default 0 is equivalent to 500.
//' @param hin const numeric scalar (default \code{0}), value of the initial step size to be attempted. Default 0.0 corresponds to \code{cvodes}' default value which is internally estimated at \eqn{t0}.
//' @param hmax const numeric scalar (default \code{0}), maximum absolute value of the time step size (>= 0.0). The default 0 value is equivalent to Inf.
//' @param hmin const numeric scalar (default \code{0}), minimum absolute step size (>= 0.0)
//' @param constraints const numeric vector (default \code{as.numeric(c())}), if non empty, defines constraint flags. \cr If constraints[i] is \cr 0 then no constraint is imposed on \eqn{y[i]}; \cr 1.0 then \eqn{y[i]} will be constrained to be \eqn{y[i]} >= 0.0; \cr -1.0 then \eqn{y[i]} will be constrained to be \eqn{y[i]} =< 0.0.; \cr 2.0 then \eqn{y[i]} will be constrained to be \eqn{y[i]} > 0.0.; \cr -2.0 then \eqn{y[i]} will be constrained to be \eqn{y[i]} < 0.0. \cr If a time step respecting hmin and avoiding constraint violation cannot be found, an error is triggered.
//' @param fjac R function of XPtr pointer to Rcpp function (default \code{NULL}), users supplied function that calculates Jacobian matrix df/dy which can be dense or sparse (cf. details for parameter list). If \code{fjac} is not supplied, a SUNDIALS internal approximation to Jacobian is used when needed.
//' @param nz const integer scalar (default \code{0}), number of non zero elements in Jacobian. If > 0, this parameter triggers a sparse Jacobian usage. In this case, the previous parameter \code{fjac}, must be a function (or a pointer to a function) with appropriate parameter list (cf. details) filling a sparse matrix. It is advised to have \code{nz} accounting not only for non zeros in the Jacobian itself but also for diagonal terms even if they are zeros. In this case, a supplementary memory allocation can be avoided. If a sparse Jacobian is used, corresponding sparse linear system is solved with the help of \pkg{rmumps} package.
//' @param rmumps_perm integer scalar (default \code{as.integer(c())}), defines permutation method that will be used by \pkg{rmumps} during symbolic analysis before solving sparse linear systems. Default value is equivalent to \code{RMUMPS_PERM_AUTO}. Possible values (defined in \pkg{rmumps} package) are \code{RMUMPS_PERM_AMF}, \code{RMUMPS_PERM_AMD}, \code{RMUMPS_PERM_AUTO}, \code{RMUMPS_PERM_QAMD}, \code{RMUMPS_PERM_PORD}, \code{RMUMPS_PERM_METIS}, \code{RMUMPS_PERM_SCOTCH}. An appropriate choice of permutation type can reduce memory requirements as well as calculation speed for sparse system solving. If a vector is supplied, only the first value is read. If \code{nz=0}, this parameter is not used.
//' @param nroot const integer scalar (default \code{0}), defines number of roots that user wishes to track.
//' @param froot R function of XPtr pointer to Rcpp function (default \code{NULL}), user defined function calculating root vector of length \code{nroot}. When at least one component of this vector crosses 0, a root finding event is triggered and corresponding root handling function is called (cf. \code{fevent}). This defines a very flexible root tracking mechanism. User can define a root based both on time values and state vector (y).
//' @param fevent R function of XPtr pointer to Rcpp function (default \code{NULL}), defines function for root proceeding. When a root is encountered, this function is called. Its return value defines what to do with this root. User's reaction (i.e. the return value of \code{fevent}) can be one of three types: 
//'   \cr \code{R2SUNDIALS_EVENT_IGNORE} - do nothing (can be helpful if, for example, 0 is crossed not in pertinent sens); \cr \code{R2SUNDIALS_EVENT_HOLD} - the time point at which the root happened as well as the corresponding \code{rootsfound} vector are added to the root matrix delivered with output attributes (cf. details).
//'   This time point is also added to the ODE solution which can lead to a new time point, originally absent in 'times' vector; \cr \code{R2SUNDIALS_EVENT_STOP} - stops the ODE solving. If it happens before reaching the last value in 'times', the ODE solution can have less time points than initially defined in 'times' parameter.
//'   \code{fevent} is allowed to modify the vector state y in an arbitrary way. It can be helpful for modeling some controlling events for example adding some compound to chemical mix or modifying speed vector if an obstacle is hitted. If such modification takes place, the ODE system is restarted to handle in appropriate way induced discontinuities.
//' @param Ns const integer scalar (default \code{0}), number of parameters for which forward sensitivity system has to be solved.
//' @param psens numeric vector (default \code{as.numeric(c())}), if not empty, defines a vector of parameters for which (or for some components of which) forward sensitivity system has to be solved. If its length > \code{Ns}, then vector \code{psens_list} of length \code{Ns} must define which are the \code{Ns} components of \code{psens} with respect to which sensitivities are to be computed. Note that parameters used in sensitivity calculations must be put in \code{psens} vector (and not somewhere in 'param' object) only if user wish to rely on SUNDIALS internal procedure for estimation of sensitivity rhs. If user provides its own function for sensitivity rhs, he is free to use either of 'param' or \code{psens} for passing and accessing sensitivity parameters.
//' @param sens_init numeric matrix (default \code{as.numeric(c())}), matrix with \code{Neq} rows and \code{Ns} columns defining initial condition for sensititivy system. Default value is equivalent to \code{Neq} x \code{Ns} zero matrix which means that \eqn{y0} is not dependent on \eqn{p}.
//' @param psens_bar numeric vector (default \code{as.numeric(c())}), a vector of \code{Ns} positive scaling factors. The values provided in this vector, describe orders of magnitude for \code{psens} components. Having such estimations can help to increase the accuracy for sensitivity vector calculations or for the right-hand side calculation for the sensitivity systems based on the internal difference-quotient function. If not defined, it is equivalent to set all \code{psens_bar} components to 1
//' @param psens_list const integer vector (default \code{as.integer(c())}), if \code{length(psens) > Ns}, this index vector can be used to indicate which are components of \code{psens} that are concerned with sensitivity calculations. Default value is valid only when \code{length(psens)==Ns}, i.e. \code{psens_list} is not used. Values in \code{psens_list} are 1-based (as regular integer indexes in R)
//' @param fsens R function of XPtr pointer to Rcpp function (default \code{NULL}), user defined sensitivity rhs for all \code{psens} at once.
//' @param fsens1 R function of XPtr pointer to Rcpp function (default \code{NULL}), user defined sensitivity rhs for \code{is}-th \code{psens} component (\code{is} is passed as input parameter to \code{fsens1}, cf. details). Only one of \code{fsens} or \code{fsens1} can be supplied. If none of them is supplied, SUNDIALS' internal difference-quotient function is used. In this case, sensitivity parameters must be passed only through \code{psens} argument.
//' @param sens_method integer scalar (default \code{as.integer(c())}), constant defining the method used for solving sensitivity system. Allowed values are \code{CV_SIMULTANEOUS} or \code{CV_STAGGERED}. Default value is equivalent to \code{CV_SIMULTANEOUS}. \itemize{
//' \item \code{CV_SIMULTANEOUS} means that the state and sensitivity variables are corrected at the same time.
//' \item \code{CV_STAGGERED} means that the correction step for the sensitivity variables takes place at the same time for all sensitivity equations, but only after the correction of the state variables has converged and the state variables have passed the local error test
//'}
//' @param errconS constant logical scalar (default \code{TRUE}), specifies whether sensitivity variables are to be included (\code{TRUE}) or not (\code{FALSE}) in the error control mechanism
//'
//' @return numeric matrix, ODE solution where each column corresponds to a state vector at a given time point. The columns (their number is refered to as \code{Nt}) are named by time points while the rows heritates the names from \code{yv}. If no names are found in \code{yv}, the rows are simply named 'V1', 'V2' and so on. After a normal execution and without root handling, column number is equal to the length of \code{times}. However, if root handling is used, it can add or remove some time points from  \code{times}. So the user must not assume that column number of output matrix is equal to \code{length(times)}. Instead, actual number of time points for which the solution was calculated can be retrieved from an attribute named "times".
//' Moreover, several attributes are defined in the returned matrix. We have mentioned "times", the others are: \describe{
//' \item{stats}{Some statistics about various events that could happen during \code{cvodes} run like the number of rhs calls, Jacobian calls, number of internal time steps, failure number and so on. Any component <name> in stats vector corresponds to SUNDIALS function pattern CVodeGet<name>, e.g. "NumRhsEvals" was obtained with CVodeGetNumRhsEvals() call. For detailed meaning of each statistics, user is invited to refer to \href{https://computation.llnl.gov/sites/default/files/public/cvs_guide.pdf}{SUNDIALS documentation};}
//' \item{roots}{matrix with row number \code{nroot+1} and column number equal to number of roots found by the \code{cvodes()} and retained by the user. Each column is a composite vector made of time point and \code{rootsfound} vector described here-before.}
//' \item{sens}{sensitivity 3D array with dimensions \code{Neq} x \code{Nt} x \code{Ns}}
//' }
//'
//' @details The package \pkg{r2sundials} was designed to avoid as much as possible memory reallocation in callback functions (\code{frhs} and others). C++ variants of these functions are fully compliant with this design principle. While R counterparts are not (as per R design). Here, we define callback function interfaces that user has to abide to. Pointers to C++ variants to be passed to \code{r2cvodes()} can be obtained with the help of \pkg{RcppXPtrUtils}. See examples for illustrations of such use.
//' \cr Right hand side function \code{frhs} provided by user calculates derivative vector \eqn{y'}. This function can be defined as classical R function or a Rcpp/RcppArmadillo function. In the first case, it must have the following list of input arguments\cr
//' \code{frhs(t, y, param, psens)}\cr
//' and return a derivative vector of length \code{Neq}. Here \code{t} is time point (numeric scalar), \code{y} current state vector (numeric vector of length \code{Neq}), \code{param} and \code{psens} are passed through from \code{r2cvodes()} arguments.
//' \cr In the C++ case, it is defined as\cr
//' \code{int (*frhs)(double t, const vec &y, vec &ydot, RObject &param, NumericVector &psens)}\cr
//' and return an integer status flag, e.g. \code{CV_SUCCESS}. For other possible status flags see the original \href{https://computation.llnl.gov/sites/default/files/public/cvs_guide.pdf}{SUNDIALS documentation}. The derivatives are stored in-place in \code{ydot} vector. See examples section for a usage sample.
//' \cr \code{fjac} is a function calculating Jacobian matrix. Its definition varies depending on 1) kind of used Jacobian: dense or sparse and 2) on programming language used: R or C++ (i.e. Rcpp/RcppArmadillo).
//' \itemize{
//' \item{For dense Jacobian calculated in R, the arguments are:\cr
//'   \code{fjac(t, y, ydot, param, psens)}\cr
//'   and the expected return value is \code{Neq} x \code{Neq} dense Jacobian matrix df/dy.}
//' \item{For dense Jacobian calculated in C++ the definition is following:\cr
//'   \code{int (*fjac)(double t, const vec &y, const vec &ydot, mat &J, RObject &param, NumericVector &psens, vec &tmp1, vec &tmp2, vec &tmp3)}\cr
//'   It must return a status flag. The resulting Jacobian is stored in-place in the the matrix \code{J}. Auxiliary vectors \code{tmp1} to \code{tmp3} are of length \code{Neq} and are made available for intermediate storage thus avoiding memory reallocation at each call to \code{fjac()}.}
//' \item{For sparse Jacobian calculated in R, the arguments are:\cr
//'   \code{fjac(t, yv, ydotv, param, psens)}\cr
//'   The return value is a list with fields \code{i} (row indices), \code{p} (column pointers) and \code{v} (matrix values) defining the content of sparse Jacobian in CSC (condensed sparse column) format. The values stored in i and p vectors are supposed to be 1-based, as it is common in R language.}
//' \item{For sparse Jacobian calculated in C++ the definition is following:\cr
//'   \code{int (*fjac)(double t, const vec &y, const vec &ydot, uvec &i, uvec &p, vec &v, int n, int nz, RObject &param, NumericVector &psens, vec &tmp1, vec &tmp2, vec &tmp3)}\cr
//'   here \code{n=Neq}, \code{nz} is passed through from r2cvodes() arguments. The resulting sparse Jacobian is stored in-place in vectors \code{i}, \code{p}, \code{v} corresponding to the CSC (Compressed Sparse Column) format. Their respective dimensions are \code{nz}, \code{n+1} and \code{nz}. The values stored in \code{i} and \code{p} must be 0 based as per usage in C++. The return value is a status flag.}
//' }
//' \code{froot} calculates a root vector, i.e. a vector whose components are tracked for 0 crossing during the time course in ODE solving. If written in R, its call follows the following pattern:\cr
//' \code{froot(t, y, param, psens)}\cr
//' and it must return a numeric vector of length 'nroot'. If written in C++, it is defined as\cr
//' \code{int (*froot)(double t, const vec &y, vec &vroot, RObject &param, NumericVector &psens)}\cr
//' The tracked values are stored in-place in \code{vroot}. The returned value is a status flag.
//' \cr \code{fevent} handles the event of root finding. If written in R, the calling pattern is\cr
//' \code{fevent(t, yvec, Ns, ySm, rootsfound, param, psens)}\cr
//' and the return value is a list with named component "ynew", "flag" and "sens_init". The last item is required only when \code{Ns > 0}. Current value of sensitivity matrix \eqn{dy/dp} can be found in parameter \code{ySm}. Integer vector 'rootsfound' of length 'nroot' provides information on \code{vroot} components that triggered the root event. If \code{rootsfound[i] != 0}, it means that \code{vroot[i]} is a root otherwise it is not. Moreover, the sign of \code{rootsfound[i]} is meaningful. If \code{rootsfound[i] > 0} then \code{vroot[i]} is increasing at 0 crossing. Respectively, if \code{rootsfound[i] < 0} then \code{vroot[i]} is decreasing. The vector 'ynew' in the output list can define a new state vector after event handling (for example, an abrupt change in velocity direction and/or magnitude after an obstacle hit). The field 'flag' in the output list is authorized to take only three values: \code{R2SUNDIALS_EVENT_IGNORE}, \code{R2SUNDIALS_EVENT_HOLD} and \code{R2SUNDIALS_EVENT_STOP} described here-before. The matrix \code{sens_init} is used for a possible restarting of sensitivity calculation. It must contain a derivative matrix \eqn{dynew/dp} of size \code{Neq x Ns}. If this field is absent (NULL) then a zero matrix is assumed.
//' \cr If written in C++, this function is defined as\cr
//' \code{int (*fevent)(double t, const vec &y, vec &ynew, int Ns, std::vector<vec> &ySv, const ivec &rootsfound, RObject &param, NumericVector &psens)}\cr
//' The new state vector can be stored in-place in \code{ynew} and the status flag indicating what to do with this event is the return value. If sensitivity calculation is going on (i.e. \code{Ns > 0}), the current value of sensitivity vectors can be found in \code{ySv} and their new values can be stored in-place in the same parameter.
//' Note that if \code{ynew} is different from the vale of \code{y} when the root was found, the ODE is restarted from this time point to handle correctly the discontinuity. As a result, there will be two columns corresponding to the same time point: one with the state vector at root finding and one with \code{ynew} values. The same is true for sensitivity output, if it is part of the problem.
//' \cr \code{fsens} calculates rhs for sensitivity system. If written in R, it must be defined as\cr
//' \code{fsens(Ns, t, y, ydot, ySm, param, psens)}\cr
//' and return a dataframe in which i-th column correspond to \eqn{s'[i]} sensitivity derivative vector. Among other parameters, it receives \code{ySm} which is a \code{Neq} x \code{Ns} matrix having the current values of sensitivity vector (i-th vector is in i-th column).
//' \cr If written in C++, it has to be defined as\cr
//' \code{int (*fsens)(int Ns, double t, const vec &y, const vec &ydot, const std::vector<vec> &ySv, std::vector<vec> &ySdotv, RObject &param, NumericVector &psens, const vec &tmp1v, const vec &tmp2v)}\cr
//' Note a slight difference in the input parameters compared with the R counterpart. Here \code{ySv} plays the role of \code{ySm} and is not a matrix but a vector of Armadillo vectors. To access m-th component of \eqn{s[i]}, one can simply do \code{ySv[i][m]} and the whole \eqn{s[i]} is selected as \code{ySv[i]}. Such data structure was retained to keep as low as possible new memory reallocation. The resulting sensitivity derivatives are to be stored in-place in \code{ySdotv} according to the same data organization scheme as in \code{ySv}. This function returns a status flag.
//' \cr \code{fsens1} does the same as \code{fsens} but provides derivatives of sensitivity vectors on one-by-one basis. This second form is provided for user's convenience as in some cases the code can become more readable if it calculates only one vector s'[i] at a time. If written in R, this function has to be defined as\cr
//' \code{fsens1(Ns, t, y, iS, ydot, yS, param, psens)}\cr
//' here \code{iS} is the index of calculated vector \eqn{s'[iS]} and \code{yS} contains the current value of \eqn{s[iS]}.
//' If written in C++, this functions has to be defined as\cr
//' \code{int (*fsens1)(int Ns, double t, const vec &yv, const vec &ydotv, int iS, const vec &ySv, vec &ySdotv, RObject &param, NumericVector &psens, const vec &tmp1v, const vec &tmp2v)}\cr
//' The result, i.e. \eqn{s'[iS]} is to be stored in-place in \code{ySdotv} vector. This function returns a status flag.
//' @examples
//' # Ex.1. Solve a scalar ODE describing exponential transition form 0 to 1
//' # y'=-a*(y-1), y(0)=0, a is a parameter that we arbitrary choose to be 2.
//' # define rhs function (here in R).
//' frhs_exp=function(t, y, p, psens) -p$a*(y-1)
//' # define parameter list
//' p=list(a=2)
//' # define time grid from 0 to 5 (arbitrary units)
//' ti=seq(0, 5, length.out=101)
//' # define initial state vector
//' y0=0
//' # we are set for a very simple r2cvodes() call
//' res_exp=r2sundials::r2cvodes(y0, ti, frhs_exp, param=p)
//' # compare the result to theoretical values: 1-exp(-a*t)
//' stopifnot(diff(range(1-exp(-p$a*ti) - res_exp)) < 1.e-6)
//' \donttest{
//' # Ex. 2. Same problem but frhs is written in C++
//' library(RcppXPtrUtils)
//' ptr_exp=cppXPtr(code='
//' int rhs_exp(double t, const vec &y, vec &ydot, RObject &param, NumericVector &psens) {
//'   NumericVector p(param);
//'   ydot[0] = -p["a"]*(y[0]-1);
//'   return(CV_SUCCESS);
//' }
//' ', depends=c("RcppArmadillo","r2sundials","rmumps"),
//'  includes="using namespace arma;\n#include <r2sundials.h>", cacheDir="lib", verbose=FALSE)
//' # For ease of use in C++, we convert param to a numeric vector instead of a list.
//' pv=c(a=p$a)
//' # new call to r2cvodes() with XPtr pointer ptr_exp.
//' res_exp2=r2sundials::r2cvodes(y0, ti, ptr_exp, param=pv)
//' stopifnot(diff(range(res_exp2 - res_exp)) < 1.e-14)
//'
//' # Ex.3. Bouncing ball simulation.
//' # A ball falls from a height y=5 m with initial vertical speed vy=0 m/s
//' # and horizontal speed vx=1 m/s. The forces acting on the ball are 1) the gravity
//' # (g=9.81 m/s^2) and 2) air resistance f=-k_r*v (k_r=0.1 N*s/m).
//' # When the ball hits the ground, it bounces instantly retaining k=0.9 part
//' # of its vertical and horizontal speed. At the bounce, the vertical speed
//' # changes its sign to the opposite while horizontal speed keeps the original sign.
//' # Simulation should stop after the 5-th bounce or at tmax=10 s which ever comes first.
//' # This example illustrates usage of root finding and handling.
//' # We give here an example of callback functions for root handling in C++.
//' yv=c(x=0, y=5, vx=1, vy=0) # initial state vector
//' pv=c(g=9.81, k_r=0.1, k=0.9, nbounce=5, tmax=10) # parameter vector
//' ti=seq(0, 20, length.out=201L) # time grid
//'
//' # rhs
//' ptr_ball=cppXPtr(code='
//' int rhs_ball(double t, const vec &y, vec &ydot, RObject &param, NumericVector &psens) {
//'   NumericVector p(param);
//'   ydot[0] = y[2]; // dx/dt=vx
//'   ydot[1] = y[3]; // dy/dt=vy
//'   ydot[2] = -p["k_r"]*y[2]; // dvx/dt= -k_r*vx
//'   ydot[3] = -p["g"] - p["k_r"]*y[3]; // dvy/dt=-g -k_r*vy
//'   return(CV_SUCCESS);
//' }
//' ', depends=c("RcppArmadillo","r2sundials","rmumps"),
//'  includes="using namespace arma;\n#include <r2sundials.h>", cacheDir="lib", verbose=FALSE)
//' 
//' # root function
//' ptr_ball_root=cppXPtr(code='
//' int root_ball(double t, const vec &y, vec &vroot, RObject &param, NumericVector &psens) {
//'   NumericVector p(param);
//'   vroot[0] = y[1]; // y==0
//'   vroot[1] = t-p["tmax"]; // t==p["tmax"]
//'   return(0);
//' }
//' ', depends=c("RcppArmadillo","r2sundials","rmumps"),
//'  includes="using namespace arma;\n#include <r2sundials.h>", cacheDir="lib", verbose=FALSE)
//' 
//' # event handler function
//' ptr_ball_event=cppXPtr(code='
//' int event_ball(double t, const vec &y, vec &ynew, int Ns, std::vector<vec> &ySv,
//'           const ivec &rootsfound, RObject &param, NumericVector &psens) {
//'   NumericVector p(param);
//'   static int nbounce=0;
//'   if (rootsfound[1] != 0) // time is out
//'      return(R2SUNDIALS_EVENT_STOP);
//'   if (rootsfound[0] > 0)
//'     // cross 0 in ascending sens, can happen when y < 0 in limits of abstol
//'     return(R2SUNDIALS_EVENT_IGNORE);
//'   ynew=y;
//'   if (++nbounce < p["nbounce"]) {
//'     // here nbounce=1:4
//'     ynew[2] *= p["k"]; // horizontal speed is lowered
//'     ynew[3] *= -p["k"]; // vertical speed is lowered and reflected
//'     return(R2SUNDIALS_EVENT_HOLD);
//'   } else {
//'     // here nbounce=5
//'     nbounce=0; // reinit counter for possible next calls to r2cvodes
//'     return(R2SUNDIALS_EVENT_STOP);
//'   }
//' }
//' ', depends=c("RcppArmadillo","r2sundials","rmumps"),
//'  includes="using namespace arma;\n#include <r2sundials.h>", cacheDir="lib", verbose=FALSE)
//'
//' # ODE solving and plotting
//' res_ball <- r2sundials::r2cvodes(yv, ti, ptr_ball, param=pv, nroot=2L,
//'   froot=ptr_ball_root, fevent=ptr_ball_event)
//' plot(res_ball["x",], res_ball["y",], xlab="X [m]", ylab="Y [m]",
//'   t="l", main="Bouncing ball simulation")
//' 
//' # Ex.4. Robertson chemical reactions
//' # This example is often used as an illustration and a benchmark for stiff ODE.
//' # We will demonstrate here:
//' # • how to use sparse Jacobian (not really meaningful for 3x3 sytem but just to give a hint);
//' # • how to make sensitivity calculations with user provided rhs.
//' #
//' # Let simulate the following chemical system of 3 compounds y1, y2 and y3
//' #  y1' = -k1*y1 + k3*y2*y3
//' #  y2' =  k1*y1 - k2*y2*y2 - k3*y2*y3
//' #  y3' =  k2*y2*y2
//' # Jacobian df/dy is
//' # 
//' # | -k1 |      k3*y3       |  k3*y2 |
//' # |-----+------------------+--------|
//' # |  k1 | -2*k2*y2 - k3*y3 | -k3*y2 |
//' # |-----+------------------+--------|
//' # |  0  |     2*k2*y2      |   0    |
//' 
//' 
//' yv <- c(y1=1, y2=0, y3=0) # initial values
//' pv <- c(k1 = 0.04, k2 = 3e7, k3 = 1e4) # parameter vector
//' ti=10^(seq(from = -5, to = 11, by = 0.1)) # exponential time grid
//' 
//' # pointer to rhs function
//' ptr_rob=cppXPtr(code='
//' int rhs_rob(double t, const vec &y, vec &ydot, RObject &param, NumericVector &psens) {
//'   NumericVector p(param);
//'   ydot[0] = -p["k1"]*y[0] + p["k3"]*y[1]*y[2];
//'   ydot[2] = p["k2"]*y[1]*y[1];
//'   ydot[1] = -ydot[0] - ydot[2];
//'   return(CV_SUCCESS);
//' }
//' ', depends=c("RcppArmadillo","r2sundials","rmumps"),
//'  includes="using namespace arma;\n#include <r2sundials.h>", cacheDir="lib", verbose=FALSE)
//' # pointer to sparse jacobian function
//' ptr_rob_jacsp=cppXPtr(code='
//' int spjac_rob(double t, const vec &y, const vec &ydot, uvec &ir, uvec &pj, vec &v, int n, int nz,
//'               RObject &param, NumericVector &psens, vec &tmp1, vec &tmp2, vec &tmp3) {
//'   if (nz < 8)
//'     stop("spjac_rob: not enough room for non zeros, must have at least 8, instead got %d", nz);
//'   NumericVector p(param);
//'   int i=0;
//'   pj[0] = 0; // init pj
//'   // first column
//'   ir[i] = 0;
//'   v[i++] = -p["k1"];
//'   ir[i] = 1;
//'   v[i++] = p["k1"]; 
//'   pj[1] = i;
//'   // second column
//'   ir[i] = 0;
//'   v[i++] = p["k3"]*y[2];
//'   ir[i] = 1;
//'   v[i++] = -p["k3"]*y[2]-2*p["k2"]*y[1];
//'   ir[i] = 2;
//'   v[i++] = 2*p["k2"]*y[1];
//'   pj[2] = i;
//'   // third column
//'   ir[i] = 0;
//'   v[i++] = p["k3"]*y[1];
//'   ir[i] = 1;
//'   v[i++] = -p["k3"]*y[1];
//'   ir[i] = 2;
//'   v[i++] = 0; // just to have the main diagonal fully in Jacobian
//'   pj[3] = i;
//'   return(0);
//' }
//' ', depends=c("RcppArmadillo","r2sundials","rmumps"),
//'  includes="using namespace arma;\n#include <r2sundials.h>", cacheDir="lib", verbose=FALSE)
//' # pointer to sensitivity rhs function
//' ptr_rob_sens1=cppXPtr(code='
//' int sens_rob1(int Ns, double t, const vec &y, const vec &ydot, int iS, const vec &yS, vec &ySdot,
//'               RObject &param, NumericVector &psens, vec &tmp1, vec &tmp2) {
//'   // calculate (df /dy)s_i(t) + (df /dp_i) for i = iS
//'   NumericVector p(param);
//'   // (df/dy)s_i(t)
//'   ySdot[0] = -p["k1"]*yS[0] + p["k3"]*y[2]*yS[1] + p["k3"]*y[1]*yS[2];
//'   ySdot[1] = p["k1"]*yS[0] - (p["k3"]*y[2]+2*p["k2"]*y[1])*yS[1] - p["k3"]*y[1]*yS[2]; 
//'   ySdot[2] = 2*p["k2"]*y[1]*yS[1];
//'   // + (df/dp_i)
//'   switch(iS) {
//'     case 0:
//'       ySdot[0] -= y[0];
//'       ySdot[1] += y[0];
//'       break;
//'     case 1:
//'       ySdot[1] -= y[1]*y[1];
//'       ySdot[2] += y[1]*y[1];
//'       break;
//'     case 2:
//'       ySdot[0] += y[1]*y[2];
//'       ySdot[1] -= y[1]*y[2];
//'   }
//'   return(CV_SUCCESS);
//' }
//' ', depends=c("RcppArmadillo","r2sundials","rmumps"),
//'  includes="using namespace arma;\n#include <r2sundials.h>", cacheDir="lib", verbose=FALSE)
//' # Note that we don't use psens param for sensitivity calculations as we provide our own fsens1.
//' res_rob <- r2sundials::r2cvodes(yv, ti, ptr_rob, param=pv, nz=8, fjac=ptr_rob_jacsp, Ns=3,
//'                               fsens1=ptr_rob_sens1)
//' # plot ODE solution
//' layout(t(1:3)) # three sublots in a row
//' for (i in 1:3)
//'    plot(ti, res_rob[i,], log="x", t="l", xlab="Time", ylab=names(yv)[i])
//' # plot sensitivities
//' layout(matrix(1:9, nrow=3)) # 9 subplots in a square
//' for (j in 1:3) # run through pv
//'    for (i in 1:3) # run through y
//'       plot(ti, attr(res_rob, "sens")[i,,j], log="x", t="l", xlab="Time",
//'            ylab=parse(text=paste0("partialdiff*y[", i, "]/partialdiff*k[", j, "]")))
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix r2cvodes(const NumericVector &yv, const vec &times, const RObject &frhs, RObject param=R_NilValue, const NumericVector tstop=NumericVector::create(), const double abstol=1.e-8, const double reltol=1.e-8, IntegerVector integrator=IntegerVector::create(), const int maxord=0, const int maxsteps=0, const double hin=0., const double hmax=0., const double hmin=0., const vec &constraints=NumericVector::create(), const RObject fjac=R_NilValue, const int nz=0, IntegerVector rmumps_perm=IntegerVector::create(), const int nroot=0, const RObject froot=R_NilValue, const RObject fevent=R_NilValue,
const int Ns=0, NumericVector psens=NumericVector::create(), NumericVector sens_init=NumericVector::create(), NumericVector psens_bar=NumericVector::create(), const IntegerVector psens_list=IntegerVector::create(), const RObject fsens=R_NilValue, const RObject fsens1=R_NilValue, IntegerVector sens_method=IntegerVector::create(), const bool errconS=true) {
  //long clk_tck = CLOCKS_PER_SEC;
  //clock_t t1, t2;
  //t1 = clock();
  Sunmem<int> mem;
  UserData udata;
  realtype t;
  N_Vector nv_y, nv_constraints, *yS;
  SUNMatrix A;
  SUNLinearSolver LS;
  int retval, retevent, i;
  ivec rootsfound, plist;
  mat res, mroots, msens_init, msens;
  vec ti;
  cube asens;
  Function rf_event(Environment::global_env()["ls"]); // just a placeholder
  rsunEventFn user_event_fn=NULL;

  if (integrator.size() == 0)
    integrator=IntegerVector(1, CV_BDF);
  if (rmumps_perm.size() == 0)
    rmumps_perm=IntegerVector(1, RMUMPS_PERM_AUTO);
  if (sens_method.size() == 0)
    sens_method=IntegerVector(1, CV_SIMULTANEOUS);
  if (frhs.sexp_type() != CLOSXP && frhs.sexp_type() != EXTPTRSXP)
    stop("r2sundials: frhs must be either R function or of classes 'XPtr' or 'externalptr' instead got '%s'", Rf_type2char(frhs.sexp_type()));
  if (!fjac.isNULL() && fjac.sexp_type() != CLOSXP && fjac.sexp_type() != EXTPTRSXP)
    stop("r2sundials: if not NULL, fjac must be either R function or of classes 'XPtr' or 'externalptr'");
  if (!fevent.isNULL() && fevent.sexp_type() != CLOSXP && fevent.sexp_type() != EXTPTRSXP)
    stop("r2sundials: if not NULL, fevent must be either R function or of classes 'XPtr' or 'externalptr'");

  if (integrator[0] != CV_ADAMS && integrator[0] != CV_BDF)
    stop("cvode: invalide integrator '%d' (must be either NULL, %d (CV_ADAMS) or %d (CV_BDF default))", integrator[0], CV_ADAMS, CV_BDF);

  udata.lp=List::create(_("param")=param, _("frhs")=frhs, _("fjac")=fjac, _("froot")=froot, _("fevent")=fevent, _("nroot")=nroot, _("fsens")=fsens, _("fsens1")=fsens1);
  udata.psens=psens;
//print(wrap(udata.lp));
  void *cvode_mem, *user_data=static_cast<void*>(&udata);
  
  int neq=yv.size(), nti=times.size();
  res.set_size(neq, nti); // can be resized if roots are added and/or event handling causes premature stop
  ti=times; // can be resized as res
  // copy times in the first row of result matrix
  asens.set_size(neq, nti, Ns);
  std::vector<vec> ySv(Ns);
  
  /* Create serial vector of length neq from yv (initial conditions)*/
  getmem(nv_y, N_VNew_Serial(neq));
  mem.add((void **) &nv_y, (funfree) N_VDestroy);
  vec yvec(NV_DATA_S(nv_y), neq, false), ynew; // vec proxy for y and y after event treatment by user function
  // copy init values
  yvec.subvec(0, neq-1)=vec(yv);
//yvec.print("yvec");
//print(yv);
//Rcout << "NV_DATA_S(nv_y)=" << NV_DATA_S(nv_y) << "\n";
//Rcout << "yvec.memptr()=" << yvec.memptr() << "\n";
  
  // create the solver memory and specify the Backward Differentiation Formula
  getmem(cvode_mem, CVodeCreate(integrator[0]));
  mem.add((void **) &cvode_mem, (funfreep) CVodeFree);
  check_retval(CVodeSetErrHandlerFn(cvode_mem, rsunerr, NULL));
  // Set cvode_mem and put different solver components
  if (tstop.size() > 0)
    check_retval(CVodeSetStopTime(cvode_mem, tstop[0]));
  if (maxord > 0)
    check_retval(CVodeSetMaxOrd(cvode_mem, maxord));
  if (maxsteps > 0)
    check_retval(CVodeSetMaxNumSteps(cvode_mem, maxsteps));
  if (hin > 0.)
    check_retval(CVodeSetInitStep(cvode_mem, hin));
  if (hmin > 0.)
    check_retval(CVodeSetMinStep(cvode_mem, hmin));
  if (hmax > 0.)
    check_retval(CVodeSetMaxStep(cvode_mem, hmax));
  check_retval(CVodeInit(cvode_mem, rhswrap, times[0], nv_y));
  check_retval(CVodeSStolerances(cvode_mem, reltol, abstol));
  check_retval(CVodeSetUserData(cvode_mem, user_data));
  if (constraints.size() != 0 && (int) constraints.size() != neq)
    stop("cvode: if not empty, constraints vector must be of the same size as yv (%d), actually its size is %d", neq, constraints.size());
  if (any(constraints)) {
    // check that only 0, ±1 and ±2 are present
    for (double it: constraints) {
      double ita=fabs(it);
      if (ita != 0. && ita != 1. && ita != 2.)
        stop("cvode: values in constraints must be 0, ±1 or ±2, instead %g found", it);
    }
    getmem(nv_constraints, N_VNewEmpty_Serial(neq));
    mem.add((void **) &nv_constraints, (funfree) N_VDestroy);
    NV_DATA_S(nv_constraints) = (realtype *) constraints.begin();
    check_retval(CVodeSetConstraints(cvode_mem, nv_constraints));
  }

  // which jacobian: dense or sparse
  if (nz > 0) {
    if (fjac.isNULL())
      stop("cvode: fjac must not be NULL if nz > 0 (%d)", nz);
    getmem(A, SUNSparseMatrix(neq, neq, nz, CSC_MAT));
    mem.add((void **) &A, (funfree) SUNMatDestroy);
    getmem(LS, SUNLinSol_RMUMPS(nv_y, A, rmumps_perm[0]));
    mem.add((void **) &LS, (funfree) SUNLinSolFree);
    check_retval(CVodeSetLinearSolver(cvode_mem, LS, A));
    check_retval(CVodeSetJacFn(cvode_mem, spjacwrap)); //jacsps)); //
  } else {
    // dense SUNMatrix for use in linear solves
    getmem(A, SUNDenseMatrix(neq, neq));
    mem.add((void **) &A, (funfree) SUNMatDestroy);
    // dense SUNLinearSolver object
    getmem(LS, SUNLinSol_Dense(nv_y, A));
    mem.add((void **) &LS, (funfree) SUNLinSolFree);
    check_retval(CVodeSetLinearSolver(cvode_mem, LS, A));
    if (!fjac.isNULL())
      check_retval(CVodeSetJacFn(cvode_mem, jacwrap)); // Jac)); //
  }
  // setup root treatment
  if (nroot > 0) {
    if (froot.isNULL())
      stop("r2sundials: froot cannot be NULL when nroot > 0 (%d)", nroot);
    if (froot.sexp_type() != CLOSXP && froot.sexp_type() != EXTPTRSXP)
      stop("r2sundials: froot must be an R function or an XPtr returned by RcppXPtrUtils::cppXPtr() or 'externalptr'");
    if (fevent.isNULL())
      stop("r2sundials: fevent cannot be NULL when nroot > 0 (%d)", nroot);
    rootsfound.resize(nroot);
    ynew.set_size(neq);
    mroots.set_size(nroot+1, 0);
    check_retval(CVodeRootInit(cvode_mem, nroot, rootwrap));
    if (fevent.sexp_type() == CLOSXP)
      rf_event=as<Function>(fevent);
  }
  // setup sensitivity
  bool sens_user=!fsens.isNULL() || !fsens1.isNULL(); // user provided fn_sens_rhs
  if (Ns > 0) {
    if (psens.size() < Ns && !sens_user) {
      stop("r2sundials: if Ns > 0 (%d) and no fsens neither fsens1 is provided then length(psesn) (%d) must be greater or equal to Ns", Ns, psens.size());
    }
    if(!fsens.isNULL() && !fsens1.isNULL()) {
      stop("r2sundials: fsens and fsens1 cannot be both provided. At least one of them or both must be NULL");
    }
    if (!sens_user) {
      // sanity check
      if (psens.size() > Ns && psens_list.size() == 0)
        stop("r2sundials: if length(psens) > Ns, psens_list must not be empty");
      // prepare param vectors
      if (psens_list.size() == 0) {
        plist=regspace<ivec>(0, Ns-1);
      } else {
        if (psens_list.size() != Ns)
          stop("r2sundials: if psens_list is not empty, its size (%s) must be equal to Ns (%d)", psens_list.size(), Ns);
        if (psens_list.size() > psens.size())
          stop("r2sundials: length(psens) must not be shorter than length(psens_list) (%d < %d)", psens.size(), psens_list.size());
        ivec u=unique(psens_list);
        if ((int) u.size() < psens_list.size())
          stop("r2sundials: psens_list has non unique entries");
        if (min(psens_list) < 1)
          stop("r2sundials: minimal value in psens_list is less than 1 (%d)", min(psens_list));
        if (max(psens_list) > psens.size())
          stop("r2sundials: maximal value in psens_list is greater than length(psens) (%d > %d)", max(psens_list), psens.size());
        plist=as<ivec>(psens_list);
        plist -= 1;
      }
    }
    // prepare yS at t0
    ivec d;
    if (sens_init.size() > 0 && !sens_init.hasAttribute("dim"))
      stop("r2sundials: if not empty, sens_init must be a dense (neq x Ns) matrix");
    if (sens_init.size() > 0)
      d=as<ivec>(sens_init.attr("dim"));
    if (sens_init.size() != 0 && d[0] != neq)
      stop("r2sundials: nrow(sens_init) (%d) must be equal to length(yv) (%d)", d[0], neq);
    if (sens_init.size() != 0 && d[1] != Ns)
      stop("r2sundials: ncol(sens_init) (%d) must be equal to Ns (%d)", d[1], Ns);
    if (sens_init.size() != 0)
      msens_init=mat(sens_init.begin(), d[0], d[1], false);
    else
      msens_init=zeros<mat>(neq, Ns);
//msens_init.print("msens_init");
    yS = N_VCloneVectorArray_Serial(Ns, nv_y);
    mem.add((void **) &yS, (funfree1<int>) N_VDestroyVectorArray, (int) Ns);
    for (int is=0; is < Ns; is++) {
      ySv[is]=vec(NV_DATA_S(yS[is]), neq, false); // vec proxies of yS "matrix"
      ySv[is]=msens_init.col(is);
    }
    // set sens problem
    if (!fsens.isNULL()) {
      check_retval(CVodeSensInit(cvode_mem, Ns, sens_method[0], sensrhswrap, yS));
      mem.add((void **) &cvode_mem, (funfree) CVodeSensFree);
    } else if (!fsens1.isNULL()) {
      check_retval(CVodeSensInit1(cvode_mem, Ns, sens_method[0], sensrhs1wrap, yS));
      mem.add((void **) &cvode_mem, (funfree) CVodeSensFree);
    } else {
      check_retval(CVodeSensInit1(cvode_mem, Ns, sens_method[0], NULL, yS));
      mem.add((void **) &cvode_mem, (funfree) CVodeSensFree);
    }
    check_retval(CVodeSetSensParams(cvode_mem, psens.begin(), psens_bar.size() == 0 ? NULL : psens_bar.begin(), plist.begin()));
    check_retval(CVodeSensEEtolerances(cvode_mem));
    check_retval(CVodeSetSensErrCon(cvode_mem, errconS));
  }

  // time course
  // copy t0
  res.col(0)=yvec;
  if (Ns > 0) {
    for (int is=0; is < Ns; is++)
      asens.slice(is).col(0)=ySv[is];
  }
  // main loop
  for (int iout=1, insroot=0; iout < nti; iout++) {
    retval=CVode(cvode_mem, times[iout], nv_y, &t, CV_NORMAL);
    if (retval == CV_ROOT_RETURN) {
      // root found proceed it
      check_retval(CVodeGetRootInfo(cvode_mem, rootsfound.begin()));
      if (fevent.sexp_type() == CLOSXP) {
        // call R fevent
//Rcout << "rf_event=" << rf_event << "\n";
//print(fevent);
        if (Ns > 0) {
          // set from ySv
          for (int is=0; is < Ns; is++)
            msens.col(is)=ySv[is];
        }
        List lres=rf_event(t, yvec, Ns, msens, rootsfound, param, psens);
//Rcout << "lres=\n";
//print(lres);
        ynew=as<vec>(lres["ynew"]);
        retevent=as<int>(lres["flag"]);
        if (Ns > 0) {
          if (lres["sens_init"] == R_NilValue)
            msens.fill(0.);
          else
            msens=as<mat>(lres["sens_init"]);
          // reinit yS
          for (int is=0; is < Ns; is++)
            ySv[is]=msens.col(is);
        }
      } else {
        // call XPtr fevent
        if (user_event_fn == NULL)
          user_event_fn=*(as<XPtr< rsunEventFn >>(fevent));
        retevent=user_event_fn(t, yvec, ynew, Ns, ySv, rootsfound, param, psens);
      }
      // treat retevent
      if (retevent == R2SUNDIALS_EVENT_IGNORE) {
        ; // do nothing
      } else {
        mroots.insert_cols(mroots.n_cols, 1, false);
        mroots.col(mroots.n_cols-1)[0]=t;
        mroots.col(mroots.n_cols-1).subvec(1, nroot)=conv_to<vec>::from(rootsfound);
        // insert this time point and current state
        if (t < times[iout]) {
          ti.insert_rows(iout+insroot, 1);
          ti[iout+insroot]=t;
          res.insert_cols(iout+insroot, yvec);
          if (Ns > 0) {
            CVodeGetSens(cvode_mem, &t, yS);
            asens.insert_cols(iout+insroot, 1, false);
            for (int is=0; is < Ns; is++)
              asens.slice(is).col(iout+insroot)=ySv[is];
          }
          insroot++;
        }
        if (any(ynew != yvec)) {
          // insert new state too
          ti.insert_rows(iout+insroot, 1);
          ti[iout+insroot]=t;
          res.insert_cols(iout+insroot, ynew);
          if (Ns > 0) {
            // reinit yS to ySv
            asens.insert_cols(iout+insroot, 1, false);
            for (int is=0; is < Ns; is++) {
              asens.slice(is).col(iout+insroot)=ySv[is];
            }
          }
          insroot++;
          // reinit cvode
          yvec=ynew;
          check_retval(CVodeReInit(cvode_mem, t, nv_y));
          if (Ns > 0) {
            // reinit sensitivities, yS is already set by fevent() via ySv
            CVodeSensFree(cvode_mem);
            if (!fsens.isNULL()) {
              check_retval(CVodeSensInit(cvode_mem, Ns, sens_method[0], sensrhswrap, yS));
            } else if (!fsens1.isNULL()) {
              check_retval(CVodeSensInit1(cvode_mem, Ns, sens_method[0], sensrhs1wrap, yS));
            } else {
              check_retval(CVodeSensInit1(cvode_mem, Ns, sens_method[0], NULL, yS));
            }
          }
        }
        if (retevent == R2SUNDIALS_EVENT_STOP) {
          // free extra cols in res and ti
          res.resize(neq, iout+insroot);
          ti.resize(iout+insroot);
          break;
        } else if (retevent != R2SUNDIALS_EVENT_HOLD) {
          stop("cvode: fevent() returned unknown flag %d", retevent);
        }
      }
      if (t < times[iout])
        iout--; // stay at the same time point for the next iteration
    } else {
      // simple storage of the current solution
      res.col(iout+insroot)=yvec;
      if (Ns > 0) {
        check_retval(CVodeGetSens(cvode_mem, &t, yS));
        for (int is=0; is < Ns; is++)
          asens.slice(is).col(iout+insroot)=ySv[is];
      }
    }
    checkUserInterrupt();
  }
  // final stats as attribute to output
  std::vector<long int> stats(15, 0L);
  // stats
  i=0;
  check_retval(CVodeGetNumSteps(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumRhsEvals(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumLinSolvSetups(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumErrTestFails(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumNonlinSolvIters(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumNonlinSolvConvFails(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumJacEvals(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumLinRhsEvals(cvode_mem, (long int*) &stats[i++]));
  check_retval(CVodeGetNumGEvals(cvode_mem, (long int*) &stats[i++]));

  if (Ns > 0) {
    check_retval(CVodeGetSensNumRhsEvals(cvode_mem, (long int*) &stats[i++]));
    check_retval(CVodeGetNumRhsEvalsSens(cvode_mem, (long int*) &stats[i++]));
    check_retval(CVodeGetSensNumLinSolvSetups(cvode_mem, (long int*) &stats[i++]));
    if (errconS) {
      check_retval(CVodeGetSensNumErrTestFails(cvode_mem, (long int*) &stats[i]));
    }
    i++;
    if (sens_method[0] == CV_STAGGERED) {
      check_retval(CVodeGetSensNumNonlinSolvIters(cvode_mem, (long int*) &stats[i++]));
      check_retval(CVodeGetSensNumNonlinSolvConvFails(cvode_mem, (long int*) &stats[i++]));
    }
    //CVodeSensFree(cvode_mem);
    //N_VDestroyVectorArray(yS, Ns);
  }
  // convert stats to int
  IntegerVector statsi(15);
  std::copy(stats.begin(), stats.end(), statsi.begin());
  statsi.attr("names")=CharacterVector::create("NumSteps", "NumRhsEvals", "NumLinSolvSetups", "NumErrTestFails", "NumNonlinSolvIters", "NumNonlinSolvConvFails", "NumJacEvals", "NumLinRhsEvals", "NumGEvals",
    "SensNumRhsEvals", "NumRhsEvalsSens", "SensNumLinSolvSetups", "SensNumErrTestFails", "SensNumNonlinSolvIters", "SensNumNonlinSolvConvFails"
  );
  //N_VDestroy(nv_y);
  //CVodeFree(&cvode_mem);
  //SUNLinSolFree(LS);
  //SUNMatDestroy(A);
  
  // set colnames in resout
  StringVector colnm(neq);
  if (yv.hasAttribute("names"))
    colnm=yv.attr("names");
  else
    for (int i=0; i < neq; i++)
      colnm[i]="V"+std::to_string(i+1);
  NumericMatrix resout=wrap(res);
  resout.attr("dimnames")=List::create(colnm, CharacterVector(ti.begin(), ti.end()));
  resout.attr("times")=ti;
  resout.attr("stats")=statsi;
  if (nroot > 0)
    resout.attr("roots")=wrap(mroots);
  if (Ns > 0)
    resout.attr("sens")=wrap(asens);
  //t2 = clock();
  //(void)printf("Temps consomme (s) : %lf \n", (double)(t2-t1)/(double)clk_tck);
  mem.freeall();
  return(resout);
}

int rhswrap(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  UserData *udata=(UserData *) user_data;
  List lp=udata->lp;
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  vec yv(NV_DATA_S(y), neq, false), ydotv(NV_DATA_S(ydot), neq, false);
//yv.print("yv in rhswrap");
//Rcout << "NV_DATA_S(y)=" << NV_DATA_S(y) << "\n";
  if (TYPEOF(lp["frhs"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["frhs"]);
    ydotv=as<vec>(user_fr(wrap(t), wrap(yv), param, udata->psens));
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunRhsFn > xuser_fn=as<XPtr< rsunRhsFn >>(lp["frhs"]);
    rsunRhsFn user_fn=*xuser_fn;
    return(user_fn(t, yv, ydotv, param, udata->psens));
  }
}
int jacwrap(realtype t, N_Vector y, N_Vector ydot, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  UserData *udata=(UserData *) user_data;
  List lp=udata->lp;
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  vec yv(NV_DATA_S(y), neq, false),
    ydotv(NV_DATA_S(ydot), neq, false),
    tmp1v(NV_DATA_S(tmp1), NV_LENGTH_S(tmp1), false),
    tmp2v(NV_DATA_S(tmp2), NV_LENGTH_S(tmp2), false),
    tmp3v(NV_DATA_S(tmp3), NV_LENGTH_S(tmp3), false);
  mat ja(SM_DATA_D(J), neq, neq, false);
  if (TYPEOF(lp["fjac"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["fjac"]);
    ja.submat(0, neq-1, 0, neq-1)=as<mat>(user_fr(wrap(t), wrap(yv), wrap(ydotv), param, udata->psens));
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunJacFn > xuser_fn=as<XPtr< rsunJacFn >>(lp["fjac"]);
    rsunJacFn user_fn=*xuser_fn;
    return(user_fn(t, yv, ydotv, ja, param, udata->psens, tmp1v, tmp2v, tmp3v));
  }
}
int spjacwrap(realtype t, N_Vector y, N_Vector ydot, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  UserData *udata=(UserData *) user_data;
  List lp=udata->lp;
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y), nz=SM_NNZ_S(J);
  vec yv(NV_DATA_S(y), neq, false), ydotv(NV_DATA_S(ydot), neq, false),
    tmp1v(NV_DATA_S(tmp1), NV_LENGTH_S(tmp1), false),
    tmp2v(NV_DATA_S(tmp2), NV_LENGTH_S(tmp2), false),
    tmp3v(NV_DATA_S(tmp3), NV_LENGTH_S(tmp3), false);
#if defined(SUNDIALS_INT32_T)
  uvec ir((unsigned int *) SM_INDEXVALS_S(J), nz, false);
  uvec p((unsigned int *) SM_INDEXPTRS_S(J), neq+1, false);
#else
  Col<sunindextype> irsun(SM_INDEXVALS_S(J), nz, false), psun(SM_INDEXPTRS_S(J), neq+1, false);
  uvec ir=conv_to<uvec>::from(irsun);
  uvec p=conv_to<uvec>::from(psun);
#endif
  vec v((double *) SM_DATA_S(J), nz, false);
  if (TYPEOF(lp["fjac"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["fjac"]);
    List res=user_fr(wrap(t), wrap(yv), wrap(ydotv), param, udata->psens);
    uvec rir=res["i"];
    uvec rp=res["p"];
    vec rv=res["v"];
#if defined(SUNDIALS_INT32_T)
    ir=rir-1;
    p=rp-1;
#else
    irsun=conv_to<Col<sunindextype>>::from(rir-1);
    psun=conv_to<Col<sunindextype>>::from(rp-1);
#endif
    v=rv;
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunSpJacFn > xuser_fn=as<XPtr< rsunSpJacFn >>(lp["fjac"]);
    rsunSpJacFn user_fn=*xuser_fn;
    int res=user_fn(t, yv, ydotv, ir, p, v, neq, nz, param, udata->psens, tmp1v, tmp2v, tmp3v);
#if !defined(SUNDIALS_INT32_T)
    irsun=conv_to<Col<sunindextype>>::from(ir);
    psun=conv_to<Col<sunindextype>>::from(p);
#endif
    return(res);
  }
}
int rootwrap(realtype t, N_Vector y, realtype *rootout, void *user_data) {
  UserData *udata=(UserData *) user_data;
  List lp=udata->lp;
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  int nroot=lp["nroot"];
  vec yv(NV_DATA_S(y), neq, false);
  vec vroot((double *) rootout, nroot, false);
  if (TYPEOF(lp["froot"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["froot"]);
    vroot=as<vec>(user_fr(wrap(t), wrap(yv), param, udata->psens));
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunRootFn > xuser_fn=as<XPtr< rsunRootFn >>(lp["froot"]);
    rsunRootFn user_fn=*xuser_fn;
    return(user_fn(t, yv, vroot, param, udata->psens));
  }
}
int sensrhswrap(int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
  UserData *udata=(UserData *) user_data;
  List lp=udata->lp;
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  vec yv(NV_DATA_S(y), neq, false);
  vec ydotv(NV_DATA_S(ydot), neq, false);
  vec tmp1v(NV_DATA_S(tmp1), NV_LENGTH_S(tmp1), false);
  vec tmp2v(NV_DATA_S(tmp2), NV_LENGTH_S(tmp2), false);
  std::vector<vec> ySv(Ns), ySdotv(Ns);
  // init stdvec
  for (int i=0; i < Ns; i++) {
    ySv[i]=vec(NV_DATA_S(yS[i]), neq, false);
    ySdotv[i]=vec(NV_DATA_S(ySdot[i]), neq, false);
  }
  if (TYPEOF(lp["fsens"]) == CLOSXP) {
    // plain R function
    mat ySm(neq, Ns);
    for (int is=0; is < Ns; is++)
      ySm.col(is)=ySv[is];
    Function user_fr=as<Function>(lp["fsens"]);
    DataFrame res=user_fr(wrap(Ns), wrap(t), wrap(yv), wrap(ydotv), wrap(ySm), param, udata->psens);
    for (int i=0; i < Ns; i++)
      ySdotv[i]=as<vec>(res[i]);
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunSensFn > xuser_fn=as<XPtr< rsunSensFn >>(lp["fsens"]);
    rsunSensFn user_fn=*xuser_fn;
    return(user_fn(Ns, t, yv, ydotv, ySv, ySdotv, param, udata->psens, tmp1v, tmp2v));
  }
}
int sensrhs1wrap(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
  UserData *udata=(UserData *) user_data;
  List lp=udata->lp;
  RObject param=lp["param"];
  int neq=NV_LENGTH_S(y);
  vec yv(NV_DATA_S(y), neq, false);
  vec ydotv(NV_DATA_S(ydot), neq, false);
  vec ySv(NV_DATA_S(yS), neq, false);
  vec ySdotv(NV_DATA_S(ySdot), neq, false);
  vec tmp1v(NV_DATA_S(tmp1), NV_LENGTH_S(tmp1), false);
  vec tmp2v(NV_DATA_S(tmp2), NV_LENGTH_S(tmp2), false);
  
  if (TYPEOF(lp["fsens1"]) == CLOSXP) {
    // plain R function
    Function user_fr=as<Function>(lp["fsens1"]);
    ySv=as<vec>(user_fr(wrap(Ns), wrap(t), wrap(yv), wrap(iS), wrap(ydotv), wrap(ySv), param, udata->psens));
    return(CV_SUCCESS);
  } else {
    // get user function pointer
    XPtr< rsunSens1Fn > xuser_fn=as<XPtr< rsunSens1Fn >>(lp["fsens1"]);
    rsunSens1Fn user_fn=*xuser_fn;
    return(user_fn(Ns, t, yv, ydotv, iS, ySv, ySdotv, param, udata->psens, tmp1v, tmp2v));
  }
}

// helpers
void rsunerr(int error_code, const char *module, const char *function, char *msg, void *eh_data) {
  throw Rcpp::exception(tfm::format("%s: %s", function, msg).c_str(), false);
  //stop("%s: %s", function, msg);
}

#define cnst_pair(c) {#c, c}

//' @title Get Predefined Constant
//' @description This function returns numerical value of predefined named constants. Constants defined for this package are: \itemize{
//' \item CV_ADAMS
//' \item CV_BDF
//' \item CV_SIMULTANEOUS
//' \item CV_STAGGERED
//' \item CV_STAGGERED1
//' \item CV_SUCCESS
//' \item R2SUNDIALS_EVENT_HOLD
//' \item R2SUNDIALS_EVENT_IGNORE
//' \item R2SUNDIALS_EVENT_STOP
//' }
//' All these constants are exported from the package and available for use in callback functions written in R.
//' @param s character scalar, name of the constant whose value is to get.
//' @return integer scalar, the value of the constant named by \code{s}
//' @examples
//' get_cnst("CV_SUCCESS")
//' # [1] 0
//' @export
// [[Rcpp::export]]
int get_cnst(std::string s) {
  static std::map<std::string, int> dict={
    cnst_pair(CV_SUCCESS),
    
    cnst_pair(CV_BDF),
    cnst_pair(CV_ADAMS),
    
    cnst_pair(R2SUNDIALS_EVENT_IGNORE),
    cnst_pair(R2SUNDIALS_EVENT_HOLD),
    cnst_pair(R2SUNDIALS_EVENT_STOP),
    cnst_pair(CV_SIMULTANEOUS),
    cnst_pair(CV_STAGGERED),
    cnst_pair(CV_STAGGERED1)
  };
  if (dict.count(s))
    return(dict[s]);
  else
    stop("get_cnst: constant '%s' is not in dictionary", s);
  return(NA_INTEGER);
}
