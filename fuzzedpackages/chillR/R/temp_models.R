#' Calculation of cumulative temperature metric according to a user-defined
#' stepwise weight function
#'
#' This function calculates heat for temperate trees according to a stepwise
#' model provided by the user.
#'
#' Temperature-based metric calculated according to the user-defined model.
#'
#' @param HourTemp Vector of hourly temperatures.
#' @param df data.frame with three columns: lower, upper and weight. lower
#' should contain the lower boundary of a chilling weight interval and upper
#' should contain the upper boundary. weight indicates the weighting to be
#' applied to the respective temperature interval.
#' @param summ Boolean parameter indicating whether calculated metrics should
#' be provided as cumulative values over the entire record (TRUE) or as the
#' actual accumulation for each hour (FALSE).
#' @return Vector of length length(HourTemp) containing the cumulative
#' temperature metric over the entire duration of HourTemp.
#' @author Eike Luedeling
#' @keywords chill and heat calculation
#' @examples
#'
#'
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#'
#' stack<-stack_hourly_temps(weather,latitude=50.4)
#'
#' df=data.frame(
#'   lower=c(-1000,1,2,3,4,5,6),
#'   upper=c(1,2,3,4,5,6,1000),
#'   weight=c(0,1,2,3,2,1,0))
#'
#' custom<-function(x) step_model(x,df)
#'
#' custom(stack$Temp)
#'
#' models<-list(Chilling_Hours=Chilling_Hours,Utah_Chill_Units=Utah_Model,
#' Chill_Portions=Dynamic_Model,GDH=GDH,custom=custom)
#'
#' tempResponse(stack,Start_JDay = 305,End_JDay = 60,models)
#'
#' @export step_model
step_model<-function(HourTemp,
                     df=data.frame(
                       lower=c(-1000,1.4,2.4,9.1,12.4,15.9,18),
                       upper=c(1.4,2.4,9.1,12.4,15.9,18,1000),
                       weight=c(0,0.5,1,0.5,0,-0.5,-1)),summ=TRUE)
{lower<-df$lower;upper<-df$upper;weight<-df$weight
if (summ==TRUE) return(cumsum(sapply(HourTemp,function(x) weight[which(x>lower&x<=upper)]))) else
  return(sapply(HourTemp,function(x) weight[which(x>lower&x<=upper)]))
}

#df=data.frame(
#  lower=c(-1000,1,2,3,4,5,6),
#  upper=c(1,2,3,4,5,6,1000),
#  weight=c(0,1,2,3,2,1,0),
#  lower_include_equal=c(F,F,T,T,T,T,T),
#  upper_include_equal=c(F,F,T,T,T,T,T))




#' Calculation of cumulative chill according to the Utah Model
#'
#' This function calculates winter chill for temperate trees according to the
#' Utah Model.
#'
#' Units of the Utah Model are calculated as suggested by Richardson et al.
#' (1974) (different weights for different temperature ranges, and negation of
#' chilling by warm temperatures).
#'
#' @param HourTemp Vector of hourly temperatures.
#' @param summ Boolean parameter indicating whether calculated metrics should
#' be provided as cumulative values over the entire record (TRUE) or as the
#' actual accumulation for each hour (FALSE).
#' @return Vector of length length(HourTemp) containing the cumulative Utah
#' Chill Units over the entire duration of HourTemp.
#' @note After doing extensive model comparisons, and reviewing a lot of
#' relevant literature, I do not recommend using the Utah Model, especially in
#' warm climates! The Dynamic Model (Chill Portions), though far from perfect,
#' seems much more reliable.
#' @author Eike Luedeling
#' @references Utah Model reference:
#'
#' Richardson EA, Seeley SD, Walker DR (1974) A model for estimating the
#' completion of rest for Redhaven and Elberta peach trees. HortScience 9(4),
#' 331-332
#' @keywords chill and heat calculation
#' @examples
#'
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#'
#' stack<-stack_hourly_temps(weather,latitude=50.4)
#'
#' Utah_Model(stack$hourtemps$Temp)
#'
#' @export Utah_Model
Utah_Model<-function(HourTemp,summ=TRUE)
  return(step_model(HourTemp,df=data.frame(lower=c(-1000,1.4,2.4,9.1,12.4,15.9,18),upper=c(1.4,2.4,9.1,12.4,15.9,18,1000),weight=c(0,0.5,1,0.5,0,-0.5,-1)),summ=summ))



#Utah_Model<-function(HourTemp)
#  {#Utah Model
#  Utah_range_0.5<-which(HourTemp<=2.4&HourTemp>1.4|
#                          HourTemp<=12.4&HourTemp>9.1)
#  Utah_range_1.0<-which(HourTemp<=9.1&HourTemp>2.4)
#  Utah_range_min0.5<-which(HourTemp<=18.0&HourTemp>15.9)
#  Utah_range_min1.0<-which(HourTemp>18.0)
#  Utah_weights<-rep(0,length(HourTemp))
#  Utah_weights[Utah_range_0.5]<-0.5
#  Utah_weights[Utah_range_1.0]<-1
#  Utah_weights[Utah_range_min0.5]<-(-0.5)
#  Utah_weights[Utah_range_min1.0]<-(-1)
#  return(cumsum(Utah_weights))
#}



#' Calculation of cumulative chill according to the Chilling Hours Model
#'
#' This function calculates winter chill for temperate trees according to the
#' Chilling Hours Model.
#'
#' Chilling Hours are calculated as suggested by Bennett (1949) (all hours with
#' temperatures between 0 and 7.2 degrees C are considered as one Chilling
#' Hour.
#'
#' @param HourTemp Vector of hourly temperatures.
#' @param summ Boolean parameter indicating whether calculated metrics should
#' be provided as cumulative values over the entire record (TRUE) or as the
#' actual accumulation for each hour (FALSE).
#' @return Vector of length length(HourTemp) containing the cumulative Chilling
#' Hours over the entire duration of HourTemp.
#' @note After doing extensive model comparisons, and reviewing a lot of
#' relevant literature, I do not recommend using the Chilling Hours, especially
#' in warm climates! The Dynamic Model (Chill Portions), though far from
#' perfect, seems much more reliable.
#' @author Eike Luedeling
#' @references Chilling Hours references:
#'
#' Weinberger JH (1950) Chilling requirements of peach varieties. Proc Am Soc
#' Hortic Sci 56, 122-128
#'
#' Bennett JP (1949) Temperature and bud rest period. Calif Agric 3 (11), 9+12
#' @keywords chill and heat calculation
#' @examples
#'
#'
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#'
#' hourtemps<-stack_hourly_temps(weather,latitude=50.4)
#'
#' Chilling_Hours(hourtemps$hourtemps$Temp)
#'
#' @export Chilling_Hours
Chilling_Hours<-function(HourTemp,summ=TRUE)
{
  CH_range<-which(HourTemp<=7.2&HourTemp>=0)
  CH_weights<-rep(0,length(HourTemp))
  CH_weights[CH_range]<-1
  if(summ==TRUE) return(cumsum(CH_weights)) else
    return(CH_weights)
}



#' @title Dynamic_Model
#'
#' @description
#' Calculation of cumulative chill according to the Dynamic Model
#'
#' This function calculates winter chill for temperate trees according to the
#' Dynamic Model.
#'
#' Chill Portions are calculated as suggested by Erez et al. (1990).
#'
#' @param HourTemp Vector of hourly temperatures in degree Celsius.
#' @param summ Boolean parameter indicating whether calculated metrics should
#' be provided as cumulative values over the entire record (TRUE) or as the
#' actual accumulation for each hour (FALSE).
#' @param A0 numeric. Parameter \eqn{A_0}{A0} of the dynamic model
#' @param A1 numeric. Parameter \eqn{A_1}{A1} of the dynamic model
#' @param E0 numeric. Parameter \eqn{E_0}{E0} of the dynamic model
#' @param E1 numeric. Parameter \eqn{E_1}{E1} of the dynamic model
#' @param slope numeric. Slope parameter for sigmoidal function
#' @param Tf numeric. Transition temperature (in degree Kelvin) for the
#' sigmoidal function.
#' @return Vector of length length(HourTemp) containing the cumulative Chill
#' Portions over the entire duration of HourTemp.
#' @author Eike Luedeling
#' @references Dynamic Model references:
#'
#' Erez A, Fishman S, Linsley-Noakes GC, Allan P (1990) The dynamic model for
#' rest completion in peach buds. Acta Hortic 276, 165-174
#'
#' Fishman S, Erez A, Couvillon GA (1987a) The temperature dependence of
#' dormancy breaking in plants - computer simulation of processes studied under
#' controlled temperatures. J Theor Biol 126(3), 309-321
#'
#' Fishman S, Erez A, Couvillon GA (1987b) The temperature dependence of
#' dormancy breaking in plants - mathematical analysis of a two-step model
#' involving a cooperative transition. J Theor Biol 124(4), 473-483
#' @keywords chill and heat calculation
#' @examples
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#'
#' hourtemps<-stack_hourly_temps(weather,latitude=50.4)
#'
#' res <- Dynamic_Model(hourtemps$hourtemps$Temp)
#'
#' @export Dynamic_Model
Dynamic_Model <- function(HourTemp,
                          summ = TRUE,
                          E0 = 4153.5,
                          E1 = 12888.8,
                          A0 = 139500,
                          A1 = 2567000000000000000,
                          slope = 1.6,
                          Tf = 277)
{
  if(missing(HourTemp)) {
    stop("HourTemp must be present! Aborting.")
  }
  ## compute temperatures in Kelvin
  TK <- HourTemp + 273
  ## pre-compute some constants
  aa <- A0/A1
  ee <- E1-E0
  sr <- exp(slope*Tf*(TK-Tf)/TK)
  xi <- sr/(1+sr)
  xs <- aa*exp(ee/TK)
  eak1 <- exp(-A1*exp(-E1/TK))
  
  ## initialise PDBF with zero
  x=0
  for (l in c(2:length(HourTemp)))  {
    S <- x[l-1]
    if(x[l-1] >= 1) {
      S <- S*(1-xi[l-2])
    }
    x[l] <- xs[l-1]-(xs[l-1]-S)*eak1[l-1]
  }
  ## chill portions
  delta <- rep(0,length(HourTemp))
  ii <- which(x >= 1)
  delta[ii] <- x[ii]*xi[ii-1]
  if (summ) return(cumsum(delta))
  return(delta)
}

## x_s(t_i) according to the first reference
xsfct  <- function(A0, A1, E0, E1, temp) {
  return ( A0/A1 * exp(-(E0 - E1)/temp) )
}
## k_1(t_i) according to the first reference
k1fct <- function(A1, E1, temp) {
  return( A1*exp(-E1/temp) )
}

#' @title DynModel_driver
#'
#' @description
#' Calculation of cumulative chill according to the Dynamic Model
#'
#' This function calculates winter chill for temperate trees according to the
#' Dynamic Model.
#'
#' Chill Portions are calculated as suggested by Erez et al. (1990).
#'
#' @details This function gives idential results as \link{Dynamic_Model} for hourly
#' temperature data, returns more details but is also a bit slower in the R code
#' version
#' @param temp Vector of temperatures.
#' @param times numeric vector. Optional times at which the temperatures where measured,
#'   if not given, hourly temperatures will be assumed
#' @param A0 numeric. Parameter \eqn{A_0}{A0} of the dynamic model
#' @param A1 numeric. Parameter \eqn{A_1}{A1} of the dynamic model
#' @param E0 numeric. Parameter \eqn{E_0}{E0} of the dynamic model
#' @param E1 numeric. Parameter \eqn{E_1}{E1} of the dynamic model
#' @param slope numeric. Slope parameter for sigmoidal function
#' @param Tf numeric. Transition temperature (in degree Kelvin) for the
#' sigmoidal function
#' @param deg_celsius boolean. whether or not the temperature vector
#' and the model temperature parameters are
#' in degree Celsius (Kelvin otherwise)
#' @return List containint four vectors of length(temp) with elements
#' \code{x} is the PDBF, \code{y} the accumulated chill, \code{delta} the
#' chill portions and \code{xs}, which is
#' \eqn{x_s=A_0/A_1\exp(-(E_0-E_1)/T)}{xs=A0/A1\exp(-(E0-E1)/T}
#' Portions over the entire duration of HourTemp.
#' @author Carsten Urbach <urbach@hiskp.uni-bonn.de>
#' @references Dynamic Model references:
#'
#' Erez A, Fishman S, Linsley-Noakes GC, Allan P (1990) The dynamic model for
#' rest completion in peach buds. Acta Hortic 276, 165-174
#'
#' Fishman S, Erez A, Couvillon GA (1987a) The temperature dependence of
#' dormancy breaking in plants - computer simulation of processes studied under
#' controlled temperatures. J Theor Biol 126(3), 309-321
#'
#' Fishman S, Erez A, Couvillon GA (1987b) The temperature dependence of
#' dormancy breaking in plants - mathematical analysis of a two-step model
#' involving a cooperative transition. J Theor Biol 124(4), 473-483
#' @keywords chill and heat calculation
#' @examples
#'
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#' hourtemps<-stack_hourly_temps(weather,latitude=50.4)
#' res2 <- DynModel_driver(temp=hourtemps$hourtemps$Temp)
#'
#' @export DynModel_driver
DynModel_driver  <- function(temp,
                             times,
                             A0 = 139500,
                             A1 = 2567000000000000000,
                             E0 = 4153.5,
                             E1 = 12888.8,
                             slope = 1.6,
                             Tf = 4, 
                             deg_celsius=TRUE) {
  if(missing(times)) {
    times <- c(1:length(temp))
  }
  stopifnot(length(temp) == length(times))
  if(deg_celsius) {
    temp <- temp + 273
    Tf <- Tf + 273
  }
  ## pre-compute x_s and k_1 for all temperatures
  ## x_s(t_i) according to Eq. 6
  xs <- A0/A1 * exp(-(E0 - E1)/temp)
  ## k_1(t_i) according to Eq. 6 multiplied by time differences and
  ## exponentiated
  ek1 <- exp(-A1*exp(-E1/temp)*c(diff(x=times, lag=1), 0))
  
  ## initialise with zero
  x <- 0.
  y  <- 0.
  ## now we compute x_i and y_i recursively for
  ## all discrete times 
  for(i in c(1:(length(times)-1))) {
    ## x at t_i+1 according to Eq. 9
    x[i+1] <- xs[i] - (xs[i] - x[i])*ek1[i]
    y[i+1] <- 0
    ## did x exceed 1?
    if(x[i+1] >= 1) {
      ## now use the modified procedure
      xtmp <- slope*Tf*(temp[i]-Tf)/temp[i]
      ## avoid instabilities in sigmoidal ratio
      delta <- x[i+1]
      if(xtmp < 17) {
        sr <- exp(xtmp)
        delta <- delta*sr/(1+sr)
      }
      y[i+1] <- delta
      x[i+1] <- x[i+1] - delta
    }
  }
  return(invisible(list(x=x, y=cumsum(y), delta=y, xs=xs)))
}



#' Calculation of cumulative heat according to the Growing Degree Hours Model
#'
#' This function calculates heat for temperate trees according to the Growing
#' Degree Hours Model.
#'
#' Growing Degree Hours are calculated as suggested by Anderson et al. (1986).
#'
#' @param HourTemp Vector of hourly temperatures.
#' @param summ Boolean parameter indicating whether calculated metrics should
#' be provided as cumulative values over the entire record (TRUE) or as the
#' actual accumulation for each hour (FALSE).
#' @return Vector of length length(HourTemp) containing the cumulative Growing
#' Degree Hours over the entire duration of HourTemp.
#' @author Eike Luedeling
#' @references Growing Degree Hours reference:
#'
#' Anderson JL, Richardson EA, Kesner CD (1986) Validation of chill unit and
#' flower bud phenology models for 'Montmorency' sour cherry. Acta Hortic 184,
#' 71-78
#' @keywords chill and heat calculation
#' @examples
#'
#'
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#'
#' hourtemps<-stack_hourly_temps(weather,latitude=50.4)
#'
#' GDH(hourtemps$hourtemps$Temp)
#'
#' @export GDH
GDH<-function(HourTemp,summ=TRUE)
{Stress<-1
Tb<-4
Tu<-25
Tc<-36

GDH_weight<-rep(0,length(HourTemp))
GDH_weight[which(HourTemp>=Tb&HourTemp<=Tu)]<-Stress*(Tu-Tb)/2*
  (1+cos(pi+pi*(HourTemp[which(HourTemp>=Tb&HourTemp<=Tu)]-Tb)/(Tu-Tb)))
GDH_weight[which(HourTemp>Tu&HourTemp<=Tc)]<-Stress*(Tu-Tb)*
  (1+cos(pi/2+pi/2*(HourTemp[which(HourTemp>Tu&HourTemp<=Tc)]-Tu)/(Tc-Tu)))
if (summ) return(cumsum(GDH_weight)) else
  return(GDH_weight)
}
#' Calculation of cumulative heat according to the Growing Degree Hours Model
#' (alternative function name)
#'
#' This function calculates heat for temperate trees according to the Growing
#' Degree Hours Model.
#'
#' Growing Degree Hours are calculated as suggested by Anderson et al. (1986).
#'
#' @param HourTemp Vector of hourly temperatures.
#' @param summ Boolean parameter indicating whether calculated metrics should
#' be provided as cumulative values over the entire record (TRUE) or as the
#' actual accumulation for each hour (FALSE).
#' @return Vector of length length(HourTemp) containing the cumulative Growing
#' Degree Hours over the entire duration of HourTemp.
#' @author Eike Luedeling
#' @references Growing Degree Hours reference:
#'
#' Anderson JL, Richardson EA, Kesner CD (1986) Validation of chill unit and
#' flower bud phenology models for 'Montmorency' sour cherry. Acta Hortic 184,
#' 71-78
#' @keywords chill and heat calculation
#' @examples
#'
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#'
#' hourtemps<-stack_hourly_temps(weather,latitude=50.4)
#'
#' GDH_model(hourtemps$hourtemps$Temp)
#'
#' @export GDH_model
GDH_model<-function(HourTemp,summ=TRUE)
{Stress<-1
Tb<-4
Tu<-25
Tc<-36

GDH_weight<-rep(0,length(HourTemp))
GDH_weight[which(HourTemp>=Tb&HourTemp<=Tu)]<-Stress*(Tu-Tb)/2*
  (1+cos(pi+pi*(HourTemp[which(HourTemp>=Tb&HourTemp<=Tu)]-Tb)/(Tu-Tb)))
GDH_weight[which(HourTemp>Tu&HourTemp<=Tc)]<-Stress*(Tu-Tb)*
  (1+cos(pi/2+pi/2*(HourTemp[which(HourTemp>Tu&HourTemp<=Tc)]-Tu)/(Tc-Tu)))
if (summ) return(cumsum(GDH_weight)) else
  return(GDH_weight)
}
#'
#'  Calculation of cumulative heat according to the Growing Degree Day Model
#'
#' This function calculates heat for temperate trees according to the Growing
#' Degree Day Model. Note that the calculuation differs slightly from the original,
#' in which it is based on daily temperature extremes only. This equation here works
#' with hourly temperatures. The normal GDD equation is GDD=(Tmax-Tmin)/2-Tbase, with
#' Tmax=30 for Tmax>30, and Tmin=10 for Tmin<10. Tbase is a species-specific base
#' temperature.
#' The first part of the equation is the arithmetic mean of daily temperature extremes.
#' In the present equation, this is replaced by Thourly/24 for each hourly temperature
#' value. If chillR was using a triangular daily temperature curve, the result would
#' be the same for both equations. Since chillR uses a sine function for daytime
#' warming and a logarithmic decay function for nighttime cooling, however, there
#' will be a slight deviation. This could be handled by defining a function the runs
#' with daily weather data. chillR doesn't currently have this capability, since
#' its primary focus is on metrics that require hourly data.
#'
#' Growing Degree Hours are calculated as suggested by Anderson et al. (1986).
#'
#' @param HourTemp Vector of hourly temperatures.
#' @param summ Boolean parameter indicating whether calculated metrics should
#' be provided as cumulative values over the entire record (TRUE) or as the
#' actual accumulation for each hour (FALSE).
#' @param Tbase Base temperature, above which Growing Degrees accrue.
#' @return Vector of length length(HourTemp) containing the cumulative Growing
#' Degree Days over the entire duration of HourTemp.
#' @author Eike Luedeling
#' @references Growing Degree Days reference:
#'
#' http://agron-www.agron.iastate.edu/Courses/agron212/Calculations/GDD.htm
#'
#' @keywords chill and heat calculation
#' @examples
#'
#'
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#'
#' hourtemps<-stack_hourly_temps(weather,latitude=50.4)
#'
#' GDD(hourtemps$hourtemps$Temp)
#'
#' @export GDD
GDD<-function(HourTemp,summ=TRUE,Tbase=5)
{
  Tlow<-10
  Tup<-30
  
  GDD_weight<-rep(0,length(HourTemp))
  GDD_weight[which(HourTemp>=Tlow&HourTemp<=Tup)]<-(HourTemp[which(HourTemp>=Tlow&HourTemp<=Tup)]-Tbase)/24
  GDD_weight[which(HourTemp<Tlow)]<-(Tlow-Tbase)/24
  GDD_weight[which(HourTemp>Tup)]<-(Tup-Tbase)/24
  
  if (summ) return(cumsum(GDD_weight)) else
    return(GDD_weight)
}
