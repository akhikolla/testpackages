#' Confidence intervals for ranks
#' 
#' This function calculates simultaneous confidence (sets) intervals (CIs) at a pre-specified level (1-alpha) for the ranks of centers mu_1,...,mu_n which are observed through a sample y using multiple testing techniques. Several possibilities are presented through a "Method" variable. There are bascially two main choices; one which uses the partitioing principle and the likelihood ratio test and the the other is based on Tukey's pairwise comparison procedure. See choices below, and for more details see the references.
#' @param y a real vector of observed data.
#' @param sigma a vector of standard deviations. If sigma is a single value, then we consider that all centers have the same standard deviation.
#' @param Method a character indicating the method used to produce the confidence intervals. 
#'  The "ExactLR" produces confidence intervals using the partitioning principle and the likelihood ratio test. 
#'  The "BoundLR" choice produces lower- or upper-bound confidence intervals (according to the "BoundChoice") for the ranks using a fast algorithm. 
#'  The "Tukey" choice produces simultaneous confidence intervals for the ranks using Tukey's HSD. 
#'  The "SeqTukey" produces simultaneous confidence intervals for the ranks using a sequential-rejective algorithm. 
#'  The "Approximate" choice provides approximate confidence intervals which are shorter than the exact ones by considering a subset of the partitions (the correctly ordered ones, see refs and below for details). 
#'  The "TukeyNoTies" choice calculates a readustement for Tukey's method under the assumption that there are no ties and then use Tukey's method again with adjusted level.
#'  The "RescaledExactLR" choice calculates a readustement for the "ExactLR" method by adjusting each and every local test.
#'  The "RescaledTukey"  choice calculates a readustement for the "Tukey" method by pluging it into a partitioning procedure and then adjusting each and every local test.
#' @param BoundChoice a character entry which is only relevant if the "Bound" choice is picked in the Method parameter. The default value is "Upper" which results in the upper-bound CIs for the ranks. If "Lower" is chosen, then the lower-bound CIs are generated.
#' @param ApproxAlgo a character entry ("Upper" by default). This parameter controls which approximation is to be used. 
#' @param alpha the significance level of the internal tests we perform (which corresponds to the FWER control of the corresponding multiple testing procedure). CIs have simultaneous significance level of 1-alpha.
#' @param control is a list of control parameters.
#' @param crit is the critical value for Tukey's HSD. If it is kept NULL, then it is calculated internally. The use of this parameter becomes handful in case the user wishes to make several simulations. By providing it, we avoid repeating a Monte-Carlo estimation of the quantile and thus we gain in execution time. In some cases (espcially when all centers have the same standard deviation), the critical value for Tukey's HSD can be found in some statistical tables.
#' @param trace a logical parameter which supresses the printing of the details of the method which was chosen. The default is TRUE (shows details).
#' @param adjustL a logical variable (default to FALSE) indicating if an adjustment on the lower bound according to the data must be considered (if possible). This choice is only relevenat if Method is chosen as "BoundLR" and BoundChoice is chosen as "Lower".
#' @param adjustU a logical variable (default to FALSE) which gives the user the choice to adjust the upper bound CIs through the parameter "n_adjust". This choice is only relevenat if Method is chosen as "BoundLR" and BoundChoice is chosen as "Upper".
#' @param n_adjust an integer-valued entry for advanced control over the lower- or upper-bound algorithms. When the "adjustL" parameter is TRUE, the new value of n_adjust is chosen automatically as the best adjustment on the lower affine bound of the chi-square quantiles according to the data. If adjustU is TRUE, then n_adjust contains the point at which the upper affine bound is tangent on the chi-square quantiles. Possible values {1,...,n-1}. If both adjustL and adjustU variables are left FALSE, then the default choice is that the lower affine bound passes between the chi-square quantiles at 1 and n-1 degrees of freedom, and the upper affine bound is tangent on n-1.
#' @param N the number of iterations used in order to calculate the Studentized range quantile for Tukey's algorithms.
#' @param MM the number of Monte-Carlo simultations required to estimate the (simultaneous) coverage. This is used in all rescaling methods.
#' @param RandPermut is the number of additional permutations to perform when using either the "ExactLR" or the "BoundLR" algorithms and only when the standard deviations are different. When the standard deviations are the same, this has no influence on the result.
#' @param SwapPerm corresponds to performing swap permutations (yes = TRUE, no = FALSE). This is used in all the methods except for "Tukey" and "ExactLR" (the latter when the standard deviations are not the same).
#' @details The vector of observations needs to be sorted. Otherwise, it is done internally. The observations are supposed to be independent realizations of Guassian distributions with unknown centers mu_1,...,mu_n and known standard deviations sigma = (sigma_1,...,sigma_n). 
#' @details The exact-partitioning confidence intervals (option "ExactLR") are calculated using an algorithm with exponential complexity. The hypotheses in each level of the partitioning are coded using the combinatorial number system.
#' @details The lower- and upper-bound CIs are calculated with a polynomial algorithm. The bracketing obtained from the lower and upper bounds is generally very narrow with a maximum gap of 1. Moreover, in regular situations, the lower and upper bounds coincide on at least 50 percent of the centers yielding the exact-partitioning result. Thus, the bracketing is an alternative for an exact-partitioning algorithm for medium- and large-size samples (n>50). When a calculus of the lower- and upper-bound CIs is required, the default choice is when no adjustment on neither the lower nor the upper bounds is taken into account. Thus, the lower affine bound of the chi-square is a line passing by the quantiles at 1 and n-1 degrees of freedom, whereas the upper affine bound is a line tangent on the chi-square quantiles at n-1 degrees of freedom. The adjustment on the lower bound CIs can in some contexts improve on the CIs and increase the number of centers where the lower and upper bounds coincide. The best option is to adjust for both the lower and upper bounds (separately).
#' @details Both "Tukey" and "SeqTukey" are based on multiple comparison testing and are superior to the LR-based CIs if the centers are far apart from each other and if the standard deviations are not significantly different from each other. The sequential rejective variant of Tukey's HSD rejects at least as much as Tukey's HSD and thus produces generally shorter confidence intervals for the ranks. 
#' @details The "TukeyNoTies" method assumes that the true vector of parameters has no ties and therefore, instead of calculating a quantile q corresponding to mu=0 with set rank [1,n] for mu_i, we calculate a quantile corresponding to mu=0 with rank \{i\} for mu_i. The method provides shorter SCI for the ranks but is still conservative.
#' @details When the standard deviations are not the same for all the means, the methods based on the partitioning principle are not guaranteed to produce the same results. The "Block" algorithm, however, is always compatible with the lower and upper CIs provided by option "BoundLR". When the number of means exceeds 10, then performing any method based on the partitioning procedure requires a long execution time since the complexity of the algorithm is super exponential of order exp(exp(n)).
#' @details When the standard deviations are not the same the approximate methods based on the LRT are not guaranteed to cover and if the standard deviations are very different, the resulting SCIs are anticonservative. If the standard deviations are close to each other, then the result is still conservative.
#' @details In terms of execution time. The Tukey method is the fastest. It can be used always. The methods based on the partitioning principle have all exponential complexity. Therefore, when the standard deviations are the same, the "ExactLR" would produce results up to 40 means. When they are not the same, no method based on the partitioning principle can be used for more than 10 means unless we limit the number of random permutations that we use which in case of great differences in the standard deviations might lead to anticonservative results. More details can be found in the references. 
#' @author Diaa Al Mohamad and Jelle J. Goeman and Erik W. van Zwet. Correspondence can be made to diaa.almohamad@@gmail.com
#' @return a list of two vectors containing the lower and upper bounds of the confidence intervals for the sorted observed centers.
#' @references Diaa Al Mohamad and Erik W. van Zwet and Jelle J. Goeman and Aldo Solari, Simultaneous confidence sets for ranks using the partitioning principle - Technical report (2017). https://arxiv.org/abs/1708.02729
#' @references Diaa Al Mohamad and Jelle J. Goeman and Erik W. van Zwet, An improvement of Tukey's HSD with application to ranking institutions (2017). https://arxiv.org/abs/1708.02428
#' @references Diaa Al Mohamad and Jelle J. Goeman and Erik W. van Zwet, Simultaneous Confidence Intervals for Ranks With Application to Ranking Institutions (2018). https://arxiv.org/abs/1812.05507
#' @examples
#' n = 5; TrueCenters = 1:n
#' alpha = 0.05; sigma = rep(0.5,n)
#' y = as.numeric(sapply(1:n, function(ll) rnorm(1,TrueCenters[ll],sd=sigma[ll])))
#' ind = sort.int(y, index.return = TRUE)$ix
#' y = y[ind]
#' sigma = sigma[ind] # The sigmas need to follow the order of the y's
#' res = ic.ranks(y, sigma, Method = "ExactLR",alpha = 0.05, control = list(trace = TRUE))
#' LowerExact = res$Lower; UpperExact = res$Upper
#' #res = ic.ranks(y, sigma, Method = "BoundLR", BoundChoice = "Lower",
#' #   control = list(adjustL = FALSE, adjustU = FALSE))
#' #LowerL = res$Lower; UpperL = res$Upper
#' #res = ic.ranks(y, sigma, Method = "BoundLR", BoundChoice = "Upper",
#' #   control = list(adjustL = FALSE, adjustU = FALSE, trace=FALSE))
#' #LowerU = res$Lower; UpperU = res$Upper
#' res = ic.ranks(y, sigma, Method = "Tukey")
#' LowerTuk = res$Lower; UpperTuk = res$Upper
#' res = ic.ranks(y, sigma, Method = "SeqTukey")
#' LowerTukSeq = res$Lower; UpperTukSeq = res$Upper
#' res = ic.ranks(y, sigma, Method = "TukeyNoTies")
#' LowerTukNoTies = res$Lower; UpperTukNoTies = res$Upper
#' resLR1 = ic.ranks(y, sigma, Method = "RescaledExactLR", alpha = alpha, 
#'   control = list(trace = TRUE, gridSize = 4, MM = 100, RandPermut=factorial(n)))
#' LowerExact 
#' #LowerL
#' #LowerU
#' LowerTuk
#' resLR1$Lower
#' resLR1$Upper 
#' @export
ic.ranks = function(y, sigma = rep(1,length(y)), Method = c("ExactLR","BoundLR","Tukey","SeqTukey","ApproximateLR", "TukeyNoTies", "RescaledExactLR", "RescaledTukey"), BoundChoice = c("Upper", "Lower"), ApproxAlgo = c("Exact","Upper"), alpha = 0.05, control = list(crit = NULL, trace = TRUE, adjustL = FALSE, adjustU = FALSE, n_adjust = length(y)-1, N = 10^4, MM = 10^3, gridSize = 5, RandPermut = 0, SwapPerm = TRUE))
{
# .................................................
# --------------------------
# Some checking before we start...
# --------------------------
# .................................................
if(length(Method) != 1) Method = "SeqTukey"
trace = control$trace
if(is.null(trace)) trace = TRUE
RandPermut = control$RandPermut
SwapPerm = control$SwapPerm
if(is.null(RandPermut)) RandPermut = 0
if(is.null(SwapPerm)) SwapPerm = TRUE
if(!(Method %in% c("ExactLR","BoundLR","Tukey","SeqTukey","ApproximateLR", "TukeyNoTies", "OnlyBlock", "RescaledExactLR", "RescaledTukey"))) {print("Error! Method not supported."); return(0)}
n = length(y)
if(length(sigma) == 1) sigma = rep(sigma,n)
if(length(sigma) != n) {print("Error: sigma and y must have the same length!"); return(0)}
if(n == 1) return(1)

# Sort the observations if it is not the case.
ind = sort.int(y, index.return = T)$ix
y = y[ind]
sigma = sigma[ind]
if(sum(ind != 1:n)) print("The sample had to be sorted in ascending way. Results are shown for the sorted sample.")
ranks = NULL
if(n <= 2 & Method == "BoundLR") {print("Upper- and Lower-bound CIs require at least three centers"); return(0)}
EqSigIndex = FALSE
if(length(unique(round(sigma,16)))==1) 
{
	#RandPermut=0
	EqSigIndex = TRUE
}
# .................................................
# --------------------------
# The Full Partitioning Principle
# --------------------------
# .................................................
if(Method == "ExactLR")
{
  crit = qchisq(1-alpha,(n-1):1)
  if(sum((y - sum(y/sigma^2)/(sum(1/sigma^2)))^2/sigma^2) < qchisq(1-alpha,n-1))
  {
	if(trace==TRUE) cat("Process ended with trivial confidence intervals.\n")
	return(list(Lower = rep(1,n), Upper = rep(n,n)))
  }

  if(EqSigIndex)
  {
	ranks = PartitioningRankingLevelEqSig(y, sigma, crit, n, trace)
  }else
  {
	ranks = PartitioningRankingLevelUneqSig(y, sigma, crit, n, trace, RandPermut, SwapPerm)
  }
  
  ranks = list(Lower = ranks[,1], Upper = ranks[,2])
}
if(Method == "RescaledExactLR")
{
  MM = control$MM
  if(is.null(MM)) MM = 1000
  gridSize = control$gridSize
  if(is.null(gridSize)) gridSize = 5
  if(EqSigIndex)
  {
	#crit = sapply((n-1):1, function(ss) qchisq(1-alpha*(1:gridSize),ss))
	crit = matrix(0, nrow = n-1, ncol = gridSize)
	alph = seq(from=alpha, to = 0.4, length = gridSize)
	for(ss in 1:gridSize)
	{
		#library(gmp)
		for(i in 3:(n+1))
		{
		  w = as.numeric(abs(Stirling1.all(i))/factorial(i))
		  crit[n-(i-2),ss] = GeneralizedInvCDF(function(x) 1-sum(w[(i-1):1]*pchisq(x,df=1:(i-1),lower.tail=FALSE)),1-alph[ss], Bsup=100,npoints=10^4)
		}
	}
	ranks = PartitioningRankingLevelEqSigRescaled(y, sigma, crit, matrix(rnorm(n*MM,sd=sigma[1]),nrow=n,ncol=MM),MM, n, RandPermut, alpha, gridSize, trace)
	ranks = list(Lower = ranks[,1], Upper = ranks[,2])
  }else
  {
	cat("\n The standard deviations are not the same.\n The rescaled partitioning procedure is not implemented with this option.\n")
  }

}
if(Method == "OnlyBlock")
{
  crit = qchisq(1-alpha,(n-1):1)
  if(sum((y - sum(y/sigma^2)/(sum(1/sigma^2)))^2/sigma^2) < qchisq(1-alpha,n-1))
  {
	if(trace==TRUE) cat("Process ended with trivial confidence intervals.\n")
	return(list(Lower = rep(1,n), Upper = rep(n,n)))
  }
  ranks = OnlyBlockRanking(y, sigma, crit, n, trace, RandPermut, SwapPerm)
  ranks = list(Lower = ranks[,1], Upper = ranks[,2])
}

# .................................................
# --------------------------
# The Bracketing algorithm for the FULL Partitioning Principle
# --------------------------
# .................................................
if(Method == "BoundLR")
{
  if(length(BoundChoice)!= 1) BoundChoice = "Upper"
  if(!(BoundChoice %in% c("Upper", "Lower"))){print("Error! Could not recognize your choice whether it is upper of lower bound."); return(0)}
  adjustL = control$adjustL
  if(is.null(adjustL)) adjustL = FALSE
  adjustU = control$adjustU
  if(is.null(adjustU)) adjustU = FALSE
  n_adjust = control$n_adjust
  if(is.null(n_adjust)) n_adjust = n-1
  n_adjust = floor(n_adjust)
  if(n_adjust > n-1 | n_adjust<1){n_adjust = n-1; cat(paste("n_adjust can take values only between 1 and ", n-1)); cat(". Default value is considered.")}
 
  if(sum((y - sum(y/sigma^2)/(sum(1/sigma^2)))^2/sigma^2) < qchisq(1-alpha,n-1))
  {
	if(trace==TRUE) cat("Process ended with trivial confidence intervals.\n")
	return(list(Lower = rep(1,n), Upper = rep(n,n)))
  }
  if(BoundChoice == "Lower")
  {
	if(trace == TRUE) cat('\n Calculate lower bounds for simultaneous confidence intervals for ranks using the partitioning principle and the LRT.\n')
	ranks = ApproximatePartition(y,sigma,"Lower", alpha, 1, trace, RandPermut)
	if(adjustL == TRUE)
	{
	  ind = which(ranks$Upper == n)[1]
	  n_adjust = min(ranks$BlockMax[1:ind])+1
	  ranks = ApproximatePartition(y,sigma,"Lower", alpha, n_adjust, trace, RandPermut, SwapPerm)
	  if(trace == TRUE) cat(paste("\n Adjustment on the lower bound. Intersection with the chi-square quantile curve at n_adjust = ",n_adjust))
	}	
  }else{
	# The upper choice is left as the alternative anyway.
	if(trace == TRUE) cat('\n Calculate upper (conservative) bounds for simultaneous confidence intervals for ranks using the partitioning principle and the LRT.\n')
	if(adjustU == FALSE) n_adjust = n-1 # keep regular bounds if no adjustment is required
	if(adjustU == TRUE & trace == TRUE) cat(paste("Adjustment on the upper bound by the user. Tangent on the chi-square quantile at n_adjust = ", n_adjust))
	ranks = ApproximatePartition(y,sigma,"Upper", alpha, n_adjust, trace, RandPermut, SwapPerm)
  }
}
# .................................................
# --------------------------
# The Partitioning Principle when only the correctly ordered hypotheses are considered
# --------------------------
# .................................................
if(Method == "ApproximateLR")
{
 if(length(ApproxAlgo)!= 1) ApproxAlgo = "Upper"
 if(!(ApproxAlgo %in% c("Exact","Upper"))) {print("Error! Approximate algorithm not supported."); return(0)}
 if(trace == TRUE) cat('\n Calculate approximate simultaneous confidence intervals for ranks using the partitioning principle and the LRT.\n')
 if(ApproxAlgo == "Upper")
 {
   if(trace == TRUE) cat('\n A fast (cubic complex) algorithm is being used.\n')
   ranks = ApproximatePartitionCorrectOrder(y,sigma,"Upper",alpha,n-1)
 }else
 {
   crit = qchisq(1-alpha,(n-1):1)
   if(trace == TRUE) cat('\n A slow (exponentially complex) algorithm is being used.\n')
   res = ApproximatePartitionCorrectOrder(y, sigma, BoundChoice = "Lower", alpha,1)
   ind = which(res$Upper == n)[1]
   n_adjust = min(res$BlockMax[1:ind])+1
   res = ApproximatePartitionCorrectOrder(y,sigma,"Lower", alpha, n_adjust)
   Lower = res$Lower-1; Upper = res$Upper-1; MinBlock = res$BlockMax
   res = ApproximatePartitionCorrectOrder(y,sigma,"Upper", alpha, floor(n/2))
   MaxBlock = res$BlockMax
   ranks = PartitioningRankingBlockCorrectOrder(y, sigma, crit, MinBlock, MaxBlock, Lower, Upper, n, trace)
   ranks = list(Lower = ranks[,1], Upper = ranks[,2])
 }

}

# .................................................
# --------------------------
# General Tukey (ties are allowed)
# --------------------------
# .................................................
if(Method == "Tukey")
{
 N = control$N
 if(is.null(N)) N = 10^4
 crit = control$crit
 if(is.null(crit))
 {
 if(length(unique(sigma)) == 1 & alpha<0.3) # a quick calculation
 {
   crit = qtukey(1-alpha,n,Inf)/sqrt(2)
 }else
 {
 x=t(mapply(rnorm,N,0,sigma))
 if(n<100){
	Cp=contrMat(rep(1,n), type ="Tukey")    # only dependence on multcomp 
	S<-Cp%*%diag(sigma^2)%*%t(Cp)             # covariance of differences 
	std=diag(S)^(-1/2)                      # 1/(standard deviations of differences)
	d=diag(std)%*%Cp%*%x                    # standardized differences
	crit=quantile(apply(abs(d),2,max),1-alpha)   
	rm(d); rm(std); rm(S); rm(Cp); rm(x)
  }else
  {
	# Quantile calculus for larger sample size
	Diff = numeric(N)
	for(i in 1:N)
	{
	  for(j in 1:(n-1))
	  {
		Diff[i] = max(Diff[i],abs(x[j,i]-x[(j+1):n,i])/sqrt(sigma[j]^2+sigma[(j+1):n]^2))
	  }
	}
  
	crit = quantile(Diff,1-alpha)
	rm(Diff)
	rm(x)
  }
 }
 }
 
 ranks = tukey(y,sigma,crit)
 if(trace == TRUE) cat(paste("\n Confidence intervals for ranks calculated using Tukey's HSD procedure at simultaneous level", 1-alpha))
}

# .................................................
# --------------------------
# Sequential Tukey
# --------------------------
# .................................................
if(Method == "SeqTukey")
{
 N = control$N
 if(is.null(N)) N = 10^4
 ranks = StepDownTukeySeqRej(y,sigma,alpha, N)
 if(trace == TRUE)
 {
   cat(paste("\n Confidence intervals for ranks calculated using a sequential-rejective variant of Tukey's HSD procedure at simultaneous level", 1-alpha))
   cat(paste("\n Number of iterations = ",ranks$NbSteps))
   cat("\n")
 }
}
# .................................................
# --------------------------
# Tukey under the assumption that there are no ties
# --------------------------
# .................................................
if(Method == "TukeyNoTies")
{
if(trace == TRUE) cat('\n Caclulating an adjusted alpha...\n')
N = control$N
MM = control$MM
if(is.null(N)) N = 10^4
if(is.null(MM)) MM = 10^3
d = numeric(N) # a vector of differences that might be useful later
sigmaTree1 = sigma
if(length(unique(sigma)) > 1 | alpha<0.3)
{
 OddInd = seq(from=1, to = n, by = 2); EvenInd = seq(from=2,to = n, by = 2)
 sigmaTree1 = sort(sigma)
 sigmaTree1 = sigmaTree1[c(EvenInd, OddInd[((n+1)/2):1])]
 x=t(mapply(rnorm,N,0,sigmaTree1))
 if(n<100)
 {
	Cp=contrMat(rep(1,n), type ="Tukey")    # only dependence on multcomp 
	S<-Cp%*%diag(sigmaTree1^2)%*%t(Cp)             # covariance of differences 
	std=diag(S)^(-1/2)                      # 1/(standard deviations of differences)
	d=diag(std)%*%Cp%*%x                    # standardized differences
	d = apply(abs(d),2,max)
	rm(std); rm(S); rm(Cp); rm(x)
  }else
  {
	# Quantile calculus for larger sample size
	d = numeric(N)
	for(i in 1:N)
	{
	  for(j in 1:(n-1))
	  {
		d[i] = max(d[i],abs(x[j,i]-x[(j+1):n,i])/sqrt(sigmaTree1[j]^2+sigmaTree1[(j+1):n]^2))
	  }
	}
  	rm(x)
  }
}

# Estimate the coverage of Tukey at alpha when mu=0
x=t(mapply(rnorm,MM,0,sigmaTree1))
TukeyCoverage = function(a)
{
 if(trace == TRUE) {cat(a);cat('\n')}
 q = 1
 if(length(unique(sigmaTree1)) > 1 | alpha<0.3)
 {
	q=quantile(d,1-a)
 }else
 {
	q = qtukey(1-a,n,Inf)/sqrt(2)
 }
 TrueLowerRank = 1:n; TrueUpperRank = 1:n
 coverageTuk = MM
 for(i in 1:MM)
 {
  y = x[,i]
  ind = sort.int(y, index.return = T)$ix
  y = y[ind]
  resTukey = tukey(y,sigmaTree1[ind], q)
  if(sum(TrueLowerRank[ind]<resTukey$Lower | TrueUpperRank[ind]>resTukey$Upper)>0) coverageTuk = coverageTuk - 1
 }
 coverageTuk/MM
}

# Calculate a modified significance level
alphaTuk = uniroot(function(a)TukeyCoverage(a)-(1-alpha), c(alpha,0.9), maxiter=15)$root
rm(x)
if(trace == TRUE) cat("Applying Tukey's HSD using the new alpha...\n")
# Use the new alpha to calculate the SCI for ranks
crit = 1
if(length(unique(sigma)) == 1 & alphaTuk<0.3) # a quick calculation
{
   crit = qtukey(1-alphaTuk,n,Inf)/sqrt(2)
}else
{
 x=t(mapply(rnorm,N,0,sigma))
 if(n<100)
 {
	Cp=contrMat(rep(1,n), type ="Tukey")    # only dependence on multcomp 
	S<-Cp%*%diag(sigma^2)%*%t(Cp)             # covariance of differences 
	std=diag(S)^(-1/2)                      # 1/(standard deviations of differences)
	d=diag(std)%*%Cp%*%x                    # standardized differences
	crit=quantile(apply(abs(d),2,max),1-alphaTuk)   
	rm(d); rm(std); rm(S); rm(Cp); rm(x)
 }else
 {
	# Quantile calculus for larger sample size
	Diff = numeric(N)
	for(i in 1:N)
	{
	  for(j in 1:(n-1))
	  {
		Diff[i] = max(Diff[i],abs(x[j,i]-x[(j+1):n,i])/sqrt(sigma[j]^2+sigma[(j+1):n]^2))
	  }
	}
  
	crit = quantile(Diff,1-alphaTuk)
	rm(Diff)
	rm(x)
 }
}
 
 ranks = tukey(y,sigma, crit)
 if(trace == TRUE) 
 {
	cat(paste("\n Confidence intervals for ranks calculated using a rescaled version \n of Tukey's HSD procedure at simultaneous level", 1-alpha))
	cat(paste("\n Rescaled significance level is ",alphaTuk)); cat(".")
 }
 
  
}


# .................................................
# --------------------------
# Tukey Rescaled
# --------------------------
# .................................................
if(Method == "RescaledTukey")
{
 N = control$N
 MM = control$MM
 gridSize = control$gridSize
 if(is.null(N)) N = 10^4
 if(is.null(MM)) MM = 10^3
 if(is.null(gridSize)) gridSize = 5
 if(EqSigIndex == 1 & alpha<0.3) # a quick calculation
 {
	crit = qtukey(1-seq(from=alpha, to = 0.4, length = gridSize),n,Inf)/sqrt(2)
 }else
 {
 	x=t(mapply(rnorm,N,0,sigma))
	Cp=contrMat(rep(1,n), type ="Tukey")    # only dependence on multcomp 
	S<-Cp%*%diag(sigma^2)%*%t(Cp)             # covariance of differences 
	std=diag(S)^(-1/2)                      # 1/(standard deviations of differences)
	d=diag(std)%*%Cp%*%x                    # standardized differences
	crit=quantile(apply(abs(d),2,max),1-seq(from=alpha, to = 0.4, length = gridSize))   
	rm(d); rm(std); rm(S); rm(Cp); rm(x)
 }
 
 if(EqSigIndex == 1)
 {
 	ranks = TukeyRankingLevelEqSigRescaled(y, sigma, as.matrix(crit), t(mapply(rnorm,MM,0,sigma)), MM, n, RandPermut, alpha, gridSize, trace)
 }else
 {
	ranks = TukeyRankingLevelUneqSigRescaled(y, sigma, as.matrix(crit), t(mapply(rnorm,MM,0,sigma)), MM, n, RandPermut, alpha, gridSize, trace)
 }
 if(trace == TRUE) cat(paste("\n Confidence intervals for ranks calculated using a rescaled version \n of Tukey's HSD procedure at simultaneous level", 1-alpha))
 ranks = list(Lower = ranks[,1], Upper = ranks[,2])
}



################
#-------------------
# END OF FUNCTION IC.RANKS
#-------------------
################

 if(trace == TRUE) cat(paste(paste("\n Number of compared centers is ",n),"\n"))
 return(list(Lower = ranks$Lower, Upper = ranks$Upper))
}

############################################
# Auxiliary functions. Internal usage only
# Generalized inverse of a cumulative distribution function
GeneralizedInvCDF = function(CdfFun, proba = 0.95, Binf = 0, Bsup = 100,npoints=1000)
{
  knots = seq(from=Binf,to=Bsup,length=npoints)
  yVal = as.numeric(sapply(1:npoints, function(ll)CdfFun(knots[ll])))
  ind = min(which(yVal>=proba))
  knots[ind]
}

# The approximatePartition function
ApproximatePartition = function(y, sigma, BoundChoice = c("Upper", "Lower"), alpha = 0.05, n_adjust, trace = FALSE, RandPermut = 0, SwapPerm = TRUE)
{
critFun = function(x)
{
 if(x<=0) return(0)
 
 slop*x
}

n = length(y)
z = qchisq(1-alpha,1:(n-1))
slop = z[n_adjust] - z[n_adjust-1]; Intercept = z[n_adjust] - slop*n_adjust

if(BoundChoice == "Lower") {slop = (z[n-1] - z[n_adjust])/(n-n_adjust-1); Intercept = z[n_adjust] - slop*n_adjust}

### Part 0: calculate the SCI when no permutation is applied
# Calculate the matrix of contributions
if(trace == TRUE) cat('\n Caclulate simultaneous confidence intervals using the correctly ordered hypotheses.\n')
EmpOrder = 1:n
res = ApproximatePartitionCorrectOrder(y, sigma, BoundChoice = BoundChoice, alpha = alpha, n_adjust=n_adjust)
Lower = res$Lower; Upper = res$Upper
if(sum(Lower==rep(1,n) & Upper==rep(n,n)) == n) return(list(Lower=Lower,Upper=Upper))
minY = min(y)
maxY = max(y)

res = ApproximatePartitionPermutations(y, sigma, Lower, Upper, n, slop, Intercept, minY, maxY, trace, SwapPerm, RandPermut)
Lower = res[,1]; Upper = res[,2]
return(list(Lower = Lower, Upper = Upper))
}

# The function without permutations and the centers are assumed ordered.
ApproximatePartitionCorrectOrder = function(y, sigma, BoundChoice = c("Upper", "Lower"), alpha = 0.05, n_adjust)
{
#library(numDeriv)
n = length(y)
if(sum((y - sum(y/sigma^2)/(sum(1/sigma^2)))^2/sigma^2) < qchisq(1-alpha,n-1))
{
	return(list(Lower = rep(1,n), Upper = rep(n,n), BlockMax = rep(n-1,n)))
}

critFun = function(x)
{
 if(x<=0) return(0)
 
 slop*x
}

z = qchisq(1-alpha,1:(n-1))
slop = z[n_adjust] - z[n_adjust-1]; Intercept = z[n_adjust] - slop*n_adjust

if(BoundChoice == "Lower") {slop = (z[n-1] - z[n_adjust])/(n-n_adjust-1); Intercept = z[n_adjust] - slop*n_adjust}

# Calculate the matrix of contributions
LogL = matrix(0, nrow = n, ncol = n)
IndividContribBlock = matrix(0, nrow = n, ncol = n)
for(j in 2:n)
{
 for(i in (j-1):1)
 {
	LogL[i,j] = sum((y[i:j] - sum(y[i:j]/sigma[i:j]^2)/sum(1/sigma[i:j]^2))^2/sigma[i:j]^2)
	
	IndividContribBlock[i,j] = min(IndividContribBlock[i,i:(j-1)]+IndividContribBlock[(i+1):j,j], LogL[i,j] - critFun(j-i))
	# The intercept needs to be substracted so that we do not count it several times, because the sum of contributions is not the contribution of the sum
 }
}

# Test the blocks by adding hypotheses to the left and to the right of it
Lower = 1:n; Upper = 1:n
BlockMax = numeric(n)
# Treat the case of mu_1=...=mu_j + hypothesis
for(j in (n-1):2)
{
	if(LogL[1,j] - critFun(j-1) + IndividContribBlock[j+1,n] - Intercept < 0) # Block [1,j] is accepted. Update the CIs.
	{
	  Lower[1:j] = 1
	  Upper[1:j] = pmax(Upper[1:j], j)
	  BlockMax[1] = j - 1
	  break
	}
}
# Treat the case of hypothesis + mu_i=...=mu_n
for(i in 2:(n-1))
{
	if(LogL[i,n] - critFun(n-i) + IndividContribBlock[1,i-1] - Intercept < 0) # Block [i,n] is accepted. Update the CIs.
	{
	  Lower[i:n] = pmin(Lower[i:n],i)
	  Upper[i:n] = n
	  break
	}
}
# Treat the case of hypothesis + mu_i=...=mu_j + hypothesis
for(i in 2:(n-2))
{#print(i)
 for(j in (n-1):(i+1))
 {
	if(LogL[i,j] - critFun(j-i) + IndividContribBlock[1,i-1] + IndividContribBlock[j+1,n] - Intercept < 0) # Block [i,j] is accepted. Update the CIs.
	{
	  Lower[i:j] = pmin(Lower[i:j],i)
	  Upper[i:j] = pmax(Upper[i:j], j)
	  BlockMax[i] = j-i
	  break
	}
 }
}
return(list(Lower = Lower, Upper = Upper, BlockMax = BlockMax))
}

####################################################
# The simple HSD
#library(multcomp)
tukey = function(y,sigma, qq) {
  n=length(y)
  ranks=matrix(0,n,2)
  for(j in 1:n)
  {
	stat = (y[j]-y)/sqrt(sigma[j]^2+sigma^2)
	ranks[j,1]=1+sum(stat>qq)
	ranks[j,2]=n-sum(stat<(-qq))
  }

  return(list(Lower = ranks[,1], Upper = ranks[,2]))
}

###########################################################
# Step-down Tukey with sequential rejection
StepDownTukeySeqRej = function(y,sigma,alpha=0.05,N = 10^4)
{
  n = length(y)
 # Calculate the critical value
 Diff = NULL
 PosPairs = matrix(0,ncol=2,nrow = n*(n-1)/2)
 NegPairs = matrix(0,ncol=2,nrow = n*(n-1)/2)
 for(j in 1:(n-1))
 {
  PosPairs[(j-1)*n-j*(j-1)/2+1:(n-j),1] = rep(j,n-j)
  PosPairs[(j-1)*n-j*(j-1)/2+1:(n-j),2] = (j+1):n  
 }
 NegPairs[,2] = PosPairs[,1] # Those pairs stay untouched because they are never rejected
 NegPairs[,1] = PosPairs[,2] 
# set.seed(16021988)
 x=t(mapply(rnorm,N,0,sigma))  
 Diff = numeric(N)
 for(k in 1:N)
 {
   for(j in 1:(n-1))
  {
	Diff[k] = max(Diff[k],abs(x[j,k]-x[(j+1):n,k])/sqrt(sigma[j]^2+sigma[(j+1):n]^2))
  }
 }

 qDown = quantile(Diff,1-alpha)
 rm(Diff)

 NbNRejOld = n*(n-1)/2; NbNRejNew = 0
 NBSteps = 0
 while(NbNRejOld>NbNRejNew)
 {
  NBSteps = NBSteps +1
  NewPairs = NULL

  for(i in 1:length(PosPairs[,1]))
  {
	T = abs(y[PosPairs[i,1]]-y[PosPairs[i,2]])/sqrt(sigma[PosPairs[i,1]]^2+sigma[PosPairs[i,2]]^2)
	if(T<qDown)
	{
	  NewPairs = rbind(NewPairs,PosPairs[i,])
	  #Diff = c(Diff,T)
	}
  }
 NbNRejOld = length(PosPairs[,1])
 PosPairs = NewPairs
 NbNRejNew = length(PosPairs[,1])
 
 # Calculate again the quantile  
 AllPairs = rbind(PosPairs, NegPairs)
 Diff = numeric(N)
 for(i in 1:N)
 {
	Diff[i] = max((x[AllPairs[,2],i]-x[AllPairs[,1],i])/sqrt(sigma[AllPairs[,1]]^2+sigma[AllPairs[,2]]^2))
 }
 
 qDown = quantile(Diff,1-alpha)
 rm(Diff)

 }
rm(x)
# Extract the ranks
ranks = matrix(0,nrow = n,ncol = 2)
ranks[,1] = 1:n
ranks[,2] = 1:n
ranks[1,2] = 1+sum(PosPairs[,1]==1)
for(i in 2:(n-1))
{
	ranks[i,1] = i - sum(PosPairs[,2] == i)
	ranks[i,2] = i + sum(PosPairs[,1]==i)
}
ranks[n,1] = n - sum(PosPairs[,2] == n)

return(list(Lower = ranks[,1], Upper = ranks[,2], NbSteps = NBSteps))
}




