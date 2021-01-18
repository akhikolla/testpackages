\name{ISOpure.step2.PPE}
\alias{ISOpure.step2.PPE}
\title{Perform second step of ISOpure purification algorithm}
\description{Performs the second step of the ISOpure purification algorithm, taking tumor data and
	normal profiles and returning the a list, ISOpureS2model, with all the updated parameters.}
\usage{
ISOpure.step2.PPE(tumordata, BB, ISOpureS1model, MIN_KAPPA, logging.level) 
}
\arguments{
	% Right now there is not full backwards compatibility with ISOLATE
	\item{tumordata}{(same as for ISOpureS1) a GxD matrix representing gene expression profiles of 
		heterogeneous (mixed) tumor samples, where G is the number of genes, D is the number of tumor
		samples.}

	\item{BB}{(same as for ISOpureS1) represents B = [b_1 ... b_(K-1)] matrix (from Genome Medicine 
		paper) a Gx(K-1) matrix, where (K-1) is the number of normal profiles 
		\eqn{(\beta_1,...,\beta_(K-1))}, G is the number of genes.  These are the normal profiles 
		representing normal cells that contaminate the tumor samples  (i.e. normal samples from the 
		same tissue location as the tumor).  The minimum element of BB must be greater than 0 -- 
		i.e. every gene/transcript must be observed on some level in each normal sample.}

	\item{ISOpureS1model}{output model list from ISOpureS1 code} %If only primary 
		%tumors from the same site of origin are represented in tumordata, then this is the same as
		%BB (default behavior).   This parameter is for backwards compatibility and replacing the 
		%original ISOLATE code from the 2009 Bioinformatics paper, and can represent potential sites 
		%of origins of the metastatic tumor (in which case tumordata represents one or more 
		%expression profiles of the secondary tumor).  Set PP=BB for default behavior, and if you 
		%need to specify MIN_KAPPA.} % *Optional* in ISOLATE!

	\item{MIN_KAPPA}{(optional) The minimum value allowed for the strength parameter kappa placed
		over the reference cancer profile m (see Quon et al, 2013).  By default, this is set to 
		1/min(BB), such that the log likelihood of the model is always finite.  However, when the 
		min(BB) is very small, this forces MIN_KAPPA to be very large, and can sometimes cause the 
		reference profile m to look too much like a 'normal profile' (and therefore you may observe 
		the tumor samples having low \% cancer content estimates).  If this is the case, you can try 
		setting MIN_KAPPA=1, or some other small value.  For reference, for the data presented in 
		Quon et al., 2013, MIN_KAPPA is on the order of 10^5.}

	\item{logging.level}{(optional) A string that gives the logging threshold for futile.logger. The
		possible options are 'TRACE', 'DEBUG', 'INFO', 'WARN', 'ERROR', 'FATAL'. Currently the 
		messages in ISOpureR are only in the categories 'INFO', 'WARN', and 'FATAL', and the default 
		setting is 'INFO'.  Setting a setting for the entire package will over-ride the setting for
		a particular function.}
}
\references{
  G Quon, S Haider, AG Deshwar, A Cui, PC Boutros, QD Morris. \emph{Computational purification of 
  individual tumor gene expression profiles}. Genome Medicine (2013) 5:29, 
  \url{http://genomemedicine.com/content/5/3/29}.
  
  G Quon, QD Morris. \emph{ISOLATE: a computational strategy for identifying the primary origin of
  cancers using high-thoroughput sequencing}. Bioinformatics 2009, 25:2882-2889 
  \url{http://bioinformatics.oxfordjournals.org/content/25/21/2882}.

}
\value{ISOpureS2model, a list with the following important fields:

	\item{theta}{a DxK matrix, giving the fractional composition of each tumor sample.  Each row 
		represents a tumor sample that was part of the input, and the first K-1 columns correspond 
		to the fractional composition with respect to the Source Panel contaminants.  The last 
		column represents the fractional composition of the pure cancer cells.  In other words,	
		each row sums to 1, and element (i,j) of the matrix denotes the fraction of tumor i 
		attributable to component j (where the last column refers to cancer cells, and the first 
		K-1 columns refer to different 'normal cell' components).  The 'cancer', or tumor purity, 
		estimate of each tumor is simply the last column of theta.}

	\item{alphapurities}{(same as ISOpureS1) tumor purities (alpha_i in paper), same as the last 
		column of the theta variable, pulled out for user convenience - not changed in step 2}

	\item{cc_cancerprofiles}{purified cancer profiles.  This matrix is of the same dimensionality as 
		tumordata, and is also on the same scale (i.e. although ISOpureS2 treats purified cancer 
		profiles as parameters of a multinomial distribution, we re-scale them to be on the same 
		scale as the input tumor profiles -- see Genome Medicine paper). Column i of
		cc_cancerprofiles corresponds to column i of tumordata.}

	\item{total_loglikelihood}{log likelihood of the model}

	\item{omega}{(internal parameter, same as ISOpureS1) prior over the reference cancer profile -
		not changed in step 2}

	\item{vv}{(internal parameter) hyper-parameters from Dirichlet distribution, representing both
		mean and strength of a Dirichlet distribution over theta}

	\item{kappa}{(internal parameter) the strength parameter over the Dirichlet distribution over 
		cc, given the reference cancer parameter, mm}

	\item{mm_weights, theta_weights, omega_weights}{(internal parameters) used in the optimization
    	of mm, theta, and omega (instead of performing constrained optimization on these 
    	positively constrained variables directly, we optimize their logs in an unconstrained 
    	fashion.)}

    \item{log_BBtranspose, PPtranspose, log_all_rates:}{ (internal parameters) used in the 
    	calculations of loglikelihood}

    \item{MIN_KAPPA}{(internal parameter) as described in the Arguments section}
    }
\author{Gerald Quon, Catalina Anghel, Francis Nguyen}
\keyword{math}
\concept{deconvolution}
\concept{cancer}
\concept{mRNA expression}
\concept{Bayesian}
\concept{statistics}
