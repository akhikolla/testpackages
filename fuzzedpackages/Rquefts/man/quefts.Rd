\name{quefts}


\alias{quefts}
\alias{run}
\alias{run,Rcpp_QueftsModel-method}

\alias{quefts_crop}
\alias{quefts_soil}
\alias{quefts_fert}
\alias{quefts_biom}
\alias{soil<-}
\alias{crop<-}
\alias{fert<-}
\alias{biom<-}

\alias{soil<-,Rcpp_QueftsModel,list-method}
\alias{crop<-,Rcpp_QueftsModel,list-method}
\alias{fert<-,Rcpp_QueftsModel,list-method}
\alias{biom<-,Rcpp_QueftsModel,list-method}


\alias{[,Rcpp_QueftsModel,character,missing-method}
\alias{[,Rcpp_QueftsModel,character,character-method}
\alias{[,Rcpp_QueftsCrop,character,missing-method}
\alias{[,Rcpp_QueftsSoil,character,missing-method}

\alias{[<-,Rcpp_QueftsModel,character,missing-method}
\alias{[<-,Rcpp_QueftsModel,character,character-method}
\alias{[<-,Rcpp_QueftsCrop,character,missing-method}
\alias{[<-,Rcpp_QueftsSoil,character,missing-method}


\title{
QUEFTS model
}

\description{
Create a QUEFTS model, set parameters, and run it to compute nutrient requirements and nutrient limited yield.

A number of default crop parameter sets are provided, as well as one example soil. You need to provide attainable crop production (in this context that is the maximum production in the absence of nutrient limitation), or target dry-matter biomass for leaves, stems and the storage organ (e.g. grain, root or tuber). Some crops are grown for the stems/leaves, in which case there is no relevant storage organ (e.g. sugarcane, jute). production yield estimates can be obtained with a crop growth model. 

For a cereal crop you can assume that 50% of the total biomass is grain, and about 25% is stem and 25% is leave biomass.
}

\usage{
quefts(soil, crop, fert, biom)
quefts_soil()
quefts_fert()
quefts_crop(name="")
quefts_biom()
crop(x) <- value
soil(x) <- value
fert(x) <- value
biom(x) <- value
run(x, ...)
}

\arguments{
 \item{soil}{list with named soil parameters. See Details. An example is returned by \code{quefts_soil()}}
 \item{crop}{list with named crop parameters. See Details. An example is returned by \code{quefts_crop()}} 
 \item{fert}{list with named fertilizer parameters (N, P and K). An example is returned by \code{quefts_fert()}} 
 \item{biom}{list with named biomass and growing season length parameters. An example is returned by \code{quefts_biom()}} 
 \item{name}{character. crop name} 
 \item{x}{QueftsModel object}
 \item{value}{list with soil, crop, fertilizer, or biomass parameters as above}
 \item{...}{additional arguments. None implemented}
}


\details{

\tabular{lll}{
	\tab \bold{Input Parameters} \tab \bold{Explanation} \cr
	\tab \bold{Soil} \tab . \cr
	\tab \code{N_base_supply, P_base_supply, K_base_supply} \tab Potential supply (kg/ha) of N, P and K of the (unfertilized) soil in a growing season of standard length (default is 120 days). See \code{\link{nutSupply2}} to estimate that. \cr
	\tab \code{N_recovery, P_recovery, K_recovery} \tab Fertilizer recovery, that is, the fraction of applied fertilizer that can be taken up by the plant. \cr
	\tab \code{UptakeAdjust} \tab Two-column matrix to compute the fraction uptake from soil supply as function of length of season. The default standard season length is 120 days. The first column is the length of the growing season, the second column is the fraction uptake. Intermediate values are computed by linear interpolation. \cr
	
	\tab \bold{Crop} \tab . \cr
	\tab \code{_minVeg, _maxVeg, _minStore, _maxStore} \tab minimum and maximum concentration of "\code{_}" (N, P, or K) in vegetative organs and in storage organs (kg/kg) \cr
	\tab \code{Yzero} \tab the maximum biomass of vegetative organs at zero yield of storage organs (kg/ha) \cr
	\tab \code{Nfix} \tab the fraction of a crop's nitrogen uptake supplied by biological fixation \cr

	\tab \bold{Management} \tab . \cr
	\tab \code{N, P, K} \tab N, P, and K fertilizer applied. \cr

	\tab \bold{Crop yield} \tab . \cr
	\tab \code{leaf_att, stem_att, store_att} \tab Attainable (in the absence of nutrient limitation), or target crop biomass (dry-matter, kg/ha) for leaves, stems and storage organs.\cr
	\tab \code{SeasonLength} \tab Length of the growing season (days)\cr

	
	\tab . \tab \cr
	\tab . \tab \cr
	\tab \bold{Output Variables} \tab \bold{Explanation} \cr
	\tab \code{N_actual_supply, P_actual_supply, K_actual_supply} \tab nutrient uptake from soil (not fertilizer) (kg/ha) \cr
	\tab \code{leaf_lim, stem_lim, store_lim} \tab nutrient limited biomass of leaves, stems, and storage organ (kg/ha) \cr
	\tab \code{N_gap, P_gap, K_gap} \tab fertilizer required to reach the specified biomass (kg/ha) \cr
}
}


\value{
vector with output variables as described in the Details
}

\references{
Janssen B.H., F.C.T. Guiking, D. van der Eijk, E.M.A. Smaling, J. Wolf and H. van Reuler, 1990. A system for the quantitative evaluation of the fertility of tropical soils (QUEFTS). Geoderma 46: 299-318

Sattari, S.Z., M.K. van Ittersum, A.F. Bouwman, A.L. Smit, and B.H. Janssen, 2014. Crop yield response to soil fertility and N, P, K inputs in different environments: Testing and improving the QUEFTS model. Field Crops Research 157: 35-46
}


\examples{

# create a QUEFTS model
# 1. get parameters
soiltype <- quefts_soil()
barley <- quefts_crop("Barley")
fertilizer <- list(N=0, P=0, K=0)
att_yield <- list(leaf_att=2200, stem_att=2700, store_att=4800, SeasonLength=110)

# 2. create a model
q <- quefts(soiltype, barley, fertilizer, att_yield)

# 3. run the model
run(q)

# change some parameters
q$SeasonLength <- 162
q$leaf_att <- 2651
q$stem_att <- 5053
q$store_att <- 8208

q$N <- 100
q$P <- 50
q$K <- 50

run(q)


## note that Rquefts uses C++ reference classes. 
## This means that if you copy a quefts model, you do not create a 
## new instance of the model, but you point to the same one!
q <- quefts()
q["N"]
k <- q
k["N"] <- 150
k["N"]
# the value of q has also changed!
q["N"]

## different ways of subsetting / replacement
q <- quefts()
q$N
q$N <- 30
q["N"]
q["N"] <- 90
q["model", "N"]
q["model", "N"] <- 60
q$N

q$soil$N_recovery
q["soil$N_recovery"]
q["soil$N_recovery"] <- .6
q["soil", "N_recovery"]
q["soil", "N_recovery"] <- .4
q$soil$N_recovery

}
