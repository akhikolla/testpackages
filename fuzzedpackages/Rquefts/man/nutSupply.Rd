\name{nutSupply}

\alias{nutSupply1}
\alias{nutSupply2}

\title{
Soil nutrients supply for QUEFTS model
}


\description{
\code{nutSupply1} computes the base (unfertilized) soil supply of N, P and K according to Janssen et al. (1990), Table 2. For use with the QUEFTS model.

\code{nutSupply2} is a modified version following Sattari et al. (2014). It has an additional variable "temperature", and P-total is required. Sattari et al suggest that, for soils that have not been fertilized with P, you can estimate P-total as \code{95 * P-Olsen}. Emperically, using soilGrids, I found \code{55 * P-Olsen}.
 
}

\usage{
nutSupply1(pH, SOC, Kex, Polsen, Ptotal=NA)
nutSupply2(temp, pH, SOC, Kex, Polsen, Ptotal)
}

\arguments{
  \item{temp}{average growing season temperature (C)}
  \item{pH}{soil pH (H2O)}
  \item{SOC}{soil organic carbon (g/kg)}
  \item{Kex}{exchangeble K in the soil (mmol/kg)}
  \item{Polsen}{soil P measured with the P-Olsen method (mg/kg)}  
  \item{Ptotal}{total soil P (mg/kg)}  
}

\value{
Matrix with three columns: Nsup, Psup and Ksup. These are the potential supply of N, P and K of the unfertilized soil (kg/ha).
}

\references{
Janssen B.H., F.C.T. Guiking, D. van der Eijk, E.M.A. Smaling, J. Wolf and H. van Reuler, 1990. A system for the quantitative evaluation of the fertility of tropical soils (QUEFTS). Geoderma 46: 299-318

Sattari, S.Z., M.K. van Ittersum, A.F. Bouwman, A.L. Smit, and B.H. Janssen, 2014. Crop yield response to soil fertility and N, P, K inputs in different environments: Testing and improving the QUEFTS model. Field Crops Research 157: 35-46
}


\examples{
s1 <- nutSupply1(6, c(23, 11, 35), 15, c(1.6, 2.6, 2.4))
s1
s2 <- nutSupply2(20, 6, c(23, 11, 35), 15, c(1.6, 2.6, 2.4), 225)
s2
}

