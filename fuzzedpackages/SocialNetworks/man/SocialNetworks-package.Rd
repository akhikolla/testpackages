\name{SocialNetworks}
\alias{SocialNetworks}
\docType{package}
\title{
Generates social networks 
}
\description{
Generates a social network based on user inputted spatial xy coordinates. 
}
\details{
\tabular{ll}{
Package: \tab SocialNetworks\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2014-08-20\cr
License: \tab GPL\cr
}

}
\author{
Glenna Nightingale

Peter Nightingale

Maintainer: Glenna Nightingale <glenna.evans@gmail.com>

}
\references{
Bradbury, J. W., and S.L. Vehrencamp (1998).  Principles of animal communication. Sinauer Associates.


Croft, D.P., James, R., and Krause, J. (2008).  Exploring animal social networks.  USA: Princeton University Press.


Funwi-Gabga, N., and Mateu, J. (2012).  Understanding the nesting spatial behaviour of gorillas in the Kagwene Sanctuary, Cameroon.  
Stochastic Environmental Research and Risk Assessment, vol 26, pp. 793-811.


Hoppitt, W. and Laland, K. N. (2013). Social Learning: An Introduction to Mechanisms, Methods, and Models.  Princeton University Press.

Illian, Moller, Waagepetersen, (2009). Hierarchical spatial
point process analysis for a plant community of high
biodiversity.Environ. Ecol. Stat. vol 16, pp 389-405



}

\examples{





#Using pairwise distances to calcualate inter-individual associations 
# generate a social network from a regular spatial point pattern ir = 0.06
#-------------------------------------------------------------------------

x = c(0.1023117, 0.1119260, 0.1625270, 0.3594291, 0.4220571, 0.4606205, 0.5927459,
 0.6847543, 0.7065195, 0.7760657, 0.9827536)

y = c(0.2525266, 0.3346728, 0.5275355, 0.2447207, 0.2765606, 0.4999600, 0.5928410,
 0.8356211, 0.2506116, 0.8994760, 0.1432255)

#plot(x,y)
irset = c(rep(0.06,11))
calculateassociations(x,y,irset)


# generate a social network from a clustered spatial point pattern
#----------------------------------------------------------------

 x =
  c(0.77302412, 0.82946034, 0.65776305, 0.62294479, 0.58577335, 0.39332654,
  0.36893684, 0.40518735, 0.53956642, 0.56596859, 0.62802969, 0.10380876,
 0.71058751, 0.65943692, 0.88056259, 0.90567566, 0.91166684, 0.89489341,
 0.92668619, 0.01544599, 0.30499431, 0.28249059, 0.30733518, 0.73165075,
 0.17712420, 0.80869511, 0.77351717, 0.75508022, 0.79445346, 0.73134413,
 0.62448310, 0.60180882, 0.66741081, 0.45884352, 0.45282315, 0.45614636,
 0.45270694, 0.44764728, 0.53259346)


 y=
 c(0.943378357, 0.933698623, 0.123641160, 0.146773076, 0.135097659, 0.978760171,
 0.981407654, 0.937111187, 0.080617391, 0.114438404, 0.061834776, 0.370322731,
 0.036576942, 0.003974257, 0.830356964, 0.837171526, 0.884801445, 0.797794654,
 0.844312417, 0.969982888, 0.672246284, 0.692111852, 0.671098280, 0.999097233,
 0.003736065, 0.255322335, 0.282689074, 0.310793806, 0.229047375, 0.266413304,
 0.324984514, 0.279652338, 0.287134158, 0.331962948, 0.365469720, 0.343868765,
0.378876999, 0.331915785, 0.368805652)


#plot(x,y)
irset = c(rep(0.05,length(x)))
calculateassociations(x,y,irset)



# generate a social network from a random spatial point pattern
#----------------------------------------------------------------

  x =
 c( 0.74905296, 0.38309725, 0.98627509, 0.02242039, 0.54703348, 0.59173730,
 0.82340399, 0.18718650, 0.49200511, 0.86098261, 0.24848640, 0.15843825,
0.72875205 )
y  =
 c(0.73521480, 0.01661629, 0.51564570, 0.61856835, 0.20815448, 0.29431260,
 0.35507188, 0.18940107, 0.98721494, 0.98129752, 0.76510267, 0.43541222,
0.04601392)


#plot(x,y)
irset = c(rep(0.1,length(x)))
calculateassociations(x,y,irset)


# Using the area of overlap between territorial zones to calcualate associations
# generate a social network for four individuals, where an interaction radius is 
# specified for each individual
# Note that the interaction radius for a group of individuals can be identical
#-------------------------------------------------------------------------------

calculate.areas(c(0.1, 0.2, 0.3,0.4), c(0.1, 0.2, 0.3, 0.4), c(0.1, 0.2, 0.3, 0.4),1000000) 
# create a social network for four individuals, with a separate interaction radius 
#for each individual

}
