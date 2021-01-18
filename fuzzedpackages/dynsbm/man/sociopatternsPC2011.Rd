\name{sociopatternsPC2011}
\docType{data}
\alias{sociopatternsPC2011}
\title{High school dynamic contact networks in the PC class in 2011 - Sociopatterns collaboration}
\description{
The dataset  consists in face-to-face encounters of high school students 
of a class from a French high school.
In this class called 'PC', interactions were recorded during  4 days (Tuesday to Friday) in Dec. 2011.
We kept only the 27 (out of 31) students that appear every day.
Interaction times were aggregated by days to form a sequence of 4 different networks.
These are undirected and weighted networks, the weight of an interaction between
two  individuals being  the  number of  interactions  between these  2
individuals divided by the number of time points for which at least two individuals interacted;
we discretize these  data   into 3  bins.

We sincerely acknowledge the SocioPatterns collaboration for providing
the dataset (http://www.sociopatterns.org)
}
\usage{data(sociopatternsPC2011)}
\format{An array of size 4x27x27.}
\source{http://www.sociopatterns.org/datasets/high-school-dynamic-contact-networks/}
\references{
Julie Fournet and Alain Barrat, \emph{Contact patterns among high school students}, PLoS ONE 9(9):e107878 (2014).

Catherine Matias and Vincent Miele, \emph{Statistical clustering of temporal networks through a dynamic stochastic block model}, Journal of the Royal Statistical Society: Series B (2016)
http://dx.doi.org/10.1111/rssb.12200
http://arxiv.org/abs/1506.07464

http://www.sociopatterns.org/
}
