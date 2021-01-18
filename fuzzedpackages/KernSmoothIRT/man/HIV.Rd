\name{HIV}
\alias{HIV}
\alias{HIVkey}
\docType{data}
\title{
The HIV Data
}
\description{
These data come from the Voluntary HIV Counseling and Testing Efficacy Study (2000a,b), conducted in 1995-1997, by the Center for AIDS Prevention Studies, University of California, San Francisco. 
This study was concerned with the effectiveness of HIV counseling and testing in reducing risk behavior for the sexual transmission of HIV. 
As part of this study, respondents were surveyed about their attitude toward condom use via a bank of 15 items. 
Respondents were asked how much they agreed with each of the statements on a 4-point response scale.  
}
\usage{
data(HIV)

}

\format{
	\describe{
	\item{\code{SITE}}{  the site the participant attended. "\code{Ken}" for Kenia, "\code{Tan}" for Tanzania and "\code{Tri}" for Trinidad.}
	\item{\code{GENDER}}{gender of participant}
	\item{\code{AGE}}{   age (at last birthday) of participant} 
	\item{\code{I1}}{		 Using condoms is a good way to protect you and your sex partner from diseases people can get from having sex.}	
	\item{\code{I2}}{ 	 Sex doesn't feel as good when you use a condom.}
	\item{\code{I3}}{ 	 It's embarrassing to buy condoms.}
	\item{\code{I4}}{ 	 Using condoms is a good way to prevent an unwanted pregnancy.}
	\item{\code{I5}}{ 	 It's embarrassing to put on a condom (put a condom on a man).}
	\item{\code{I6}}{ 	 Most of your friends think you should use condoms whenever you have sex, including with your spouse.}
	\item{\code{I7}}{ 	 Condoms often break or slip off during sex.}
	\item{\code{I8}}{ 	 If your sex partner wants to use a condom, I'd suspect he/she may be having sex with someone else.}
	\item{\code{I9}}{ 	 Most of your friends use condoms when they have sex with a new partner.}
	\item{\code{I10}}{	 It is easy to buy condoms in the area where I now live.}
	\item{\code{I11}}{	 Most of your friends think that condoms are uncomfortable.}
	\item{\code{I12}}{	 Most of your friends think you should always use condoms when having sex with a new person.}
	\item{\code{I13}}{ 	 It is easy to get condoms in my home town without having to pay for them.}
	\item{\code{I14}}{	 Condoms cost too much.}
	\item{\code{I15}}{	 Most people your age are using condoms these days.}
	}

}

\details{
%%  ~~ If necessary, more details than the __description__ above ~~
Given the "negative" wording of the items \code{I2}, \code{I3}, \code{I5}, \code{I7}, \code{I8}, \code{I11}, and \code{I14}, in order to uniform the data, their scores has been reversed.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Voluntary HIV Counseling and Testing Efficacy Study (2000a,b), conducted in 1995-1997, by the Center for AIDS Prevention Studies, University of California, San Francisco.

}


\references{
Bertoli-Barsotti, L., Muschitiello, C. and Punzo A (2010). 
"Item Analysis of a Selected Bank from the Voluntary HIV-1 Counseling and Testing Efficacy Study Group." 
Technical Report 1, Dipartimento di Matematica, Statistica, Informatica e Applicazioni (Lorenzo Mascheroni), 
Universita degli Studi di Bergamo.

The Voluntary HIV-1 Counseling and Testing Efficacy Study Group (2000a). 
Efficacy of voluntary HIV-1 counselling and testing in individual and couples in Kenya, Tanzania, and Trinidad: A randomised trial. 
\emph{Lancet}, \bold{356}, 103-112.

The Voluntary HIV-1 Counseling and Testing Efficacy Study Group (2000b). 
The Voluntary HIV-1 Counseling and Testing Efficacy Study: Design and Methods. 
\emph{AIDS and Behavior}, \bold{4}(1), 5-14. 
}

\keyword{datasets}