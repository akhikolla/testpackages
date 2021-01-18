\name{high_school}
\alias{high_school}
\docType{data}

\title{Interactions between high school students}

\description{
The data concern face to face interactions among 327 high school students in Marseilles,
France, and were collected by means of wearable sensors over a period of 5 days in December 2013.
Students wore a sensor badge on their chest and the instrument recorded when they were facing each
other with a time resolution of 20 seconds. Thus, any pair of students was considered interacting
face-to-face when the sensors of the two were exchanging data packets at any given time during the
20 seconds interval.
Additional information on the students is available from the same dataset. Students may have 4
different main specializations: biology (BIO), mathematics and physics (MP), physics and chemistry
(PC), and engineering studies (PSI).
}

\format{
The list contains:
\describe{
\item{adj}{An adjacency list indicating whether any pair of students had at least one interaction during the 5 days of the study.}
\item{edgelist}{An edgelist in a format that can be handled by this package.}
\item{program}{Clustering variable indicating the program each student is registered to.}
\item{program_levels}{Names of the different programs.}
\item{program_aggr}{Aggregated version of the previous clustering variable, where programs are aggregated into 4 areas.}
\item{program_levels_aggr}{Names of the different areas, these correspond to biology (BIO), mathematics and physics (MP), physics and chemistry (PC), and engineering studies (PSI)}
\item{sex}{Clustering given by the sex of the students.}
\item{sex_levels}{Labels for each of the sex classes.}
}
}

\usage{data(high_school)}

\references{
R. Mastrandrea, J. Fournet, and A. Barrat (2015). "Contact patterns in a high school: A comparison between data collected using wearable sensors, contact diaries and friendship surveys".
\emph{PLOS ONE 10.9, pp. 1-26}.
}

\keyword{datasets}
