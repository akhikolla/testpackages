\name{cognitive}
\alias{cognitive}
\docType{data}
\title{Kenya School Lunch Intervention Cognitive Dataset}
\description{
In the Kenya school lunch intervention, children were given one of four school lunch interventions: meat, milk, calorie, or control. The first three groups were fed a school lunch of a stew called githeri supplemented with either meat, milk, or oil to create a lunch with a given caloric level, while the control group did not receive a lunch. Three schools were randomized to each group and the lunch program is the same for all children within a school. The data is available in \cite{modeling-longitudinal-data-rob-weiss} and is broken up into sub data sets from four domains: Anthropometry, Cognitive, Morbidity, and Nutrition. We will be using the cognitive dataset for analyzing how the cognition level of the school children change over time and how the change is associated with other variables. The main cognitive measures is Raven's colored progressive matrices (Raven's), a measure of cognitive ability. There are three additional response variables: arithmetic score (arithmetic), verbal meaning (vmeaning), and total digit span score (dstotal) where digit span is a test of memory while others are considered measures of intelligence or education. The cognitive measurement baseline was taken prior to the lunch program onset and measurements were assessed at up to five times, called rounds, for each subject. More information about this dataset please see the reference:

Robert E Weiss.Modeling longitudinal data.  Springer Science & Business Media, 2005.

}
\usage{data(cognitive)}
\format{
  \describe{
   A data frame of 1562 observations and 26 variables.
   \item{id}{Grouping variable. Unique ID for each subject.}
   \item{schoolid}{School id 1-12.}
   \item{treatment}{Calorie, meat, milk, control}
   \item{rn}{round. }
   \item{year}{Time in years from baseline.}
   \item{revans}{Raven's colored matrices score.}
   \item{arithmetic}{Arithmetic score.}
   \item{vmeaning}{Verbal meaning.}
   \item{dstotal}{Total digit span score.}
   \item{sex}{Girl or Boy.}
   \item{age_at_time0}{age at baseline.}
   \item{height}{height at baseline. }
   \item{weight}{weight at baseline. }
   \item{head_circ}{Head circumference at baseline. }
   \item{ses}{Socio-Economic Status score. }
   \item{mom_read}{Mother's reading test. }
   \item{mom_write}{Mother's writing test. }
   \item{mom_edu}{Mother's years of educations. }
   \item{morbscore}{Morbidity score: none/mild/severe. }
   \item{complete}{Logical variable specifying whether the subject has all five rounds. 1-Yes, 0-No. }
   \item{rnone}{Logical variable specifying whether the observation is the baseline. 1-round one (baseline), 0-not round one. }
   \item{relmonth}{Time in months from baseline.}
  }
}


\examples{
data(cognitive)
}
\keyword{datasets}
