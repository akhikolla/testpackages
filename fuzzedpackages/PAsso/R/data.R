#' The 2016 ANES Time Series Study with pre-election interview
#'
#' A subset of 2,188 participants of the 2016 American National Election Time Series Study,
#' which was to track the enduring social trend and record the political moment of 2016 (DeBell, 2018).
#' This study consisted of two surveys with same population. The pre-election interview
#' was during the weeks before the 2016 general election, including 4,271 respondents in
#' total. The post-election interview is the re-interview during the weeks after the
#' election, including 3,649 respondents (662 respondents did not complete post-interviews).
#'
#' The Pre-election preference is recorded as "PreVote" and the "PreVote.num"
#' is the numeric of it. Observations with missing values, or "No thought"
#' responses have been removed. Respondents expressing a voting preference
#' other than Clinton or Trump have been removed.
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 2188 rows and 10 variables.
#' \itemize{
#'   \item \code{age} Respondent's age in years.
#'
#'   \item \code{edu.year} Respondent's education year, which is mapped from \code{education}
#'   level: \code{MS}=8, \code{HSdrop}=11, \code{HS}=12, \code{Coll}=14,
#'   \code{CCdeg}=15, \code{BAdeg}=17, \code{MAdeg}=19.
#'
#'   \item \code{education} Respondent's education level.
#'
#'   \item \code{income.num} Respondent's family income in thousands: an numerical variable.
#'   It is median value of the range of each \code{income} level
#'
#'   \item \code{income} Respondent's family income level:
#'
#'   '(01) 01. Under $5,000' = 5,
#'
#'   '(02) 02. $5,000-$9,999',
#'
#'   '(03) 03. $10,000-$12,499',
#'
#'   '(04) 04. $12,500-$14,999',
#'
#'   '(05) 05. $15,000-$17,499',
#'
#'   '(06) 06. $17,500-$19,999',
#'
#'   '(07) 07. $20,000-$22,499',
#'
#'   '(08) 08. $22,500-$24,999',
#'
#'   '(09) 09. $25,000-$27,499',
#'
#'   '(10) 10. $27,500-$29,999',
#'
#'   '(11) 11. $30,000-$34,999',
#'
#'   '(12) 12. $35,000-$39,999',
#'
#'   '(13) 13. $40,000-$44,999',
#'
#'   '(14) 14. $45,000-$49,999',
#'
#'   '(15) 15. $50,000-$54,999',
#'
#'   '(16) 16. $55,000-$59,999',
#'
#'   '(17) 17. $60,000-$64,999',
#'
#'   '(18) 18. $65,000-$69,999',
#'
#'   '(19) 19. $70,000-$74,999',
#'
#'   '(20) 20. $75,000-$79,999',
#'
#'   '(21) 21. $80,000-$89,999',
#'
#'   '(22) 22. $90,000-$99,999',
#'
#'   '(23) 23. $100,000-$109,999',
#'
#'   '(24) 24. $110,000-$124,999',
#'
#'   '(25) 25. $125,000-$149,999',
#'
#'   '(26) 26. $150,000-$174,999',
#'
#'   '(27) 27. $175,000-$249,999',
#'
#'   '(28) 28. $250,000 or more' = 250.
#'
#'   \item \code{PID} Party identification: a numeric variable with value from 1 to 7
#'   representing strong Democrat, \code{strDem} < weak Democrat, \code{weakDem} <
#'   independent Democrat, \code{indDem} < independent independent \code{indind} <
#'   independent Republican, \code{indRep} < weak Republican, \code{weakRep} <
#'   strong Republican, \code{strRep}.

#'   \item \code{selfLR} The respondent' self-placement about own left-right in 7 ordinal levels
#'   (from extremely liberal to extremely conservative). \code{extLib}: extremely liberal, \code{Lib}: liberal,
#'   \code{sliLib}: slightly liberal, \code{Mod}: moderate, \code{sliCon}:
#'   slightly conservative, \code{Con}: conservative, \code{extCon}: extremely
#'   conservative. \code{extLib} < \code{Lib} < \code{sliLib} <
#'   \code{Mod} < \code{sliCon} < \code{Con} < \code{extCon}.
#'
#'   \item \code{TrumpLR} The respondent's opinion about Donald Trump's left-right
#'   placement (same scale as selfLR).
#'
#'   \item \code{ClinLR} The respondent's opinion about Hilary Clinton's left-right
#'   placement (same scale as selfLR).
#'
#'   \item \code{PreVote} The respondent's voting preference between
#'   Donald Trump and Hilary Clinton two months preceeding the
#'   November election (Pre-election interview). It is a factor with levels
#'   \code{HillaryClinton} and \code{DonaldTrump}.
#'
#'   \item \code{PreVote.num} Recode the PreVote to numeric values,
#'   'HillaryClinton'=0, 'DonaldTrump'=1.
#'
#'   \item \code{WeightforPreVote} Pre-election weight of a respondent.
#'
#' }
#'
#' @references
#' DeBell, Matthew, Jon A. Krosnick, Katie Gera, David S. Yeager, and Michael P. McDonald.
#' The turnout gap in surveys: Explanations and solutions. \emph{Sociological Methods & Research}, 2018:
#' 0049124118769085.
#'
#' Enamorado, Ted, Benjamin Fifield, and Kosuke Imai. \emph{User’s guide and codebook for the ANES
#' 2016 time series voter validation supplemental data}. Tech. rep.,
#' American National Election Studies, 2018.
#'
#' @name ANES2016
#'
#' @usage
#' data(ANES2016)
#'
#' @examples
#' head(ANES2016)
NULL

#' Simulated quadratic data
#'
#' Data simulated from a probit model with a quadratic trend. The data are
#' described in Example 2 of Liu and Zhang (2017).
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 2000 rows and 2 variables.
#' \itemize{
#'   \item \code{x} The predictor variable.
#'   \item \code{y} The response variable; an ordered factor.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach. \emph{Journal of the American Statistical Association}
#'
#' @name df1
#'
#' @usage
#' data(df1)
#'
#' @examples
#' head(df1)
#' @keywords internal
#'
NULL
