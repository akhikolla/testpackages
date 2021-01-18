

#' Model Response Styles in Partial Credit Models
#' 
#' Performs PCMRS, a method to model response styles in Partial Credit Models
#' 
#' 
#' @name PCMRS-package
#' @docType package
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}\cr
#' \url{https://www.sg.tum.de/epidemiologie/team/schauberger/}
#' @seealso \code{\link{PCMRS}}, \code{\link{person.posterior}}, \code{\link{tenseness}}, \code{\link{emotion}}
#' @references Tutz, Gerhard, Schauberger, Gunther and Berger, Moritz (2018): 
#' Response Styles in the Partial Credit Model, \emph{Applied Psychological Measurement}, \url{https://journals.sagepub.com/doi/10.1177/0146621617748322}
#' @keywords package Partial Credit Response Style PCMRS
#' @examples
#' \dontshow{
#' k <- 4; n <- 80; I <- 4
#' set.seed(1860)
#' Y <- as.data.frame(matrix(sample(1:k, I*n, TRUE),nrow = n))
#' Y <- data.frame(lapply(Y, as.ordered))
#' 
#' mini.ex <- PCMRS(Y, cores = 2)
#' mini.ex
#' }
#' \dontrun{
#' ################################################
#' ## Small example to illustrate model and person estimation
#' ################################################
#' 
#' data(tenseness)
#' 
#' set.seed(5)
#' samples <- sample(1:nrow(tenseness), 100)
#' tense_small <- tenseness[samples,1:4]
#' 

#' m_small <- PCMRS(tense_small, cores = 2)
#' m_small
#' plot(m_small)
#' 
#' persons <- person.posterior(m_small, cores = 2)
#' plot(jitter(persons, 100))
#' 
#' ################################################
#' ## Example from Tutz et al. 2017:
#' ################################################
#' 
#' data(emotion)
#' m.emotion <- PCMRS(emotion)
#' m.emotion
#' 
#' plot(m.emotion)
#' }
NULL


#' Tenseness data from the Freiburg Complaint Checklist (tenseness)
#' 
#' Data from the Freiburg Complaint Checklist. 
#' The data contain all 8 items corresponding to the scale \emph{Tenseness} for 2042 participants of the 
#' standardization sample of the Freiburg Complaint Checklist. 
#' 
#' @name tenseness
#' @docType data
#' @format A data frame containing data from the Freiburg Complaint Checklist with 2042 observations. 
#' All items refer to the scale \emph{Tenseness} and are measured on a 5-point Likert scale where low numbers 
#' correspond to low frequencies or low intensitites of the respective complaint and vice versa. 
#' \describe{ 
#' \item{Clammy hands}{Do you have clammy hands?}
#' \item{Sweat attacks}{Do you have sudden attacks of sweating?}
#' \item{Clumsiness}{Do you notice that you behave clumsy?}
#' \item{Wavering hands}{Are your hands wavering frequently, e.g. when lightning a cigarette or when holding a cup?}
#' \item{Restless hands}{Do you notice that your hands are restless?}
#' \item{Restless feet}{Do you notice that your feet are restless?}
#' \item{Twitching eyes}{Do you notice unvoluntary twitching of your eyes?}
#' \item{Twitching mouth}{Do you notice unvoluntary twitching of your mouth?}
#'  }
#' @references Tutz, Gerhard, Schauberger, Gunther and Berger, Moritz (2018): 
#' Response Styles in the Partial Credit Model, \emph{Applied Psychological Measurement}, \url{https://journals.sagepub.com/doi/10.1177/0146621617748322}
#' @source 
#' ZPID (2013). PsychData of the Leibniz Institute for Psychology Information ZPID. Trier: Center for Research Data in Psychology.
#' 
#' Fahrenberg, J. (2010). Freiburg Complaint Checklist [Freiburger Beschwerdenliste (FBL)]. Goettingen, Hogrefe.
#' @keywords datasets
#' @examples
#' \dontrun{
#' data(tenseness)
#' 
#' set.seed(1860)
#' samples <- sample(1:nrow(tenseness), 300)
#' tense_small <- tenseness[samples,]
#' 
#' m_small <- PCMRS(tense_small, cores = 25)
#' m_small
#' plot(m_small)
#' 
#' persons <- person.posterior(m_small, cores = 25)
#' plot(jitter(persons,100))
#' }
NULL



#' Emotional reactivity data from the Freiburg Complaint Checklist (emotion)
#' 
#' Data from the Freiburg Complaint Checklist. 
#' The data contain all 8 items corresponding to the scale \emph{Emotional reactivity} for 2032 participants of the 
#' standardization sample of the Freiburg Complaint Checklist. 
#' 
#' @name emotion
#' @docType data
#' @format A data frame containing data from the Freiburg Complaint Checklist with 2032 observations. 
#' All items refer to the scale \emph{Emotional reactivity} and are measured on a 5-point Likert scale where low numbers 
#' correspond to low frequencies or low intensitites of the respective complaint and vice versa. 
#' \describe{ 
#' \item{Feel upset in whole body}{Do you feel it in the whole body when you get upset about something?}
#' \item{Eyes well up with tears}{Do your eyes well up with tears in certain situations?}
#' \item{Stammer}{Do you sometimes start stammering in certain situations?}
#' \item{Blush}{Do you blush?}
#' \item{Gasp for air}{Do you have to gasp for air in exciting situations, so that you have to take a deep breath?}
#' \item{Rapid heartbeat in excitement}{Do you feel a rapid heartbeat in excitement?}
#' \item{Urge to defecate in excitement}{Do you feel the urge to defecate in excitement?}
#' \item{Trembling knees}{Do you start trembling in excitement or do you get trembling knees?}
#'  }
#' @references Tutz, Gerhard, Schauberger, Gunther and Berger, Moritz (2018): 
#' Response Styles in the Partial Credit Model, \emph{Applied Psychological Measurement}, \url{https://journals.sagepub.com/doi/10.1177/0146621617748322}
#' @source 
#' ZPID (2013). PsychData of the Leibniz Institute for Psychology Information ZPID. Trier: Center for Research Data in Psychology.
#' 
#' Fahrenberg, J. (2010). Freiburg Complaint Checklist [Freiburger Beschwerdenliste (FBL)]. Goettingen, Hogrefe.
#' @keywords datasets
#' @examples
#' \dontrun{
#' data(emotion)
#' m.emotion <- PCMRS(emotion)
#' m.emotion
#' 
#' plot(m.emotion)
#' }
NULL
