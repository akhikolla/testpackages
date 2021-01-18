#' Mouse-tracking experiment of a memory task
#' 
#' @description This dataset contains a subset of data originally presented in Coco & Duran (2016). In this task 
#' participants see sentence and scene pairs that varied in plausibility and are requested 
#' to classify the pairs as congruent or incongruent. The experimental variables are \emph{congruency} with two categorical levels (i.e., congruent, incongruent) and
#' \emph{plausibility} with two categorical levels (i.e., plausible, implausible). Participants have to classify each stimulus as belonging to one of these four levels. 
#' 
#' The dataset contains two participants (I=2), each measured along three trials, two categorical variables (Q=2) each with two levels (K=2). The total number of trials is J=12.
#' Mouse-tracking trajectories are raw-data, i.e. they have not been previously pre-processed.
#' 
#' @format A long-format dataframe of 728 observations containing 
#'   information on the following variables.
#' \describe{
#'  \item{sbj}{The ID number of participants}
#'  \item{trial}{The ID number of trials}
#'  \item{congruency}{A factor of levels \code{congruent}, \code{incongruent}}
#'  \item{plausibility}{A factor of levels \code{plausible}, \code{implausible}}
#'  \item{timestep}{The ID number of the recorded x-y trajectories} 
#'  \item{x}{The recorded x-trajectories}
#'  \item{y}{The recorded y-trajectories}
#' }
#' 
#' @source Coco, M. I., & Duran, N. D. (2016). When expectancies collide: 
#' Action dynamics reveal the interaction between stimulus plausibility and congruency. 
#' \emph{Psychonomic bulletin & review}, 23(6), 1920-1931.
"congruency"

#' Mouse-tracking experiment of a lexical decision task
#' 
#' @description This dataset contains a subset of data originally presented in Barca & Pezzullo (2012). In this task 
#' participants see a printed stimulus on the screen (e.g., water) and are requested 
#' to perform a dichotomous choice task where the stimulus can be classified as word or non-word. The experimental variable is the \emph{stimulus type} with four categorical levels (i.e., high-frequency word, low-frequency word, pseudowords, and strings of letters). 
#' Participants have to classify each stimulus as belonging to word or non-word categories. 
#' 
#' The dataset contains five participants (I=5), each measured along three trials, one categorical variable (Q=1) with four levels (K=4). The total number of trials is J=12. 
#' Mouse-tracking trajectories have previously been pre-processed with N=101 timesteps, translated into the first quadrant,
#' and rotated so that the Target point (y_T) is always on the right-side. 
#' 
#' 
#' @format A long-format dataframe of 6060 observations containing 
#'   information on the following variables.
#' \describe{
#'  \item{sbj}{The ID number of participants}
#'  \item{condition}{A factor of levels \code{HF}, \code{LF}, \code{PW}, \code{NW} indicating the type of stimulus}
#'  \item{timestep}{The ID number of the recorded x-y trajectories} 
#'  \item{x}{The recorded x-trajectories}
#'  \item{y}{The recorded y-trajectories}
#'  \item{trial}{The ID number of trials}
#' }

#' 
#' @source Barca, L., & Pezzulo, G. (2012). 
#' Unfolding visual lexical decision in time. 
#' \emph{PloS one}, 7(4), e35932.
"language"