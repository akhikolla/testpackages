#' Trial Matrix Data
#' 
#' This data set contains the subject's responses to items. 
#' Correct answers are denoted by 1 and incorrect answers are denoted by 0.
#' 
#' @format A data frame with 252 observations on the following 30 variables.
#' \describe{
#'   \item{\code{t1}}{Subject's Response to Item 1.}
#'   \item{\code{t2}}{Subject's Response to Item 2.}
#'   \item{\code{t3}}{Subject's Response to Item 3.}
#'   \item{\code{t4}}{Subject's Response to Item 4.}
#'   \item{\code{t5}}{Subject's Response to Item 5.}
#'   \item{\code{t6}}{Subject's Response to Item 6.}
#'   \item{\code{t7}}{Subject's Response to Item 7.}
#'   \item{\code{t8}}{Subject's Response to Item 8.}
#'   \item{\code{t9}}{Subject's Response to Item 9.}
#'   \item{\code{t10}}{Subject's Response to Item 10.}
#'   \item{\code{t11}}{Subject's Response to Item 11.}
#'   \item{\code{t12}}{Subject's Response to Item 12.}
#'   \item{\code{t13}}{Subject's Response to Item 13.}
#'   \item{\code{t14}}{Subject's Response to Item 14.}
#'   \item{\code{t15}}{Subject's Response to Item 15.}
#'   \item{\code{t16}}{Subject's Response to Item 16.}
#'   \item{\code{t17}}{Subject's Response to Item 17.}
#'   \item{\code{t18}}{Subject's Response to Item 18.}
#'   \item{\code{t19}}{Subject's Response to Item 19.}
#'   \item{\code{t20}}{Subject's Response to Item 20.}
#'   \item{\code{t21}}{Subject's Response to Item 21.}
#'   \item{\code{t22}}{Subject's Response to Item 22.}
#'   \item{\code{t23}}{Subject's Response to Item 23.}
#'   \item{\code{t24}}{Subject's Response to Item 24.}
#'   \item{\code{t25}}{Subject's Response to Item 25.}
#'   \item{\code{t26}}{Subject's Response to Item 26.}
#'   \item{\code{t27}}{Subject's Response to Item 27.}
#'   \item{\code{t28}}{Subject's Response to Item 28.}
#'   \item{\code{t29}}{Subject's Response to Item 29.}
#'   \item{\code{t30}}{Subject's Response to Item 30.}
#' }
#' 
#' @source 
#' Choice38 Experiment at UIUC during Spring 2014 - Fall 2014
#' 
#' @author 
#' Steven Andrew Culpepper and James Joseph Balamuta
"trial_matrix"


#' Payout Matrix Data
#' 
#' This data set contains the payout information for each subject. 
#' 
#' @format A data frame with 252 observations on the following 4 variables.
#' \describe{
#'   \item{\code{Participant}}{Subject ID}
#'   \item{\code{cum_sum}}{Sum of all payouts}
#'   \item{\code{num_correct_choices}}{Total number of correct choices (out of 15)}
#'   \item{\code{num_correct_trials}}{Total number of correct trials (out of 30)}
#' }
#' 
#' @source 
#' Choice38 Experiment at UIUC during Spring 2014 - Fall 2014
#' 
#' @author 
#' Steven Andrew Culpepper and James Joseph Balamuta
#' 
"payout_matrix"

#' Choice Matrix Data
#' 
#' This data set contains the subject's choices and point values for the 
#' difficult questions.
#' 
#' @format A data frame with 3780 observations on the following 5 variables.
#' \describe{
#'   \item{\code{subject_id}}{Research Participant Subject ID. There are 102 IDs and each ID has 15 observations.}
#'   \item{\code{hard_q_id}}{The item ID of the hard question assigned to the student (16-30)}
#'   \item{\code{easy_q_id}}{The item ID of the easy question assigned to the student (1-15)}
#'   \item{\code{choose_hard_q}}{Selected either: Difficult Question (1) or Easy Question (0)}
#'   \item{\code{high_value}}{Range of values associated with Difficult Question that span from 12 to 16, repeated three times per subject}
#'   \item{\code{low_value}}{Range of values associated with Easy Question that span from 4 to 6, repeated five times per subject}
#'   \item{\code{is_correct_choice}}{Did the user select an item that was answered correctly?}
#' }
#' 
#' @source 
#' Choice38 Experiment at UIUC during Spring 2014 - Fall 2014
#' 
#' @author 
#' Steven Andrew Culpepper and James Joseph Balamuta
"choice_matrix"

#' Survey Data
#' 
#' This data set contains the subject's responses survey questions administered
#' using Choice38.
#'  
#' @format A data frame with 102 observations on the following 2 variables.
#' \describe{
#'   \item{\code{id}}{Subject's Assigned Research ID}
#'   \item{\code{sex}}{Subject's sex: 
#'     \itemize{
#'        \item Male
#'        \item Female
#'     }
#'   }
#' }
#' 
#' @source 
#' Choice38 Experiment at UIUC during Spring 2014 - Fall 2014
#' 
#' @author 
#' Steven Andrew Culpepper and James Joseph Balamuta
"survey_data"
