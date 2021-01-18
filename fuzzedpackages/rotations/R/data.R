#' Drill data set
#'
#' The \code{drill} data set was collected to assess variation in human movement
#' while performing a task (Rancourt, 1995). Eight subjects drilled into a metal
#' plate while being monitored by infared cameras. Quaternions are used to
#' represent the orientation of each subjects' wrist, elbow and shoulder in one
#' of six positions. For some subjects several replicates are available. See
#' Rancourt et al. (2000) for one approach to analyzing these data.
#'
#' @format A data frame with 720 observations on the following 8 variables:
#' \describe{
#'   \item{\code{Subject}}{Subject number (1-8)}
#'   \item{\code{Joint}}{Joint name (Wrist, elbow, shoulder)}
#'   \item{\code{Position}}{Drilling position (1-6)}
#'   \item{\code{Replicate}}{Replicate number (1-5)}
#'   \item{\code{Q1}}{First element of orientation (quaternion)}
#'   \item{\code{Q2}}{Second element of orientation (quaternion)}
#'   \item{\code{Q3}}{Third element of orientation (quaternion)}
#'   \item{\code{Q4}}{Fourth element of orientation (quaternion)}
#' }
#'
#' @source
#'   \url{https://www.mat.ulaval.ca/lrivest/louis-paul-rivest/publications/}
#'
#' @references \enumerate{
#'   \item Rancourt, D. (1995). "Arm posture and hand
#'   mechanical impedance in the control of a hand-held power drill." Ph.D.
#'   Thesis, MIT.
#'   \item Rancourt, D., Rivest, L. & Asselin, J. (2000). "Using
#'   orientation statistics to investigate variations in human kinematics."
#'   Journal of the Royal Statistical Society: Series C (Applied Statistics),
#'   49(1), pp. 81-94.
#' }
#'
#' @examples
#' # Estimate central orientation of the first subject's wrist
#' Subject1Wrist <- subset(drill, Subject == 1 & Joint == "Wrist")
#' Qs <- as.Q4(Subject1Wrist[, 5:8])
#' mean(Qs)
#'
#' \donttest{
#'   # Plot Subject 1's wrist measurements using the connection to rotation matrices
#'   plot(Qs, col = c(1, 2, 3))
#' }
#'
#' # Translate the quaternion measurements into rotations and
#' # estimate the central orientation in terms of rotations
#' Rs <- as.SO3(Qs)
#' mean(Rs)
"drill"

#' Nickel electron backscatter diffraction data set
#'
#' This data set consists of electron backscatter diffraction (EBSD) data
#' obtained by scanning a fixed 12.5 \eqn{\mu}m-by-10 \eqn{\mu}m nickel surface
#' at individual locations spaced 0.2 \eqn{\mu}m apart. This scan was repeated
#' 14 times for each of the 3,449 locations yielding a total of 48,286
#' observations. Every observation corresponds to the orientation, expressed as
#' a rotation matrix, of a cubic crystal on the metal surface at a particular
#' location. Be aware that there are missing values and erroneous scans at some
#' locations and scans. See Bingham et al. (2009) and Bingham et al. (2010) for
#' more details and analysis.
#'
#' @format A data frame with 48,286 rows and the following 13 columns:
#' \describe{
#'   \item{\code{xpos}}{location x position}
#'   \item{\code{ypos}}{location y position}
#'   \item{\code{location}}{Location number for easy reference}
#'   \item{\code{rep}}{Replicate scan identifier}
#'   \item{\code{V1}}{First element of x-axis describing crystal orientation at corresponding location}
#'   \item{\code{V2}}{Second element of x-axis describing crystal orientation at corresponding location}
#'   \item{\code{V3}}{Third element of x-axis describing crystal orientation at corresponding location}
#'   \item{\code{V4}}{First element of y-axis describing crystal orientation at corresponding location}
#'   \item{\code{V5}}{Second element of y-axis describing crystal orientation at corresponding location}
#'   \item{\code{V6}}{Third element of y-axis describing crystal orientation at corresponding location}
#'   \item{\code{V7}}{First element of z-axis describing crystal orientation at corresponding location}
#'   \item{\code{V8}}{Second element of z-axis describing crystal orientation at corresponding location}
#'   \item{\code{V9}}{Third element of z-axis describing crystal orientation at corresponding location}
#' }
#'
#' @source The data set was collected by the Ames Lab located in Ames, IA.
#'
#' @references \enumerate{
#'   \item Bingham, M. A., Nordman, D., & Vardeman, S. (2009). "Modeling and
#'   inference for measured crystal orientations and a tractable class of
#'   symmetric distributions for rotations in three dimensions." Journal of the
#'   American Statistical Association, 104(488), pp. 1385-1397.
#'   \item Bingham, M. A., Lograsso, B. K., & Laabs, F. C. (2010). "A
#'   statistical analysis of the variation in measured crystal orientations
#'   obtained through electron backscatter diffraction." Ultramicroscopy,
#'   110(10), pp. 1312-1319.
#'   \item Stanfill, B., Genschel, U., & Heike, H. (2013). "Point estimation of
#'   the central orientation of random rotations". Technometrics, 55(4), pp.
#'   524-535.
#' }
#'
#' @examples
#' # Subset the data to include only the first scan
#' Rep1 <- subset(nickel, rep == 1)
#'
#' # Get a rough idea of how the grain map looks by plotting the first
#' # element of the rotation matrix at each location
#' ggplot2::qplot(xpos, ypos, data = Rep1, colour = V1, size = I(2))
#'
#' # Focus in on a particular location, for example location 698
#' Rs <- subset(nickel, location == 698)
#'
#' # Translate the Rs data.frame into an object of class 'SO3'
#' Rs <- as.SO3(Rs[,5:13])
#'
#' # Some observations are not rotations, remove them
#' Rs <- Rs[is.SO3(Rs),]
#'
#' # Estimate the central orientation with the average
#' mean(Rs)
#'
#' # Re-estimate central orientation robustly
#' median(Rs)
#'
#' \donttest{
#'   # Visualize the location, there appears to be two groups
#'   plot(Rs, col = c(1, 2, 3))
#' }
"nickel"
