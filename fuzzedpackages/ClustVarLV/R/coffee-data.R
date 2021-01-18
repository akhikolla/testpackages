#' coffee data
#'
#' Case study pertaining to consumer emotions associations for a
#' variety of 12 coffee aromas. The participants were asked to complete each rating (i.e., rating the odor of 12 aromas on 15 emotion terms) on a 5-point rating scale.
#'
#'
#' @docType data
#'
#' @usage data(coffee)
#'
#' @format An object of class \code{"array"} with 12 odors (mode 1), 84 subjects (mode 2) and 15 emotions (mode 3):
#' \describe{
#' \item{odors}{Vanilla, B.Rice, Lemon, Coffee.Flower, Cedar, Hazelnut,
#' Coriander.Seed, Honey, Medicine, Apricot, Earth, Hay}
#' \item{subjects}{persons from Oniris}
#' \item{emotions}{Amused, Angry, Calm, Disappointed, Disgusted,
#' Energetic, Excited, Free, Happy, Irritated,
#' Nostalgic, Surprised, Unique, Unpleasant and Well}
#' }
#'
#' @keywords datasets
#'
#' @references Cariou, V., & Wilderjans, T. F. (2018). Consumer segmentation in multi-attribute product evaluation by means of non-negatively constrained CLV3W. Food Quality and Preference, 67, 18-26.
#' (\href{https://doi.org/10.1016/j.foodqual.2017.01.006}{ScienceDirect})
#'
#' @examples
#' data(coffee)
#' str(coffee)
"coffee"
