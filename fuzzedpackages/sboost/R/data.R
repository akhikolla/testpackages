#' Malware System Calls
#'
#' System call data for apps identified as malware and not malware.
#'
#' Experimental data generated in this research paper:
#'
#' M. Dimjašević, S. Atzeni, I. Ugrina, and Z. Rakamarić,
#' "Evaluation of Android Malware Detection Based on System Calls,"
#' in Proceedings of the International Workshop on Security and Privacy Analytics (IWSPA), 2016.
#'
#' Data used for kaggle competition: \url{https://www.kaggle.com/c/ml-fall2016-android-malware}
#'
#' @format A data frame with 7597 rows and 361 variables: \emph{outcomes} 1 if malware, 0 if not.
#' \emph{X1... X360} system calls.
#'
#' @source \url{https://zenodo.org/record/154737#.WtoA1IjwaUl}
"malware"



#' Mushroom Classification
#'
#' A classic machine learning data set describing hypothetical samples from the Agaricus and Lepiota family.
#'
#' Data gathered from:
#'
#' Mushroom records drawn from The Audubon Society Field Guide to North American Mushrooms (1981).
#' G. H. Lincoff (Pres.), New York: Alfred A. Knopf
#'
#' @format A data frame with 7597 rows and 361 variables:
#' \describe{
#'   \item{outcomes}{p=poisonous, e=edible}
#'   \item{cap_shape}{bell=b, conical=c, convex=x, flat=f, knobbed=k, sunken=s}
#'   \item{cap_surface}{fibrous=f, grooves=g, scaly=y, smooth=s}
#'   \item{cap_color}{brown=n, buff=b, cinnamon=c, gray=g, green=r, pink=p, purple=u, red=e, white=w, yellow=y}
#'   \item{bruises}{bruises=t, no=f}
#'   \item{odor}{almond=a, anise=l, creosote=c, fishy=y, foul=f, musty=m, none=n, pungent=p, spicy=s}
#'   \item{gill_attachment}{attached=a, descending=d, free=f, notched=n}
#'   \item{gill_spacing}{close=c, crowded=w, distant=d}
#'   \item{gill_size}{broad=b, narrow=n}
#'   \item{gill_color}{black=k, brown=n, buff=b, chocolate=h, gray=g, green=r, orange=o, pink=p, purple=u, red=e, white=w, yellow=y}
#'   \item{stalk_shape}{enlarging=e, tapering=t}
#'   \item{stalk_root}{bulbous=b, club=c, cup=u, equal=e, rhizomorphs=z, rooted=r, missing=?}
#'   \item{stalk_surface_above_ring}{fibrous=f, scaly=y, silky=k, smooth=s}
#'   \item{stalk_surface_below_ring}{fibrous=f, scaly=y, silky=k, smooth=s}
#'   \item{stalk_color_above_ring}{brown=n, buff=b, cinnamon=c, gray=g, orange=o, pink=p, red=e, white=w, yellow=y}
#'   \item{stalk_color_below_ring}{brown=n, buff=b, cinnamon=c, gray=g, orange=o, pink=p, red=e, white=w, yellow=y}
#'   \item{veil_type}{partial=p, universal=u}
#'   \item{veil_color}{brown=n, orange=o, white=w, yellow=y}
#'   \item{ring_number}{none=n, one=o, two=t}
#'   \item{ring_type}{cobwebby=c, evanescent=e, flaring=f, large=l, none=n, pendant=p, sheathing=s, zone=z}
#'   \item{spore_print_color}{black=k, brown=n, buff=b, chocolate=h, green=r, orange=o, purple=u, white=w, yellow=y}
#'   \item{population}{abundant=a, clustered=c, numerous=n, scattered=s, several=v, solitary=y}
#'   \item{habitat}{grasses=g, leaves=l, meadows=m, paths=p, urban=u, waste=w, woods=d}
#' }
#'
#' @source \url{https://archive.ics.uci.edu/ml/datasets/mushroom}
"mushrooms"
