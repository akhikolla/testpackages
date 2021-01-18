#' Transpose a color into the string format used in html attributes
#'
#' @param col The name of the color
#' @param alpha Optionally, the alpha (transparency), with 0 being fully
#'              transparent and 1 being fully colorized.
#'
#' @return The string used to specify a color in an html tag attribute
#' @export
#'
#' @examples
#' set_col('red')
#' set_col('red', alpha=0.5)
set_col <- function(col, alpha=1) {
  ucol = unique(col)
  rgb_col = grDevices::col2rgb(ucol)
  rgb_substring = sprintf("%s, %s, %s", rgb_col[1,], rgb_col[2,], rgb_col[3,])
  rgb_substring = rgb_substring[match(col, ucol)]
  sprintf("rgba(%s, %s)", rgb_substring, alpha)
}

#' Create a highlight color for a html style attribute
#'
#' Designed to be used together with the attr_style function.
#' The return value can directly be used to set the color in
#' an html tag attribute (e.g., color, background-color)
#'
#' @param value Either a logical vector or a numeric vector with
#'              values between 0 and 1. If a logical vector is used, then tokens
#'              with TRUE will be highlighted (with the color specified in pos_col).
#'              If a numeric vector is used, the value determines the alpha (transparency),
#'              with 0 being fully transparent and 1 being fully colored.
#' @param col The color used to highlight
#'
#' @return The string used to specify a color in an html tag attribute
#' @export
#'
#' @examples
#' highlight_col(c(NA, 0, 0.1,0.5, 1))
#'
#' ## used in combination with attr_style()
#' attr_style(color = highlight_col(c(NA, 0, 0.1,0.5, 1)))
#'
#' ## note that for background-color you need inversed quotes to deal
#' ## with the hyphen in an argument name
#' attr_style(`background-color` = highlight_col(c(NA, 0, 0.1,0.5, 1)))
#'
#' tag_attr(class = c(1, 2),
#'          style = attr_style(`background-color` = highlight_col(c(FALSE,TRUE))))
highlight_col <- function(value, col='yellow') {
  value = as.numeric(value)
  if (any(value < 0 | value > 1, na.rm = T)) stop('highlight value has to be logical (TRUE/FALSE) or a number between 0 and 1')

  ifelse(value > 0 & !is.na(value),  ## also consider zero value as NA, because it would be fully transparent anyway (alpha = 0)
         yes = set_col(col, value),
         no  = NA)
}

#' Create a scale color for a html style attribute
#'
#' Designed to be used together with the attr_style function. The return value
#' can directly be used to set the color in an html tag attribute (e.g., color, background-color)
#'
#' @param value A numeric vector with values between -1 and 1. Determines the color mixture
#'              of the scale colors specified in col_range
#' @param alpha Optionally, the alpha (transparency) can be specified, with 0 being fully
#'              transparent and 1 being fully colored. This can be a vector to specify a
#'              different alpha for each value.
#' @param col_range The colors used in the scale.
#'
#' @return The string used to specify a color in a html tag attribute
#' @export
#'
#' @examples
#' scale_col(c(NA, -1, 0, 0.5, 1))
#'
#' ## used in combination with attr_style()
#' attr_style(color = scale_col(c(NA, -1, 0, 0.5, 1)))
#'
#' ## note that for background-color you need inversed
#' ## quotes to deal with the hyphen in an argument name
#' attr_style(`background-color` = scale_col(c(NA, -1, 0, 0.5, 1)))
#'
#' tag_attr(class = c(1, 2),
#'          style = attr_style(`background-color` = scale_col(c(-1,1))))
scale_col <- function(value, alpha=1, col_range=c('red','blue')) {
  if (any(abs(value) > 1, na.rm = T)) stop('scale value has to be a number between -1 and 1')
  if (any(alpha < 0 | alpha > 1, na.rm = T)) stop('alpha value has to be a number between 0 and 1')

  value = (value + 1) / 2 ## colscale wants values between 0 and 1
  ifelse(!is.na(value),
         yes = colscale_to_attr_rgb(value, col_range, alpha = alpha),
         no  = NA)
}

colscale_to_attr_rgb <- function(value, colors=c('red','blue'), alpha=1) {
  cramp = grDevices::colorRamp(colors)
  rgb_col = cramp(value)
  rgb_col = round(rgb_col)
  sprintf("rgba(%s, %s, %s, %s)", rgb_col[,1], rgb_col[,2], rgb_col[,3], alpha)
}

