#' Wrap values in an HTML tag
#'
#' @param x a vector of values to be wrapped in a tag
#' @param tag A character vector of length 1, specifying the html tag (e.g., "div", "h1", "span")
#' @param attr_str A character string of the same length as x (or of length 1).
#' @param ignore_na If TRUE, do not add tag if value is NA
#' @param span_adjacent If TRUE, include adjacent tokens with identical attr_str within the same tag
#' @param doc_id  If span_adjacent is TRUE, The document ids are required to ensure that tags do not span from one document to another.
#'
#' @export
#'
#' @return a character vector
#' @examples
#' x = c("Obama","Bush")
#' add_tag(x, 'span')
#'
#' ## add attributes with the tag_attr function
#' add_tag(x, 'span',
#'         tag_attr(class = "president"))
#'
#' ## add style attributes with the attr_style function within tag_attr
#' add_tag(x, 'span',
#'         tag_attr(class = "president",
#'                  style = attr_style(`background-color` = 'rgba(255, 255, 0, 1)')))
add_tag <- function(x, tag, attr_str=NULL, ignore_na=F, span_adjacent=F, doc_id=NULL) {
  if (!is.null(attr_str)) {
    attr_str = ifelse(is.na(attr_str), yes = if (ignore_na) NA else '',
                                       no = stringi::stri_paste(' ', attr_str, sep=''))
  } else attr_str = if (ignore_na) NA else ''


  if (span_adjacent) {

    if (is.null(doc_id)) stop('If span_adjacent is used, document ids must be given')
    same_as_next = str_is_lead(attr_str) & str_is_lead(doc_id)
    same_as_prev = str_is_lag(attr_str) & str_is_lag(doc_id)

    x = ifelse(!same_as_prev, stringi::stri_paste('<',tag, attr_str,'>',  x, sep=''), x)
    x = ifelse(!same_as_next, stringi::stri_paste(x, '</',tag,'>', sep=''), x)
  } else x = stringi::stri_paste('<',tag, attr_str,'>',  x,  '</',tag,'>',  sep='')
  x
}

str_is_lag <- function(x) {
  x = as.character(x)
  c(F, x[-length(x)] == x[-1])
}

str_is_lead <- function(x) {
  x = as.character(x)
  c(x[-1] == x[-length(x)], F)
}


#' HTML tables for meta data per document
#'
#' Each row of the data.frame is transformed into a html table with two columns: name and value.
#' The columnnames of meta are used as names.
#'
#' @param meta a data.frame where each row represents the meta data for a document
#' @param ignore_col optionally, a character vector with names of metadata columns to ignore
#'
#' @return a character vector where each value contains a string for an html table.
#' @export
#' @examples
#' tabs = create_meta_tables(sotu_data$meta)
#' tabs[1]
create_meta_tables <- function(meta, ignore_col=NULL) {
  if (ncol(meta) > 0) {
    html_table = ''
    for (col in colnames(meta)) {
      if (col == 'NAVIGATION_STRING') next
      if (col %in% ignore_col) next
      colval = as.character(meta[[col]])
      colval[is.na(colval)] = ''
      html_table = stringi::stri_paste('\n', html_table, '<tr><th>', col, '</th><td>', colval, '</td></tr>')
    }
    add_tag(html_table, 'table', tag_attr(class='meta_table', style=attr_style(`border-spacing`= '0px')))
  } else ''
}

