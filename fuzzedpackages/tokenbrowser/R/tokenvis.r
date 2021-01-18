#' Convert tokens into full texts in an HTML file
#'
#' @param tokens    A data.frame with a column for document ids (doc_col)
#'                  and a column for tokens (token_col)
#' @param meta      A data.frame with a column for document_ids (doc_col). All other columns are added
#'                  to the browser as document meta
#' @param doc_col   The name of the document id column
#' @param token_col The name of the token column
#' @param space_col  Optionally, a column with space indications (" ", "\\n", etc.) per token (which is how some NLP parsers indicate spaces)
#' @param doc_nav   The name of a column (factor or character) in meta, used to create a navigation bar for selecting document groups.
#' @param token_nav Alternative to doc_nav, a column in the tokens. Navigation filters will then be used to select documents in which
#'                  the value occurs at least once.
#' @param filename  Name of the output file. Default is temp file
#' @param css_str   A character string, to be directly added to the css style header
#' @param header    Optionally, specify the header
#' @param subheader Optionally, specify a subheader
#' @param n         If TRUE, report N in header
#' @param navfilter If TRUE (default) enable filtering with nav(igation) bar.
#' @param top_nav   A number. If token_nav is used, navigation filters will only apply to the top x values with highest token occurence in a document
#' @param thres_nav Like top_nav, but specifying a threshold for the minimum number of tokens.
#' @param colors    Optionally, a vector with color names for the navigation bar. Length has to be identical to
#'                  unique non-NA items in the navigation.
#' @param style_col1 Color of the browser header
#' @param style_col2 Color of the browser background
#'
#' @return The name of the file where the browser is saved. Can be opened conveniently from within R using browseUrl()
#' @export
#'
#' @examples
#' url = create_browser(sotu_data$tokens, sotu_data$meta, token_col = 'token', header = 'Speeches')
#'
#' \donttest{
#' view_browser(url)   ## view browser in the Viewer
#' }
#' if (interactive()) {
#' browseURL(url)     ## view in default webbrowser
#' }
create_browser <- function(tokens, meta=NULL, doc_col='doc_id', token_col='token', space_col=NULL, doc_nav=NULL, token_nav=NULL, filename=NULL, css_str=NULL, header='', subheader='', n=TRUE, navfilter=TRUE, top_nav=NULL, thres_nav=1, colors=NULL, style_col1="#7D1935", style_col2="#F5F3EE"){
  tokens[[doc_col]] = factor(as.character(tokens[[doc_col]]), levels=unique(tokens[[doc_col]]))

  docs = wrap_documents(tokens, meta, doc_col, token_col, space_col, nav=doc_nav, token_nav = token_nav, top_nav=top_nav, thres_nav=thres_nav)
  docstring = stringi::stri_paste(docs, collapse='\n\n')

  doc_ids = unique(tokens[[doc_col]])
  n_doc = length(doc_ids)

  nav = NULL
  if (!is.null(doc_nav)) {
    radio = TRUE
    nav = if (methods::is(meta[[doc_nav]], 'factor')) levels(meta[[doc_nav]]) else unique(meta[[doc_nav]])
  }
  if (!is.null(token_nav)) {
    radio = if (!is.null(top_nav) && top_nav == 1) T else F
    nav = if (methods::is(tokens[[token_nav]], 'factor')) levels(tokens[[token_nav]]) else unique(tokens[[token_nav]])
  }
  if (!is.null(nav)) nav = id_nav(nav, colors, navfilter, radio=radio)

  template = html_template('browser', css_str=css_str, col1=style_col1, col2=style_col2)
  template$header = gsub('$NAVIGATION$', if (is.null(nav)) '' else nav, template$header, fixed = T)

  if (!is.null(top_nav) | thres_nav > 1) {
    navmeta = nav_meta_label(top_nav,thres_nav)
    template$header = gsub('$NAVIGATION_META$', navmeta, template$header, fixed = T)
  } else {
    template$header = gsub('$NAVIGATION_META$', '', template$header, fixed = T)
  }

  if (!n && subheader != '') subheader = sprintf('<i>%s</i>', subheader)
  if (n && subheader != '') subheader = sprintf('<i>%s, N = <ndoc>%s</ndoc></i>', subheader, n_doc)
  if (n && subheader == '') subheader = sprintf('<i>N = <ndoc>%s</ndoc></i>', n_doc)
  template$header = gsub('$SUBHEADER$', subheader, template$header, fixed = T)

  template$header = gsub('$HEADER$', header, template$header, fixed = T)


  save_html(docstring, template, filename)
}

#' Convert tokens into full texts in an HTML file with highlighted tokens
#'
#' @param tokens    A data.frame with a column for document ids (doc_col)
#'                  and a column for tokens (token_col)
#' @param value     Either a logical vector or a numeric vector with
#'                  values between 0 and 1. If a logical vector is used, then tokens
#'                  with TRUE will be highlighted (with the color specified in pos_col).
#'                  If a numeric vector is used, the value determines the alpha (transparency),
#'                  with 0 being fully transparent and 1 being fully colored.
#' @param meta      A data.frame with a column for document_ids (doc_col). All other columns are added
#'                  to the browser as document meta
#' @param col       The color used to highlight
#' @param doc_col   The name of the document id column
#' @param token_col The name of the token column
#' @param doc_nav   The name of a column in meta, used to set a navigation tag
#' @param token_nav Alternative to doc_nav, a column in the tokens, used to set a navigation tag
#' @param filename  Name of the output file. Default is temp file
#' @param unfold    Either a character vector or a named list of vectors of the same length as tokens. If given, all tokens with a tag can be clicked on to unfold the given text. If a list of vectors is given,
#'                  the values of the columns are concatenated with the column name. E.g. list(doc_id = 1, sentence = 1) will be [doc_id = 1, sentence = 2].
#' @param span_adjacent If TRUE, include adjacent tokens with identical attributes within the same tag
#' @param ...       Additional formatting arguments passed to create_browser()
#'
#' @return The name of the file where the browser is saved. Can be opened conveniently from within R using browseUrl()
#' @export
#'
#' @examples
#' ## as an example, highlight words based on word length
#' highlight = nchar(as.character(sotu_data$tokens$token))
#' highlight = highlight / max(highlight)
#' highlight[highlight < 0.3] = NA
#' url = highlighted_browser(sotu_data$tokens, value = highlight, sotu_data$meta)
#'
#' \donttest{
#' view_browser(url)   ## view browser in the Viewer
#' }
#' if (interactive()) {
#' browseURL(url)     ## view in default webbrowser
#' }
highlighted_browser <- function(tokens, value, meta=NULL, col='yellow', doc_col='doc_id', token_col='token', doc_nav=NULL, token_nav=NULL, filename=NULL, unfold=NULL, span_adjacent=T, ...){
  tokens[[token_col]] = highlight_tokens(tokens[[token_col]], value=value, col = col, unfold=unfold, span_adjacent = span_adjacent, doc_id=tokens[[doc_col]])
  create_browser(tokens, meta, doc_col, token_col, doc_nav=doc_nav, token_nav=token_nav, filename=filename, ...)
}

#' Convert tokens into full texts in an HTML file with color ramp highlighting
#'
#' @param tokens    A data.frame with a column for document ids (doc_col)
#'                  and a column for tokens (token_col)
#' @param value     A numeric vector with values between -1 and 1. Determines the color
#'                  mixture of the scale colors specified in col_range
#' @param alpha     Optionally, the alpha (transparency) can be specified, with 0 being fully
#'                  transparent and 1 being fully colored. This can be a vector to specify a
#'                  different alpha for each value.
#' @param meta      A data.frame with a column for document_ids (doc_col). All other columns are added
#'                  to the browser as document meta
#' @param col_range The color used to highlight
#' @param doc_col   The name of the document id column
#' @param token_col The name of the token column
#' @param doc_nav   The name of a column in meta, used to set a navigation tag
#' @param token_nav Alternative to doc_nav, a column in the tokens, used to set a navigation tag
#' @param filename  Name of the output file. Default is temp file
#' @param unfold  Either a character vector or a named list of vectors of the same length as tokens. If given, all tokens with a tag can be clicked on to unfold the given text. If a list of vectors is given,
#'                the values of the columns are concatenated with the column name. E.g. list(doc_id = 1, sentence = 1) will be [doc_id = 1, sentence = 2].
#' @param span_adjacent If TRUE, include adjacent tokens with identical attributes within the same tag
#' @param ...       Additional formatting arguments passed to create_browser()
#'
#' @return The name of the file where the browser is saved. Can be opened conveniently from within R using browseUrl()
#' @export
#'
#' @examples
#' ## as an example, scale word colors based on number of characters
#' scale = nchar(as.character(sotu_data$tokens$token))
#' scale[scale>6] = scale[scale>6] +20
#' scale = rescale_var(sqrt(scale), -1, 1)
#' scale[abs(scale) < 0.5] = NA
#' url = colorscaled_browser(sotu_data$tokens, value = scale, meta=sotu_data$meta)
#'
#' \donttest{
#' view_browser(url)   ## view browser in the Viewer
#' }
#' if (interactive()) {
#' browseURL(url)     ## view in default webbrowser
#' }
colorscaled_browser <- function(tokens, value, alpha=0.4, meta=NULL, col_range=c('red','blue'), doc_col='doc_id', token_col='token', doc_nav=NULL, token_nav=NULL, filename=NULL, unfold=NULL, span_adjacent=T, ...){
  tokens[[token_col]] = colorscale_tokens(tokens=tokens[[token_col]], value=value, col_range = col_range, alpha=alpha, unfold=unfold, span_adjacent = span_adjacent, doc_id=tokens[[doc_col]])
  create_browser(tokens, meta, doc_col, token_col, doc_nav=doc_nav, token_nav=token_nav, filename=filename, ...)
}

#' Convert tokens into full texts in an HTML file with category highlighting
#'
#' @param tokens    A data.frame with a column for document ids (doc_col)
#'                  and a column for tokens (token_col)
#' @param category  Either a numeric vector with values representing categories, or a factor vector, in which case
#'                  the values are used as labels. If a numeric vector is used, the labels can also be specified in the labels argument
#' @param alpha     Optionally, the alpha (transparency) can be specified, with 0 being fully
#'                  transparent and 1 being fully colored. This can be a vector to specify a
#'                  different alpha for each value.
#' @param labels    A character vector giving names to the unique category values. If category is a factor vector, the factor levels are
#'                  used.
#' @param meta      A data.frame with a column for document_ids (doc_col). All other columns are added
#'                  to the browser as document meta.
#' @param colors    A character vector with color names for unique values of the category argument. Has to be the same length
#'                  as unique(na.omit(category))
#' @param doc_col   The name of the document id column
#' @param token_col The name of the token column
#' @param filename  Name of the output file. Default is temp file
#' @param span_adjacent If TRUE, include adjacent tokens with identical attributes within the same tag
#' @param unfold  Either a character vector or a named list of vectors of the same length as tokens. If given, all tokens with a tag can be clicked on to unfold the given text. If a list of vectors is given,
#'                the values of the columns are concatenated with the column name. E.g. list(doc_id = 1, sentence = 1) will be [doc_id = 1, sentence = 2].
#' @param ...       Additional formatting arguments passed to create_browser()
#'
#' @return The name of the file where the browser is saved. Can be opened conveniently from within R using browseUrl()
#' @export
#'
#' @examples
#' ## as an example, use simple grep to code tokens
#' code = rep(NA, nrow(sotu_data$tokens))
#' code[grep('war', sotu_data$tokens$token)] = 'War'
#' code[grep('mother|father|child', sotu_data$tokens$token)] = 'Family'
#' code = as.factor(code)
#' url = categorical_browser(sotu_data$tokens, category=code, meta=sotu_data$meta)
#'
#' \donttest{
#' view_browser(url)   ## view browser in the Viewer
#' }
#' if (interactive()) {
#' browseURL(url)     ## view in default webbrowser
#' }
categorical_browser <- function(tokens, category, alpha=0.3, labels=NULL, meta=NULL, colors=NULL, doc_col='doc_id', token_col='token', filename=NULL, unfold=NULL, span_adjacent=T, ...){
  if (methods::is(category, 'character')) category = as.factor(category)
  if (methods::is(category, 'numeric') && is.null(labels)) labels = stats::na.omit(unique(category))
  if (methods::is(category, 'factor')) {
    category = droplevels(category)
    if (is.null(labels)) labels = stats::na.omit(levels(category))
    category = as.numeric(category)
  }

  if (is.null(meta)) {
    meta = data.frame(doc_id = unique(tokens[[doc_col]]))
    colnames(meta) = doc_col
  }
  if (is.null(colors)) colors = grDevices::rainbow(length(unique(stats::na.omit(category))))

  tokens[[token_col]] = category_highlight_tokens(tokens[[token_col]], category=category, labels=labels, alpha=alpha, colors = colors, unfold=unfold, span_adjacent = span_adjacent, doc_id=tokens[[doc_col]])
  tokens[['multi_cat']] = factor(category, labels=labels)
  create_browser(tokens, meta, doc_col, token_col, token_nav='multi_cat', filename= filename, colors=colors, ...)
}

