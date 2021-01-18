token_nav_string <- function(tokens, meta, doc_col, token_nav, top_nav=NULL, thres_nav=NULL) {
  if (!is.null(top_nav) | !is.null(thres_nav)) {
    utok = stats::aggregate(tokens$doc_id, by=list(doc_id=tokens[[doc_col]], nav=tokens[[token_nav]]), FUN='length')
    if (!is.null(thres_nav)) utok = utok[utok$x >= thres_nav,]
    if (!is.null(top_nav)) {
      utok = utok[order(-utok$x),]
      utok = split(as.character(utok$nav), utok[[doc_col]])
      utok = lapply(utok, utils::head, n=top_nav)
    } else {
      utok = split(as.character(utok$nav), utok[[doc_col]])
    }
  } else {
    utok = stats::na.omit(unique(tokens[,c(doc_col,token_nav)]))
    utok = split(utok[[token_nav]], utok[[doc_col]])
  }

  #browser()
  ids = names(utok)

  #ids
  #navstring

  navstring = sapply(utok, function(x) paste(stats::na.omit(sprintf('<tag>%s</tag>', x)), collapse=', '), simplify = F)
  #navstring = stringi::stri_paste_list(navstring, sep = ', ')
  #navstring

  out = rep('', nrow(meta))
  out[match(ids, meta[[doc_col]])] = navstring
  out
}

nav_meta_label <- function(top_nav, thres_nav) {
  if (!is.null(top_nav)) {
    if (thres_nav <= 1)
      navmeta = sprintf('*filter applies to top %s', if (top_nav == 1) 'label' else paste0(top_nav, ' labels'))
    else
      navmeta = sprintf('*filter applies to top %s with at least %s tokens', if (top_nav == 1) 'label' else paste0(top_nav, ' labels'), thres_nav)
  } else {
    if (thres_nav > 1) navmeta = sprintf('*filter applies to labels with at least %s tokens', thres_nav)
  }
  navmeta
}
