#' create the html template
#'
#' @param template The name of the template to be used
#' @param css_str A character string, to be directly added to the css style header
#' @param col1    The first style color (top bar color)
#' @param col2    The second style color (background color)
#'
#' @return A list with the html header and footer
html_template <- function(template, css_str=NULL, col1="#7D1935", col2="#F5F3EE") {
  TEMPLATE = system.file(sprintf("template/%s.html", template), package="tokenbrowser", mustWork=T)
  html = readChar(TEMPLATE, file.info(TEMPLATE)$size)
  html = stringi::stri_split(str = html, fixed = '$CONTENT$')[[1]]
  html = list(header = html[1], footer = html[2])

  css = create_css(template, css_str=css_str, col1=col1, col2=col2)

  html$header = gsub("$CSS$", css, html$header, fixed = T)
  html
}

create_css <- function(template, css_str=NULL, col1="#7D1935", col2="#F5F3EE") {
  CSS_TEMPLATE = system.file(sprintf("template/%s.css", template), package="tokenbrowser", mustWork=T)
  css = readLines(CSS_TEMPLATE, warn=F)

  ## add custom settings
  if (!is.null(css_str)) {
    css = c(css, css_str)
  }

  ## replace custom settings
  css = stringi::stri_paste(css, collapse='\n')
  css = gsub('$col1$', col1, css, fixed = T)
  css = gsub('$col2$', col2, css, fixed = T)

  css
}

