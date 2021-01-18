tabDocu = fluidPage(align="left",
                    column(
                      box(
                        h1(paste("Genepop version", config$version_genepop)),
                        #div(withMathJax(includeMarkdown(getFile("GenepopS.Rmd")))),
                        width = 0
                      ),
                      width = 12)
)
