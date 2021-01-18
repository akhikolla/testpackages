shinyUI(
    fluidPage(
        # CSS
        tags$head(
            tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
            tags$title("Phase I: Vaccine Dose Selection Design and Analysis"),
            tags$style(".shiny-file-input-progress {display: none}"),
        tags$script(src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML", type="text/javascript")
        ),


        withTags({
            div(class = "cheader",
                h5("Phase I: Vaccine Dose Selection Design and Analysis"),
                tags$button(
                    id = 'close',
                    type = "button",
                    class = "btn action-button",
                    onclick = "setTimeout(function(){window.close();},500);",
                    "Exit",
                    style="float: right;"
                )
            )
        }),

        # Main Page
        uiOutput("mainpage")
    )
)
 
