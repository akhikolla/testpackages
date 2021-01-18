#'   @import ggplot2
#'   @import ggthemes


## Shiny GUI for Eagle
## Developer:  Andrew W. George
## Version: 2.1.0



rootdir <-  c('Home' = Sys.getenv("HOME"))
if(.Platform$OS.type == "windows") {
     rootdir <-  c('Home' = paste0(rootdir, "\\..\\"))
}



##---------------------------
## Analyse Page Functions
##~~~~~~~~~~~~~~~~~~~~~~~~~~

bannerAnal <- function()
{
   page =  fluidPage(
 fluidRow( column(12, { tags$div(img(src = "images/analyse_banner.jpg", style="width: 100% ; height: 100%")) }

                      ) ## end column(12, )
              ), ## end fluidRow
              br(),
              fluidRow(column(12, 
                       shinyBS::bsButton(inputId="dummy4", label="Hover here for details",
                       style="warning", size="large", type="action", block=TRUE,
                       icon=icon("question-circle-o"))
                      ) ## end column
             ) ## end fluidRow
          ) ## end fluidPage
   return(page)
}

row1Anal <- function(){
  page <-  fluidRow(column(12,  
                    wellPanel(
                       uiOutput("analyse_names"),

                       shinyBS::bsTooltip("analyse_names",
                       title='<font size="3" > Select a single variable to be treated as the trait for the analysis  </font>',
                       placement="right", trigger="hover", options=list(container="body"))

                       ) ## end wellPanel
                     ) ## end column
                   ) ## end fluidRow                             
  return(page)
}


row2Anal <- function()
{
 page =   fluidRow(column(12,
          wellPanel(
            uiOutput("analyse_fnames"),
              shinyBS::bsTooltip("analyse_fnames",
              title='<font size="3" > Select the variables, if any, to be used as fixed effects in the analysis. If no variables are selected, then only an overall mean will be fitted. </font>',
                placement="right", trigger="hover", options=list(container="body")),
              textOutput("fmodel")
           ) ## end wellPanel
                ) ## end column
          ) ## end fluidRow
return(page)
}



row3Anal <- function()
{
page =  fluidRow(column(12,  wellPanel(
           numericInput(inputId="analyse_cpu", label=h4("Step 3: Specify number of cpu"), value=1),
                    style="padding: 1px",
                    shinyBS::bsTooltip("analyse_cpu",
           title='<font size="3" > set to the number of cpu available for distributed computing. </font>',
           placement="right", trigger="hover",
               options=list(container="body"))
           ) ## end wellpanel
        ))  ## end column and fluidRow
  return(page)
}


row4Anal <- function()
{
page =  fluidRow( column(12, 
         wellPanel(
          radioButtons(inputId="analyse_lambda", label=h4("Step 4: Specify lambda value (controls the false positive rate)"), 
        choiceNames = list(
        tags$span(style = "font-size:18px", "Set manually"),
        tags$span(style = "font-size:18px", "Set automatically (via permutation)")),
                                                    choiceValues=c("manual","auto")),
                        style="padding: 1px",
 shinyBS::bsTooltip("analyse_lambda", title='<font size="3" > Select the manual option if you want a quick analysis. If you leave the lambda value at 1, its default value, this will be a conservative analysis. Select auto if you want to perform an analysis with a specified false positive rate. This analysis will take about 5 times as long as a permutation step is performed to fine-tune the lambda value for the desired false positive rate.  </font>',
                               placement="right", trigger="hover", options=list(container="body"))

                     )  ## column 12
          ) ,  ## wellPanel 
                                          

                conditionalPanel(
                   condition = "input.analyse_lambda == 'manual'",
                        wellPanel(
                        fluidRow(column(12,
                           sliderInput(inputId="analyse_setlambda", label=h4("Specify lambda value. "),
                               value=1, min = 0, max = 1, step = 0.01),
                           style="padding: 1px",
                           shinyBS::bsTooltip("analyse_setlambda", title='<font size="3" >The lambda parameter controls the conservativeness of the model building process. Values closer to 1 (0) decrease (increase) the false positive rate. The default value is 1 - its most conservative setting. </font>',
                               placement="right", trigger="hover", options=list(container="body"))
                               ) ## colunn12,
                           ) ## fluidRow
                        ) ## wellPanel
                ), ## conditionalPanel

                conditionalPanel(
                   condition = "input.analyse_lambda == 'auto'",
                        wellPanel(
                       fluidRow(column(12, 
                           sliderInput(inputId="analyse_fpr", label=h4("Specify desired false positive rate."),
                               value=0.05, min = 0.01, max = 0.5, step = 0.01),
                           style="padding: 1px",
                           shinyBS::bsTooltip("analyse_fpr", title='<font size="3" > Set the slider to the desired false positive rate for the analysis. The default value is 0.05.   </font>',
                               placement="right", trigger="hover", options=list(container="body"))
 
                       ) ## end column 12
                       ),  ## end fluidRow



                        fluidRow(column(12,
                   sliderInput(inputId="analyse_numreps", label=h4(" Specify number of replicates."),
                                       value=200, min = 30, max = 1000, step = 5),
                                        style="padding: 1px",
                   shinyBS::bsTooltip("analyse_numreps", title='<font size="3" > To set the number of replicates, start with 200 replicates and increase in 50 replicate increments. Stop when the lambda value stabilizes.  </font>',
                          placement="right", trigger="hover", options=list(container="body"))



                               ) ## colunn12,
                           ) ## fluidRow
                        ) ## wellPanel
                ) ## conditionalPanel

 
                                           
    )  ## end fluidRow 
                                           
                                           
                                           
                        
return(page)                                         
}  ## end function row4Analyse 







row5Anal <- function()
{
   page = fluidRow(column(12,
             wellPanel(
                h4("Step 5: Perform genome-wide analysis"),
                            actionButton(inputId="analyse_go",label="", width='35%', style='padding:5px 5px 5px 5px; font-size:180%',
                                         icon=icon("upload", lib="glyphicon")),
                                          style='padding: 1px',
                                        shinyBS::bsTooltip("analyse_go",
                     title='<font size="3" >  Click here to find the set of snp in strongest association with the trait. Manually setting the lambda value is much faster than having the lambda value set automatically. However, automatically is the preferred option.  </font>',
                          placement="right", trigger="hover",
                          options=list(container="body"))
             ) ## wellPanel
            )   ## column12
         ) ## end fluidRow
  return(page)
}


##-----------------------------------
## Rows for Plotting Page 
##-----------------------------------
banner1Plot <- function(){
  page <- fluidRow(
             column(12, {
                  tags$div(img(src = "images/plot_banner.jpg", style="width: 100% ; height: 100%"))
                        }
                   ) ## end column(12, )
                  ) ## end fluidRow
     return(page)
}

banner2Plot <- function(){
   page <- fluidRow(column(12, 
                       shinyBS::bsButton(inputId="plot_overview", label="Hover here for details",
                       style="warning", size="large", type="action", block=TRUE,
                       icon=icon("question-circle-o"))
                           ) ## end column
                   ) ## end fluidRow
           return(page)             
}

col1Plot <-   function(){ 
   page <-  column(3,
                   wellPanel( 
                      uiOutput("plot_choice"), 

                      shinyBS::bsTooltip("plot_choice",
                     title='<font size="3" > Manhattan plot type -  score statistics are converted into  p-values and their -log values plotted on the y-axis. Score statistic plot type - the score statistics are plotted on the y-axis. The score statistics are used by Eagle to identify the SNP in strongest association with the trait.  A new set of score statistics are generated at each iteration of the model buidling process. </font>',
                      placement="right", trigger="hover", options=list(container="body"))

                 ) ## end wellPanel
             )  # column
  return(page)
}



col2Plot <- function(){
  page <-  column(3,
                  wellPanel(
                    uiOutput("plot_modelits"),
                    shinyBS::bsTooltip("plot_modelits",
                       title='<font size="3" > SNP-trait associations are found by building the best model iteratively.  </font>',
                       placement="right", trigger="hover", options=list(container="body"))

                  ) ## end wellPanel
           ) ## end column
  return(page)
}


col3Plot <- function(){

  page <-  column(3,
             wellPanel(
               uiOutput("plot_chromosomes"),

               shinyBS::bsTooltip("plot_chromosommes",
               title='<font size="3" > Select a chromosome or the entire genome (All)  </font>',
               placement="right", trigger="hover", options=list(container="body"))

            ) ## end wellPanel
          ) ## end column
  return(page)
}


col4Plot <- function()
{
   page = column(3,
             wellPanel(
                shinyjs::useShinyjs(),
                            actionButton(inputId="plot_go",label=HTML('<font size="4">Generate Plot</font>'), width='100%', style='padding:1px 1px 1px 1px; font-size:100%'),
                                          style='padding: 1px',
                                        shinyBS::bsTooltip("plot_go",
                     title='<font size="3" >   Press this button to generate the plot. If you change your selections after the plot has been created, you will need to click this button again to recreate the plot.  </font>',
                          placement="left", trigger="hover",
                          options=list(container="body"))
            ) ## wellPanel
         )   ## column3
  return(page)
}

##---------------------------
## Rows for Help page
##---------------------------

row1Help <- function(){
  page <- fluidRow(
             column(12, 
               includeHTML("help.html")
             )  ## end column
          ) ## end fluidRow
  return(page)
 }





#------------------------------
# Misc functions
#------------------------------


home_intro <- function(){
  txt <- "
strong(Eagle)  is a software package for genome-wide association mapping.
It differs from most other association mapping packages in that it fits all marker-trait associations simultaneously, an 
returns the 'best' set of snp loci in strongest association with a trait as its findings. 
Eagle can handle data collected from populations of arbitrary structure. The populations can contain inbred or outbred individuals. 
br()
An analysis is performed by reading in the marker data (Read Genotypes), reading in the phenotypic data (Read Phenotypes), reading in the 
marker map if known (Read Map), reading in the Z matrix if needed, and performing the genome-wide analysis (Analyse).  
br()
Help is available by hovering over the widgets or by clicking on the help tab at the top of the screen.  "
    return(txt)
}

read_geno_intro <- function(){
  txt <- "
  Eagle can handle two different types of marker data; genotype data in a plain text space separated file 
  "
  return(txt)
}
read_pheno_intro <- function(){
  txt <- "adfadf"
  
  return(txt)
}  
  
    
  

##-------------------
## ui Body
##--------------------

## ui.R ##
#library(shiny)
#library(shinythemes)
#library(shinyBS)
#library(shinyjs)
#library(shinyFiles)
#library(ggthemes)

FullPage <- navbarPage(title="Eagle: Genome-wide association mapping",  
                      theme = shinythemes::shinytheme("flatly"),
                #       theme = shinytheme("paper"),
                #      theme = shinytheme("united"),

                      
                       ##----------------------------##
                       ##   Home Page                ##
                       ##----------------------------##
                            tabPanel("Home", icon=icon("fa-home", class="fa fa-home  fa-lg"),
                            tags$head(includeCSS("css.css")),
                            fluidPage(
                              fluidRow(
                                column(12,
                                tags$div(img(src = "images/homescreen.jpg", 
                                             style="width: 100% ; height: 100%"))
                                )
                              ) ## end fluidRow
                            
                                
                              
                            ) ## end fluidPage
                                ), ## end tabPanel("Home") 

                       ##-----------------------##
                       ##     Read Genotypes    ##
                       ##-----------------------##
                                 
                      tabPanel("Read Genotypes",  icon=icon("file-o"), 
                               tags$head(tags$style(HTML('

                                                         .popover {
                                                         max-width: 80%;
                                                         
                                                         }
                                                         '))
                               ),




                            fluidPage(
                              fluidRow(
                                column(12,  {
                                       tags$div(img(src = "images/marker_banner.jpg", 
                                                    style="width: 100% ; height: 100%;"))
                               
                                }
                                       ) ## end column(12, )
                              ), ## end fluidRow
                              br(),
                              fluidRow(column(12, 
                                            
                                              shinyBS::bsButton(inputId="dummy1", label="Hover here for details", 
                                                    style="warning", size="large", type="action", block=TRUE, 
                                                    icon=icon("question-circle-o")
                                                    )
                                           
                                       ) ## end column
                              ), ## end fluidRow
                              
                              
                              br(),


                              fluidRow(
                                column(5, 
                                       fluidPage(
                                         fluidRow(
                                           column(12,
                                                  wellPanel(
                                                  radioButtons(inputId="filetype", label=h4("Step 1: Choose file type"), 
                                                    choiceNames=list(

        tags$span(style = "font-size:18px", "vcf"), 
        tags$span(style = "font-size:18px", "PLINK"), 
        tags$span(style = "font-size:18px", "Text/ASCII")), 
                                                    choiceValues=c("vcf", "plink","text")),
                                                  style="padding: 1px",
                                                  shinyBS::bsTooltip("filetype",
title='<font size="3" > click on file type </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  ),  ## wellPanel
                                           
                                           conditionalPanel(
                                             condition = "input.filetype == 'text'",
                                             wellPanel(
                                               fluidPage(
                                                 wellPanel(
                                                   h5("Assign marker genotypes to snp genotypes AA, AB, BB, and missing"),
                                                 fluidRow(
                                                  column(4, textInput(inputId="AA",label="AA", value="") ),
                                                  column(4, textInput(inputId="AB",label="AB", value="") ),
                                                  column(4, textInput(inputId="BB",label="BB", value="") ),
                                                  column(4, textInput(inputId="missing",label="missing", value="") ) ,
shinyBS::bsTooltip("AB", 
title='<font size="3" > Only a single value can be entered. If inbreds, leave blank  </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")),
shinyBS::bsTooltip("missing",
title='<font size="3" > Enter genotype code used in file that is to be missing. Leave blank if data contains no missing marker genotypes </font>' , 
placement="right", trigger="hover",
                                                            options=list(container="body"))





                                                ) ## end fluidRow
                                               ) ## end inner wellPanel
                                               )  ## end fluidPage

                                                 ) ## end Wellpanel
                                            
                                           ) ## end conditionalPanel
                                           
                                           
                                           
                                           
                                                  ) ## end column
                                           
                                         ), ## end fluidRow choose file type
                        
                                         fluidRow(
                                           column(12, wellPanel(
                                                  numericInput(inputId="memsize", label=h4("Step 2: Specify available memory in Gbytes"), 
                                                               value=8, min = 2, max = NA, step = NA),
                                                  style="padding: 1px",
                                                  shinyBS::bsTooltip("memsize",
title='<font size="3" > set to maximum available memory in Gbytes  </font>',
placement="right", trigger="hover",
                                                          options=list(container="body"))
                                                  )) ## end column
                                           
                                           
                                         ), ## end fluidRow specify amout of memory
                                        
                                         
                                      fluidRow(column(12, 
                                        wellPanel(
                                        h4("Step 3: Select marker file"),
                                        shinyFilesButton('choose_marker_file', 'Select File', 'Please select file', FALSE),
                                        textInput("choose_marker_file_text", label = h5("or enter file name (including full path)"))
                                         
                                           
                                          )  ## end wellPannel
                                         )
                                         ), ## end fluidRow


                                      






 
                                         
                                         fluidRow(column(12, 
                                                       wellPanel(
                                                          shinyjs::useShinyjs(),
                                                          h4("Step 4: Upload file"),




                                                          actionButton(inputId="marker_go",label="", width='35%', style='padding:5px 5px 5px 5px; font-size:180%',
                                                                       icon=icon("upload", lib="glyphicon")),



                                                          style='padding: 1px',
                                                          shinyBS::bsTooltip("marker_go", 
title='<font size="3" > Upload marke data file. <br> This may take some time if the file is large.  </font>',
placement="right", trigger="hover",
                                                                     options=list(container="body"))




 
                                                        )
                                                  )
                                         ) ## end fluidRow
                                         
                                         
                                         
                                       ) ## end fluidPage -- widgets on left hand side
                                       
                                       
                                       
                                       
                                       ), ## end column(6,  )  -- left half of page
                                          ## for input widgets
                                column(7, 
                                        verbatimTextOutput("ReadMarker", placeholder=TRUE) 




                                )  ## end column(6, ) -- right half of page
                                
                              ) ## end fluidRow
                              
                              
                              
                              
                              
                            ) ## end fluidPage    
                                
                                
                                
                                
                                
                      ),  ## end tabPanel("Read Genotypes")
                       
                       
                       ##----------------------##
                       ## Read Phenotypes      ##
                       ##----------------------##
                      
                                    
                      tabPanel("Read Phenotypes",  icon=icon("file-o"), 
                               tags$head(tags$style(HTML('
                                                         .popover {
                                                         max-width: 80%;
                                                         
                                                         }
                                                         '))
                               ),


                            fluidPage(
                              fluidRow(
                                column(12, {
                                       tags$div(img(src = "images/pheno_banner.jpg", 
                                                    style="width: 100% ; height: 100%; "))
                               
                                }
                                       ) ## end column(12, )
                              ), ## end fluidRow
                              br(),
                              fluidRow(column(12, 
                                             shinyBS::bsButton(inputId="dummy2", label="Hover here for details",
                                                    style="warning", size="large", type="action", block=TRUE, 
                                                    icon=icon("question-circle-o")
                                                    )
                                           
                                       ) ## end column
                              ), ## end fluidRow
                              
                              
                              br(),
                              fluidRow(
                                column(5, 
                                       fluidPage(
                                         fluidRow(
                                           column(12,
                                                  wellPanel(
                                                  radioButtons(inputId="pheno_header", label=h4("Step 1: Select if file contains column names"), 
                                                               choices=c("yes"="yes","no"="no" ), selected="yes"),
                                                  style="padding: 1px",
                                                  shinyBS::bsTooltip("pheno_header",
title='<font size="3" > click on yes if the first row of the file contains the column names. Generic names will be assigned if no is clicked.  </font>',
placement="right", 
trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  )  ## wellPanel
                                           
                                           
                                           
                                             ) ## end column
                                           
                                         ), ## end fluidRow choose file type
                        
                                         fluidRow(
                                           column(12, wellPanel(
                                                  radioButtons(inputId="pheno_csv", label=h4("Step 2: Is the file comma separated"),
                                                               choices=c("yes"="yes","no"="no" ), selected="no" ),
                                                  style="padding: 1px",
                                                  shinyBS::bsTooltip("pheno_csv",
title='<font size="3" > click on yes if the file is a csv file. </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  )  ## wellPanel
                                             ) ## end column
                                           
                                         ), ## end fluidRow 
                                        

                                         fluidRow(
                                           column(12, wellPanel(
                                                   textInput(inputId="pheno_missing", label=h4("Step 3: Code for missing value", value="") ),
                                                  style="padding: 1px",
                                                  shinyBS::bsTooltip("pheno_missing",
title='<font size="3" > Assign value that denotes a missing value. Leave blank if file does not contain missing data. </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  )  ## wellPanel
                                             ) ## end column

                                         ), ## end fluidRow 

                                      fluidRow(column(12,
                                        wellPanel(
                                        h4("Step 4: Select phenotypic file"),
                                        shinyFilesButton('choose_pheno_file', 'Select File', 'Please select file', FALSE),
                                        textInput("choose_pheno_file_text", label = h5("or enter file name (including full path)"))
                                          )  ## end wellPannel
                                         )
                                         ), ## end fluidRow




                                         fluidRow(column(12, 
                                                       wellPanel(
                                                          shinyjs::useShinyjs(),
                                                          h4("Step 5: Upload file"),




                                                          actionButton(inputId="pheno_go",label="", width='35%', style='padding:5px 5px 5px 5px; font-size:180%',
                                                                       icon=icon("upload", lib="glyphicon")),



                                                          style='padding: 1px',
                                                          shinyBS::bsTooltip("pheno_go",
title='<font size="3" > Upload phenotypic file.  This file can contain multiple trait and covariate data.</font>',
placement="right", trigger="hover",
                                                                     options=list(container="body"))




 
                                                        )
                                                  )
                                         ) ## end fluidRow
                                         
                                         
                                       ) ## end fluidPage -- widgets on left hand side
                                       
                                       
                                       
                                       
                                       ), ## end column(6,  )  -- left half of page
                                          ## for input widgets

                               column(7,
                                        verbatimTextOutput("ReadPheno", placeholder=TRUE)
                                )  ## end column(6, ) -- right half of page
                              ) ## end fluidRow
                              
                              
                              
                              
                              
                            ) ## end fluidPage    
                                
                                
                                
                                
                                
                      ),  ## end tabPanel("Read Phenotypes")





#              #-----------------------TESTING ======================================= 
#                                    
#                     tabPanel("Read Phenotypes",  icon=icon("file-o"), 
#                               tags$head(tags$style(HTML('
#                                                         .popover {
#                                                         max-width: 80%;
#                                                         
#                                                         }
#                                                         '))
#                               ),
#
#
#                            fluidPage(
#                              fluidRow(
#                                column(12, {
#                                       tags$div(img(src = "images/pheno_banner.jpg", 
#                                                    style="width: 100% ; height: 100%; "))
#                               
#                                }
#                                       ) ## end column(12, )
#                              ), ## end fluidRow
#                              br(),
#                              fluidRow(column(12, 
#                                             shinyBS::bsButton(inputId="dummy2", label="Hover here for details",
#                                                    style="warning", size="large", type="action", block=TRUE, 
#                                                    icon=icon("question-circle-o")
#                                                    )
#                                           
#                                       ) ## end column
#                              ), ## end fluidRow
#                              
#                              
#                              br(),
#                              fluidRow(
#                                column(5, 
#                                       fluidPage(
#                                         fluidRow(
#                                           column(12,
#                                                  wellPanel(
#                                                  radioButtons(inputId="pheno_header", label=h4("Step 1: Select if file contains column names"), 
#                                                               choices=c("yes"="yes","no"="no" )),
#                                                  style="padding: 1px",
#                                                  shinyBS::bsTooltip("pheno_header",
#title='<font size="3" > click on yes if the first row of the file contains the column names. Generic names will be assigned if no is clicked.  </font>',
#placement="right", 
#trigger="hover",
#                                                            options=list(container="body")
#                                                      )
#                                                  )  ## wellPanel
#                                           
#                                           
#                                           
#                                             ) ## end column
#                                           
#                                         ), ## end fluidRow choose file type
#
#
# fluidRow(column(12,
#                                          wellPanel(
#                                            h4("Step 3: Select marker file"),
#
#                                           actionButton(inputId="choose_marker_file", h6("Choose File")), br(),
#                                           textOutput("choose_marker_file"),
#                                           style='padding: 1px',
#                                           shinyBS::bsTooltip("choose_marker_file",
#title='<font size="3" >WARNING! File browser window may open behind web browser  </font>',
#placement="right",
#trigger="hover",
#                                                     options=list(container="body"))
#
#
#                                          )
#                                         )
#                                         ) ## end fluidRow
#
#
#
#
#
#                        
#)
#)
#))),
#




                      ##-----------------------------------##
                      ## Read Z matrix (if needed)         ##
                      ##-----------------------------------##
                      
                       tabPanel("Read Z matrix (if needed)", icon=icon("file-o"), 
                               tags$head(tags$style(HTML('

                                                         .popover {
                                                         max-width: 80%;
                                                         
                                                         }
                                                         '))
                               ),


                            fluidPage(
                              fluidRow(
                                column(12, {
                                       tags$div(img(src = "images/Zmat_banner.jpg", 
                                                    style="width: 100% ; height: 100%"))
                               
                                }
                                       ) ## end column(12, )
                              ), ## end fluidRow
                              br(),
                              fluidRow(column(12, 
                                          shinyBS::bsButton(inputId="Zmat1", label="Hover here for details",
                                          style="warning", size="large", type="action", block=TRUE,
                                          icon=icon("question-circle-o")
                                          )

  
                                           
                                       ) ## end column
                              ), ## end fluidRow
                              
                              
                              br(),
                              fluidRow(
                                column(5, 
                                       fluidPage(
                        



                                      fluidRow(column(12,
                                        wellPanel(
                                        h4("Step 1: Select Z matrix file"),
                                        shinyFilesButton('choose_Zmat_file', 'Select File', 'Please select file', FALSE),
                                        textInput("choose_Zmat_file_text", label = h5("or enter file name (including full path)"))


                                          )  ## end wellPannel
                                         )
                                         ), ## end fluidRow









 
                                         
                                         fluidRow(column(12, 
                                                       wellPanel(
                                                          shinyjs::useShinyjs(),
                                                          h4("Step 2: Upload file"),




                                                          actionButton(inputId="Zmat_go",label="", width='35%', style='padding:5px 5px 5px 5px; font-size:180%',
                                                                       icon=icon("upload", lib="glyphicon")),



                                                          style='padding: 1px',
                                                          shinyBS::bsTooltip("Zmat_go",
title='<font size="3" > Upload Z matrix file.    </font>',
placement="right", trigger="hover",
                                                                     options=list(container="body"))
 
                                                        )
                                                  )
                                         ) ## end fluidRow
                                       ) ## end fluidPage -- widgets on left hand side
                                       
                                       
                                       
                                       
                                       ), ## end column(6,  )  -- left half of page
                                          ## for input widgets

                               column(7,
                                        verbatimTextOutput("ReadZmat", placeholder=TRUE)
                                )  ## end column(6, ) -- right half of page

                                
                              ) ## end fluidRow
                              
                            ) ## end fluidPage    
                                
                                
                      ),  ## end tabPanel("Read Z matrix")











                      ##-------------------------##
                      ## Read Marker map         ##
                      ##-------------------------##
                      
                       tabPanel("Read Map (optional)", icon=icon("file-o"), 
                               tags$head(tags$style(HTML('

                                                         .popover {
                                                         max-width: 80%;
                                                         
                                                         }
                                                         '))
                               ),


                            fluidPage(
                              fluidRow(
                                column(12, {
                                       tags$div(img(src = "images/map_banner.jpg", 
                                                    style="width: 100% ; height: 100%"))
                               
                                }
                                       ) ## end column(12, )
                              ), ## end fluidRow
                              br(),
                              fluidRow(column(12, 
                                          shinyBS::bsButton(inputId="dummy3", label="Hover here for details",
                                          style="warning", size="large", type="action", block=TRUE,
                                          icon=icon("question-circle-o")
                                          )

  
                                           
                                       ) ## end column
                              ), ## end fluidRow
                              
                              
                              br(),
                              fluidRow(
                                column(5, 
                                       fluidPage(
                        

                                         fluidRow(
                                           column(12,
                                                  wellPanel(
                                                  radioButtons(inputId="map_header", label=h4("Step 1: Select if file contains column names"),
                                                               choices=c("yes"="yes","no"="no" )),
                                                  style="padding: 1px",
                                                  shinyBS::bsTooltip("map_header",
title='<font size="3" > click on yes if the first row of the file contains the column names. Generic names will be assigned if no is clicked. </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  )  ## wellPanel



                                             ) ## end column

                                         ), ## end fluidRow choose file type


                                         fluidRow(
                                           column(12, wellPanel(
                                                  radioButtons(inputId="map_csv", label=h4("Step 2: Is the file comma separated"),
                                                               choices=c("yes"="yes","no"="no" ), selected="no"),
                                                  style="padding: 1px",
                                                  shinyBS::bsTooltip("map_csv",
title='<font size="3" > click on yes/no  </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  )  ## wellPanel
                                             ) ## end column
                                           
                                         ), ## end fluidRow 
                                        


                                      fluidRow(column(12,
                                        wellPanel(
                                        h4("Step 3: Select map file"),
                                        shinyFilesButton('choose_map_file', 'Select File', 'Please select file', FALSE),
                                        textInput("choose_map_file_text", label = h5("or enter file name (including full path)"))
                                          )  ## end wellPannel
                                         )
                                         ), ## end fluidRow



                                       
                                         
                                         fluidRow(column(12, 
                                                       wellPanel(
                                                          shinyjs::useShinyjs(),
                                                          h4("Step 4: Upload file"),




                                                          actionButton(inputId="map_go",label="", width='35%', style='padding:5px 5px 5px 5px; font-size:180%',
                                                                       icon=icon("upload", lib="glyphicon")),



                                                          style='padding: 1px',
                                                          shinyBS::bsTooltip("map_go",
title='<font size="3" > Upload marker map file.    </font>',
placement="right", trigger="hover",
                                                                     options=list(container="body"))
 
                                                        )
                                                  )
                                         ) ## end fluidRow
                                       ) ## end fluidPage -- widgets on left hand side
                                       
                                       
                                       
                                       
                                       ), ## end column(6,  )  -- left half of page
                                          ## for input widgets

                               column(7,
                                        verbatimTextOutput("ReadMap", placeholder=TRUE)
                                )  ## end column(6, ) -- right half of page


                                
                              ) ## end fluidRow
                              
                            ) ## end fluidPage    
                                
                                
                      ),  ## end tabPanel("Read Map")



                  ##-------------------------------------##
                  ## Analysis                            ##
                  ##-------------------------------------##

                      
                  tabPanel("Analyse", icon=icon("fa-area-chart", class = "fa fa-area-chart fa-lg", lib = "font-awesome"), 
                           tags$head(tags$style(HTML(' 
                                                         .popover {
                                                         max-width: 80%;
                                                         
                                                         }
                                                    '))
                                   ),
                           bannerAnal(),

                           br(),

                           fluidRow(
                               # left half of page
                               column(6,
                                   fluidPage(

                                     row1Anal(),

                                     row2Anal(), 

                                     row3Anal(),


                                     row4Anal(),


                                     row5Anal()

                    

                                   ) # end fluidPage
                              ),  ## end column 

                               # right half of page
                               column(6, 
                                   verbatimTextOutput("AM", placeholder=TRUE)
                                )  ## end column(6, ) -- right half of page





                           ) ## end fluidRow

                                      ), ## end tablPanel










                      navbarMenu("More", 


 
                       tabPanel("Findings", icon=icon("fa-puzzle-piece", class="fa fa-puzzle-piece fa-lg"), 
                               tags$head(tags$style(HTML('

                                                         .popover {
                                                         max-width: 80%;

                                                         }
                                                         '))
                               ),


                            fluidPage(
                              fluidRow(
                                column(12, {
                                       tags$div(img(src = "images/findings_banner.jpg",
                                                    style="width: 100% ; height: 100%"))

                                }
                                       ) ## end column(12, )
                              ), ## end fluidRow
                              br(),
                              fluidRow(column(12,
                                          shinyBS::bsButton(inputId="dummy5", label="Hover here for details",
                                          style="warning", size="large", type="action", block=TRUE,
                                          icon=icon("question-circle-o")
                                          )



                                       ) ## end column
                              ), ## end fluidRow
                           br(), 
                              fluidRow(
                                column(4,
                                       fluidPage(

                                         fluidRow(
                                           column(12, wellPanel(


                                                  h4("Click to calculate additional summary information"),
                                                  actionButton(inputId="pvalue_go",label="", width='35%', style='padding:5px 5px 5px 5px; font-size:180%',
                                                  icon=icon("glyphicon glyphicon-stats", lib="glyphicon")),

                                                  style="padding: 1px",
                                                  shinyBS::bsTooltip("summary",
title='<font size="5" > click on yes/no </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  )  ## wellPanel
                                             ) ## end column

                                         ) ## end fluidRow 
                                ) ## ene fluidPage
                         ),  ## end column
                          column(8, 
                              fluidPage(





                               fluidRow(
                                   column(12, 
tags$div(
         HTML(paste( tags$span(style="color: #ad1d28; font-size: 22px", "Marker-trait Associations"), sep = ""))),
                                       tableOutput("findings")
                                  ) ## end column
                               ), ## end fluidRow

                             fluidRow(
                                column(12, 
                                    conditionalPanel(condition="input.pvalue_go > 0", 
tags$div(
         HTML(paste( tags$span(style="color: #ad1d28; font-size: 22px", "Size and Significance of Effects"), sep = ""))),
                                    tableOutput("size")

                                    ) ## end conditionalPanel
                               ) ## end column
                           ), ## end fluidRow



                            fluidRow(
                                column(12, 
                                    conditionalPanel(condition="input.pvalue_go > 0", 
tags$div(
         HTML(paste( tags$span(style="color: #ad1d28; font-size: 22px", "Extra Summary Information"), sep = ""))),
                                    tableOutput("R")
                                    )
                               ) ## end column
                           ) ## end fluidRow


                         ) ## ene fluidPage 
                       ) # end column
                       ) # end fluidRow
                       ) # end fluidPage

),  ## end tabPAnel pvalue

##----------------------------##
## Plotting                   ##
##----------------------------##
  tabPanel("Plots", icon=icon("fa-pie-chart", class="fa fa-pie-chart fa-lg"),
     tags$head(tags$style(HTML('.popover { max-width: 80%; } '))),
   fluidPage(
      banner1Plot(),
      br(),
      banner2Plot(),
      br(),
      fluidRow(
         col1Plot(),
         col2Plot(),
         col3Plot(),
         col4Plot()
      ), ## end fluidRow
      fluidRow(
        textOutput("caption")
      ),
      fluidRow(
         ## plotly plot
         plotOutput("plot")
      ) ## end fluidRow
     ) ## end fluidPage    
   ),  ## end tabPanel("Plots")


##------------------------------##
## Help                         ##
##------------------------------##

   tabPanel("Help",  icon=icon("question-circle-o", class="fa fa-question-circle-o fa-lg  "),
   fluidPage(
     row1Help()


     ) ## end fluidPage
    ) ## end tabPanel("Help") 

       ) ## end navbarMenu







                       ) ## end navbarPage


## displays eagle log on navbar
FullPage[[3]][[1]]$children[[1]]$children[[1]]$children[[1]] <- 
  tags$img(src = 'images/logo.jpg', width = 80, height = 60)
ui <- FullPage

get_path <- function (defaultpath="/R/library/Eagle/shiny_app/shinydata/genoDemo.dat") {
            path_to_file_res <- tryCatch({
                if(.Platform$OS.type=="unix"){
                    path_to_file_res <- tk_choose.files()
                    #print(path_to_file_res)
                } else {
                    path_to_file_res <- file.choose()                   
                }                
                }, warning = function(war) {
                    print(paste("Eagle::get_path() Warning: ",war))
                    path_to_file_res<-defaultpath
                    return (path_to_file_res)
                }, error = function(err) {
                    print(paste("Eagle::get_path() Error: ",err))
                    path_to_file_res<-defaultpath
                    return (path_to_file_res)
                }, finally = {
                   # path_to_file_res<-"/R/library/Eagle/shiny_app/shinydata/genoDemo.dat"
                  #  return (path_to_file_res)
                }) # END tryCatch
    
            return (path_to_file_res)
      }


server <- function(input, output, session){
  library("Eagle")

  readM <- reactiveValues(path_to_marker_file=NA)
  readP <- reactiveValues(path_to_pheno_file=NA)

  ##------------------------------------------
  ## Intros to pages
  ##-------------------------------------------
 
  output$read_geno_intro <- renderText(read_geno_intro())
  output$read_pheno_intro <- renderText(read_pheno_intro())
  
  ##----------------------------------------
  ##  Read marker path and file name
  ##---------------------------------------- 
  ## upload path and file name
    shinyFileChoose(input=input, id='choose_marker_file', session=session, roots=rootdir)
    observeEvent(input$choose_marker_file, {
           inFile <- parseFilePaths(roots=rootdir, input$choose_marker_file)
           updateTextInput(session, "choose_marker_file_text", value =  as.character(inFile$datapath))
           # path_to_marker_file  <<- as.character(inFile$datapath)
           readM$path_to_marker_file <- as.character(inFile$datapath)
           readM$path_to_marker_file  <- "hmmm"
    })

    observeEvent(input$choose_marker_file_text, {
           #path_to_marker_file  <<- as.character(input$choose_marker_file_text)
           readM$path_to_marker_file  <- as.character(input$choose_marker_file_text)
    })

 



   ## Read marker information
   ##~~~~~~~~~~~~~~~~~~~~~~~~~

#
#  observeEvent({input$choose_marker_file 
#   observeEvent(input$marker_go,  {
#   withProgress(message = 'Loading marker data', value = 1, {
#   
#     if(input$filetype == "vcf"){
#       withCallingHandlers({
#                 shinyjs::html("ReadMarker", "")
#                 if (file.exists(readM$path_to_marker_file) == TRUE) {
#                   geno <<- ReadMarker(filename = readM$path_to_marker_file,  type="vcf", availmemGb = input$memsize, quiet = TRUE)
#                 } else {
#                    shinyjs::html(id = "ReadMarker", html = paste0("ReadMarker", "  File does not exist:", readM$path_to_marker_file))
#              }
#          }, ## end withCallingHandlers
#              message = function(m) {
#                 shinyjs::html(id = "ReadMarker", html = m$message, add = TRUE)
#             })
#
#     }
#
#
# 
#     if(input$filetype == "plink"){
#       withCallingHandlers({
#                 shinyjs::html("ReadMarker", "")
#                 if (file.exists(readM$path_to_marker_file) == TRUE) {
#                   geno <<- ReadMarker(filename = readM$path_to_marker_file, type = "PLINK", availmemGb = input$memsize, quiet = TRUE)
#                 } else {
#                    shinyjs::html(id = "ReadMarker", html = paste0("ReadMarker", "  File does not exist:", readM$path_to_marker_file))
#              }
#          }, ## end withCallingHandlers
#              message = function(m) {
#                 shinyjs::html(id = "ReadMarker", html = m$message, add = TRUE)
#             })
#
#     }
#
#     if(input$filetype == "text"){
#             withCallingHandlers({
#                 shinyjs::html("ReadMarker", "")
#                 aa <- input$AA
#                 ab <- input$AB
#                 bb <- input$BB
#                 missing <- input$missing
#                 if(input$AA=="")
#                     aa <- NULL
#                 if(input$AB=="")
#                     ab <- NULL
#                 if(input$BB=="")
#                     bb <- NULL
#                 if(input$missing=="")
#                     missing <- NULL
#                 if (file.exists(readM$path_to_marker_file) == TRUE) {
#                 geno <<- ReadMarker(filename = readM$path_to_marker_file, type = "text", AA = aa, 
#                            AB = ab  , BB = bb, availmemGb = input$memsize,  quiet = TRUE , missing=missing) 
#                } else {
#                    shinyjs::html(id = "ReadMarker", html = paste0("ReadMarker", "  File does not exist:", readM$path_to_marker_file))
#                } 
#
#              },  ## end withCallingHandlers
#              message = function(m) {
#                 shinyjs::html(id = "ReadMarker", html = m$message, add = TRUE)
#             })
#
#
#     }  ## end if(input$filetype == "text")
#
#
#  })  ## withProgress
#
#  })
#
#
#  })
#


#observeEvent(input$choose_marker_file, {

observeEvent(input$marker_go,  {
   withProgress(message = 'Loading marker data', value = 1, {
   
     if(input$filetype == "vcf"){
       withCallingHandlers({
                 shinyjs::html("ReadMarker", "")
                 if (file.exists(readM$path_to_marker_file) == TRUE) {
                   geno <<- ReadMarker(filename = readM$path_to_marker_file,  type="vcf", availmemGb = input$memsize, quiet = TRUE)
                 } else {
                    shinyjs::html(id = "ReadMarker", html = paste0("ReadMarker", "  File does not exist:", readM$path_to_marker_file))
              }
          }, ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadMarker", html = m$message, add = TRUE)
             })

     }


 
     if(input$filetype == "plink"){
       withCallingHandlers({
                 shinyjs::html("ReadMarker", "")
                 if (file.exists(readM$path_to_marker_file) == TRUE) {
                   geno <<- ReadMarker(filename = readM$path_to_marker_file, type = "PLINK", availmemGb = input$memsize, quiet = TRUE)
                 } else {
                    shinyjs::html(id = "ReadMarker", html = paste0("ReadMarker", "  File does not exist:", readM$path_to_marker_file))
              }
          }, ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadMarker", html = m$message, add = TRUE)
             })

     }

     if(input$filetype == "text"){
             withCallingHandlers({
                 shinyjs::html("ReadMarker", "")
                 aa <- input$AA
                 ab <- input$AB
                 bb <- input$BB
                 missing <- input$missing
                 if(input$AA=="")
                     aa <- NULL
                 if(input$AB=="")
                     ab <- NULL
                 if(input$BB=="")
                     bb <- NULL
                 if(input$missing=="")
                     missing <- NULL
                 if (file.exists(readM$path_to_marker_file) == TRUE) {
                 geno <<- ReadMarker(filename = readM$path_to_marker_file, type = "text", AA = aa, 
                            AB = ab  , BB = bb, availmemGb = input$memsize,  quiet = TRUE , missing=missing) 


                } else {
                    shinyjs::html(id = "ReadMarker", html = paste0("ReadMarker", "  File does not exist:", readM$path_to_marker_file))
                } 

              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadMarker", html = m$message, add = TRUE)
             })


     }  ## end if(input$filetype == "text")


  })  ## withProgress

  }
   ) 
# })  ## end observeEvent






  ##----------------------------------------
  ##  Read phenotypic path and file name
  ##---------------------------------------- 
  ## upload path and file name



        shinyFileChoose(input=input, id='choose_pheno_file', session=session, roots=rootdir )

        observeEvent(input$choose_pheno_file, {
           inFile <- parseFilePaths(roots=rootdir, input$choose_pheno_file)
           updateTextInput(session, "choose_pheno_file_text", value =  as.character(inFile$datapath))
           readP$path_to_pheno_file  <- as.character(inFile$datapath)
           readP$path_to_pheno_file  <- "pathplaceholder"
       })

        observeEvent(input$choose_pheno_file_text, {
           readP$path_to_pheno_file  <<- as.character(input$choose_pheno_file_text)
         })



   ## Read phenotypic  information
   ##~~~~~~~~~~~~~~~~~~~~~~~~~
   pheno <- NULL
   observeEvent(input$pheno_go, {
   withProgress(message = 'Loading phenotypic data', value = 1, {

   header_flag <- FALSE
   if(input$pheno_header == "yes")
      header_flag <- TRUE
   csv_flag <- FALSE
   if(input$pheno_csv == "yes")
      csv_flag <- TRUE

   pheno_missing <- input$pheno_missing
   if(input$pheno_missing=="")
      pheno_missing <- "NA"



   withCallingHandlers({
                shinyjs::html("ReadPheno", "")
                 if (file.exists(readP$path_to_pheno_file) == TRUE) {
                 pheno  <<- ReadPheno(filename = readP$path_to_pheno_file, header=header_flag, csv=csv_flag, missing= pheno_missing)
                 } else {
                    shinyjs::html(id = "ReadPheno", html = paste0("ReadPheno", "File does not exist:", readP$path_to_pheno_file))
                 }
              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadPheno", html = m$message, add = TRUE)
       })



  })  ## end withProgress

  })  ## end observeEvent





  ##------------------------------##
  ## Read Z matrix                ## 
  ##------------------------------ss


  ##----------------------------------------
  ##  Read Z matrix path and file name
  ##---------------------------------------- 
  ## upload path and file name
        shinyFileChoose(input=input, id='choose_Zmat_file', session=session, roots=rootdir )

        observeEvent(input$choose_Zmat_file, {
           inFile <- parseFilePaths(roots=rootdir, input$choose_Zmat_file)
           updateTextInput(session, "choose_Zmat_file_text", value =  as.character(inFile$datapath))
           path_to_Zmat_file  <- as.character(inFile$datapath)
       })
        observeEvent(input$choose_Zmat_file_text, {
           path_to_Zmat_file  <<- as.character(input$choose_Zmat_file_text)
         })



   ## Read Zmat  information
   ##~~~~~~~~~~~~~~~~~~~~~~~~~
   Zmat <- NULL
   observeEvent(input$Zmat_go, {
     withProgress(message = 'Loading Z matrix file', value = 1, {

       withCallingHandlers({
                 shinyjs::html("ReadZmat", "")
         
                 if (file.exists(path_to_Zmat_file) == TRUE) {
                 Zmat  <<- ReadZmat(filename = path_to_Zmat_file)
                 } else {
                    shinyjs::html(id = "ReadZmat", html = paste0("ReadZmat", "File does not exist:", path_to_Zmat_file))
                 }
         
                 # Zmat  <<- ReadZmat(filename = path_to_Zmat_file
              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadZmat", html = m$message, add = TRUE)
       })

  })

  })  ## end observeEvent







  ##----------------------------------------
  ##  Read map path and file name
  ##---------------------------------------- 
  ## upload path and file name
        shinyFileChoose(input=input, id='choose_map_file', session=session, roots=rootdir )

        observeEvent(input$choose_map_file, {
           inFile <- parseFilePaths(roots=rootdir, input$choose_map_file)
           updateTextInput(session, "choose_map_file_text", value =  as.character(inFile$datapath))
           path_to_map_file  <- as.character(inFile$datapath)
       })

        observeEvent(input$choose_map_file_text, {
           path_to_map_file  <<- as.character(input$choose_map_file_text)
         })
   





   ## Read map  information
   ##~~~~~~~~~~~~~~~~~~~~~~~~~
   map <- NULL
   observeEvent(input$map_go, {
    withProgress(message = 'Loading marker map file', value = 1, {


     csv_flag <- FALSE
     if(input$map_csv == "yes")
        csv_flag <- TRUE


   map_header_flag <- FALSE
   if(input$map_header == "yes")
      map_header_flag <- TRUE

       withCallingHandlers({
                 shinyjs::html("ReadMap", "")
         
                 if (file.exists(path_to_map_file) == TRUE) {
                    if(csv_flag){
                       map  <<- ReadMap(filename = path_to_map_file, sep="," , header= map_header_flag)            
                    } else {
                       map  <<- ReadMap(filename = path_to_map_file, header= map_header_flag)            
                    }   ## ed if csv_flag 
                 }  else {
                    shinyjs::html(id = "ReadMap", html = paste0("ReadMap", "File does not exist:", path_to_map_file))
                 }
         
               #   map  <<- ReadMap(filename = path_to_map_file, csv=csv_flag, header= map_header_flag)
         
              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadMap", html = m$message, add = TRUE)
       })




  })
  })  ## end observeEvent




  ##-------------------
  ## Analyse Data
  ##-------------------
  #traitn <- NULL
  ## gets column names of pheno file
  nms <- reactive({
     if(input$pheno_go && input$pheno_header == "yes")
        return(names(pheno))
     if(input$pheno_go && input$pheno_header == "no")
     {  ## pheno file is not named
        nms <- paste("V", 1:ncol(pheno), sep="")
        return(nms)
     } 


     })

  ## get column names for fixed effects after trait has been selected
  fnms <- reactive({
     indx <- NULL
     fixednames <- NULL
     indx <- which(nms()==input$nmst)
      fixednames  <- nms()[-indx]
      return(fixednames )

  })



  ## gets column names length of  pheno file
  sz <- reactive({
     if(input$pheno_go && input$pheno_header == "yes")
     {
        return(length(names(pheno)))
     }
     if(input$pheno_go && input$pheno_header == "no")
     {  ## pheno file is not named
        nms <- paste("V", 1:ncol(pheno), sep="")
        return(length(nms))
     }


     })


sz <- 0
  output$analyse_names <- renderUI({
   if (length(nms()) < 5){
       sz <- length(nms())
   } else {
       sz <- 5
   }



  selectInput(inputId="nmst", label=h4("Step 1: Choose trait"), choices=nms(), size = sz   , selectize=FALSE )     
  })  ## end renderUI








  output$analyse_fnames <- renderUI({
      checkboxGroupInput("nmsf", h4("Step 2: Choose fixed effects"), fnms() , inline=TRUE)
    })  ## end renderUI

  fform <- NULL
   output$fmodel <- renderText({
           fform <<- paste(input$nmsf, collapse="+")
 
   })






 ##  AM analysis for calculation of FPR
# res <- NULL
setlambda <- 1

   #observeEvent((input$analyse_go & input$pheno_go & input$marker_go) , {
 observeEvent(input$marker_go, {
  observeEvent(input$pheno_go, {
   observeEvent(input$analyse_go  , {
   withProgress(message = 'Analysing data', value = 1, {
       fform <<- paste(input$nmsf, collapse="+")

     if(input$analyse_lambda=="manual"){

        withCallingHandlers({
                  shinyjs::html("AM", "")
                  res <<- AM(trait=input$nmst , fformula=fform , 
                             lambda=input$analyse_setlambda,
                             ncpu = input$analyse_cpu,  pheno = pheno, geno=geno, map=map, Zmat=Zmat)
                  setlambda <<- input$analyse_setlambda 
               },  ## end withCallingHandlers
               message = function(m) {
                  shinyjs::html(id = "AM", html = m$message, add = TRUE)
        })



     }


     if(input$analyse_lambda=="auto"){
       
           withCallingHandlers({
                 shinyjs::html("AM", "")
                 print(" testing FPR4AM inputs ")
                 print(  input$analyse_cpu )

 
                 res <<- FPR4AM(numreps = input$analyse_numreps,  falseposrate=input$analyse_fpr,
                            trait=input$nmst , fformula=fform , 
                            ncpu = input$analyse_cpu,  pheno = pheno, geno=geno, map=map, Zmat = Zmat) 
                  if(!is.null(res)){ 
                    setlambda <<- res$setlambda
                    res <<- AM(trait=input$nmst , fformula=fform , 
                             lambda=res$setlambda,
                             ncpu = input$analyse_cpu,  pheno = pheno, geno=geno, map=map, Zmat=Zmat)
                  }

                 },  ## end withCallingHandlers
                    message = function(m) {
                    shinyjs::html(id = "AM", html = m$message  , add = TRUE)
       })  ## withCallingHandlers

     }
  })  ## end withProgress
  })  ## end observeEvent
 })
 })


 ##----------------------------------
 ## Plotting  - server code -
 ##----------------------------------

  #----------------------------------------------------
  ## get chromosome labels

  chromnames <- reactive({
             if (input$analyse_go){
                 # map has been entered
                  if (input$map_go){
                     return(c(unique(map[,2]), "All") )
                  } else {
                    # no map 
                     return(1)
                  }
             }
           })



  sz <- 0
  output$plot_chromosomes <- renderUI({

   if (length(chromnames) < 5){
       sz <- length(chromnames)
   } else {
       sz <- 5
   }
  
  selectInput(inputId="chosenchrm", label=h4("Choose chromosome"), choices=chromnames() , size = sz   , selectize=FALSE )     
  })  ## end renderUI
  #----------------------------------------------------


  #----------------------------------------------------
  ## get number of iterations of model building 
  modelits <- reactive({
                if(input$analyse_go){
                  return(c(1:length(res$Mrk)))
                }
  })


  sz <- 0
  output$plot_modelits <- renderUI({

   if (length(modelits) < 5){
       sz <- length(modelits)
   } else {
       sz <- 5
   }

  selectInput(inputId="chosenits", label=h4("Choose iteration"), choices=modelits() , size = sz   , selectize=FALSE )    
  })  ## end renderUI

  #----------------------------------------------------



  #----------------------------------------------------

   output$plot_choice <- renderUI({
       radioButtons(inputId="plotchoice", label=h4("Choose plot type"),
                       choiceNames=c("Manhattan (-log(p))","Score statistics"),
                       choiceValues=c("Manhattan", "Score") )
   })

  #----------------------------------------------------

   IsItBigger <- function(vals, itnum, xindx=NULL){
   ## calculates percentage change in size (increase or decrease)
   bigger <- NULL
   percentagechange <- NULL
   if(as.numeric(itnum) > 1){
       # entire genome
       bigger <- rep("" , length(vals[[as.numeric(itnum)]]) ) ## initialize 
       percentagechange <- rep(0, length(vals[[as.numeric(itnum)]]) ) ## initialize
       a <-  vals[[as.numeric(itnum)]]
       b <- vals[[as.numeric(itnum) - 1 ]]
       # a > b       
       indx <- which( a  >  b )
       bigger[indx] <- "Increased value"
       percentagechange[indx] <- ( (a  -    b ) )[indx] 
       # a < b
       indx <- which( a <=  b  )
       bigger[indx] <- "Decreased value"
       percentagechange[indx] <- ( (b  -    a ) )[indx] 
      
       if(!is.null(xindx)){
          # reduce to chromosome
          bigger <- bigger[xindx]
          percentagechange <- percentagechange[xindx]
       }
     }
     res <- list(bigger=bigger, percentagechange=percentagechange)
     return(res)
  }


  #---------------------------------------------------
  # ggplot plotting function
 observeEvent(input$plot_go , {
     output$plot <- renderPlot( PlotAM(res, itnum=as.numeric(input$chosenits), chr=input$chosenchrm, type=input$plotchoice    ,  interactive=FALSE))


#     # we do not have a map
#     xindx <- 1:length( res$outlierstat[[as.numeric(input$chosenits)]] )
#     xvals <- xindx
#
#     yvals <- res$outlierstat[[as.numeric(input$chosenits)]]
#     isit  <- IsItBigger(vals=res$outlierstat, itnum=input$chosenits )
#     bigger <- isit$bigger
#     percentagechange <- isit$percentagechange
#     chrm <- rep(1, length(xindx))
#     pos  <- xvals 
# 
#     # map exisits
#     if(!is.null(map)){
#        if(input$chosenchrm != "All"){
#          # picking a single chrm to plot
#          xindx <- which(as.character(map[,2]) == input$chosenchrm)
#          xvals <- map[xindx, ncol(map)]
#          chrm <-  map[xindx, 2]
#          pos <- xvals 
#          yvals <-  res$outlierstat[[as.numeric(input$chosenits)]][xindx]
#          isit <- IsItBigger(vals=res$outlierstat, itnum=input$chosenits, xindx=xindx )
#          bigger <- isit$bigger
#          percentagechange <- isit$percentagechange
#       } else {
#          # plotting all the chromosomes - more difficult
#          # reordering based on chrm then map position
#          oindx <- order(map[,2], map[, ncol(map)])
#          yvals <- res$outlierstat[[as.numeric(input$chosenits)]][oindx]  ## reordering yvals
#
#          if( as.numeric(input$chosenits)  > 1){
#             bigger <- rep(""    , length(  res$outlierstat[[as.numeric(input$chosenits)]][oindx] ) )
#             percentagechange <- rep(0, length(  res$outlierstat[[as.numeric(input$chosenits)]][oindx] ) )
#
#
#            a <-  res$outlierstat[[as.numeric(input$chosenits)]][oindx]
#            b <- res$outlierstat[[as.numeric(input$chosenits) - 1 ]][oindx]
#
#             indx <- which(  a >  b )
#             bigger[indx] <- "Increased value"
#             percentagechange[indx] <-  (( b - a)) [indx]
#
#             indx <- which(  a <=  b )
#             bigger[indx] <- "Decreased value"
#             percentagechange[indx] <-  (( a - b)) [indx]
#
#          }
#
#          mapordered <- map[oindx,]
#          # map position is within chrm, need cumulative postion. 
#          chrms <- unique(mapordered[,2])
#          xvals <- mapordered[, ncol(mapordered)]
#          if (length(chrms) > 1){
#            xvals <- rep(0, nrow(mapordered))
#            indx <- which(mapordered[,2] == chrms[1])
#            xvals[indx] <- mapordered[indx,ncol(mapordered)]
#            genometot <- max(xvals)
#             for(ii in chrms[-1]){
#               indx <- which(mapordered[,2] == ii)
#               xvals[indx] <- mapordered[indx, ncol(mapordered)] + genometot
#               genometot <- max(xvals)
#              }  ## end for
#          } ## end if length(chrms)
#           chrm <- mapordered[,2]
#           pos <- mapordered[, ncol(mapordered)]
#       }  # if else 
#       
#
#     }  ##  if(!is.null(map))
#
#    xlabel <- "Map Position (bp)"
#    if(is.null(map))
#      xlabel <- "Column Position of SNP"
#
#    ylabel <- "Score Statistic"
#    if(input$plotchoice=="manhattan")
#       ylabel <- "-log10(p value)"
#    
#
#     # addition on SNP-trait positions on map
#     if(length(res$Chr) > 1){  ## first entry of list is always NA
#         # found associations 
#          found.chr <- res$Chr[!is.na(res$Chr)]
#          found.pos <- res$Pos[!is.na(res$Pos)]
#          found.label <- 1:length(found.chr)  ## used for annotation in plot
#      }
#
#
#     # place on -lgo10 scale if manhattan selected
#     if(input$plotchoice=="manhattan"){
#       yvals[is.nan(yvals)] <- 0
#       yvals[yvals < 0] <- 0  ## rounding error - very close to 0 when negative
#       ts <- sqrt(yvals)
#       pval <- 1 - pnorm(ts)
#       logp <- -1*log10(pval)
#       yvals <- logp
#     }
#
#
#     # create data frame for plotting 
#     df <- data.frame(xvals=xvals, yvals=yvals, chrm=chrm, pos=pos, foundchr=FALSE, foundpos=FALSE, foundlabel=0 )
#
#     # check for SNP findings from AM
#     if (length(res$Chr)>1){
#          for(ii in 1:length(found.chr)){
#             indx <- which(df$chrm == found.chr[ii])
#             if(!is.null(indx))
#                   df$foundchr[indx] <- TRUE
#             indx <- which(df$pos == found.pos[ii])
#             if(!is.null(indx))
#                   df$foundpos[indx] <- TRUE
#              indx <- which(df$foundchr & df$foundpos & df$foundlabel==0)
#              if(!is.null(indx))
#                   df$foundlabel[indx] <- ii
#       }  ## end for ii
#     }  ## end if length()       
#     geomX <- with(df, xvals[foundchr&foundpos])
#     geomLabels <- with(df, foundlabel[foundchr&foundpos])
#     geomX <- geomX[order(geomLabels)]
#     geomLabels <- geomLabels[order(geomLabels)]
#
#     
#     if(is.null(bigger)){
#       p  <- ggplot(data=df, aes(x=xvals, y=yvals )) + geom_point()
#
#     } else {
# # scale percentagechange to be between 1 and 2. 
#       percentagechange <- 0.25 + ( percentagechange -  min(percentagechange) )/(max(percentagechange) - min(percentagechange))
#
#       p  <- ggplot(data=df, aes(x=xvals, y=yvals , color=bigger ) )  + geom_point( size = percentagechange ) +  scale_color_manual(values=c("#3b5998","#cae1ff" ))
#    }
#
#  p <- p + theme_hc()
#            p <- p + ylab(ylabel) + xlab(xlabel)
#            p <- p +  theme(legend.title=element_blank())  ## no legend title
#            p <- p + theme(legend.position="right")
#
#            if(!is.null(geomX)){
#             if (itnum >1){
#              for(ii in geomX[1:(itnum-1) ] ){
#               yadj <- sample(seq(0.5,0.9,0.1), 1)
#               p <- p + geom_vline(xintercept = ii, linetype="solid", color="#FFE4B5", size=0.5)
#               p <- p + annotate("text", size=8, label=geomLabels[which(geomX==ii)] , 
#                       x=(ii  - ( diff(range(df$xvals))*0.02) )   , 
#                        y = max(df$yvals)*yadj )
#              }  ## end for ii
#
#              if (itnum < (length(geomX)+1) ){
#                ii  <- geomX[itnum]
#                p <- p + geom_vline(xintercept = ii, linetype="solid", color="red", size=0.5)
#              }
#
#             } else {
#
#         ii  <- geomX[1]    
#                p <- p + geom_vline(xintercept = ii, linetype="solid", color="red", size=0.5)
#
#
#
#            }   ## end if itnum > 1
#              p <- p + scale_size(guide='none')
#            } ## if !is.null
#           p <- p + theme( axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), 
#                           axis.title.x=element_text(size=16), axis.title.y=element_text(size=16), 
#                           legend.text=element_text(size=14) )
#           p <- p  + guides(colour = guide_legend(override.aes = list(size=8)))  ## changing point size in legend
#          
#           p <- p + theme(legend.position = 'bottom', legend.spacing.x = unit(0.5, 'cm')) 
#           output$plot <- renderPlot(p)
#
#


           if(input$chosenchrm=="All"){
              # entire genome
               txt2 <- "across all chromosomes"
               if(input$plotchoice=="manhattan"){
                  txt1 <- " -log p value of the score statistic"
                  txt3 <- "-log p value"
               } else {
                 txt1 <- "score statistic"
                  txt3 <- txt1
               } ## inner else
            } else{
               # chrm selected
               txt2 <- paste("on chromosome", input$chosenchrm)
               if(input$plotchoice=="manhattan"){
                  txt1 <- " -log p value of the score statistic"
                  txt3 <- "-log p value"
               } else {
                  txt1 <- "score statistic"
                  txt3 <- txt1

               } ## inner else
            }  ## end outer else
         
           if (input$chosenits==1){
           txt <- paste("Figure: ", 
                        "A plot of the", txt1, 
                        "verse map position", txt2 , 
                        "at the first iteration of the model building process. If there is a red horizontal line, it denotes the position of a new SNP-trait association.") 

          } else { 
           txt <- paste("Figure: ", "A plot of the", txt1, "verse map position", txt2 , "at iteration", input$chosenits, "of the model building process.. The orange horizontal lines denote the position of SNP found by Eagle to be in association with the trait, ", input$nmst, ". A red horizontal line denotes the position of a new SNP-trait association. The numbers are the order in which the SNP-trait associations were found. Purple (green) points denote a", txt3, "that has decreased (increased) in size from the previous iteration. This is useful for inspecting how different parts of the chromosome gain or lose importance as SNP-trait associations are found. The size of the point is proportional to the size of the increase/decrease.")

          } 
           output$caption <- renderText(txt)

  })  ## end observeEvent


 #--------------------------------------------------

 ##--------------------------
 ## Print findings .... 
 ##--------------------------

 ## form data frame of results 
 observeEvent(input$marker_go, {
 observeEvent(input$pheno_go, {
 observeEvent(input$analyse_go, {
 dfparams <- NULL

 if(!is.null(fform))
    dfparams <- data.frame(Parameters=c("Trait", "Fixed effects", "Working memory", "Number CPU", "Gamma"), Settings=c(input$nmst, as.character(fform), input$memsize, input$analyse_cpu, round(setlambda,3)))

 if(is.null(fform))
   dfparams <- data.frame(Parameters=c("Trait", "Fixed effects", "Working memory", "Number CPU", "Gamma"), Settings=c(input$nmst, "overall mean" , input$memsize, input$analyse_cpu, round(setlambda,3)))   

 output$parameters <- renderTable(dfparams) 



 dfres <- NULL
 if(length(res$Mrk)>1){
    dfres <- data.frame(snps=res$Mrk[-1], chrm=res$Chr[-1], position=res$Pos[-1])

 }

 if(is.null(dfres)){
    ## no associations found
    output$findings <- renderText("No significant associations between snp and trait were found.")
  }
  if(!is.null(dfres)){
    output$findings <- renderTable(dfres)
  }

     observeEvent(input$pvalue_go, {
   withProgress(message = ' Calculating additional summary measures', value = 1, {



      withCallingHandlers({
                 shinyjs::html("summary", "")
                  sumres <- SummaryAM(AMobj=res )
                  output$pvalue <- renderTable(sumres[["pvalue"]], digits=-1, hover=TRUE, bordered=TRUE)
                  output$size <- renderTable(sumres[["size"]], digits=-1, hover=TRUE, bordered=TRUE)
                  output$R <- renderTable(sumres[["summarylist"]],  hover=TRUE, bordered=TRUE)


              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "summary", html = m$message, add = TRUE)
       })



  })

  })  ## end observeEvent




}) }) }) ## end observeEvent


##----------------------
## Help - Docs
##---------------------



## FAQ html 
#   output$faq <- renderUI({
#        HTML(markdown::markdownToHTML(knit('faq.rmd', 
#              quiet = TRUE),options=c('toc'), fragment.only=TRUE))
#    })

# observeEvent(input$quickstart, {
#      RShowDoc("QuickStart", package="Eagle")
#})

  ##---------------------------------
  ## Help files   - addPopover
  ##--------------------------------

 
  shinyBS::addPopover(session, "dummy1", "Details", content = HTML("
Eagle can handle three types of marker genotype file; a vcf file, a space separated plain text file and PLINK
ped file. We assume the marker loci are snps. 
Missing marker genotypes are allowed but the 
proportion of missing genotypes is assumed to be low. 
<br><br>
For the vcf file, version 4.0 is assumed where data have been collected on snp genotypes. 
Since vcf files also contain map information, there is no need to load a separate map file as the 
map is extracted from the vcf file.  
<br><br>


The marker genotype file should not contain column names. 
We also assume that each row of data in the file corresponds to data on a different 
individual. The ordering of the rows, by individual, must be the same for the marker genotype file and phenotypic file.<br><br> 


If the file is a plain text file, then character or numeric genotypes can be used to 
denote the snp genotypes. However, Eagle needs 
to map these genotypes to its internal snp genotypes. We do this by asking the user to assign their snp genotype codes to our AA, 
AB and BB codes. If data are collected on inbred individuals, only AA and BB need be specified. 
<br> <br>

To load the marker genotype data into Eagle, follow the four steps.  Upon cliking Upload, Eagle checks the genotype file for errors, 
and recodes the genoytpes for later analysis. If the marker genotype is large (many Gbytes), this step can take several minutes. <br><br>
Output from reading in the marker genotype file will appear in the right hand-side panel. 
                                                           "), trigger = 'hover')

  
  shinyBS::addPopover(session, "dummy2", "Details", content = HTML("  
  Eagle assumes the phenotypic file is either a space separated or comma separated file. The rows correspond to data 
  collected on the individuals. The first row of the file can contain column headings or not.  
  The number of rows of data in the phenotypic file must equal the number of rows in the marker genotype file otherwise an error occurs.  
  Also, Eagle assumes the phenotypic data is row ordered by individual in the same way as the marker genotype data. 
 <br> <br>
Data on multiple traits and fixed effects that may or may not be used in the analysis can be included in this file.  <br> <br>
Missing values are allowed.  <br> <br>
Output from reading in the phenotypic file will appear in the right hand-side panel. 
  "), trigger = 'hover')


shinyBS::addPopover(session, "Zmat1", "Details", content = HTML("
The Z matrix contains only zeros and ones. The number of rows must be greater than the number of columns. 
It is used for those situations where multiple observations of the same trait have been recorded for an individual. 
"), trigger = "hover") 





  
shinyBS::addPopover(session, "dummy3", "Details", content = HTML("
    Eagle does not
     require a known marker map in order to analyse the data.  
     If a map file is read into Eagle, then the
     marker names are used when results are reported in  'Findings'. If a
     map file is not supplied, generic names M1, M2, ..., are
     given to the marker loci. 
      <br> <br>
     The map file can have three or four columns. If the
     map file has three columns, then it is assmed that the three
     columns are the marker locus names, the chromosome number, and the
     map position (in any units). If the map file has four columns as
     with a PLINK map file, then the columns are assumed to be the
     marker locus names, the chromosome number, the map position in
     centimorgans, and the map position in base pairs.
      <br> <br>
     Missing values are not allowed.
      <br> <br>
    The order of the marker loci in this file is assumed to be in the
     same order as the loci (or columns) in the marker data file.
  "), trigger = "hover") 



shinyBS::addPopover(session, "dummy4", "Details", content = HTML(paste("
This page goes through the steps that are needed to analyse the data. 
In the first step, the column in the phenotype file containing the trait data is specified. 
In the second step, any fixed effects are specified. If no fixed effects are selected, the fixed effects part of the model only contains an overall mean. 
In the thrid step, the number of available CPU is set. The default is 1 but if more are available, increasing this number will improve performance significantly. 
The fourth step is to set the lambda parameter. The lambda parameter controls the conservativeness (or false positive rate) of the model building process. For a quick preliminary analysis of the 
data, choose the manual option and leave the parameter at its default setting. For a more detailed analysis of the data where the false positive rate is prespecified, 
choose the auto option. 
<br><br>
To perform the analysis, click on the button in step 5.   Output will start appearing in the right hand panel.  A table of results is given in 'Findings'. 
<br><br>", tags$span(style="color:red",
"Once an analysis has been completed, a new analysis can be performed  
by changing any of the choices in steps 1 to 4 and clicking the 'Perform genome-wide analysis' button.", sep=""))
), trigger = "hover")



shinyBS::addPopover(session, "dummy5", "Details", content = HTML(paste("
By default, the 'best' set of snp in strongest association with the trait are reported. 
These results are given in table form. 
<br><br>
By clicking the 'Additional Summary' button on the left, two additional tables of results are shown; a table on the significance of the snp 
in the model and a table for the amount of phenotypic variance explained as they are added one at a time to the model.
<br><br>
There is additional computation needed to produce these extra tables. It may take a few minutes before these tables appear. 
     ", sep="")


), trigger = "hover")



shinyBS::addPopover(session, "plot_overview", "Details", content = HTML(paste("
Eagle finds SNP-trait associations by building a model iteratively. At each iteration of the model building process, the next 'best' SNP is found. This is done by identifying the SNP with the largest score statistic. A new score statitic is calculated at each iteration of the model building process.<br><br>
Here, the score statistics or their -log p-values can be plotted. A user can see how these score statistics change as the model is built. Red (blue) points mean the score statistic has increased (decreased) from the previous iteration.   
<br><br>
The vertical dotted lines mark the location of the SNP-trait findings. The number is the order in which the SNP-trait associations were found by Eagle. 
     ", sep="") ), trigger = "hover")




 
session$onSessionEnded(stopApp)
  
}



shinyApp(ui=ui, server=server)


