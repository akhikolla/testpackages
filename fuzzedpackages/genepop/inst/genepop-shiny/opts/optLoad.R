GenepopFile <- reactive({
  input$file1
})

output$contents <- DT::renderDataTable({
  inFile <- input$file1
  if (is.null(inFile))
    return(matrix())
  donnees <- read.genepop(inFile$datapath, ncode = 4L, quiet = TRUE)
  DT::datatable(donnees$datas, options = list(pageLength = 10 ))
})


output$distPops <- renderPlot({
  inFile = input$file1
  if (! is.null(inFile)) {
    donnees <- read.genepop(inFile$datapath, ncode = 4L, quiet = TRUE)
    if (length(donnees) > 0)
    {
      df=data.frame(PopSize=donnees$nind.bypop, popnames=donnees$pops)
      p<-ggplot(data=df, aes(x=popnames, y=PopSize)) +  geom_bar(stat="identity",fill="steelblue")+  theme_minimal()
      print(p + ggtitle(inFile$name) + xlab("population") + ylab("sample size"))
    }
  }
})

output$InputFileData <- renderText({
  filein <- GenepopFile()$datapath
  filePath <- toString(filein)
  conn <- file(filePath, open="r")
  linn <-readLines(conn)
  Res <- ""
  for (i in 1:length(linn)){
    tmp <- Res
    Res <- paste(tmp,linn[i],sep = "\n")
  }
  close(conn)
  Res
})