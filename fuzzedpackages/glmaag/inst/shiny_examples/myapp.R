#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

ui <- fluidPage(

  # App title ----
  titlePanel("Uploading Files"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      fileInput('xtrain', 'X train',
                multiple = T, accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
      fileInput(inputId = 'ytrain', label = 'Y train',
                multiple = T, accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
      fileInput(inputId = 'xtest', label = 'X test',
                multiple = T, accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
      fileInput(inputId = 'ytest', label = 'Y test',
                multiple = T, accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
      fileInput(inputId = 'L', label = 'Laplacian matrix (optional)',
                multiple = T, accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
      selectInput(inputId = 'fam', label = 'Phenotype',
                  choices = c('Gaussian', 'Logistic', 'Cox'), selected = 'Gaussian'),
      selectInput(inputId = 'est', label = 'Estimated network or Elastic Net?',
                  choices = c('Estimated', 'Elastic Net'), selected = 'Estimated'),
      selectInput(inputId = 'cvss', label = 'Cross validation or stability selection?',
                  choices = c('cv', 'ss'), selected = 'cross validation'),
      selectInput(inputId = 'dev', label = 'Cross validation by deviance?',
                  choices = c('deviance', 'robust'), selected = 'robust'),
      selectInput(inputId = 'max', label = 'Cross validation by max or 1se?',
                  choices = c('max', '1se'), selected = '1se'),
      checkboxInput('tune', 'Tune the network (input network required)?', F),
      checkboxInput('al1', 'Adapt l1?', T),
      checkboxInput('al2', 'Adapt quadratic?', T),
      checkboxInput('parallel', 'parallel', F),
      sliderInput(inputId = 'nfolds', label = 'Number of folds for cross validation',
                  min = 3, max = 10, value = 5),
      sliderInput(inputId = 'nsam', label = 'Number of samples for stability selection',
                  min = 100, max = 500, value = 100),
      sliderInput(inputId = 'beta', label = 'Cut off for stability selection',
                  min = 0, max = 1, value = .2),
      actionButton('submit', 'Ready? Go!')
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Data file ----
      h4('Plot'),
      plotOutput('predplot'),
      h4('Goodness of prediction'),
      tableOutput('predtab'),
      h4('Download'),
      downloadButton('coefficients.csv', 'Download coefficients including intercept'),
      downloadButton('features.csv', 'Download selected features')
    )

  )
)

# Define server logic to read selected file ----
server <- function(input, output) {

  model <- eventReactive(input$submit, {
    req(input$xtrain)
    req(input$xtest)
    req(input$ytrain)
    req(input$ytest)
    x_train <- as.matrix(fread(input$xtrain$datapath))
    x_test <- as.matrix(fread(input$xtest$datapath))
    varnam <- colnames(x_train)
    if (input$fam == 'Cox') {
      y_train <- Surv(as.matrix(fread(input$ytrain$datapath)))
      y_test <- as.matrix(fread(input$ytest$datapath))
    } else {
      y_train <- fread(input$ytrain$datapath, data.table = F)[, 1]
      y_test <- fread(input$ytest$datapath, data.table = F)[, 1]
    }

    if(!is.null(input$L)) {
      L <- as.matrix(fread(input$L$datapath))
    } else {
      if (input$est == 'Estimated') {
        L <- getS(x_train)
      } else {
        L <- diag(ncol(x_train))
      }
    }
    meas <- input$dev == 'deviance'
    type1se <- input$max == '1se'

    mod <- switch (input$cvss,
                   cv = cv_glmaag(y_train, x_train, L, fam = input$fam, nfolds = input$nfolds, measdev = meas, type1se = type1se, adaptl1 = input$al1, adaptl2 = input$al2, parallel = input$parallel),
                   ss = ss_glmaag(y_train, x_train, L, fam = input$fam, nfolds = input$nfolds, nsam = input$nsam, beta = input$beta, type1se = type1se, adaptl1 = input$al1, adaptl2 = input$al2, parallel = input$parallel)
    )

    ypre <- predict(mod, x_test, type = 'response')

    if (input$fam == 'Gaussian') {
      plotobj <- ggplot2::qplot(y_test, ypre) + ggplot2::xlab('true value') + ggplot2::ylab('predicted value') + ggplot2::geom_abline(intercept = 0, slope = 1, color = 'grey') + ggplot2::geom_smooth(method = 'lm') + ggplot2::theme_bw()
    } else if(input$fam == 'Logistic') {
      aucp <- pROC::auc(y_test ~ ypre)
      datroc <- data.frame(test = y_test, pre = ypre)
      plotobj <- ggplot2::ggplot(datroc, ggplot2::aes(d = test, m = pre)) + plotROC::geom_roc(n.cuts = 0) + plotROC::style_roc() + ggplot2::annotate("text", x = .5, y = .5, label = paste("AUC", round(aucp, 3), sep = ' = '))
    } else {
      datsurv <- data.frame(train = y_train, pre = predict(mod, x_train, type = 'response'))
      cutp <- maxstat::maxstat.test(train ~ pre, smethod = 'LogRank', datsurv)$estimate
      test1 <- y_test[, 1]
      test2 <- y_test[, 2]
      group <- as.numeric(ypre <= cutp)
      datsurvtest <- data.table(test1 = test1, test2 = test2, group = group)
      plotobj <- survminer::ggsurvplot(survminer::surv_fit(Surv(test1, test2) ~ group, data = datsurvtest), risk.table = T, conf.int = T, legend.labs = c('high risk', 'low risk'), pval = T, ggtheme = ggplot2::theme_bw())
    }
    datab <- switch (input$fam,
                     Gaussian = data.frame(t(evaluate(ypre, y_test))),
                     Logistic = data.frame(t(evaluate(ypre, y_test, getcut(predict(mod, type = 'response'), y_train, 'Logistic'), 'Logistic'))),
                     Cox <- data.frame(Concodimrdance = survConcordance(Surv(y_test[, 1], y_test[, 2]) ~ ypre)$concordance)
    )
    if (input$cvss == 'ss') {
      coefmod <- coef(mod)
    } else {
      coefmod <- coef(mod, type1se = type1se)
    }
    if (input$fam == 'Cox') {
      feasel <- coefmod != 0
    } else {
      feasel <- coefmod[-1] != 0
    }
    datab <- cbind(datab, n_para = sum(feasel))
    feature <- varnam[feasel]
    if (input$fam == 'Cox') {
      names(coefmod) <- varnam
    } else {
      names(coefmod) <- c('Intercept', varnam)
    }
    list(plotobj = plotobj, datab = datab, coefs = data.table(t(coefmod)), selvar = data.table(feature))
  })

  output$predplot <- renderPlot({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    model()$plotobj
  })

  output$predtab <- renderTable({
    model()$datab
  })

  output$coefficients.csv <- downloadHandler(
    filename = function() {
      paste('coefficients.csv')
    },
    content = function(filename) {
      fwrite(model()$coefs, filename)
    }
  )

  output$features.csv <- downloadHandler(
    filename = function() {
      paste('features.csv')
    },
    content = function(filename) {
      fwrite(model()$selvar, filename)
    }
  )
}
# Run the application
shinyApp(ui = ui, server = server)

