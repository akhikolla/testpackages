# Libraries are included to keep intellisense quiet
library(shinyjs)
library(shinyAce)
library(shinyBS)

shinyServer(function(input, output, session) {

  getData = reactive({
    # Read and pre-process editor data
    data = input$data
    # Replace multiple spaces or tabs by single tab
    data = str_replace_all(data,"([\t ]+)","\t")
    data = str_replace_all(data,",",".")
    if (nchar(data) < 10) return(NULL)
    tc = textConnection(data)
    d = na.omit(read.table(tc, sep = "\t", header = TRUE))
    close(tc)
    d$record = as.factor(d$record)
    validate(
      need(input$method_a == "stan" || nrow(d) > 20,
           "At least 20 data values required. Try Bayesian method instead"),
      need(input$method_a == "stan" || nlevels(d$record) > 3,
           "At least 3 records required. Try Bayesian method instead.")
    )
    comment = paste(unlist(str_extract_all(data, "^#.*\\n")), collapse = "\n")
    comment = str_replace_all(comment,"\\t", " ")
    comment(d) = comment
    d
  })


  pc = reactive({
    # Compute fit
    d = getData();
    if (is.null(d)) return(NULL)
    if (input$method_a == "nlme") {
      model = eval(parse(text = input$fit_model))
      variant = input$variant
      ng = nlme_gastempt(d, model = model, variant = variant )
      comment(ng) = comment(d)
    } else {
      model_name = stan_models[input$cov_model, input$fit_model]
      ng = stan_gastempt(d, model_name = model_name, lkj = input$lkj,
                         student_df = input$student_df, chains = 1)
      comment(ng) = comment(d)
    }
    ng
  })

  popit = function(session, show, id, title, placement = "right" ){
    if (show) {
      content = pop_content[id]
      if (is.na(pop_content[id]))
        content = ""
      addPopover(session, id, title, content, placement)
    } else {
      removePopover(session, id)
    }
  }

  observe({
    show = input$show_pop
    popit(session, show, "method_a", "Fitting method")
    popit(session, show, "fit_model",  "Modelling curves")
    popit(session, show, "variant",  "Available Variants")
    popit(session, show, "cov_model",  "Covariance of Stan Model")
    popit(session, show, "seed",  "Randomization Seed")
    popit(session, show, "model_s",  "Curves created as")
    popit(session, show, "manual", "Manual set parameter of simulated data")
    popit(session, show, "lkj", "LKJ parameter for covariance")
    popit(session, show, "student_df", "Residual error outliers")
    popit(session, show, "tempt_mean",  "Mean of emptying time constant")
    popit(session, show, "tempt_std_perc",
              "Between-record standard deviation of emptying time constant")
    popit(session, show, "v0_mean", "Mean of initial volume")
    popit(session, show, "v0_std_perc",
          "Between-record standard deviation of initial volume")
    popit(session, show, "kappa_beta_mean",  "Mean of kappa or beta")
    popit(session, show, "kappa_beta_std_perc",
          "Between-record standard deviation of kappa or beta")
    popit(session, show, "student_t_df",  "Type of noise")
    popit(session, show, "noise_perc", "Amplitude of noise")
    popit(session, show, "missing", "Fraction of data missing")
    popit(session, show, "data",  "Entering data", "bottom")
  })

  observe({
    # Create dependency on all numeric fields
    preset  = input$preset
    if (is.null(preset)) return(NULL)
    ss = presets %>% filter(id == preset)
    num_presets = ss[,numcols]
    lapply(seq_along(num_presets), function(i){
      name = names(num_presets)[i]
      updateNumericInput(session, name, value = num_presets[[name]] )
    })
    updateSelectInput(session, "model_s", selected =  ss$model_s)
  })

  observe({
    # Clear ace editor
    if (input$clearButton == 0)
      return(NULL)
    updateAceEditor(session, "data",value = 1)
  })

  observe({
    # Create simulated data
    n_records = input$n_records
    v0_mean = input$v0_mean
    v0_std = input$v0_std_perc*input$v0_mean/100
    tempt_mean = input$tempt_mean
    tempt_std = input$tempt_std_perc*tempt_mean/100
    kappa_mean = input$kappa_beta_mean
    kappa_std = input$kappa_beta_std_perc*kappa_mean/100
    beta_mean = kappa_mean
    beta_std = kappa_std
    noise = input$noise_perc*v0_mean/100.
    student_t_df = as.integer(input$student_t_df)
    missing = as.double(input$missing)/100.
    model_name = input$model_s
    model = eval(parse(text = model_name))
    # Compute simulated data
    d = simulate_gastempt(n_records, v0_mean, v0_std, tempt_mean,
      tempt_std, kappa_mean, kappa_std, beta_mean, beta_std, noise,
      student_t_df, missing, model, seed = input$seed)
    # Copy simulated data to editor
    tc = textConnection("dt","w")
    comment = str_replace_all(comment(d$data),"\\n", " ")
    writeLines(paste0("# ", comment), con = tc)
    suppressWarnings(write.table(d$data, file = tc, append = TRUE,
                row.names = FALSE, sep = "\t", quote = FALSE))
    updateAceEditor(session, "data", value = paste(dt, collapse = "\n") )
    close(tc)
  })

  observe({
    # Update preset popover TODO: does not work reliably
    removePopover(session, "preset")
    addPopover(session, "preset",  "Simulated Sample Data",
               preset_description(input$preset), "right")

  })

  output$residual_plot = renderPlot({
    p = pc()
    if (is.null(p)) return(NULL)
    # Todo: residuals for Bayes
    if (class(p) == "nlme_gastempt") {
      aic = AIC(p$nlme_result)
      max_resid = max(abs(summary(p$nlme)$residuals))+0.2
      plot(p$nlme_result, pch = 16, id = 0.05,
           main = paste("Standardized residuals of fit; AIC =", round(aic)),
           ylim = c(-max_resid, max_resid),
           xlab = "fitted volumes (ml)"
           )
    } else {
      plot(x = 1, main = "Residuals for Bayesian fit not yet implemented",
           ylab = "Just nothing yet", xlab = "Just nothing yet")
    }
  }, height =  500, width = 700) # Make height variable

  output$table = DT::renderDataTable({
    p1 = pc()
    if (is.null(p1)) return(NULL)
    p = coef(p1, signif = 3)
    if (is.null(p)) return(NULL)
    DT::datatable(p, rownames = FALSE, caption = comment(p1),
      options = list(
        paging = FALSE,
        searching = FALSE
      ))
  })

  output$fit_plot = renderPlot({
    p1 = pc()
    if (is.null(p1)) return(NULL)
    p1$plot
  })

  output$download_coef = downloadHandler(
    filename = function() {
      paste('gastempt_', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      p = pc()
      if (is.null(p)) {
        writeLines(paste0("# No valid data"), con = con)
        return(NULL)
      }
      cf = coef(p, signif = 3)
      comment = comment(p)
      if (!is.null(comment) || comment != ""){
        comment = str_replace_all(comment,"\\n", " ")
        print(str(comment))
        writeLines(paste0("# ", comment), con = con)

      }
      suppressWarnings(write.table(cf, file = con, append = TRUE,
               row.names = FALSE, sep = ",", quote = FALSE))
    }
  )

})