## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PandemicLP)
library(stats)

## -----------------------------------------------------------------------------
regions <- c("PR", "SC","RS")
last_date <- "2020-10-01"
case_type <- "deaths"

## -----------------------------------------------------------------------------
data <- list()
outputs <- list()
preds <- list()
states <- state_list()
# Load precomputed MCMC
download.file("http://github.com/CovidLP/PandemicLP/raw/master/temp/PandemicLP_SumRegions.rda","./PandemicLP_SumRegions.rda")
load("./PandemicLP_SumRegions.rda")
for(i in 1:length(regions)) { 
  if (is.na(match(regions[i],states$state_abb))){
    data[[i]] <- load_covid(country_name=regions[i],last_date=last_date)
  } else {
    data[[i]] <- load_covid(country_name="Brazil", state_name = regions[i],last_date=last_date)
  }
  # Commenting to reduce vignette runtime
  # outputs[[i]] <- pandemic_model(data[[i]],case_type = case_type, covidLPconfig = TRUE)
  preds[[i]] <- posterior_predict(outputs[[i]])
}
# Load precomputed MCMC
download.file("http://github.com/CovidLP/PandemicLP/raw/master/temp/PandemicLP_SumRegions.rda","./PandemicLP_SumRegions.rda")
load("./PandemicLP_SumRegions.rda")

## -----------------------------------------------------------------------------
pred_final <- preds[[1]] # Just to get the class and format
mu_t <- t(preds[[1]]$pastMu)
mu_final <- data.frame(data = preds[[1]]$data$date,mu_t)
names_mu <- names(preds[[1]]$pastMu)

bind_regions <- regions[1]
for (i in 2:length(regions)) {
bind_regions <- paste(bind_regions, "and", regions[i])
}
pred_final$location <- bind_regions

## -----------------------------------------------------------------------------
for (i in 2:length(preds)){
  pred_final$predictive_Long <- pred_final$predictive_Long + preds[[i]]$predictive_Long
  pred_final$predictive_Short <- pred_final$predictive_Short + preds[[i]]$predictive_Short
  pred_final$futMu <- pred_final$futMu + preds[[i]]$futMu
  
  mu_t <- t(preds[[i]]$pastMu)
  mu_2 <- data.frame(data = preds[[i]]$data$date,mu_t)
  names_mu <- c(names_mu,names(preds[[i]]$pastMu))
  mu_final <- rbind(mu_final,mu_2)
  
  # Merging the data -> can't sum directly because the dates may be different
  # Use merge to avoid sum different dates
  pred_final$data <- merge(pred_final$data,preds[[i]]$data, by = "date", all = TRUE)
  pred_final$data[is.na(pred_final$data)] = 0
  pred_final$data$cases.x = pred_final$data$cases.x + pred_final$data$cases.y
  pred_final$data$deaths.x = pred_final$data$deaths.x + pred_final$data$deaths.y
  pred_final$data$new_cases.x = pred_final$data$new_cases.x + pred_final$data$new_cases.y
  pred_final$data$new_deaths.x = pred_final$data$new_deaths.x + pred_final$data$new_deaths.y
  pred_final$data <- pred_final$data[,-c(6:9)]
  names(pred_final$data) <- c("date","cases","deaths","new_cases","new_deaths")
}

# Aggregate the mu
mu_final <- aggregate(. ~ data, data=mu_final, FUN=sum)
mu_final <- mu_final[,-1]
mu_final <- t(mu_final)
names_mu <- unique(names_mu)
colnames(mu_final) <- names_mu
pred_final$pastMu <- mu_final

## ---- fig.asp=2/(sqrt(5)+1), fig.width=8, fig.align='center', warnings = FALSE----
plots <- plot(pred_final,term = "both",summary = FALSE)
plots$long
plots$short

