forecast <img src="man/figures/logo.png" align="right" />
======================

![R build status](https://github.com/robjhyndman/forecast/workflows/R-CMD-check/badge.svg)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/forecast)](https://cran.r-project.org/package=forecast)
[![cran checks](https://cranchecks.info/badges/worst/forecast)](https://cran.r-project.org/web/checks/check_results_forecast.html)
[![Lifecycle: retired](https://img.shields.io/badge/lifecycle-retired-orange.svg)](https://www.tidyverse.org/lifecycle/#retired)
[![Downloads](https://cranlogs.r-pkg.org/badges/forecast)](https://cran.r-project.org/package=forecast)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)


The R package *forecast* provides methods and tools for displaying and analysing univariate time series forecasts including exponential smoothing via state space models and automatic ARIMA modelling.

This package is now retired in favour of the [fable](http://fable.tidyverts.org/) package. The forecast package will remain in its current state, and maintained with bug fixes only. For the latest features and development, we recommend forecasting with the fable package.

## Installation
You can install the **stable** version from
[CRAN](https://cran.r-project.org/package=forecast).

```s
install.packages('forecast', dependencies = TRUE)
```

You can install the **development** version from
[Github](https://github.com/robjhyndman/forecast)

```s
# install.packages("remotes")
remotes::install_github("robjhyndman/forecast")
```

## Usage

```s
library(forecast)
library(ggplot2)

# ETS forecasts
USAccDeaths %>%
  ets() %>%
  forecast() %>%
  autoplot()

# Automatic ARIMA forecasts
WWWusage %>%
  auto.arima() %>%
  forecast(h=20) %>%
  autoplot()

# ARFIMA forecasts
library(fracdiff)
x <- fracdiff.sim( 100, ma=-.4, d=.3)$series
arfima(x) %>%
  forecast(h=30) %>%
  autoplot()

# Forecasting with STL
USAccDeaths %>%
  stlm(modelfunction=ar) %>%
  forecast(h=36) %>%
  autoplot()

AirPassengers %>%
  stlf(lambda=0) %>%
  autoplot()

USAccDeaths %>%
  stl(s.window='periodic') %>%
  forecast() %>%
  autoplot()

# TBATS forecasts
USAccDeaths %>%
  tbats() %>%
  forecast() %>%
  autoplot()

taylor %>%
  tbats() %>%
  forecast() %>%
  autoplot()
```

## For more information

  * Get started in forecasting with the online textbook at http://OTexts.org/fpp2/
  * Read the Hyndsight blog at https://robjhyndman.com/hyndsight/
  * Ask forecasting questions on http://stats.stackexchange.com/tags/forecasting
  * Ask R questions on http://stackoverflow.com/tags/forecasting+r
  * Join the International Institute of Forecasters: http://forecasters.org/

## License

This package is free and open source software, licensed under GPL-3.
