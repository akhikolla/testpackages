## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE, warning = FALSE-------------------------------------------
#  install.packages("Rdca")

## ---- fig.width = 6, fig.height= 8, fig.align = "center", warning = TRUE------
library(Rdca)
library(magrittr)
library(ggplot2)
library(ggpubr)

dcl_param_exp <- decline_param(input_unit = "Field", output_unit = "Field", fluid = "oil", 
                               model = "exponential", qi = 1000, Di = 0.0015, b = 0, q_abnd = NULL)

dcl_param_exp

decline_time_exp <- decline_time(c(1:7300), unit = "day")   

str(decline_time_exp)

decline_predict_exp <- decline_predict(dcl_param_exp, decline_time_exp)

head(decline_predict_exp, 10)

p1 <- decline_predict_exp %>% ggplot(aes(x = `Time_(day)`, y = `q_(bbl/day)`)) + 
  geom_point(color = "green4") +  
  theme_bw()

p2 <- decline_predict_exp %>% ggplot(aes(x = `Time_(day)`, y = `Q_(bbl)`)) + 
  geom_point(color = "green4") +
  theme_bw()

exp_plots <- ggarrange(p1, p2, ncol = 1, nrow = 2, align = "v")

exp_plots


## ---- fig.width = 6, fig.height= 8, fig.align = "center", warning = TRUE------
library(Rdca)
library(ggplot2)
library(ggpubr)

dcl_param_exp <- decline_param(input_unit = "Field", output_unit = "SI", fluid = "gas", 
                               model = "exponential", qi = 75000, Di = 0.03, b = 0, q_abnd = 1000)

dcl_param_exp

decline_time_exp <- decline_time(c(1:360), unit = "month")   

str(decline_time_exp)

decline_predict_exp <- decline_predict(dcl_param_exp, decline_time_exp)

head(decline_predict_exp, 10)

time_abnd <- decline_predict_exp$`time_abnd_(months)`[1]
EUR <- decline_predict_exp$`EUR_(m3)`[1]

p1 <- decline_predict_exp %>% ggplot(aes(x = `Time_(month)`, y = `q_(m3/month)`)) +
  geom_point(color = "red") +
  geom_vline(aes(xintercept = `time_abnd_(months)`), linetype = 2) + 
  annotate(geom = "text", x = time_abnd + 10, y = 5e5, label = "time_abnd", angle = 90, 
           color = "blue") +
  geom_hline(aes(yintercept = 1e6 / 35.3147), linetype = 2) + 
  annotate(geom = "text", x = 50, y = 3e6 / 35.3147, label = "rate_abnd", color = "blue") +
  theme_bw()

p2 <- decline_predict_exp %>% ggplot(aes(x = `Time_(month)`, y = `Q_(m3)`)) +
  geom_point(color = "red") +
  geom_vline(aes(xintercept = `time_abnd_(months)`), linetype = 2) + 
  annotate(geom = "text", x = time_abnd + 10, y = 6e7, label = "time_abnd", angle = 90, 
           color = "blue") +
  geom_hline(aes(yintercept = `EUR_(m3)`), linetype = 2) + 
  annotate(geom = "text", x = 50, y = 1.02 * EUR, label = "EUR", color = "blue") +
  theme_bw()

exp_plots <- ggarrange(p1, p2, ncol = 1, nrow = 2, align = "v")

exp_plots


## ---- fig.width = 6, fig.height= 8, fig.align = "center", warning = TRUE------
library(Rdca)
library(ggplot2)

dcl_time_exp <- decline_time(1:1000, unit = "day")

set.seed(123)

prod_data <- 3000 * exp(-0.00234 * dcl_time_exp$t) + 50 * rnorm(1000)

field_data <- data.frame(time = dcl_time_exp$t, q = prod_data)

dcl_fit_param_exp <- decline_fit_param(input_unit = "Field", output_unit = "Field", 
                                       fluid = "oil", model = "exponential", 
                                       fit_data = "rate", prod_data = prod_data, 
                                       initial_param = c(3000, 0.001, 0), lower = NULL, 
                                       upper = NULL)

tibble::glimpse(dcl_fit_param_exp)

dcl_fit_exp <- decline_fit(dcl_fit_param_exp, dcl_time_exp)

tibble::glimpse(dcl_fit_exp)

names(attr(dcl_fit_exp, which = "nls.out"))

attr(dcl_fit_exp, which = "nls.out")$par

attr(dcl_fit_exp, which = "nls.out")$info

attr(dcl_fit_exp, which = "nls.out")$niter

attr(dcl_fit_exp, which = "nls.out")$deviance

dcl_predict_exp <- decline_predict(dcl_fit_exp, dcl_time_exp)

field_data %>% ggplot(aes(x = time, y = q)) +
  geom_point(color = "blue", shape = 21, size = 3) +
  geom_line(aes(x = `Time_(day)`, y = `q_(bbl/day)`), data = dcl_predict_exp,
            color = "red", size = 1) +
  labs(x = "Time (days)", y = "Rate (bbl/day)") + 
  theme_bw()


## ---- fig.width = 6, fig.height= 8, fig.align = "center", warning = TRUE------
library(Rdca)
library(ggplot2)
library(ggpubr)

dcl_param_harm <- decline_param(input_unit = "SI", output_unit = "SI", fluid = "oil", 
                               model = "harmonic", qi = 1000, Di = 0.075, b = 1, q_abnd = 50)

dcl_param_harm

decline_time_harm <- decline_time(c(1:360), unit = "month")   

str(decline_time_harm)

decline_predict_harm <- decline_predict(dcl_param_harm, decline_time_harm)

head(decline_predict_harm, 10)

time_abnd <- decline_predict_harm$`time_abnd_(months)`[1]

EUR <- decline_predict_harm$`EUR_(m3)`[1]

p1 <- decline_predict_harm %>% ggplot(aes(x = `Time_(month)`, y = `q_(m3/month)`)) + 
   geom_point(color = "green4") +
  geom_vline(aes(xintercept = `time_abnd_(months)`), linetype = 2) + 
  annotate(geom = "text", x = time_abnd + 10, y = 250, label = "time_abnd", angle = 90, 
           color = "blue") +
  geom_hline(aes(yintercept = 50), linetype = 2) + 
  annotate(geom = "text", x = 50, y = 80, label = "rate_abnd", color = "blue") +
  theme_bw()

p2 <- decline_predict_harm %>% ggplot(aes(x = `Time_(month)`, y = `Q_(m3)`)) + 
  geom_point(color = "green4") +
  geom_vline(aes(xintercept = `time_abnd_(months)`), linetype = 2) + 
  annotate(geom = "text", x = time_abnd + 10, y = 35000, label = "time_abnd", angle = 90, 
           color = "blue") +
  geom_hline(aes(yintercept = `EUR_(m3)`), linetype = 2) + 
  annotate(geom = "text", x = 50, y = 1.04 * EUR, label = "EUR", color = "blue") +
  theme_bw()

harm_plots <- ggarrange(p1, p2, ncol = 1, nrow = 2, align = "v")

harm_plots


## ---- fig.width = 6, fig.height= 8, fig.align = "center", warning = TRUE------
library(Rdca)
library(ggplot2)

dcl_time_harm <- decline_time(1:360, unit = "month")

set.seed(1234)

prod_data <- 30000 / (1 + 0.02 * dcl_time_harm$t) + 500 * rnorm(360)   # rate

field_data <- data.frame(time = dcl_time_harm$t, q = prod_data)

field_data$Q <- cumsum(field_data$q)                              # cumulative

dcl_fit_param_harm <- decline_fit_param(input_unit = "SI", output_unit = "SI", fluid = "gas", 
                              model = "harmonic", fit_data = "cum", prod_data = field_data$Q, 
                              initial_param = c(40000, 0.01, 1), lower = NULL, upper = NULL,
                              control = list(maxiter = 100))

tibble::glimpse(dcl_fit_param_harm)

dcl_fit_harm <- decline_fit(dcl_fit_param_harm, dcl_time_harm)

tibble::glimpse(dcl_fit_harm)

names(attr(dcl_fit_harm, which = "nls.out"))

attr(dcl_fit_harm, which = "nls.out")$par

attr(dcl_fit_harm, which = "nls.out")$info

attr(dcl_fit_harm, which = "nls.out")$niter

attr(dcl_fit_harm, which = "nls.out")$deviance

dcl_predict_harm <- decline_predict(dcl_fit_harm, dcl_time_harm)

p_cum <- field_data %>% ggplot(aes(x = time, y = Q)) +
  geom_point(color = "blue", shape = 21, size = 3) +
  geom_line(aes(x = `Time_(month)`, y = `Q_(m3)`), data = dcl_predict_harm, 
            color = "red", size = 1) +
  labs(x = "Time (months)", y = "Cumulative Production (m3)") + 
  theme_bw()

p_cum

p_rate <- field_data %>% ggplot(aes(x = time, y = q)) +
  geom_point(color = "blue", shape = 21, size = 3) +
  geom_line(aes(x = `Time_(month)`, y = `q_(m3/month)`), data = dcl_predict_harm, 
            color = "red", size = 1) +
  labs(x = "Time (months)", y = "Production Rate (m3/month)") + 
  theme_bw()

p_rate


## ---- fig.width = 6, fig.height= 8, fig.align = "center", warning = TRUE------
library(Rdca)
library(ggplot2)
library(ggpubr)

dcl_param_hyp <- decline_param(input_unit = "Field", output_unit = "Field", fluid = "gas", 
                               model = "hyperbolic", qi = 100000, Di = 0.0055, b = 0.85, 
                               q_abnd = 2000)

dcl_param_hyp

decline_time_hyp <- decline_time(seq(as.Date("2000/1/1"), as.Date("2030/12/31"), "days"), 
                                 unit = "date")   

str(decline_time_hyp)

decline_predict_hyp <- decline_predict(dcl_param_hyp, decline_time_hyp)

head(decline_predict_hyp, 10)

time_abnd <- decline_predict_hyp$`time_abnd_(days)`[1]

EUR <- decline_predict_hyp$`EUR_(MMSCF)`[1]

p1 <- decline_predict_hyp %>% ggplot(aes(x = `Time_(day)`, y = `q_(MSCF/day)`)) + 
   geom_point(color = "red") +
  geom_vline(aes(xintercept = `time_abnd_(days)`), linetype = 2) + 
  annotate(geom = "text", x = time_abnd + 200, y = 25000, label = "time_abnd", angle = 90, 
           color = "blue") +
  geom_hline(aes(yintercept = 2000), linetype = 2) + 
  annotate(geom = "text", x = 7500, y = 5000, label = "rate_abnd", color = "blue") +
  theme_bw()

p2 <- decline_predict_hyp %>% ggplot(aes(x = `Time_(day)`, y = `Q_(MMSCF)`)) + 
  geom_point(color = "red") +
  geom_vline(aes(xintercept = `time_abnd_(days)`), linetype = 2) + 
  annotate(geom = "text", x = time_abnd + 200, y = 4e4, label = "time_abnd", angle = 90, 
           color = "blue") +
  geom_hline(aes(yintercept = `EUR_(MMSCF)`), linetype = 2) + 
  annotate(geom = "text", x = 8000, y = 1.03 * EUR, label = "EUR", color = "blue") +
  theme_bw()

hyp_plots <- ggarrange(p1, p2, ncol = 1, nrow = 2, align = "v")

hyp_plots


## ---- fig.width = 6, fig.height= 8, fig.align = "center", warning = TRUE------
library(Rdca)
library(ggplot2)

dcl_time_hyp <- decline_time(1:10000, unit = "day")

set.seed(321)

prod_data <- 4500 / (1 + 0.002 * 0.834 * dcl_time_hyp$t) ^ (1 / 0.834) + 
  25 * rnorm(10000)   # rate

field_data <- data.frame(time = dcl_time_hyp$t, q = prod_data)

dcl_fit_param_hyp <- decline_fit_param(input_unit = "Field", output_unit = "Field", fluid = "gas", 
                              model = "hyperbolic", fit_data = "rate", prod_data = prod_data, 
                              initial_param = c(1000, 0.01, 1.0), lower = c(0, 1e-6, 1e-6),
                              upper = NULL, control = list(maxiter = 100))

tibble::glimpse(dcl_fit_param_hyp)

dcl_fit_hyp <- decline_fit(dcl_fit_param_hyp, dcl_time_hyp)

tibble::glimpse(dcl_fit_hyp)

names(attr(dcl_fit_hyp, which = "nls.out"))

attr(dcl_fit_hyp, which = "nls.out")$par

attr(dcl_fit_hyp, which = "nls.out")$info

attr(dcl_fit_hyp, which = "nls.out")$niter

attr(dcl_fit_hyp, which = "nls.out")$deviance

dcl_predict_hyp <- decline_predict(dcl_fit_hyp, dcl_time_hyp)

field_data %>% ggplot(aes(x = time, y = q)) +
  geom_point(color = "blue", shape = 21, size = 3) +
  geom_line(aes(x = `Time_(day)`, y = `q_(MSCF/day)`), data = dcl_predict_hyp,
            color = "red", size = 1) +
  labs(x = "Time (days)", y = "Rate (MSCF/day)") + 
  theme_bw()


## ---- fig.width = 6, fig.height= 8, fig.align = "center", warning = TRUE------
library(Rdca)
library(ggplot2)
library(ggpubr)

dcl_param_mod_hyp <- decline_param(input_unit = "SI", output_unit = "Field", fluid = "oil", 
                               model = "modified_hyperbolic", qi = 1000, Di = 0.0055, b = 0.85, 
                               Dt = 0.0005, q_abnd = 5)
dcl_param_mod_hyp

decline_time_mod_hyp <- decline_time(seq(as.Date("2000/1/1"), as.Date("2030/12/31"), "days"), 
                                 unit = "date")   

str(decline_time_mod_hyp)

decline_predict_mod_hyp <- decline_predict(dcl_param_mod_hyp, decline_time_mod_hyp)

head(decline_predict_mod_hyp, 10)

time_abnd <- decline_predict_mod_hyp$`time_abnd_(days)`[1]

EUR <- decline_predict_mod_hyp$`EUR_(bbl)`[1]

p1 <- decline_predict_mod_hyp %>% ggplot(aes(x = `Time_(day)`, y = `q_(bbl/day)`)) + 
   geom_point(color = "red") +
  geom_vline(aes(xintercept = `time_abnd_(days)`), linetype = 2) + 
  annotate(geom = "text", x = time_abnd + 200, y = 1000, label = "time_abnd", angle = 90, 
           color = "blue") +
  geom_hline(aes(yintercept = 5 * 6.289814), linetype = 2) + 
  annotate(geom = "text", x = 8000, y = 200, label = "rate_abnd", color = "blue") +
  theme_bw()

p2 <- decline_predict_mod_hyp %>% ggplot(aes(x = `Time_(day)`, y = `Q_(bbl)`)) + 
  geom_point(color = "red") +
  geom_vline(aes(xintercept = `time_abnd_(days)`), linetype = 2) + 
  annotate(geom = "text", x = time_abnd + 200, y = 3e6, label = "time_abnd", angle = 90, 
           color = "blue") +
  geom_hline(aes(yintercept = `EUR_(bbl)`), linetype = 2) + 
  annotate(geom = "text", x = 8000, y = 1.03 * EUR, label = "EUR", color = "blue") +
  theme_bw()

hyp_plots <- ggarrange(p1, p2, ncol = 1, nrow = 2, align = "v")

hyp_plots


## ---- fig.width = 6, fig.height= 8, fig.align = "center", warning = TRUE------
library(Rdca)
library(ggplot2)
library(magrittr)

dcl_time_mod_hyp <- decline_time(1:300, unit = "month")

dcl_param_mod_hyp <- decline_param(input_unit = "Field", output_unit = "Field", fluid = "oil", 
                               model = "modified_hyperbolic", qi = 8000, Di = 0.04, b = 0.75, 
                               Dt = 0.01, q_abnd = 10)

set.seed(4321)

prod_data <- decline_predict(dcl_param_mod_hyp, dcl_time_mod_hyp)$`q_(bbl/month)` +
  rnorm(300, mean = 100, sd = 50)   # rate

field_data <- data.frame(time = dcl_time_mod_hyp$t, q = prod_data, Q = cumsum(prod_data))

dcl_fit_param_mod_hyp <- decline_fit_param(input_unit = "Field", output_unit = "Field", fluid = "oil", 
                              model = "modified_hyperbolic", fit_data = "cum", prod_data = field_data$Q, 
                              initial_param = c(10000, 0.1, 1.0, 0.01), lower = NULL,
                              upper = NULL, control = list(maxiter = 100))

tibble::glimpse(dcl_fit_param_mod_hyp)

dcl_fit_mod_hyp <- decline_fit(dcl_fit_param_mod_hyp, dcl_time_mod_hyp)

tibble::glimpse(dcl_fit_mod_hyp)

names(attr(dcl_fit_mod_hyp, which = "nls.out"))

attr(dcl_fit_mod_hyp, which = "nls.out")$par

attr(dcl_fit_mod_hyp, which = "nls.out")$info

attr(dcl_fit_mod_hyp, which = "nls.out")$niter

attr(dcl_fit_mod_hyp, which = "nls.out")$deviance

dcl_predict_mod_hyp <- decline_predict(dcl_fit_mod_hyp, dcl_time_mod_hyp)

p_cum <- field_data %>% ggplot(aes(x = time, y = Q)) +
  geom_point(color = "blue", shape = 21, size = 3) +
  geom_line(aes(x = `Time_(month)`, y = `Q_(bbl)`), data = dcl_predict_mod_hyp, 
            color = "red", size = 1) +
  labs(x = "Time (months)", y = "Cumulative Production (bbl)") + 
  theme_bw()

p_cum

p_rate <- field_data %>% ggplot(aes(x = time, y = q)) +
  geom_point(color = "blue", shape = 21, size = 3) +
  geom_line(aes(x = `Time_(month)`, y = `q_(bbl/month)`), data = dcl_predict_mod_hyp, 
            color = "red", size = 1) +
  labs(x = "Time (months)", y = "Production Rate (bbl/month)") + 
  theme_bw()

p_rate


