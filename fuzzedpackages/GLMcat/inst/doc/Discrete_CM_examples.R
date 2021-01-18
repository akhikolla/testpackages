## -----------------------------------------------------------------------------
# devtools::load_all()
library(GLMcat)

## -----------------------------------------------------------------------------
data("TravelChoice")
head(TravelChoice)
str(TravelChoice)

## -----------------------------------------------------------------------------
exp_8.4 <- Discrete_CM(
  formula = choice ~ hinc + gc + invt,
  case_id = "indv",
  alternatives = "mode",
  reference = "air",
  data = TravelChoice,
  alternative_specific = c("gc", "invt"),
  distribution = "logistic")
summary(exp_8.4)

## -----------------------------------------------------------------------------
(constant_model <- Discrete_CM(
  formula = choice ~ 1 ,
  case_id = "indv",
  alternatives = "mode",
  reference = c("air", "train", "bus", "car"),
  data = TravelChoice,
  distribution = "logistic"
))

(car_0 <- Discrete_CM(
  formula = choice ~ hinc[air] + psize[air] + gc + ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = c("air", "train", "bus", "car"),
  alternative_specific = c("gc", "ttme"),
  data = TravelChoice,
  distribution = "logistic"
))

## -----------------------------------------------------------------------------
mod_1 <- Discrete_CM(
  formula = choice ~ hinc[air] + psize[air] + gc + ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "air",
  alternative_specific = c("gc", "ttme"),
  data = TravelChoice,
  distribution = "logistic"
)
logLik(mod_1)

## -----------------------------------------------------------------------------
mod_2 <- Discrete_CM(
  formula = choice ~ hinc[air] + psize[air] + gc + ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "bus",
  alternative_specific = c("gc", "ttme"),
  data = TravelChoice,
  distribution = "student",
  freedom_degrees = 30
)
logLik(mod_2)

## -----------------------------------------------------------------------------
mod_3 <- Discrete_CM(
  formula = choice ~ hinc[air] + psize[air] + gc + ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "car",
  alternative_specific = c("gc", "ttme"),
  data = TravelChoice,
  distribution = "student",
  freedom_degrees = 0.2
)
logLik(mod_3)

## -----------------------------------------------------------------------------
mod_4 <- Discrete_CM(
  formula = choice ~ hinc[air] + psize[air] + gc + ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "train",
  alternative_specific = c("gc", "ttme"),
  data = TravelChoice,
  distribution = "student",
  freedom_degrees = 1.35
)
logLik(mod_4)

