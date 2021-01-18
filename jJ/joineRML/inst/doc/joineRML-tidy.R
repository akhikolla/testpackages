## ----load-joineRML, include = FALSE-------------------------------------------
library(joineRML)
library(knitr)

## ----vignette, eval = FALSE---------------------------------------------------
#  vignette("joineRML", package = "joineRML")
#  help("heart.valve", package = "joineRML")

## ----hvd_data-----------------------------------------------------------------
data(heart.valve)
hvd <- heart.valve[!is.na(heart.valve$grad) & !is.na(heart.valve$lvmi), ]

## ----hvd_data_small-----------------------------------------------------------
hvd <- hvd[hvd$num <= 50, ]

## ----hvd_model_fit------------------------------------------------------------
set.seed(12345)
fit <- mjoint(
  formLongFixed = list(
    "grad" = log.grad ~ time + sex + hs,
    "lvmi" = log.lvmi ~ time + sex
  ),
  formLongRandom = list(
    "grad" = ~ 1 | num,
    "lvmi" = ~ time | num
  ),
  formSurv = Surv(fuyrs, status) ~ age,
  data = list(hvd, hvd),
  timeVar = "time"
)

## ----tidy---------------------------------------------------------------------
tidy(fit)

## ----tidy-long----------------------------------------------------------------
tidy(fit, component = "longitudinal")

## ----tidy-ci------------------------------------------------------------------
tidy(fit, ci = TRUE)
tidy(fit, ci = TRUE, conf.level = 0.99)

## ----tidy-boot, eval = FALSE--------------------------------------------------
#  bSE <- bootSE(fit, nboot = 100, safe.boot = TRUE, progress = FALSE)
#  tidy(fit, boot_se = bSE, conf.int = TRUE)

## ----tidy-plotting, fig.height = 5, fig.width = 5, fig.align = "center"-------
library(ggplot2)
out <- tidy(fit, conf.int = TRUE)
ggplot(out, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_point() +
  geom_errorbar()

## ----augment-nd-grad----------------------------------------------------------
preds <- augment(fit)
head(preds[, c("num", "log.grad", ".fitted_grad_0", ".fitted_grad_1", ".resid_grad_0", ".resid_grad_1")])

## ----augment-nd-lvmi----------------------------------------------------------
head(preds[, c("num", "log.lvmi", ".fitted_lvmi_0", ".fitted_lvmi_1", ".resid_lvmi_0", ".resid_lvmi_1")])

## ----augment-plot-------------------------------------------------------------
out <- preds[preds$num %in% c(26, 36, 227, 244), ]
ggplot(out, aes(x = time, colour = num)) +
  geom_line(aes(y = log.grad, linetype = "Measured")) +
  geom_line(aes(y = .fitted_grad_1, linetype = "Fitted")) +
  labs(linetype = "Type", colour = "ID", y = "Aortic gradient")

## ----glance-fit---------------------------------------------------------------
glance(fit)

## ----glance-fit2--------------------------------------------------------------
set.seed(67890)
fit2 <- mjoint(
  formLongFixed = list(
    "grad" = log.grad ~ time + sex + hs,
    "lvmi" = log.lvmi ~ time + sex
  ),
  formLongRandom = list(
    "grad" = ~ 1 | num,
    "lvmi" = ~ 1 | num
  ),
  formSurv = Surv(fuyrs, status) ~ age,
  data = list(hvd, hvd),
  timeVar = "time"
)

## ----glance-comparison--------------------------------------------------------
glance(fit)
glance(fit2)

## ----vignette-broom, eval = FALSE---------------------------------------------
#  vignette(topic = "broom", package = "broom")

