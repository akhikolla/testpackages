## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=8, fig.height=6
)

## ---- echo = FALSE, results = "asis", message = FALSE-------------------------
library(knitr)

tex2markdown <- function(texstring) {
  writeLines(text = texstring,
             con = myfile <- tempfile(fileext = ".tex"))
  texfile <- pandoc(input = myfile, format = "html")
  cat(readLines(texfile), sep = "\n")
  unlink(c(myfile, texfile))
}

textable <- "
\\begin{table}[h!]
\\centering
\\caption{. Overview of inputs, including input parameters, to be specified by the user of the model. k indicates c, n or p. VBGC = von Bertalanffy growth curve.}
\\begin{tabular}{l l l}
\\hline
Symbol & Description & Unit\\\\
\\hline
$a_\\textrm{k}$  & Element-specific assimilation efficiency & \\_ \\\\
$l_\\textrm{t}$  & Total length of individual & cm \\\\
$linf$  & Asymptotic adult length (VBGC) & cm \\\\
$\\kappa$  & Growth rate parameter (VBGC) & $\\textrm{yr}^{-1}$ \\\\
$t_0$  & Age at settlement (VBGC) & $\\textrm{yr}$ \\\\
$lw_a$  & Parameter length-weight relationship & $\\textrm{g cm}^{-1}$ \\\\
$lw_b$  & Parameter length-weight relationship & \\_\\\\
$Q_\\textrm{k}$  & Element-specific body content percentage & \\% \\\\
$f_\\textrm{0}$  & Metabolic normalisation constant independent of body mass & $\\textrm{g C} \\textrm{g}^{-\\alpha} \\textrm{d}^{-1}$ \\\\
$alpha$  & Mass-scaling exponent &  \\_\\\\
$theta$  & Activity scope & \\_ \\\\
$v$  & Environmental temperature & \\textdegree C \\\\
$h$  & trophic level & \\_ \\\\
$r$  & Aspect ratio of caudal fin & \\_ \\\\
$F0nz$  & Mass-specific turnover rate of N & $\\textrm{g N} \\textrm{g}^{-1} \\textrm{d}^{-1}$ \\\\
$F0pz$  & Mass-specific turnover rate of P & $\\textrm{g P} \\textrm{g}^{-1} \\textrm{d}^{-1}$ \\\\
$mdw$ & Ratio of dry mass and wet mass of fish & \\_ \\\\
$D_\\textrm{k}$  & Elemental stoichiometry of diet & \\% \\\\
\\hline
\\end{tabular}
\\end{table}
"

tex2markdown(textable)

## ---- message=TRUE, echo=TRUE-------------------------------------------------
# example
fishflux::name_errors("Zebrazoma scopas")

## ----message=TRUE,echo=TRUE---------------------------------------------------
# example
fishflux::find_lw("Zebrasoma scopas", mirror = "se")

## ----message=FALSE,echo=TRUE--------------------------------------------------
# example
# The option otolith=TRUE filters out sources that used otoliths for the estimation of growth parameters
fishflux::growth_params("Sargocentron microstoma", otolith = FALSE)

## ----cache=TRUE,results='hide',message=FALSE----------------------------------
# example
zebsco <- fishflux::model_parameters("Zebrasoma scopas", family = "Acanthuridae", temp = 27, mirror = "se")
## Here we set the temperature at 27 degrees as an example, this the average sea temperature in Moorea, French Polynesia

## -----------------------------------------------------------------------------
print(zebsco)

## ---- message=FALSE-----------------------------------------------------------
## load the example parameters for Zebrasoma scopas, a list
param_zebsco <- fishflux::param_zebsco
## Run the model, specifying the target length(s) and the parameter list
model <- fishflux::cnp_model_mcmc(TL = 5:20, param = param_zebsco)

## -----------------------------------------------------------------------------
fishflux::extract(model, c("Fn","Fp"))

## ---- message=FALSE, warning=FALSE--------------------------------------------
## limitation
fishflux::limitation(model)
## Plot one variable:
fishflux::plot_cnp(model,  y = "Fp", x = "tl", probs = c(0.5, 0.8, 0.95))
## Plot multiple variables:
fishflux::plot_cnp(model,  y = c("Fp", "Gp", "Ip", "Wp"), x = "tl", probs = 0.5)

## ---- message=FALSE-----------------------------------------------------------
## General overview:
fishflux::sensitivity(TL = 10, param = param_zebsco, par = c("Dn_sd", "Dp_sd", "Qn_sd", "Qp_sd", "k_sd"), out = c("Fn", "Fp", "Ic"))

