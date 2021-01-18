library(stringr)
library(ggplot2)
library(gastempt)
suppressPackageStartupMessages(library(dplyr))
library(readxl)
suppressPackageStartupMessages(library(shinyBS))
github_repo = "https://github.com/dmenne/gastempt/blob/master/exec/"
presets = na.omit(
  suppressWarnings(read_excel("gastempt_presets.xlsx", sheet = "gastempt_samples"))) %>%
  mutate(
    id = as.character(id)
  )
numcols = which(sapply(presets, is.numeric))

pop_content = c(
  model_a = "<code>linexp</code>, <b>vol = v0 * (1 + kappa * t / tempt) * exp(-t / tempt):</b><br>Recommended for gastric emptying curves with an initial volume overshoot from secretion. With parameter kappa &gt; 1, there is a maximum after t=0. When all emptying curves start with a steep drop, this model can be difficult to fit.<hr><code>powexp</code>, <b>vol = v0 * exp(-(t / tempt) ^ beta):</b><br>The power exponential function introduced by Elashof et. al. to fit scintigraphic emptying data; this type of data does not have an initial overshoot by design. Compared to the linexp model, fitting powexp is more reliable and rarely fails to converge in the presence of noise and outliers. The power exponential can be useful with MRI data when there is an unusual late phase in emptying.",

  variant = "<b>Variant 1:</b> The most generic assumptions to estimate the per-record parameters is also the most likely to fail to converge. If this variant works, use it. Otherwise, try one of the other variants.<br><b>Variant 2</b>: A slightly more restricted version; sometimes converges when variant 1 fails.<br><b>Variant 3</b>: Parameters beta or kappa, which are most difficult to estimate for each curve individually, are computed once only for all records, so these parameters cannot be tested for between-treatment differences. If you are only interested in a good estimate of t50 and v0 and the other variants do not work, use this method.",

  cov_model = 'Try the simple model "Without covariance" first, it runs faster. Compare with the model "With covariance" when the fit is not satisfactory. The preset "rough powexp scint" provides example data where the more complex fit is markedly superior.',

  data = "Enter data from Excel-clipboard or other tab-separated data here. Column headers must be present and named <code>record, minute, vol</code>.<br>Lines starting with <code>#</code> are comments that will be shown in plots and output files.<br>Avoid editing details in this table because curves are recalculated on every key press; use the Clear button, and paste in the full edited data set from source instead.",

seed = "The seed used to initialize the random generator that produces noise and between-record variance of parameters. With the same parameters and the same seed, the full analyzed data set is exactly reproduced. Try to change the seed to see how with different record sets the fit toggles between convergence failure and success.",

model_s = "The model used to <b>generate</b> the simulated data set; this is different from the model used to <b>analyze</b> the data, which is selected in the top box. When you add your own data from the clipboard, only the top box is used. Normally, one should manually select the same model for simulation and analysis, but using different settings can help understanding the effect of fitting different models to your experimental data. Presets do not automatically switch to the correct role - this is by design to force the user to think what she is doing.",

kappa_beta = "<b>kappa</b> is used for the linexp model; values of kappa &gt; 1 indicate that a maximum after t=0 is present, i.e. an overshoot from secretion. <br><b>beta</b> is used for the powexp model only; values above 1 indicate that the curve starts with an horizontal slope at t=0, but the volume never rised above the initial volume v0.",

student_t_df = "Gaussian noise is the usual nice-behaved noise that never occurs in the real world of medical research. With three variants of Student-t distributed noise, outliers can be generated to test the methods for robustness. ",
fit_model2 = '<b>linexp or powexp?</b> The linexp function can have an initial overshoot to fit gastric emptying curves with secretion.<br><code>vol(t) = v0 * (1 + kappa * t / tempt) * exp(-t / tempt)</code>.<img src="powlinexp.png">',

fit_model = '<b>linexp or powexp?</b><img src="linexp.png"><br>The <b>linexp</b> function introduced by the author of this app can have an initial overshoot to fit gastric emptying curves with secretion.<br><code>vol(t) = v0 * (1 + kappa * t / tempt) * exp(-t / tempt)</code><br>. Fatty meals analyzed by MRI are best represented by this function. The method does not work well when gastric emptying is almost linear or exponential in time; try the Bayesian version with covariance estimation to get stable results.<br><img src="powexp.png"><br>The <b>powexp</b> function historically used by Elashof et al. is strictly montonously decreasing but has more freedom to model details in the function tail. It is the model of choice for scintigraphic data which by definition cannot overshoot.<br><code>vol(t) = v0 * exp(-(t / tempt) ^ beta)</code><br>',

method_a = "<b>Methods to fit curves</b> For gastric emptying time series without too many missing data and few outliers, the <code>nlme</code> population fit provide quick and reliable fits. If <code>nlme</code> fails, the Bayesian method will work; it needs much more time and may bring the free online account on ShinyApps to their knee; no problem if you run this on your own computer. If you have a data set that does not give a fit with the Bayesian method, please post the example on the github issues page; link see on the bottom left of the app.",

noise_perc = "Noise amplitude, measured in % of the initial volume v0",
missing = "Fraction of randomly missing data to test the method for robustness",
manual = "Expert settings for parameters of simulated emptying curves and Bayesian fit.",
lkj = 'LKJ prior for kappa/tempt correlation, only required for model with covariance. Values from 1.5 (strong correlation) to 50 (almost independent) are useful. See http://www.psychstatistics.com/2014/12/27/d-lkj-priors/ for examples.',
student_df = "Student-t degrees of freedom for residual error; default 5. Use 3 when there are heavy outliers in the data set; values above 10 are adequate for almost gaussian residuals."
)

stan_models = ####################### add powexp_gastro_1d ++++
  matrix(c("linexp_gastro_1b", "linexp_gastro_2b",
           "powexp_gastro_1b", "powexp_gastro_2c"),
         nrow = 2, dimnames = list(c("nocov", "withcov"), c("linexp", "powexp")))

preset_description = function(id){
  as.character(presets[presets$id == id, "description"])
}
