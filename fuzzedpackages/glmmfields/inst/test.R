library(ggplot2)
library(doParallel)
registerDoParallel(cores = 8L)

set.seed(1)
out <- plyr::ldply(seq_len(40), function(i) {
  N <- 200 # number of data points
  temperature <- rnorm(N, 0, 1) # simulated temperature data
  X <- cbind(1, temperature) # our design matrix
  s <- sim_glmmfields(
    n_draws = 1, gp_theta = 1.5, n_data_points = N,
    gp_sigma = 0.3, sd_obs = 0.2, n_knots = 10, obs_error = "gamma",
    covariance = "squared-exponential", X = X,
    B = c(0.5, 0.2) # B represents our intercept and slope
  )
  d <- s$dat
  d$temperature <- temperature

  m_spatial <- glmmfields(y ~ temperature,
    data = d, family = Gamma(link = "log"),
    lat = "lat", lon = "lon", nknots = 10, iter = 600, chains = 1, cores = 1,
    prior_intercept = student_t(3, 0, 10),
    prior_beta = student_t(3, 0, 3),
    prior_sigma = half_t(3, 0, 3),
    prior_gp_theta = half_t(3, 0, 10),
    prior_gp_sigma = half_t(3, 0, 3),
    seed = 123
  )

  o1 <- data.frame(t(quantile(
    rstan::extract(m_spatial$model)$gp_sigma, probs = c(0.1, 0.5, 0.9))))
  names(o1) <- paste0("gp_sigma_", names(o1))
  o2 <- data.frame(t(quantile(
    rstan::extract(m_spatial$model)$gp_theta, probs = c(0.1, 0.5, 0.9))))
  names(o2) <- paste0("gp_theta_", names(o2))
  data.frame(o1, o2)
}, .parallel = TRUE)

out$i <- seq_len(nrow(out))

median(out$gp_sigma_X50.)
median(out$gp_theta_X50.)

ggplot(out, aes(i, gp_sigma_X50., ymin = gp_sigma_X10., ymax = gp_sigma_X90.)) +
  geom_pointrange() +
  geom_hline(yintercept = 0.3) +
  scale_y_log10()

library(dplyr)
0.2 * 40
n <- filter(out, gp_sigma_X10. > 0.3 | gp_sigma_X90. < 0.3) %>% count()
n
round(1 - n / nrow(out), 2)

ggplot(out, aes(i, gp_theta_X50., ymin = gp_theta_X10., ymax = gp_theta_X90.)) +
  geom_pointrange() +
  geom_hline(yintercept = 1.5) +
  scale_y_log10()

0.2 * 40
n <- filter(out, gp_theta_X10. > 1.5 | gp_theta_X90. < 1.5) %>% count()
n
round(1 - n / nrow(out), 2)

plot(out$gp_sigma_X50., out$gp_theta_X50., log = "xy")
