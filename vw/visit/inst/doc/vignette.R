## ---- eval = TRUE, echo = FALSE, message = FALSE-------------------------
require(visit);
set.seed(10000);

## ---- eval = FALSE, echo = TRUE------------------------------------------
#  install.packages("visit");
#  require(visit);

## ---- eval = TRUE, echo = TRUE-------------------------------------------
tox       <-  c(0.07, 0.23, 0.66);
res       <-  c(0.50, 0.17, 0.59);
rho       <-  c(0.98, 0.40, 0.46);
scenario  <-  vtScenario(tox = tox, res = res, rho = rho);
summary(scenario);

## ---- eval = TRUE--------------------------------------------------------
summary(scenario);

## ---- eval = TRUE, echo = TRUE, fig.height = 4, fig.width = 8------------
oldpar <- par(mfrow = c(1,2));
plot(scenario, draw.curves = 1:2, main = "Marginal DLT Risk and Response Rates");
plot(scenario, draw.curves = 3:6, main = "Joint DLT Risk and Response Rates");
par(oldpar);

## ---- eval = TRUE, echo = TRUE-------------------------------------------
tau   <- c(0.39, 0.87, 0.49);
prior <- vtPriorPar(tau = tau, sdalpha = 10, sdrho = 10);

## ---- eval = TRUE, results = 'hide'--------------------------------------
simu <- vtSimu(n.rep = 100, trueps = scenario,
               size.cohort = 5, size.level = 10,
               etas = c(0.3, 0.7), dec.cut = c(0.45, 0.55, 0.75),
               prob.mdl = "NONPARA+");

## ---- eval = TRUE, echo = TRUE-------------------------------------------
sum.1 <- summary(simu);
print(sum.1);

## ---- eval = TRUE, echo = TRUE-------------------------------------------
sum.2 <- summary2(simu);
print(sum.2);

## ---- eval = TRUE, echo = TRUE, fig.height = 5, fig.width = 5------------
etas       <- c(0.1, 0.3)
dec.cut    <- c(0.6,0.6,0.6)
cur.obs.y  <- c(3, 2, 1, 1)
prev.obs.y <- c(5, 2, 0, 0)
rst.inter  <- vtInterim(cur.obs.y,  prev.obs.y = prev.obs.y,
                        prob.mdl = "NONPARA", etas = etas, dec.cut = dec.cut,
                        nsmp = 2000);

plot(rst.inter);

## ---- eval = TRUE, echo = TRUE, fig.height = 4, fig.width = 7------------
obs <- rbind(c(1, 6, 4, 3, 6), c(2, 4, 9, 3, 3), c(3, 2, 6, 6, 5));
vtTrack(obs, end.width = 0.8);

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  vtShiny();

