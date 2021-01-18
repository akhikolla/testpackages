## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE---------------------------------------------------------
library(LeMaRns)

## ------------------------------------------------------------------------
NS_par

## ------------------------------------------------------------------------
nfish <- nrow(NS_par)
nsc <- 32
maxsize <- max(NS_par$Linf)*1.01 # the maximum size is 1% bigger than the largest Linf.
l_bound <- seq(0, maxsize, maxsize/nsc); l_bound <- l_bound[-length(l_bound)]
u_bound <- seq(maxsize/nsc, maxsize, maxsize/nsc)
mid <- l_bound+(u_bound-l_bound)/2

## ------------------------------------------------------------------------
Linf <- NS_par$Linf # the von-Bertalanffy asymptotic length of each species (cm)
W_a <- NS_par$W_a # length-weight conversion parameter
W_b <- NS_par$W_b # length-weight conversion parameter
k <- NS_par$k # the von-Bertalanffy growth parameter
Lmat <- NS_par$Lmat # the length at which 50\% of individuals are mature (cm)


## ------------------------------------------------------------------------
tmp <- calc_phi(k, Linf, nsc, nfish, u_bound, l_bound, calc_phi_min=FALSE, phi_min=0.1)
phi <- tmp$phi
phi_min <- tmp$phi_min

## ------------------------------------------------------------------------
tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b
                             , phi_min)
ration <- tmp$ration
sc_Linf <- tmp$sc_Linf
wgt <- tmp$wgt
g_eff <- tmp$g_eff

## ------------------------------------------------------------------------
mature <- calc_mature(Lmat, nfish, mid, kappa=rep(10, nfish), sc_Linf)

## ------------------------------------------------------------------------
other <- NS_other

## ---- echo=FALSE, fig.width=4, fig.height=3.5, fig.align="center"--------
x <- seq(1, 1e10, length.out=10000)
SSB <- x/1e9 # SSB in tonnes x 10^3
y <- 1.5*SSB*exp(-0.5*SSB)*1e6
plot(x/1e6, y, ylab="Recruits", xlab="SSB (tonnes)", type="l", main="Ricker", ylim=c(0, 1.1e6))

## ---- echo=FALSE, fig.height=6.5,fig.width=7, fig.align="center"---------
x <- seq(1, 1e10, length.out=10000)
SSB <- x/1e9 # SSB in tonnes x 10^3
par(mfrow=c(2, 2))
y <- (1.5*SSB/(1+1.5*SSB))*1e6
plot(x/1e6, y, ylab="Recruits", xlab="SSB (tonnes)", type="l", main="Beverton-Holt", ylim=c(0, 1.1e6))
tmp <- 0.25*SSB; hs_b <- 1
y <- ifelse(tmp<hs_b, tmp, hs_b)*1e6
plot(x/1e6, y, ylab="Recruits", xlab="SSB (tonnes)", type="l", main="hockey-stick", ylim=c(0, 1.1e6))
plot(x/1e6, tmp*1e6, ylab="Recruits", xlab="SSB (tonnes)", type="l", main="linear")
plot(x/1e6, rep(hs_b, length(x))*1e6, ylab="Recruits", xlab="SSB (tonnes)", type="l", main="constant", ylim=c(0, 1.1e6))

## ------------------------------------------------------------------------
stored_rec_funs <- get_rec_fun(rep("hockey-stick", nfish))
recruit_params <- do.call("Map", c(c, list(a=NS_par$a, b=NS_par$b)))

## ---- echo=FALSE, fig.height=3, fig.width=7, fig.align="center"----------
x <- seq(0, 125, 1)
y <- ifelse(x<0.75*125, 0, 0.8)
par(mfrow=c(1, 3))
plot(x, y, main="std_RNM", xlab="Length", ylab="M1", ylim=c(0, 1), type="l")
plot(x, rep(0.8, length(x)), main="constant", xlab="Length", ylab="M1", ylim=c(0, 1), type="l")
plot(x, 0.8/125*x, main="linear", xlab="Length", ylab="M1", ylim=c(0, 1), type="l")

## ------------------------------------------------------------------------
M1 <- calc_M1(nsc, sc_Linf, phi_min, 
              natmort_opt=rep("std_RNM", length(sc_Linf)), 
              Nmort=rep(0.8, length(sc_Linf)), 
              prop=rep(0.75, length(sc_Linf)))

## ------------------------------------------------------------------------
prefs <- calc_prefs(pred_mu=-2.25, pred_sigma=0.5, wgt, sc_Linf)

## ------------------------------------------------------------------------
suit_M2 <- calc_suit_vect(nsc, nfish, sc_Linf, prefs, NS_tau)

## ---- echo=FALSE, fig.width=4, fig.height=3.5, fig.align="center"--------
x <- seq(0, 125, 0.1)
y <- 1/(1+exp(-0.25*(x-50)))
plot(x, y, type="l", xlab="Length", ylab="Catchability", main="logistic")

## ---- echo=FALSE, fig.width=4, fig.height=3.5, fig.align="center"--------
x <- seq(0, 125, 0.1)
y <- dlnorm(x, log(50), 1)
plot(x, y/max(y), type="l", xlab="Length", ylab="Catchability", main="log_gaussian")

## ---- echo=FALSE, fig.width=4, fig.height=3.5, fig.align="center"--------
x <- seq(0, 125, 0.1)
y <- x>50
plot(x, y, type="l", xlab="Length", ylab="Catchability", main="knife-edge")

## ------------------------------------------------------------------------
Qs <- calc_Q(nsc=nsc, nfish=nfish, mid=mid, l_bound=l_bound, u_bound=u_bound, 
             species_names=NS_par$species_names,curve=rep("logistic", nfish), 
             species=NS_par$species_names, max_catchability=rep(1, nfish), 
             gear_name=NS_par$species_names,eta=rep(0.25, nfish), L50=Lmat)

## ------------------------------------------------------------------------
NS_params <- new("LeMans_param",
                  nsc=nsc,
                  nfish=nfish,
                  phi_min=phi_min,
                  l_bound=l_bound,
                  u_bound=u_bound,
                  mid=mid,
                  species_names=NS_par$species_names,
                  Linf=Linf,
                  W_a=W_a,
                  W_b=W_b,
                  k=k,
                  Lmat=Lmat,
                  mature=mature,
                  sc_Linf=sc_Linf,
                  wgt=wgt,
                  phi=phi,
                  ration=ration,
                  other=other,
                  M1=M1,
                  suit_M2=suit_M2,
                  Qs=Qs,
                  stored_rec_funs=stored_rec_funs,
                  eps=1e-05,
                  recruit_params=recruit_params)

## ------------------------------------------------------------------------
NS_params <- LeMansParam(species_names=NS_par$species_names, 
                         Linf=Linf, k=k, W_a=W_a, W_b=W_b,
                         Lmat=Lmat, tau=NS_tau, 
                         recruit_params=list(a=NS_par$a, b=NS_par$b),
                         eta=rep(0.25, 21), L50=Lmat)

## ---- eval=FALSE---------------------------------------------------------
#  NS_params <- LeMansParam(NS_par, tau=NS_tau, eta=rep(0.25, 21), L50=NS_par$Lmat, other=NS_other)

## ---- echo=FALSE---------------------------------------------------------
NS_params <- LeMansParam(NS_par, tau=NS_tau, eta=rep(0.25, 21), L50=NS_par$Lmat, other=NS_other)

## ---- eval=FALSE---------------------------------------------------------
#  N <- get_N0(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5)

## ---- echo=FALSE,fig.height=6.5,fig.width=7,fig.align="center"-----------
N <- get_N0(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5)
par(mfrow=c(1,1), mar=c(5,5,5,0))
layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2), nrow=3, byrow=TRUE))
rainbowcols <- rainbow(ncol(N), s=0.75)
plot(c(min(mid), max(mid)), c(min(N+0.1), max(N)), type="n", xlab="Length (cm)", ylab="N", log="xy", 
     font.lab=2, cex.lab=1.5, cex.axis=1.5)
lines(mid, rowSums(N), col="black", lty=2)
for (i in 1:ncol(N)) {
  lines(mid, N[, i], col=rainbowcols[i])
}
par(mar=c(0,0,0,0))
plot(0, 0, axes=F, type="n", xlab="", ylab="")
legend("center", legend=c("Community", NS_params@species_names), col=c("black", rainbowcols), 
       lty=c(2, rep(1, length(rainbowcols))), cex=1.2, bty="n")

## ---- eval=FALSE---------------------------------------------------------
#  effort <- rep(0.5, dim(Qs)[3])
#  Fs <- matrix(0, nsc, nfish)
#  for (j in 1:length(effort)) {
#     Fs <- Fs+effort[j]*Qs[,,j]
#  }

## ---- eval=FALSE---------------------------------------------------------
#  SSB <- calc_SSB(mature, N, wgt)
#  R <- calc_recruits(SSB, stored_rec_funs, recruit_params)
#  N[1, ] <- N[1,]+R

## ---- eval=FALSE---------------------------------------------------------
#  M2 <- calc_M2(N, ration, wgt, nfish, nsc, other, sc_Linf, suit_M2)
#  Z <- Fs+M1+M2+1e-05
#  Catch <- (Fs/Z)*N*(1-exp(-Z))*wgt
#  N <- N*exp(-Z)

## ---- eval=FALSE---------------------------------------------------------
#  N <- calc_growth(N, phi, nfish, nsc)

## ------------------------------------------------------------------------
N0 <- get_N0(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5) # initialise N
years <- 50 # run for 10 years
tot_time <- years*phi_min # calculate the total number of time steps

## ------------------------------------------------------------------------
effort <- matrix(0.5, tot_time, dim(Qs)[3])
# Calculate F
Fs <- array(0, dim=c(nsc, nfish, tot_time))
for (j in 1:ncol(effort)) {
  for (ts in 1:tot_time) {
    Fs[,,ts] <- Fs[,,ts]+effort[ts, j]*Qs[,,j]
  }
}

## ------------------------------------------------------------------------
model_run <- run_LeMans(N0=N0, tot_time=tot_time, Fs=Fs, nsc=nsc, nfish=nfish, phi_min=phi_min, 
                        mature=mature, sc_Linf=sc_Linf, wgt=wgt, phi=phi, ration=ration, 
                        other=other, M1=M1, suit_M2=suit_M2, stored_rec_funs=stored_rec_funs, 
                        recruit_params=recruit_params, eps=1e-05)

## ------------------------------------------------------------------------
model_run <- run_LeMans(NS_params, N0=N0, Fs=Fs, tot_time=tot_time)

## ------------------------------------------------------------------------
effort_mat <- matrix(0.5, years, dim(Qs)[3])
model_run <- run_LeMans(NS_params, years=50, effort=effort_mat)

## ---- eval=FALSE---------------------------------------------------------
#  biomass <- get_biomass(N=N, wgt=wgt)

## ---- eval=FALSE---------------------------------------------------------
#  biomass <- get_biomass(inputs=NS_params, outputs=model_run)

## ---- eval=FALSE---------------------------------------------------------
#  biomass <- get_biomass(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:4)

## ---- eval=FALSE---------------------------------------------------------
#  biomass <- get_biomass(inputs=NS_params, outputs=model_run, time_steps=1:500,
#                         species=c("Herring", "Whiting", "Haddock", "Cod"))

## ---- eval=FALSE---------------------------------------------------------
#  plot_biomass(N=N, wgt=wgt, time_steps=1:500, species=1:4)

## ---- eval=FALSE---------------------------------------------------------
#  plot_biomass(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21)

## ---- echo=FALSE,fig.height=5,fig.width=7,fig.align="center"-------------
plot_biomass(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21)

## ---- eval=FALSE---------------------------------------------------------
#  biomass <- get_biomass(inputs=NS_params, outputs=model_run)
#  plot_biomass(biomass=biomass, time_steps=1:500, species=1:4)

## ---- eval=TRUE,fig.height=5,fig.width=7,fig.align="center"--------------
plot_biomass(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:4, 
             full_plot_only=FALSE)

## ----fig.height=5,fig.width=7,fig.align="center"-------------------------
SSB <- get_SSB(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21)
plot_SSB(SSB=SSB, time_steps=1:500, species=1:21, species_names=NS_params@species_names)

## ---- eval=FALSE---------------------------------------------------------
#  LFI <- get_LFI(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21,
#                 length_LFI=c(30, 40))
#  MML <- get_MML(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21)
#  TyL <- get_TyL(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21)
#  LQ <- get_LQ(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21,
#               prob=c(0.5, 0.95))

## ---- eval=FALSE---------------------------------------------------------
#  plot_LFI(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21,
#                 length_LFI=c(30, 40))
#  plot_MML(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21)
#  plot_TyL(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21)
#  plot_LQ(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21,
#               prob=c(0.5, 0.95))

## ---- eval=FALSE---------------------------------------------------------
#  plot_LFI(LFI, length_LFI=c(30, 40))
#  plot_MML(MML)
#  plot_TyL(TyL)
#  plot_LQ(LQ, prob=c(0.5, 0.95))

## ---- eval=TRUE, echo=FALSE, fig.width=5.5, fig.height=4, fig.align="center"----
plot_LFI(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21, 
               length_LFI=c(30, 40))
plot_MML(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21)
plot_TyL(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21)
plot_LQ(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21, 
             prob=c(0.5, 0.95))

## ---- eval=FALSE---------------------------------------------------------
#  indicators <- get_indicators(inputs=NS_params, outputs=model_run, time_steps=1:500,
#                               species=1:21,length_LFI=c(30, 40), prob=c(0.5, 0.95))

## ---- eval=FALSE---------------------------------------------------------
#  plot_indicators(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21,
#                  length_LFI=c(30, 40), prob=c(0.5, 0.95))

## ---- eval=FALSE---------------------------------------------------------
#  plot_indicators(LFI=indicators[['LFI']], MML=indicators[['MML']],
#                  TyL=indicators[['TyL']], LQ=indicators[['LQ']],
#                  time_steps=1:500, species=1:21,
#                  length_LFI=c(30, 40), prob=c(0.5, 0.95))

## ---- eval=TRUE, echo=FALSE,fig.height=5,fig.width=7,fig.align="center"----
plot_indicators(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21, 
                length_LFI=c(30, 40), prob=c(0.5, 0.95))

## ---- eval=FALSE---------------------------------------------------------
#  CPUE <- get_CPUE(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21)
#  CPG <- get_CPG(inputs=NS_params, outputs=model_run, time_steps=1:500, species=1:21)

## ------------------------------------------------------------------------
NS_mixed_fish

## ------------------------------------------------------------------------
NS_params <- LeMansParam(df=NS_par, gdf = NS_mixed_fish, tau=NS_tau, eta=NS_eta, L50=NS_L50, other=NS_other)

## ------------------------------------------------------------------------
effort_mat <- matrix(0, 50, dim(NS_params@Qs)[3])
colnames(effort_mat) <- c("Industrial", "Otter", "Beam", "Pelagic")
model_run <- run_LeMans(NS_params, years=50, effort=effort_mat)

## ----fig.height=5,fig.width=7,fig.align="center"-------------------------
plot_SSB(inputs=NS_params, outputs=model_run, full_plot_only=TRUE)

## ----fig.height=5,fig.width=7,fig.align="center"-------------------------
plot_indicators(inputs=NS_params, outputs=model_run)

## ------------------------------------------------------------------------
v_SSB <- colMeans(get_SSB(inputs=NS_params, outputs=model_run, time_steps=492:501))/1e6

## ------------------------------------------------------------------------
ef_lvl <- c(0, 0.5, 1, 1.5, 2)
efs <- expand.grid(Industrial=ef_lvl, Otter=ef_lvl, Beam=ef_lvl, Pelagic=ef_lvl)

## ------------------------------------------------------------------------
run_the_model<- function(ef){
  effort_mat <- matrix(ef, 50, dim(NS_params@Qs)[3], byrow=T)
  colnames(effort_mat) <- c("Industrial", "Otter", "Beam", "Pelagic")
  model_run <- run_LeMans(params=NS_params, N0=model_run@N[,,501], 
                          years=50, effort=effort_mat)
  return(model_run)
}

## ---- eval=F-------------------------------------------------------------
#  sce <- apply(efs, 1, run_the_model)

## ----echo=F--------------------------------------------------------------
load("sce.rda")

## ---- eval=F-------------------------------------------------------------
#  LFI <- unlist(lapply(lapply(sce, FUN=get_LFI, time_steps=492:501, inputs=NS_params),
#                       mean))
#  TyL <- unlist(lapply(lapply(sce, FUN=get_TyL, time_steps=492:501, inputs=NS_params),
#                       mean))
#  MML <- unlist(lapply(lapply(sce, FUN=get_MML, time_steps=492:501, inputs=NS_params),
#                       mean))

## ---- eval=T,fig.height=5,fig.width=7,fig.align="center"-----------------
par(mfrow=c(2, 2))
boxplot(LFI~efs[, 1], main="Industrial", xlab="Effort", ylab="LFI")
boxplot(LFI~efs[, 2], main="Otter", xlab="Effort", ylab="LFI")
boxplot(LFI~efs[, 3], main="Beam", xlab="Effort", ylab="LFI")
boxplot(LFI~efs[, 4], main="Pelagic", xlab="Effort", ylab="LFI")

par(mfrow=c(2, 2))
boxplot(TyL~efs[, 1], main="Industrial", xlab="Effort", ylab="TyL")
boxplot(TyL~efs[, 2], main="Otter", xlab="Effort", ylab="TyL")
boxplot(TyL~efs[, 3], main="Beam", xlab="Effort", ylab="TyL")
boxplot(TyL~efs[, 4], main="Pelagic", xlab="Effort", ylab="TyL")

par(mfrow=c(2, 2))
boxplot(MML~efs[, 1], main="Industrial", xlab="Effort", ylab="MML")
boxplot(MML~efs[, 2], main="Otter", xlab="Effort", ylab="MML")
boxplot(MML~efs[, 3], main="Beam", xlab="Effort", ylab="MML")
boxplot(MML~efs[, 4], main="Pelagic", xlab="Effort", ylab="MML")

## ---- eval=F-------------------------------------------------------------
#  # SSB of the final years (tonnes)
#  new_SSB <- do.call(rbind, lapply(lapply(sce, FUN=get_SSB, time_steps=492:501,
#                                          inputs=NS_params), colMeans))/1e6
#  # Relative SSB
#  rel_SSB <- t(t(new_SSB)/v_SSB)
#  colnames(rel_SSB) <- NS_params@species_names

## ---- eval=T,fig.height=5,fig.width=7,fig.align="center"-----------------
par(mfrow=c(2,2))
boxplot(rel_SSB[, "Saithe"]~efs[, 1], main="Industrial", xlab="Effort",
        ylab="Relative SSB")
boxplot(rel_SSB[, "Saithe"]~efs[, 2], main="Otter", xlab="Effort",
        ylab="Relative SSB")
boxplot(rel_SSB[, "Saithe"]~efs[, 3], main="Beam", xlab="Effort",
        ylab="Relative SSB")
boxplot(rel_SSB[, "Saithe"]~efs[, 4], main="Pelagic", xlab="Effort",
        ylab="Relative SSB")

## ---- eval=T ,fig.height=5,fig.width=7,fig.align="center"----------------
par(mfrow=c(2,2))
boxplot(rel_SSB[, "Horse mackerel"]~efs[, 1], main="Industrial", xlab="Effort",
        ylab="Relative SSB")
boxplot(rel_SSB[, "Horse mackerel"]~efs[, 2], main="Otter", xlab="Effort",
        ylab="Relative SSB")
boxplot(rel_SSB[, "Horse mackerel"]~efs[, 3], main="Beam", xlab="Effort",
        ylab="Relative SSB")
boxplot(rel_SSB[, "Horse mackerel"]~efs[, 4], main="Pelagic", xlab="Effort",
        ylab="Relative SSB")

## ---- ,fig.height=5,fig.width=7,fig.align="center",eval=F----------------
#  risk <- apply(rel_SSB, 1, function(x){return(sum(x<0.1))})
#  
#  par(mfrow=c(2, 2))
#  boxplot(risk~efs[, 1], main="Industrial", xlab="Effort", ylab="Stocks at risk")
#  boxplot(risk~efs[, 2], main="Otter", xlab="Effort", ylab="Stocks at risk")
#  boxplot(risk~efs[, 3], main="Beam", xlab="Effort", ylab="Stocks at risk")
#  boxplot(risk~efs[, 4], main="Pelagic", xlab="Effort", ylab="Stocks at risk")

## ---- eval=T, fig.height=5,fig.width=7,fig.align="center"----------------
z_mat <- outer(ef_lvl, ef_lvl, FUN=function(x, y, efs) {
  mapply(function(x, y, efs) {
  mean(risk[intersect(which(efs[, 2]==x), which(efs[, 3]==y))])}, x=x, y=y,
  MoreArgs=list(efs=efs))
}, efs=efs)
layout(matrix(c(1,1,1,1,1,2), nrow=1, byrow=TRUE))
image(z=-z_mat, x=ef_lvl, y=ef_lvl, axes=T, cex.lab=1.5, xlab="Otter effort", ylab="Beam effort")
axis(1); axis(2); box()
image(z=-matrix(sort(unique(as.numeric(z_mat))), nrow=1), 
      y=sort(unique(as.numeric(z_mat))), axes=F, cex.lab=1.5, 
      ylab="Number of stocks at risk")
axis(2); box()
box()

## ----fig.height=5,fig.width=7,fig.align="center"-------------------------
Industrial <- rep(1.5, 20)
Otter <- -1/100*1:20*(1:20-20)
Beam <- 1:20*1/20+0.25
Pelagic <- 1+c(1:5*1/5, 5:1*1/5, 1:5*1/5, 5:1*1/5)
par(mfrow=c(2, 2))
plot(1:20, Industrial, ylab="Effort", main="Industrial", xlab="Year", 
     ylim=c(0, 2), type="s")
plot(1:20, Otter, ylab="Effort", main="Otter", xlab="Year", 
     ylim=c(0, 2), type="s")
plot(1:20, Beam, ylab="Effort", main="Beam", xlab="Year", 
     ylim=c(0, 2), type="s")
plot(1:20, Pelagic, ylab="Effort", main="Pelagic", xlab="Year", 
     ylim=c(0, 2), type="s")

## ------------------------------------------------------------------------
# Set up effort for the model run
effort_mat <- cbind(Industrial, Otter, Beam, Pelagic)
colnames(effort_mat) <- c("Industrial", "Otter", "Beam", "Pelagic")
model_run_dyn <- run_LeMans(params=NS_params, N0=model_run@N[,,501],
                            years=20, effort=effort_mat)
catches <- get_annual_catch(inputs=NS_params, outputs=model_run_dyn)/1e6 # in tonnes
colnames(catches) <- NS_params@species_names

## ----fig.height=5,fig.width=7,fig.align="center"-------------------------
par(mfrow=c(2, 2))
plot(1:20, catches[, "Sprat"], type="l", main="Sprat", xlab="Year", 
     ylab="Catch (tonnes)", ylim=c(0, max(catches[, "Sprat"])))
plot(1:20, catches[, "Cod"], type="l", main="Cod", xlab="Year", 
     ylab="Catch (tonnes)", ylim=c(0, max(catches[, "Cod"])))
plot(1:20, catches[, "Sole"], type="l", main="Sole", xlab="Year", 
     ylab="Catch (tonnes)", ylim=c(0,max(catches[, "Sole"])))
plot(1:20, catches[, "Herring"], type="l", main="Herring", xlab="Year", 
     ylab="Catch (tonnes)", ylim=c(0, max(catches[, "Herring"])))

## ---- eval = T,fig.height=5,fig.width=7,fig.align="center"---------------
plot_SSB(inputs=NS_params, outputs=model_run_dyn, 
         species=c("Sprat", "Cod", "Sole", "Herring"), 
         full_plot_only=FALSE)

## ----fig.height=5,fig.width=7,fig.align="center"-------------------------
catch_per_gear <- get_CPG(inputs=NS_params, outputs=model_run_dyn, effort=effort_mat)
tot_pg <- apply(catch_per_gear, c(2, 3), sum)/1e6
year <- rep(1:20, each=10)
# Total catch per gear per year
tot_pgpy <- t(sapply(1:20, function(x, year){
  tele <- which(year==x)
  if (length(tele)>1){
    return(rowSums(tot_pg[, tele]))
  }
  return(tot_pg[, tele])
}, year=year))

par(mfrow=c(2, 2))
plot(1:20, tot_pgpy[, 1], type="l", ylim=c(0, max(tot_pgpy[, 1])), xlab="Year", 
     main="Industrial", ylab="Catch (tonnes)")
plot(1:20, tot_pgpy[, 2], type="l", ylim=c(0, max(tot_pgpy[, 2])), xlab="Year", 
     main="Otter", ylab="Catch (tonnes)")
plot(1:20, tot_pgpy[, 3], type="l", ylim=c(0, max(tot_pgpy[, 3])), xlab="Year", 
     main="Beam", ylab="Catch (tonnes)")
plot(1:20, tot_pgpy[, 4], type="l", ylim=c(0, max(tot_pgpy[, 4])), xlab="Year", 
     main="Pelagic", ylab="Catch (tonnes)")

## ----fig.height=5,fig.width=7,fig.align="center"-------------------------
CPUE <- get_CPUE(inputs=NS_params, outputs=model_run_dyn, effort=effort_mat)/1e6
cpue_py <- t(sapply(1:20, function(x, year){
  tele <- which(year==x)
  if (length(tele)>1){
    return(colMeans(CPUE[tele, ]))
  }
  return(CPUE[tele, ])
}, year=year))

colnames(cpue_py) <- NS_params@species_names

par(mfrow=c(2, 2))
plot(1:20, cpue_py[, "Sprat"], type="l", main="Sprat", xlab="Year", 
     ylab="CPUE (tonnes)", ylim=c(0, max(cpue_py[, "Sprat"])))
plot(1:20, cpue_py[, "Cod"], type="l", main="Cod", xlab="Year", 
     ylab="CPUE (tonnes)", ylim=c(0, max(cpue_py[, "Cod"])))
plot(1:20, cpue_py[, "Sole"], type="l", main="Sole", xlab="Year", 
     ylab="CPUE (tonnes)", ylim=c(0, max(cpue_py[, "Sole"])))
plot(1:20, cpue_py[, "Herring"], type="l", main="Herring", xlab="Year", 
     ylab="CPUE (tonnes)", ylim=c(0, max(cpue_py[, "Herring"])))

## ------------------------------------------------------------------------
rec_fish <- data.frame(catch_species=c("Cod", "Haddock", "Herring", "Horse mackerel", 
                                       "Mackerel", "Plaice", "Saithe", "Whiting"), 
                       curve=rep("knife-edge"), gear_name="Recreational")
rec_fish$max_catchability <- c(0.01, 0.01, 0.005, 0.05, 0.05, 0.01, 0.01, 0.02)
Lmin <- c(35, 30, 20, 15, 30, 27, 35, 27)

## ------------------------------------------------------------------------
gdf <- rbind(NS_mixed_fish, rec_fish)

## ------------------------------------------------------------------------
eta1 <- c(NS_eta, rep(0, 8))
L501 <- c(NS_L50, rep(0, 8))
Lmin1 <- c(rep(0, 21), Lmin)

## ------------------------------------------------------------------------
NS_params_rec <- LeMansParam(df=NS_par, gdf = gdf, tau=NS_tau, eta=eta1, L50=L501, 
                             other=NS_other, Lmin=Lmin1)
effort_mat1 <- cbind(effort_mat, 0.1+1:20*0.05/20)
colnames(effort_mat1)[5] <- "Recreational"

## ---- fig.width=5.5, fig.height=4, fig.align="center"--------------------
plot(1:20, effort_mat1[, "Recreational"], xlab="Years", ylab="Effort", type="l", 
     ylim=c(0, 2), main="Recreational")

## ------------------------------------------------------------------------
model_run_rec <- run_LeMans(params=NS_params_rec, N0=model_run@N[,,501], 
                            years=20, effort=effort_mat1)

## ---- fig.height=5,fig.width=7,fig.align="center"------------------------
catch_per_gear_rec <- get_CPG(inputs=NS_params_rec, outputs=model_run_rec, effort=effort_mat1)
rec_py <- t(sapply(1:20, function(x, year){
  tele <- which(year==x)
  if (length(tele)>1){
    return(rowSums(catch_per_gear_rec[, 5, tele]))
  }
  return(catch_per_gear_rec[, 5, tele])
}, year=year))/1e6
colnames(rec_py) <- NS_params_rec@species_names

par(mfrow=c(2, 2))
plot(1:20, rec_py[, "Herring"], type="l", main="Herring", xlab="Year", 
     ylab="Catch (tonnes)")
plot(1:20, rec_py[, "Mackerel"], type="l", main="Horse mackerel", xlab="Year", 
     ylab="Catch (tonnes)")
plot(1:20, rec_py[, "Cod"], type="l", main="Cod", xlab="Year", 
     ylab="Catch (tonnes)")
plot(1:20 ,rec_py[, "Saithe"], type="l", main="Saithe", xlab="Year", 
     ylab="Catch (tonnes)")

## ----eval=T--------------------------------------------------------------
NS_params_n <- LeMansParam(NS_par, tau=NS_tau, eta=NS_eta, L50=NS_L50, other=NS_other)
model_run_n_init <- run_LeMans(NS_params_n, years=50)
# Set the initial state
N0 <- model_run_n_init@N[,,501]

## ----eval=T--------------------------------------------------------------
calc_catch<-function(x,i,eff){
  eff[i] <- x
  tmp <- run_LeMans(NS_params_n, N0=N0, years=20, effort=matrix(eff, nrow=20, ncol=21, byrow=T))
  return(mean(tail(colSums(tmp@Catch[,i,]), 10)))
}

## ----eval=F--------------------------------------------------------------
#  eff <- fmsy <-c(1.3, 0.35, 0.35, 0.72, 0.6, 0.41, 0.25, 0.5, 0.33, 0.22, 0.32, 0.21, 0.27,
#                  0.27, 0.25, 0.15, 0.30, 0.11, 0.1, 0.19, 0.3)
#  stat <- rep(FALSE, 21)
#  while(any(stat==FALSE)){
#    for (i in 1:21){
#      opts <- optim(eff[i], calc_catch, i=i, eff=eff, method="Brent", lower=0, upper=2,control = list(fnscale = -1))
#      stat[i] <- abs(eff[i]-opts$par)<0.01
#      eff[i] <- opts$par
#    }
#  }

## ----echo=F--------------------------------------------------------------
load("nash.rda")
msy <- run_LeMans(NS_params_n, N0=N0, years=20, effort=matrix(fmsy_lm, nrow=20, ncol=21, byrow=T))

## ----eval=T,fig.height=5,fig.width=7,fig.align="center"------------------
nash <- run_LeMans(NS_params_n, N0=N0, years=20, effort=matrix(eff, nrow=20, ncol=21, byrow=T))
plot_SSB(inputs=NS_params_n, outputs=nash)

## ----eval=F--------------------------------------------------------------
#  fmsy_lm <- sapply(1:21, function(i, eff){optim(eff[i], calc_catch, i=i, eff=eff,
#                                             method="Brent", lower=0, upper=2,control = list(fnscale = -1))$par},
#                  eff=fmsy)
#  msy <- run_LeMans(NS_params_n, N0=N0, years=20, effort=matrix(fmsy_lm, nrow=20, ncol=21, byrow=T))

## ----eval=T,fig.width=5.5, fig.height=4, fig.align="center"--------------
plot(fmsy_lm, eff, ylab=expression(F[Nash]), xlab=expression(F[MSY]))
abline(a=0, b=1)

