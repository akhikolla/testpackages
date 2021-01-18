## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(fig.width = 12, fig.height=8, fig.align="default")
knitr::opts_chunk$set(error = TRUE)

## ----hide=TRUE----------------------------------------------------------------
oldpar <- par(no.readonly=TRUE)
set.seed(1)
library(BGVAR)

## ----"eerData"----------------------------------------------------------------
data(eerData)

## ----"eerData2"---------------------------------------------------------------
names(eerData)

## ----"eerData3"---------------------------------------------------------------
colnames(eerData$UK)

## ----"US",echo=TRUE-----------------------------------------------------------
head(eerData$US)

## ----"convert",echo=TRUE------------------------------------------------------
bigX<-list_to_matrix(eerData)

## ----"convert2",echo=TRUE-----------------------------------------------------
colnames(bigX)[1:10]

## ----"tradeW",echo=TRUE-------------------------------------------------------
head(W.trade0012)

## ----"rownames.W"-------------------------------------------------------------
all(colnames(W.trade0012)==names(eerData))

## ----"rowSums.W"--------------------------------------------------------------
rowSums(W.trade0012)
diag(W.trade0012)

## ---- "eerDatasmall", hide=TRUE-----------------------------------------------
cN<-c("EA","US","RU")
eerData<-eerData[cN]
W.trade0012<-W.trade0012[cN,cN]
W.trade0012<-apply(W.trade0012,2,function(x)x/rowSums(W.trade0012))
W.list<-lapply(W.list,function(l){l<-apply(l[cN,cN],2,function(x)x/rowSums(l[cN,cN]))})

## ----"model.1",results="hide"-------------------------------------------------
 model.1<-bgvar(Data=eerData,
                W=W.trade0012,
                draws=100,
                burnin=100,
                plag=1,
                prior="NG",
                hyperpara=NULL, 
                SV=TRUE,
                thin=1,
                trend=TRUE,
                h=0,
                save.country.store=FALSE,
                eigen=1.05
                )


## ----"SV",results="hide"------------------------------------------------------
model.1$cc.results$sig$EA[,"EA.y","EA.y"]

## ----"ng.eigen",echo=TRUE-----------------------------------------------------
model.1$stacked.results$F.eigen[1:10]

## ----"print.model",echo=TRUE--------------------------------------------------
print(model.1)

## ----"summary.model"----------------------------------------------------------
 summary(model.1)

## ----"stats",echo=TRUE, results="hide"----------------------------------------
Fmat <- coef(model.1)
Smat <- vcov(model.1)
lik  <- logLik(model.1)

## ----"insample",fig.margin=TRUE,fig.width=6,fig.height=8,fig.cap="In-sample fit for euro area variables"----
yfit <- fitted(model.1)
plot(model.1, global=FALSE, resp="EA")

## ----"ssvs.1",echo=TRUE, results="hide"---------------------------------------
 model.ssvs.1<-bgvar(Data=eerData,
                     W=W.trade0012,
                     draws=100,
                     burnin=100,
                     plag=1,
                     prior="SSVS",
                     hyperpara=NULL, 
                     SV=TRUE,
                     thin=1,
                     trend=TRUE,
                     h=0,
                     save.country.store=FALSE,
                     eigen=1.05
                     )


## ----"Pips"-------------------------------------------------------------------
model.ssvs.1$cc.results$PIP$PIP.cc$EA

## ----"pips.avg"---------------------------------------------------------------
model.ssvs.1$cc.results$PIP$PIP.avg

## ----"var.weight"-------------------------------------------------------------
eerData2<-eerData
variable.list<-list();variable.list$real<-c("y","Dp","tb");variable.list$fin<-c("stir","ltir","rer")

## ----results="hide"-----------------------------------------------------------
# weights for first variable set tradeW.0012, for second finW0711
model.ssvs.2<-bgvar(Data=eerData2,
                    W=W.list[c("tradeW.0012","finW0711")],
                    plag=1,
                    draws=100,
                    burnin=100,
                    prior="SSVS",
                    SV=TRUE,
                    thin=1,
                    hyperpara=NULL, 
                    eigen=TRUE,
                    variable.list=variable.list,
                    OE.weights=NULL, 
                    Wex.restr=NULL,
                    trend=TRUE,
                    save.country.store=FALSE
                    )

## ----"ltir.estimate", results="hide"------------------------------------------
# does include ltir* only when ltir is missing domestically
model.ssvs.3<-bgvar(Data=eerData,
                    W=W.trade0012,
                    plag=1,
                    draws=100,
                    burnin=100,
                    prior="SSVS",
                    SV=TRUE,
                    thin=1,
                    hyperpara=NULL, 
                    eigen=TRUE,
                    variable.list=NULL,
                    OE.weights=NULL, 
                    Wex.restr="ltir",
                    trend=TRUE,
                    save.country.store=FALSE
                    )

## ----"print.model.ssvs.3"-----------------------------------------------------
 print(model.ssvs.3)

## ----"OC"---------------------------------------------------------------------
eerData2$OC<-eerData$US[,c("poil"),drop=FALSE] # move oil prices into own slot
eerData2$US<-eerData$US[,c("y","Dp", "rer" , "stir", "ltir","tb")] # exclude it from US m odel

## ----"OC.weights"-------------------------------------------------------------
OC.weights<-list()
OC.weights$weights<-rep(1/3, 3)
names(OC.weights$weights)<-names(eerData2)[1:3] # last one is OC model, hence only until 3
OC.weights$variables<-c(colnames(eerData2$OC),"y") # first entry, endog. variables, second entry weighted average of y from the other countries to proxy demand
OC.weights$exo<-"poil"

## ----"OC.weights2"------------------------------------------------------------
# other entities weights with same name as new oil country
OE.weights <- list(OC=OC.weights)

## ----"estimate.OC",results="hide"---------------------------------------------
model.ssvs.4<-bgvar(Data=eerData2,
                    W=W.trade0012,
                    plag=1,
                    draws=100,
                    burnin=100,
                    prior="SSVS",
                    SV=TRUE,
                    thin=1,
                    hyperpara=NULL, 
                    variable.list=NULL,
                    OE.weights=OE.weights, 
                    trend=TRUE,
                    save.country.store=FALSE
                    )

## ----"aux"--------------------------------------------------------------------
aux1<-model.ssvs.1$cc.results$PIP$PIP.avg;aux1<-aux1[-nrow(aux1),1:6]
aux2<-model.ssvs.2$cc.results$PIP$PIP.avg;aux2<-aux2[-nrow(aux2),1:6]
aux3<-model.ssvs.3$cc.results$PIP$PIP.avg;aux3<-aux3[-nrow(aux3),1:6]
aux4<-model.ssvs.4$cc.results$PIP$PIP.avg;aux4<-aux4[-nrow(aux4),1:6]

## ----"heat1", fig.show="hold",out.width="25%",fig.cap="Heatmaps of PIPs."-----
heatmap(aux1,Rowv=NA,Colv=NA, main="Model 1")
heatmap(aux2,Rowv=NA,Colv=NA, main="Model 2")
heatmap(aux3,Rowv=NA,Colv=NA, main="Model 3")
heatmap(aux4,Rowv=NA,Colv=NA, main="Model 4")

## ----"us.mp", results="hide"--------------------------------------------------
  # US monetary policy shock
 shocks<-list();shocks$var="stir";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-100
 irf.chol.us.mp<-irf(model.ssvs.1,shock=shocks,n.ahead=24,save.store=TRUE)

## ---- "us.mp2"----------------------------------------------------------------
names(irf.chol.us.mp)

## ----"us.mp4", fig.margin=TRUE,out.width="80%",fig.cap="Responses of US country model"----
plot(irf.chol.us.mp,resp="US")

## ---- "us.gdp", results="hide"------------------------------------------------
# Recursive US GDP
shocks<-list();shocks$var="y";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-1
irf.chol.us.y<-irf(model.ssvs.1,shock=shocks,n.ahead=24)
# GIRF US GDP
shocks<-list();shocks$var="y";shocks$cN<-"US";shocks$ident="girf";shocks$scal=-1
irf.girf.us.y<-irf(model.ssvs.1,shock=shocks,n.ahead=24)

## ---- "us.gdp.plots",fig.cap="Comparison of responses Cholesky (left) and GIRF (right) to a negative GDP shock.",fig.show="hold",out.width="25%"----
plot(irf.chol.us.y,resp="US.y")
plot(irf.girf.us.y,resp="US.y")

plot(irf.chol.us.y,resp="US.rer")
plot(irf.girf.us.y,resp="US.rer")

## ---- "global.gdp",results="hide"---------------------------------------------
shocks<-list();shocks$var="y";shocks$cN<-c("EA", "US", "RU");shocks$ident="girf";shocks$scal=-1
irf.chol.ssvs<-irf(model.ssvs.1,shock=shocks,n.ahead=24)

## ---- hide=TRUE---------------------------------------------------------------
data("eerDataspf")
eerDataspf<-eerDataspf[cN]
W.trade0012.spf<-W.trade0012.spf[cN,cN]
W.trade0012.spf<-apply(W.trade0012.spf,2,function(x)x/rowSums(W.trade0012.spf))

## ---- "us.spf", results="hide"------------------------------------------------
model.ssvs.eer<-bgvar(Data=eerDataspf,
                      W=W.trade0012.spf,
                      plag=1,
                      draws=100,
                      burnin=100,
                      prior="SSVS",
                      SV=TRUE,
                      thin=1
                      )

## ---- "us.spf.sign.spec"------------------------------------------------------
 sign.constr.eer<-list()
 sign.constr.eer$shock1$shock<-"US.y" # Positive AD Shock, gdp goes up,
 sign.constr.eer$shock1$restrictions$res1<-"US.Dp" #inflation up and interest rates as well
 sign.constr.eer$shock1$sign<-c(">",">")
 sign.constr.eer$shock1$rest.horz<-c(1,1)
 sign.constr.eer$shock1$constr<-c(1,1)# no cross-country restrictions, set constr. to 1
 sign.constr.eer$shock1$scal<-1 #+1% increase
 sign.constr.eer$shock2$shock<-"US.Dp" # Negative AS shock, inflation up
 sign.constr.eer$shock2$restrictions$res1<-"US.y"
 sign.constr.eer$shock2$sign<-c(">","<")
 sign.constr.eer$shock2$rest.horz<-c(1,1)
 sign.constr.eer$shock2$constr<-c(1,1) # no cross-country restrictions, set constr. to 1
 sign.constr.eer$shock2$scal<-1 #+1% increase
 names(sign.constr.eer)<-c("AD","AS")
 sign.constr.eer$MaxTries<-10000 # maximum number of rotations per draw

## ---- "us.spf.sign",message=FALSE, results="hide"-----------------------------
irf.sign<-irf(model.ssvs.eer,n.ahead=24,sign.constr=sign.constr.eer,save.store=TRUE)

## ---- "us.spf.sign2"----------------------------------------------------------
irf.sign$rot.nr

## ----"us.spf.plots",fig.cap="Responses to AS (left) and AD (right) shock.",fig.show="hold",out.width="25%"----
plot(irf.sign,resp="US.y",shock.nr=1)
plot(irf.sign,resp="US.y",shock.nr=2)

plot(irf.sign,resp="US.Dp",shock.nr=1)
plot(irf.sign,resp="US.Dp",shock.nr=2)

plot(irf.sign,resp="US.rer",shock.nr=1)
plot(irf.sign,resp="US.rer",shock.nr=2)

## ---- "us.spf.sign3",results="hide"-------------------------------------------
sign.constr<-list()
# shock to 4-step ahead expectation of US.stir
sign.constr$shock1$shock<-"US.stir_t+4"
sign.constr$shock1$restrictions$res1<-"US.Dp_t+4"
# zero restriction on US.stir
sign.constr$shock1$restrictions$res2<-"US.stir"
sign.constr$shock1$restrictions$res3<-"US.y_t+4"
# rationality condition: US.stir_t+4 on impact is equal to average of IRF of 
# US.stir between horizon 1 and 4 (defined with rest.horz, but as period 5!)
sign.constr$shock1$restrictions$res4<-"US.stir_t+4" 
# rationality condition: US.Dp_t+4 on impact is equal to H-step ahead IRF 
# of US.Dp in horizon 4 (defined with rest.horz, but as period 5!)
sign.constr$shock1$restrictions$res5<-"US.Dp_t+4" 
# rationality condition: US.y_t+4 on impact is equal to H-step ahead IRF 
# of US.y in horizon 4 (defined with rest.horz, but as period 5!)
sign.constr$shock1$restrictions$res6<-"US.y_t+4"
sign.constr$shock1$sign<-c(">","<","0","<","ratio.avg","ratio.H","ratio.H")
sign.constr$shock1$rest.horz<-c(1,1,1,1,5,5,5)
sign.constr$shock1$constr<-c(1,1,1,1,1,1,1)
sign.constr$shock1$scal=1
sign.constr$MaxTries<-200
irf.sign.zero<-irf(model.ssvs.eer,n.ahead=20,sign.constr=sign.constr,save.store=TRUE)

## ---- "eer.spf.plots",fig.cap="Rationality conditions I.",out.width="50%",fig.show="hold"----
# rationality condition: US.stir_t+4 on impact is equal to average of IRF of 
# US.stir between horizon 2 and 5
matplot(cbind(irf.sign$posterior["US.stir_t+4",,1,"median"],
              irf.sign$posterior["US.stir",,1,"median"]),
        type="l",ylab="",main="stir",lwd=2); legend("topright",lty=c(1,2),c("expected","actual"),lwd=2,bty="n",col=c("black","red"))
abline(h=mean(irf.sign$posterior["US.stir",2:5,1,"median"]),lwd=2)
abline(v=c(2,5),lty=3,col="grey",lwd=2)

# rationality condition: US.y_t+4 on impact is equal to H-step ahead IRF 
# of US.y in horizon 5
matplot(cbind(irf.sign$posterior["US.y_t+4",,1,"median"],
              irf.sign$posterior["US.y",,1,"median"]),
        type="l",ylab="",main="y",lwd=2)
legend("topright",lty=c(1,2),c("expected","actual"),lwd=2,bty="n",col=c("black","red"))
abline(h=irf.sign$posterior["US.y_t+4",1,1,"median"],lwd=2)
abline(v=5,lty=3,col="grey",lwd=2)

## ---- "ea.data"---------------------------------------------------------------
data(monthlyData);monthlyData$OC<-NULL
names(monthlyData)
# list of weights of other entities with same name as additional country model
OE.weights = list(EB=EA.weights)
EA_countries <- c("AT", "BE", "DE","ES", "FI","FR")
                  # "IE", "IT", "NL", "PT","GR","SK","MT","CY","EE","LT","LV")

## ---- "restrict_sample", hide=TRUE--------------------------------------------
monthlyData <- monthlyData[c(EA_countries,"EB")]
W<-W[EA_countries,EA_countries]
W<-apply(W,2,function(x)x/rowSums(W))
OE.weights$EB$weights <- OE.weights$EB$weights[names(OE.weights$EB$weights)%in%EA_countries]

## ----"ea.estimate"------------------------------------------------------------
# estimates the model
model.ssvs<-bgvar(Data=monthlyData,
                  W=W,
                  draws=100,
                  burnin=100,
                  plag=1,
                  prior="SSVS",
                  thin=1,
                  eigen=1.05,
                  OE.weights=OE.weights)

## ----"ea.sign"----------------------------------------------------------------
# imposes sign restrictions on the cross-section and for a global shock
# (long-term interest rates)
sign.constr<-list()
# the variable to shock, can be imposed for more than one country
# but should be the same variable for all of them
sign.constr$shock1$shock = c(paste0(EA_countries[-c(3,12)],".ltir")) 
# restrictions (industrial production should decrease for selected countries)                         
sign.constr$shock1$restrictions$res1<-paste0(EA_countries,".y")
# another set of restrictions (inflation  should decrease for selected countries)
sign.constr$shock1$restrictions$res2<-paste0(EA_countries,".p") 
# first entry is for the shock, following entries for the restrictions 
# (ltir should go up, y and p go down)
sign.constr$shock1$sign<-c(">","<","<") 
# nr. of time periods restrictions are imposed, first entry is for the shock, 
# following entries for the restrictions
sign.constr$shock1$rest.horz<-c(1,1,1)
# are constraints binding for all (1) countries specified or for at least 50% of 
# the countries (0.5), or 75% (0.75)
sign.constr$shock1$constr<-c(1,0.5,0.5) 
# a minus 100 bp shock to long-term interest rates (on average)
sign.constr$shock1$scal=-100
# rotation specification
sign.constr$MaxTries<-200

## ----"global.restrictions"----------------------------------------------------
sign.constr$shock1$shock
sign.constr$shock1$restrictions

## ----"global.shock.irf",echo=TRUE,results="hide"------------------------------
irf.sign.ssvs<-irf(model.ssvs,             # gvar object
                   shock=NULL,             # shockc is not needed 
                   n.ahead=24,             # impulse response horizon
                   sign.constr=sign.constr # list with sign-restrictions
                   )

## ----"ea.sign.verify"---------------------------------------------------------
irf.sign.ssvs$posterior[paste0(EA_countries[-c(3,12)],".ltir"),1,1,"median"]
irf.sign.ssvs$posterior[paste0(EA_countries,".y"),1,1,"median"]
irf.sign.ssvs$posterior[paste0(EA_countries,".p"),1,1,"median"]

## ---- "ea.sign.plots",fig.show="hold",out.width="25%",fig.cap="Output responses of selected euro area countries."----
plot(irf.sign.ssvs,resp="AT.y")
plot(irf.sign.ssvs,resp="BE.y")

plot(irf.sign.ssvs,resp="DE.y")
plot(irf.sign.ssvs,resp="ES.y")

## ---- "fevd"------------------------------------------------------------------
#calculates the LN GFEVD 
gfevd.us.mp=gfevd(model.ssvs.eer,n.ahead=24,running=TRUE,cores=4)$FEVD

# get position of EA 
idx<-which(grepl("EA.",dimnames(gfevd.us.mp)[[2]]))
own<-colSums(gfevd.us.mp["EA.y",idx,])
foreign<-colSums(gfevd.us.mp["EA.y",-idx,])

## ---- "fevd.plot",fig.cap="FEVD of EA GDP.",out.width="50%"-------------------
barplot(t(cbind(own,foreign)),legend.text =c("own","foreign"))

## ---- "fevd.struc"------------------------------------------------------------
# calculates FEVD for variables US.y
fevd.us.y=fevd(irf.sign,var.slct=c("US.y"))$FEVD
idx<-which(grepl("US.",rownames(fevd.us.y)))


## ---- "fevd.struc.plot",fig.cap="FEVD of US GDP.",out.width="50%"-------------
barplot(fevd.us.y[idx,,])

## ----"hd"---------------------------------------------------------------------
 HD<-hd(irf.chol.us.mp)
 # summing them up should get you back the original time series
 org.ts<-apply(HD$hd_array,c(1,2),sum) # this sums up the contributions of all shocks + constant, initial conditions and residual component (last three entries in the third dimension of the array)

## ---- "hd.plot",fig.cap="Historical decomposition of euro area GDP.",out.width="50%"----
  matplot(cbind(HD$x[,1],org.ts[,1]),type="l",ylab="",lwd=2)
  legend("bottomright",c("hd series","original"),col=c("black","red"),lty=c(1,2),bty="n",cex=2)

## ----"fcast.est", results="hide"----------------------------------------------
model.ssvs.h8<-bgvar(Data=eerData,W=W.trade0012,plag=1,draws=100,burnin=100,prior="SSVS",
                    SV=TRUE,thin=1,eigen=1.05,trend=TRUE,h=8)

## ----"fcast.predict", results="hide"------------------------------------------
fcast <- predict(model.ssvs.h8, n.ahead=8,save.store=TRUE)

## ----"lps"--------------------------------------------------------------------
lps.h8 <- lps(fcast)
rmse.h8 <- rmse(fcast)

## ---- "fcast.plot",fig.cap="Forecast plot.",out.width="50%"-------------------
plot(fcast, resp="US.Dp", Cut=8)

## ----"cond.predict",results="hide"--------------------------------------------
# matrix with constraints
constr <- matrix(NA,nrow=fcast$n.ahead,ncol=ncol(model.ssvs.1$xglobal))
colnames(constr) <- colnames(model.ssvs.h8$xglobal)
# set "US.Dp" for five periods on its last value
constr[1:5,"US.Dp"] <-model.ssvs.h8$xglobal[nrow(model.ssvs.h8$xglobal),"US.Dp"]
# compute conditional forecast (hard restriction)
cond_fcast <- cond.predict(constr=constr, bgvar.obj=model.ssvs.h8, pred.obj=fcast, constr_sd=NULL)

## ----"cond.predict.sd",results="hide"-----------------------------------------
# add uncertainty to conditional forecasts
constr_sd <- matrix(NA,nrow=fcast$n.ahead,ncol=ncol(model.ssvs.h8$xglobal))
colnames(constr_sd) <- colnames(model.ssvs.h8$xglobal)
constr_sd[1:5,"US.Dp"] <- 0.001

# compute conditional forecast with soft restrictions
cond_fcast2 <- cond.predict(constr=constr, bgvar.obj=model.ssvs.h8, pred.obj=fcast, constr_sd=constr_sd)

## ---- "cond.plot.1",out.width="50%",fig.show="hold",fig.cap="Conditional forecast of US Inflation, top panel without uncertainty during the conditioning, bottom panel with uncertainty."----
plot(cond_fcast, resp="US.Dp", Cut=10)
plot(cond_fcast2, resp="US.Dp", Cut=10)

## ---- eval=FALSE--------------------------------------------------------------
#  # load dataset
#  data(eerData)
#  
#  # Minnesota prior and two different weight matrices and no SV
#  # weights for first variable set tradeW.0012, for second finW0711
#  variable.list      <- list()
#  variable.list$real <- c("y","Dp","tb")
#  variable.list$fin  <- c("stir","ltir","rer")
#  Hyperparm.MN <- list(a_i = 0.01, # prior for the shape parameter of the IG
#                         b_i = 0.01  # prior for the scale parameter of the IG
#                         )
#  model.MN<-bgvar(Data=eerData,
#                    W=W.list[c("tradeW.0012","finW0711")],
#                    draws=200,
#                    burnin=200,
#                    plag=1,
#                    hyperpara=Hyperparm.MN,
#                    prior="MN",
#                    thin=1,
#                    eigen=TRUE,
#                    SV=TRUE,
#                    variable.list=variable.list)
#  
#  # SSVS prior
#  Hyperparm.ssvs <- list(tau0   = 0.1,  # coefficients: prior variance for the spike
#                                        # (tau0 << tau1)
#                         tau1   = 3,    # coefficients: prior variance for the slab
#                                        # (tau0 << tau1)
#                         kappa0 = 0.1,  # covariances: prior variance for the spike
#                                        # (kappa0 << kappa1)
#                         kappa1 = 7,    # covariances: prior variance for the slab
#                                        # (kappa0 << kappa1)
#                         a_i    = 0.01, # prior for the shape parameter of the IG
#                         b_i    = 0.01, # prior for the scale parameter of the IG
#                         p_i    = 0.5,  # prior inclusion probability of coefficients
#                         q_ij   = 0.5   # prior inclusion probability of covariances
#                         )
#  model.ssvs<-bgvar(Data=eerData,
#                    W=W.trade0012,
#                    draws=100,
#                    burnin=100,
#                    plag=1,
#                    hyperpara=Hyperparm.ssvs,
#                    prior="SSVS",
#                    thin=1,
#                    eigen=TRUE)
#  
#  # Normal Gamma prior
#  data(monthlyData)
#  monthlyData$OC<-NULL
#  Hyperparm.ng<-list(c_tau    = 0.01, # covariances: prior hyperparameter for the NG-prior
#                     d_tau    = 0.01, # covariances: prior hyperparameter for the NG-prior
#                     e_lambda = 1.5,  # coefficients: prior hyperparameter for the NG-prior
#                     d_lambda = 1,    # coefficients: prior hyperparameter for the NG-prior
#                     prmean   = 0,    # prior mean for the first lag of the AR coefficients
#                     a_i      = 0.01, # prior for the shape parameter of the IG
#                     b_i      = 0.01, # prior for the scale parameter of the IG
#                     a_start  = .6,   # (hyper-)parameter for the NG
#                     sample_A = FALSE # estimate a?
#                     )
#  model.ng<-bgvar(Data=monthlyData,
#                  W=W,
#                  draws=200,
#                  burnin=100,
#                  plag=1,
#                  hyperpara=Hyperparm.ng,
#                  prior="NG",
#                  thin=2,
#                  eigen=TRUE,
#                  SV=TRUE,
#                  OE.weights=list(EB=EA.weights))

## ---- eval=FALSE--------------------------------------------------------------
#    # First example, a US monetary policy shock, quarterly data
#    library(BGVAR)
#    data(eerData)
#    model.ssvs.eer<-bgvar(Data=eerData,W=W.trade0012,draws=500,burnin=500,plag=1,prior="SSVS",thin=10,eigen=TRUE,trend=TRUE)
#    # US monetary policy shock
#    shocks<-list();shocks$var="stir";shocks$cN<-"US";shocks$ident="chol";shocks$scal=-100
#    irf.chol.us.mp<-irf(model.ssvs.eer,shock=shocks,n.ahead=24)
#  
#    # plots an impulse response function
#    plot(irf.chol.us.mp,resp="US.y")
#  
#   # calculates generalized impulse response functions for the same shock as above
#    shocks$ident="girf"
#    irf.girf.ssvs<-irf(model.ssvs.eer,shock=shocks,n.ahead=24)
#    plot(irf.girf.ssvs,resp="US.y")
#  
#  # Shock to first ordered variable yields same responses of cholesky and GIRF
#    shocks<-list();shocks$var="y";shocks$cN<-"US";shocks$ident="chol";shocks$scal=+1
#    irf.chol<-irf(model.ssvs.eer,shock=shocks,n.ahead=24)
#    shocks$ident<-"girf"
#    irf.girf<-irf(model.ssvs.eer,shock=shocks,n.ahead=24)
#    matplot(cbind(irf.chol$posterior["US.y",,1,"median"],irf.girf$posterior["US.y",,1,"median"]),type="l",ylab="")
#    matplot(cbind(irf.chol$posterior["US.Dp",,1,"median"],irf.girf$posterior["US.Dp",,1,"median"]),type="l",ylab="")
#    matplot(cbind(irf.chol$posterior["EA.y",,1,"median"],irf.girf$posterior["EA.y",,1,"median"]),type="l",ylab="")
#  
#  
#  
#    # second example, cross-country restrictions, multiple shocks and ECB country model
#    library(BGVAR)
#    data(monthlyData);monthlyData$OC<-NULL
#    OE.weights = list(EB=EA.weights)
#    # estimates the model
#    model.ssvs<-bgvar(Data=monthlyData,W=W,draws=500,burnin=500,plag=1,prior="SSVS",thin=10 ,eigen=TRUE,OE.weights=OE.weights)
#  
#    EA_countries <- c("AT", "BE", "DE","ES", "FI","FR", "IE", "IT", "NL", "PT","GR","SK") #,"MT","CY","EE","LT","LV")
#  
#    # A simultaneous Cholesky shock to long-term interest rates in the euro area countries, scaled to amount to -100 basis points (on average over the EA countries).
#    # Note that the ordering of the variables influences the response, the ordering is exactly as in the country models, to use a different order you have re-estimate
#    # the model (by bgvar)
#    shocks<-list();shocks$var="ltir";shocks$cN<-EA_countries;shocks$ident="chol";shocks$scal=-100
#    irf.chol.ssvs<-irf(model.ssvs,shock=shocks,n.ahead=24)
#  
#    # imposes sign restrictions on the cross-section and for a global shock (long-term interest rates)
#    sign.constr<-list()
#    sign.constr$shock1$shock<-c(paste0(EA_countries[-c(3,12)],".ltir")) # the variable to shock, can be imposed for more than one country
#                                                                        #but should be the same variable   for all of them
#    sign.constr$shock1$restrictions$res1<-paste0(EA_countries,".y") # restrictions (industrial production should decrease for selected countries)
#    sign.constr$shock1$restrictions$res2<-paste0(EA_countries,".p") # another set of restrictions (inflation  should decrease for selected countries)
#    sign.constr$shock1$sign<-c(">","<","<") # first entry is for the shock, following entries for the restrictions (ltir should go up, y and p go down)
#    sign.constr$shock1$rest.horz<-c(1,1,1) # nr. of time periods restrictions are imposed, first entry is for the shock, following entries for the restrictions
#    sign.constr$shock1$constr<-c(1,0.5,0.5) # are constraints binding for all (1) countries specified or for at least 50% of the countries (0.5), or 75% (0.75)
#    sign.constr$shock1$scal=-100 # a minus 100 bp shock to long-term interest rates (on average)
#    sign.constr$MaxTries<-200
#    irf.sign.ssvs<-irf(model.ssvs,shock=NULL,n.ahead=24,sign.constr=sign.constr)
#  
#  
#   # Same example but using a local (German) shock and cross-country restrictions.
#    # Note that the ordering of the variables influences the response, the ordering is exactly as in the country models, to use a different order you have re-estimate
#    # the model (by bgvar)
#    shocks<-list();shocks$var="ltir";shocks$cN<-EA_countries;shocks$ident="chol";shocks$scal=-100
#    irf.chol.ssvs<-irf(model.ssvs,shock=shocks,n.ahead=24)
#  
#    # imposes sign restrictions on the cross-section and for a global shock (long-term interest rates)
#    sign.constr<-list()
#    sign.constr$shock1$shock<-c("DE.ltir") # the variable to shock, can be imposed for more than one country
#                                                                        #but should be the same variable   for all of them
#    sign.constr$shock1$restrictions$res1<-paste0(EA_countries,".y") # restrictions (industrial production should decrease for selected countries)
#    sign.constr$shock1$restrictions$res2<-paste0(EA_countries,".p") # another set of restrictions (inflation  should decrease for selected countries)
#    sign.constr$shock1$sign<-c(">","<","<") # first entry is for the shock, following entries for the restrictions (ltir should go up, y and p go down)
#    sign.constr$shock1$rest.horz<-c(2,2,1) # nr. of time periods restrictions are imposed, first entry is for the shock, following entries for the restrictions
#    sign.constr$shock1$constr<-c(1,0.5,0.5) # are constraints binding for all (1) countries specified or for at least 50% of the countries (0.5), or 75% (0.75)
#    sign.constr$shock1$scal=-100 # a minus 100 bp shock to long-term interest rates (on average)
#    sign.constr$MaxTries<-200
#    irf.sign.ssvs<-irf(model.ssvs,shock=NULL,n.ahead=24,sign.constr=sign.constr)
#  
#   # Example with zero restriction (Arias et al., 2018) and rationality conditions (D'Amico and King, 2017).
#    data("eerDataspf")
#    model.ssvs.eer<-bgvar(Data=eerDataspf,
#                      W=W.trade0012.spf,
#                      plag=1,
#                      draws=500,
#                      burnin=500,
#                      prior="SSVS",
#                      SV=TRUE,
#                      thin=10)
#  
#    sign.constr<-list()
#    # shock to 4-step ahead expectation of US.stir
#    sign.constr$shock1$shock<-"US.stir_t+4"
#    sign.constr$shock1$restrictions$res1<-"US.Dp_t+4"
#    # zero restriction on US.stir
#    sign.constr$shock1$restrictions$res2<-"US.stir"
#    sign.constr$shock1$restrictions$res3<-"US.y_t+4"
#    # rationality condition: US.stir_t+4 on impact is equal to average of IRF of
#    # US.stir between horizon 1 and 4 (defined with rest.horz, but as period 5!)
#    sign.constr$shock1$restrictions$res4<-"US.stir_t+4"
#    # rationality condition: US.Dp_t+4 on impact is equal to H-step ahead IRF
#    # of US.Dp in horizon 4 (defined with rest.horz, but as period 5!)
#    sign.constr$shock1$restrictions$res5<-"US.Dp_t+4"
#    # rationality condition: US.y_t+4 on impact is equal to H-step ahead IRF
#    # of US.y in horizon 4 (defined with rest.horz, but as period 5!)
#    sign.constr$shock1$restrictions$res6<-"US.y_t+4"
#    sign.constr$shock1$sign<-c(">","<","0","<","ratio.avg","ratio.H","ratio.H")
#    sign.constr$shock1$rest.horz<-c(1,1,1,1,5,5,5)
#    sign.constr$shock1$constr<-c(1,1,1,1,1,1,1)
#    sign.constr$shock1$scal=1
#    sign.constr$MaxTries<-200
#    irf.sign<-irf(model.ssvs.eer,
#                  n.ahead=20,
#                  sign.constr=sign.constr,
#                  save.store=TRUE)
#    par(mfrow=c(4,1))
#   # rationality condition: US.stir_t+4 on impact is equal to average of IRF of   US.stir between horizon 1 and 4
#   matplot(cbind(irf.sign$posterior["US.stir_t+4",,1,"median"],irf.sign$posterior["US.stir",,1,"median"]),type="l",ylab="",main="stir")
#    abline(h=mean(irf.sign$posterior["US.stir",2:5,1,"median"]));abline(v=c(1,5),lty=3,col="grey")
#  
#   # rationality condition: US.y_t+4 on impact is equal to H-step ahead IRF of US.y in horizon 4
#    matplot(cbind(irf.sign$posterior["US.y_t+4",,1,"median"],irf.sign$posterior["US.y",,1,"median"]),type="l",ylab="",main="y")
#    abline(h=irf.sign$posterior["US.y_t+4",1,1,"median"]);abline(v=5,lty=3,col="grey")
#  
#   # rationality condition: US.Dp_t+4 on impact is equal to H-step ahead IRF of US.Dp in horizon 4
#    matplot(cbind(irf.sign$posterior["US.Dp_t+4",,1,"median"],irf.sign$posterior["US.Dp",,1,"median"]),type="l",ylab="",main="Dp")
#    abline(h=irf.sign$posterior["US.Dp_t+4",1,1,"median"]);abline(v=5,lty=3,col="grey")
#    par(mar=rep(0,4))
#    plot("1",type="n",axes=FALSE)
#    legend("center",c("expectation","actual"),lty=1:2,col=c("black","red"),bty="n",ncol=2)
#  
#   par(mar=c(5.1, 4.1, 4.1, 2.1))

## ---- hide=TRUE---------------------------------------------------------------
par(oldpar)

