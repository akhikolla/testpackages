## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.align = "center",
  fig.height = 4
)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(data.frame(Job=c("Craftsman","Technician","Executive","Office employee"),Habitation = c("Owner","Renter","Starter","By family"),Time_in_job = c(10,20,5,2), Children = c(0,1,2,3), Family_status=  c("Divorced","Widower","Single","Married"),Default = c("No","No","Yes","No")))

## ---- echo=TRUE, results='asis'-----------------------------------------------
x = matrix(runif(1000), nrow = 1000, ncol = 1)
p = 1/(1+exp(-3*x^5))
y = rbinom(1000,1,p)
modele_lin <- glm(y ~ x, family = binomial(link="logit"))
pred_lin <- predict(modele_lin,as.data.frame(x),type="response")
pred_lin_logit <- predict(modele_lin,as.data.frame(x))

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(data.frame(True_prob = p,Pred_lin = pred_lin)))

## ---- echo=TRUE, results='asis'-----------------------------------------------
x_disc <- factor(cut(x,c(-Inf,0.5,0.7,0.8,0.9,+Inf)),labels = c(1,2,3,4,5))
modele_disc <- glm(y ~ x_disc, family = binomial(link="logit"))
pred_disc <- predict(modele_disc,as.data.frame(x_disc),type="response")
pred_disc_logit <- predict(modele_disc,as.data.frame(x_disc))


## ---- echo=FALSE--------------------------------------------------------------

knitr::kable(head(data.frame(True_prob = p,Pred_lin = pred_lin,Pred_disc = pred_disc)))
plot(x,3*x^5,main = "Estimated logit transform of p(Y|X)", ylab = "p(Y|X) under different models")
lines(x,pred_lin_logit,type="p",col="red")
lines(x,pred_disc_logit,type="p",col="blue")


## ---- echo=TRUE, results='asis'-----------------------------------------------
x_disc_bad_idea <- factor(cut(x,c(-Inf,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,+Inf)),labels = c(1,2,3,4,5,6,7,8,9,10))

## ---- echo=FALSE, results='asis'----------------------------------------------
liste_coef <- list()

for (k in 1:10) {
     x_part <- factor(x_disc_bad_idea[((k-1)*nrow(x)/10 +1) : (k/10*nrow(x))])
     y_part <- y[((k-1)*length(y)/10 +1) : (k/10*length(y))]
     modele_part <- glm(y_part ~ x_part, family=binomial(link = "logit"))
     liste_coef[[k]] <- (modele_part$coefficients)
}

estim_coef <- matrix(NA, nrow = nlevels(x_disc_bad_idea), ncol = 10)

for (i in 1:nlevels(x_disc_bad_idea)) {
     estim_coef[i,] <- unlist(lapply(liste_coef,function(batch) batch[paste0("x_part",levels(factor(x_disc_bad_idea))[i])]))
}

stats_coef <- matrix(NA, nrow = nlevels(x_disc_bad_idea), ncol = 3)

for (i in 1:nlevels(x_disc_bad_idea)) {
     stats_coef[i,1] <- mean(estim_coef[i,], na.rm = TRUE)
     stats_coef[i,2] <- sd(estim_coef[i,], na.rm = TRUE)
     stats_coef[i,3] <- sum(is.na(estim_coef[i,]))
}

stats_coef <- stats_coef[-1,] 
row.names(stats_coef) <- levels(x_disc_bad_idea)[2:nlevels(x_disc_bad_idea)]

plot (row.names(stats_coef), stats_coef[,1],ylab="Estimated coefficient",xlab="Factor value of x", ylim = c(-1,8))
segments(as.numeric(row.names(stats_coef)), stats_coef[,1]-stats_coef[,2],as.numeric(row.names(stats_coef)),stats_coef[,1]+stats_coef[,2])
lines(row.names(stats_coef),rep(0,length(row.names(stats_coef))),col="red")

## ---- echo=TRUE, results='asis'-----------------------------------------------
x = matrix(runif(300), nrow = 100, ncol = 3)
cuts = seq(0,1,length.out= 4)
xd = apply(x,2, function(col) as.numeric(cut(col,cuts)))
theta = t(matrix(c(0,0,0,2,2,2,-2,-2,-2),ncol=3,nrow=3))
log_odd = rowSums(t(sapply(seq_along(xd[,1]), function(row_id) sapply(seq_along(xd[row_id,]),
function(element) theta[xd[row_id,element],element]))))
y = rbinom(100,1,1/(1+exp(-log_odd)))

## ---- echo=TRUE,warning=FALSE, message=FALSE, results='hide',eval=FALSE-------
#  library(glmdisc)
#  set.seed(123)
#  discretization <- glmdisc(x,y,iter=50,m_start=5,test=FALSE,validation=FALSE,criterion="aic",interact=FALSE)

## ---- echo=FALSE,warning=FALSE, message=FALSE, results='hide',eval=TRUE-------
library(glmdisc)
set.seed(1)
discretization <- glmdisc(x,y,iter=50,m_start=5,test=FALSE,validation=FALSE,criterion="aic",interact=FALSE)

## ---- echo=TRUE,warning=FALSE, message=FALSE, results='hide',eval=TRUE--------
all_formula <- list()

for (l in 1:10) {
     set.seed(l)
     x = matrix(runif(300), nrow = 100, ncol = 3)
     cuts = seq(0,1,length.out= 4)
     xd = apply(x,2, function(col) as.numeric(cut(col,cuts)))
     theta = t(matrix(c(0,0,0,2,2,2,-2,-2,-2),ncol=3,nrow=3))
     log_odd = rowSums(t(sapply(seq_along(xd[,1]), function(row_id) sapply(seq_along(xd[row_id,]),
     function(element) theta[xd[row_id,element],element]))))
     y = rbinom(100,1,1/(1+exp(-log_odd)))
     
     discretization <- glmdisc(x,y,iter=50,m_start=5,test=FALSE,validation=FALSE,criterion="aic",interact=TRUE)
     all_formula[[l]] <- discretization@best.disc$formulaOfBbestestLogisticRegression
}

#barplot(table(grepl(":",all_formula)),names.arg=c("No interaction","One interaction"),xlim = 2)


## ---- echo=TRUE,warning=FALSE, message=FALSE, results='hide'------------------

all_formula <- list()

for (l in 1:10) {
     set.seed(100+l)
     x = matrix(runif(3000), nrow = 1000, ncol = 3)
     cuts = seq(0,1,length.out= 4)
     xd = apply(x,2, function(col) as.numeric(cut(col,cuts)))
     theta = t(matrix(c(0,0,0,2,2,2,-2,-2,-2),ncol=3,nrow=3))
     log_odd = matrix(0,1000,1)
       for (i in 1:1000) {
          log_odd[i] = 2*(xd[i,1]==1)+
               (-2)*(xd[i,1]==2)+
               2*(xd[i,1]==1)*(xd[i,2]==1)+
               4*(xd[i,1]==2)*(xd[i,2]==2)+
               (-2)*(xd[i,2]==1)*(xd[i,3]==2)+
               (-4)*(xd[i,2]==2)*(xd[i,3]==1)+
               1*(xd[i,2]==3)*(xd[i,3]==1)+
               (-1)*(xd[i,1]==1)*(xd[i,3]==3)+
               3*(xd[i,1]==3)*(xd[i,3]==2)+
               (-1)*(xd[i,1]==2)*(xd[i,3]==3)
       }
     
     y = rbinom(1000,1,1/(1+exp(-log_odd)))
     
     discretization <- glmdisc(x,y,iter=50,m_start=5,test=FALSE,validation=FALSE,criterion="aic",interact=TRUE)
     all_formula[[l]] <- discretization@best.disc$formulaOfBbestestLogisticRegression
}

#barplot(table(grepl(":",all_formula)),names.arg=c("No interaction","One interaction"),xlim = 2)


## ---- echo=TRUE---------------------------------------------------------------
discretization@parameters

## ---- echo=TRUE---------------------------------------------------------------
discretization@best.disc[[1]]

# The first link function is:
discretization@best.disc[[2]][[1]]

## ---- echo=TRUE---------------------------------------------------------------
discretization@performance[[1]]

## ---- echo=TRUE---------------------------------------------------------------
discretization@performance[[2]][1:5]

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  discretization@disc.data

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(discretization@disc.data))

## ---- echo=TRUE,eval=FALSE----------------------------------------------------
#  discretization@cont.data

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(discretization@cont.data))

## ---- echo=FALSE--------------------------------------------------------------
plot(x[,1],xd[,1])
plot(discretization@cont.data[,1],discretization@disc.data[,1])

## ---- echo=TRUE---------------------------------------------------------------
print(discretization)

## ---- echo=TRUE---------------------------------------------------------------
show(discretization)

## ---- echo=TRUE, results='asis'-----------------------------------------------
x_new <- discretize(discretization,x)

## ---- echo=FALSE, warning=FALSE-----------------------------------------------
knitr::kable(head(x_new))

## ---- echo=TRUE---------------------------------------------------------------
# pred_new <- predict(discretization, data.frame(x))

## ---- echo=FALSE, warning=FALSE-----------------------------------------------
# knitr::kable(head(pred_new))

