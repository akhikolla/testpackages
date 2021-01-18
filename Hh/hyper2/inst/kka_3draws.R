## Analysis of three chess players (Kasparov, Karpov, Anand).  See
## file karpov_kasparov_anand.R for details on the origin of the data.

## This is a more sophisticated version (with more parameters) than
## inst/karpov_kasparov.R: it allows for each player to have a
## distinct draw monster.

## This file creates likelihood function 'H' 


library("hyper2")
H <- hyper2(pnames=c("Karpov","Kasparov","Anand","white","Karpov_draw","Kasparov_draw","Anand_draw"))

results <- as.list(kka)
attach(results)


karpov_vs_kasparov <- c("Karpov","Kasparov","Karpov_draw","Kasparov_draw","white")
draw1 <- c("Karpov_draw","Kasparov_draw")

H %<>% trial(c("Karpov"  ,"white"), karpov_vs_kasparov, karpov_plays_white_beats_kasparov)
H %<>% trial(c("Kasparov","white"), karpov_vs_kasparov, kasparov_plays_white_beats_karpov)
H %<>% trial(c("Karpov"  )        , karpov_vs_kasparov, kasparov_plays_white_losesto_karpov)
H %<>% trial(c("Kasparov")        , karpov_vs_kasparov, karpov_plays_white_losesto_kasparov)
H %<>% trial(draw1                , karpov_vs_kasparov, karpov_plays_white_draws_kasparov)
H %<>% trial(draw1                , karpov_vs_kasparov, kasparov_plays_white_draws_karpov)

## Kasparov vs Anand
kasparov_vs_anand <- c("Kasparov","Anand","Kasparov_draw","Anand_draw","white")
draw2 <- c("Kasparov_draw","Anand_draw")
H %<>% trial(c("Kasparov","white"), kasparov_vs_anand, kasparov_plays_white_beats_anand)
H %<>% trial(c("Anand"   ,"white"), kasparov_vs_anand, anand_plays_white_beats_kasparov)
H %<>% trial(c("Kasparov"  )      , kasparov_vs_anand, anand_plays_white_losesto_kasparov)
H %<>% trial(c("Anand")           , kasparov_vs_anand, kasparov_plays_white_losesto_anand)
H %<>% trial(draw2                , kasparov_vs_anand, kasparov_plays_white_draws_anand)
H %<>% trial(draw2                , kasparov_vs_anand, anand_plays_white_draws_kasparov)


## Karpov vs Anand
karpov_vs_anand <- c("Karpov","Anand","Karpov_draw","Anand_draw","white")
draw3 <- c("Karpov_draw","Anand_draw")

H %<>% trial(c("Karpov","white"), karpov_vs_anand, karpov_plays_white_beats_anand)
H %<>% trial(c("Anand" ,"white"), karpov_vs_anand, anand_plays_white_beats_karpov)
H %<>% trial(c("Karpov"        ), karpov_vs_anand, anand_plays_white_losesto_karpov)
H %<>% trial(c("Anand"         ), karpov_vs_anand, karpov_plays_white_losesto_anand)
H %<>% trial(draw3              , karpov_vs_anand, karpov_plays_white_draws_anand) 
H %<>% trial(draw3              , karpov_vs_anand, anand_plays_white_draws_karpov) 

detach(results)

stopifnot(H == kka_3draws)


## Test the hypothesis that all three players have the same strength

## First do the free optimization:
max_support_free <- maxp(H,give=TRUE)$value
ml_p_free    <- maxp(H)

## Now the constrained optimization.  We enforce that
## Karpov==Anand==Kasparov==p but allow the white ghost and the draw
## monster to range freely, subject to the unit sum constraint,


objective <- function(x){
  loglik(H,c(x[1],x[1],x[1],x[2:4]))
}

constrained_optimization <-
  constrOptim(
      theta = rep(0.1,4),
      f     = objective,
      grad  =  NULL,
      ui    = rbind(diag(4),-c(3,1,1,1)),
      ci    = c(0,0,0,0,-1),
      control=list(fnscale= -1)
    )

max_support_constrained <- constrained_optimization$value

jj <- constrained_optimization$par
jj <- fillup(c(jj[c(1,1,1,2:4)]))
names(jj) <- pnames(H)
ml_p_constrained <- jj

support_different_strengths <- max_support_free - max_support_constrained

cat(paste("support for different strengths = ", support_different_strengths, "\n",sep=""))
cat(paste("likelihood ratio = ", exp(support_different_strengths), "\n",sep=""))
if(support_different_strengths > 2){
  cat("two units of support criterion exceeded: strong evidence that the three players have different strengths\n")
} else {
  cat("less than two units of support: no evidence for differing players' strengths\n")
}

cat(paste("p-value = ",pchisq(2*support_different_strengths,df=1,lower.tail=FALSE)))
cat("\n\n")

## Now test the hypothesis that playing white confers no advantage:


small <- 1e-4

objective <- function(x){
  loglik(H,c(x[1:3],white=small,x[4:5])) # Anand_draw is the fillup, white advantage set to 'small'
}

constrained_optimization_nowhite <-
  constrOptim(
      theta = rep(0.1,5),
      f       = objective,
      grad    =  NULL,
      ui      = rbind(diag(5),-1),
      ci      = c(0,0,0,0,0,small-1),
      control = list(fnscale= -1)
  )

jj <- constrained_optimization_nowhite$par
jj <- fillup(c(jj[1:3],small,jj[4:5]))
names(jj) <- pnames(H)
constrained_optimization_nowhite$par <- jj

max_support_nowhite <- constrained_optimization_nowhite$value
support_no_white_advantage <- max_support_free - max_support_nowhite

cat(paste("support for white advantage = ", support_no_white_advantage, "\n",sep=""))
cat(paste("likelihood ratio = ", exp(support_no_white_advantage), "\n",sep=""))
if(support_no_white_advantage > 2){
cat("two units of support criterion exceeded: strong evidence that playing white is an advantage\n")
} else {
  cat("less than two units of support: no evidence for white being an advantage\n")
}

cat(paste("p-value = ",pchisq(2*support_no_white_advantage,df=1,lower.tail=FALSE)))
cat("\n\n")
