## Analysis of three chess players (Kasparov, Karpov, Anand).  See
## file karpov_kasparov_anand.R for details on the origin of the data.

## This is a more sophisticated version (with more parameters) than
## inst/karpov_kasparov.R: it allows each player to have a distinct
## white strength.

## This file creates likelihood function 'H', which is identical to
## dataset kka_3whites in the package.


library("hyper2")
H <- hyper2(pnames=c("Karpov","Kasparov","Anand","Karpov_white","Kasparov_white","Anand_white","draw"))

results <- as.list(kka)
attach(results)

D <- "draw"

## First: Karpov vs Kasparov
karpov_plays_white_vs_kasparov <- c("Karpov","Kasparov","Karpov_white","draw"  )  # "players" 
kasparov_plays_white_vs_karpov <- c("Karpov","Kasparov","Kasparov_white","draw")
H %<>% trial(c("Karpov"  ,"Karpov_white"  ), karpov_plays_white_vs_kasparov, karpov_plays_white_beats_kasparov  ) # Karpov wins playing white
H %<>% trial(c("Kasparov","Kasparov_white"), kasparov_plays_white_vs_karpov, kasparov_plays_white_beats_karpov  ) # Kasparov wins playing white
H %<>% trial(c("Kasparov")                 , karpov_plays_white_vs_kasparov, karpov_plays_white_losesto_kasparov) # Kasparov wins playing black
H %<>% trial(c("Karpov"  )                 , kasparov_plays_white_vs_karpov, kasparov_plays_white_losesto_karpov) # Karpov wins playing black
H %<>% trial(D                             , karpov_plays_white_vs_kasparov, karpov_plays_white_draws_kasparov  ) # Karpov white, draws
H %<>% trial(D                             , kasparov_plays_white_vs_karpov, kasparov_plays_white_draws_karpov  ) # Kasparov white, draws

## Second: Karpov vs Anand
karpov_plays_white_vs_anand <- c("Karpov","Anand","Karpov_white","draw"  )
anand_plays_white_vs_karpov <- c("Karpov","Anand","Anand_white","draw")
H %<>% trial(c("Karpov","Karpov_white"), karpov_plays_white_vs_anand, karpov_plays_white_beats_anand  ) # Karpov wins playing white
H %<>% trial(c("Anand" ,"Anand_white" ), anand_plays_white_vs_karpov, anand_plays_white_beats_karpov  ) # Anand wins playing white
H %<>% trial(c("Anand"                ), karpov_plays_white_vs_anand, karpov_plays_white_losesto_anand) # Anand wins playing black
H %<>% trial(c("Karpov"               ), anand_plays_white_vs_karpov, anand_plays_white_losesto_karpov) # Karpov wins playing black
H %<>% trial(D                         , karpov_plays_white_vs_anand, karpov_plays_white_draws_anand  ) # Karpov white, draws
H %<>% trial(D                         , anand_plays_white_vs_karpov, anand_plays_white_draws_karpov  ) # Anand white, draws

## Third: Kasparov vs Anand
anand_plays_white_vs_kasparov <- c("Anand","Kasparov","Anand_white","draw"   )
kasparov_plays_white_vs_anand <- c("Anand","Kasparov","Kasparov_white","draw")
H %<>% trial(c("Kasparov","Kasparov_white"), kasparov_plays_white_vs_anand, kasparov_plays_white_beats_anand  ) # Kasparov wins playing white
H %<>% trial(c("Anand"   ,"Anand_white"   ), anand_plays_white_vs_kasparov, anand_plays_white_beats_kasparov  ) # Anand wins playing white
H %<>% trial(c("Anand"                    ), kasparov_plays_white_vs_anand, kasparov_plays_white_losesto_anand) # Anand wins playing black
H %<>% trial(c("Kasparov"                 ), anand_plays_white_vs_kasparov, anand_plays_white_losesto_kasparov) # Kasparov wins playing black
H %<>% trial(D                             , kasparov_plays_white_vs_anand, kasparov_plays_white_draws_anand  ) # Kasparov white, draws
H %<>% trial(D                             , anand_plays_white_vs_kasparov, anand_plays_white_draws_kasparov  ) # Anand white, draws

detach(results)

stopifnot(H == kka_3whites)


## Test the hypothesis that all three players have the same strength

## First do the free optimization:
max_support_free <- maxp(H,give=TRUE)$value
ml_p_free    <- maxp(H)

## Now the constrained optimization.  

objective <- function(x){ # p=Kr, Ks, An, KrW, KsW, AnW (draw is the fillup)
  loglik(H,x[c(1,2,3,4,4,4)])
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
