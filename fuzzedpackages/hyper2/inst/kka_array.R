## This file creates a 3x3x3 array of results for the
## karpov/kasparov/anand dataset.

library("hyper2")
library("abind")
attach(as.list(kka))
        
players <- c("Anand","Karpov","Kasparov")

plays_white_wins <- matrix(NA,3,3)
dimnames(plays_white_wins) <- list(plays_white_wins=players,plays_black_loses=players)
plays_white_wins["Anand"   ,"Karpov"  ] <- anand_plays_white_beats_karpov
plays_white_wins["Anand"   ,"Kasparov"] <- anand_plays_white_beats_kasparov
plays_white_wins["Karpov"  ,"Anand"   ] <- karpov_plays_white_beats_anand
plays_white_wins["Karpov"  ,"Kasparov"] <- karpov_plays_white_beats_kasparov
plays_white_wins["Kasparov","Anand"   ] <- kasparov_plays_white_beats_anand
plays_white_wins["Kasparov","Karpov"  ] <- kasparov_plays_white_beats_karpov

plays_white_draws  <- matrix(NA,3,3)
dimnames(plays_white_draws) <- list(plays_white_draws=players,plays_black_draws=players)
plays_white_draws["Anand"   ,"Karpov"  ] <- anand_plays_white_draws_karpov
plays_white_draws["Anand"   ,"Kasparov"] <- anand_plays_white_draws_kasparov
plays_white_draws["Karpov"  ,"Anand"   ] <- karpov_plays_white_draws_anand
plays_white_draws["Karpov"  ,"Kasparov"] <- karpov_plays_white_draws_kasparov
plays_white_draws["Kasparov","Anand"   ] <- kasparov_plays_white_draws_anand
plays_white_draws["Kasparov","Karpov"  ] <- kasparov_plays_white_draws_karpov


plays_white_loses <- matrix(NA,3,3)
dimnames(plays_white_loses) <- list(plays_white_loses=players,plays_black_wins=players)
plays_white_loses["Karpov"   ,"Anand"   ] <- karpov_plays_white_losesto_anand
plays_white_loses["Kasparov" ,"Anand"   ] <- kasparov_plays_white_losesto_anand
plays_white_loses["Anand"    ,"Karpov"  ] <- anand_plays_white_losesto_karpov
plays_white_loses["Kasparov" ,"Karpov"  ] <- kasparov_plays_white_losesto_karpov
plays_white_loses["Anand"    ,"Kasparov"] <- anand_plays_white_losesto_kasparov
plays_white_loses["Karpov"   ,"Kasparov"] <- karpov_plays_white_losesto_kasparov

kka_array <- abind(
    plays_white_wins,
    plays_white_draws,
    plays_white_loses,
    along=3)
detach(as.list(kka))
