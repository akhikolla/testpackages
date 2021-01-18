## https://www.youtube.com/watch?v=XKWzlG4jDnI

library(hyper2)
library(partitions) ## needed for perms()
library(magrittr) 

team1  <- c("autimatic","tarik", "Skadoodle","Stewie2k","RUSH")   #Cloud9
team2 <- c("NiKo","olofmeister","karrigan","GuardiaN","rain")  # FaZe Clan

## This script creates a hyper2 object H which is a likelihood function
## for the strengths of the players in team1 and team2 above.  It defines
## a function counterstrike_maker() which is called in a loop at the bottom
## of the script.   The dataset is a list called zachslist, defined in
## this file below. 

## Function counterstrike_maker() needs the identities of all players
## on each team.  It assumes that players are always killed by the
## opposing team; if a player is killed by his own team,
## counterstrike_maker() may return an error.

## In the function, the 'deathorder' argument specifies the order in
## which players were killed (the first element is the first player to
## be killed and so on).  Note that the identity of the killer is not
## needed as it is assumed that a death is due to the combined
## strength of the team, rather than the individual shooter who
## actually fired the shot.

## Note that it is not required for all the players to be killed;
## short rounds are OK (but are not as statistically informative).

## In the function, the main loop iterates through the deathorder
## vector.  Suppose team1 = (a,b,c,d) and team2 = (x,y,z); suppose the
## deathorder were (a,b,y,c).  Then the likelihood function for a's
## death would be (x+y+z)/(a+b+c+d+x+y+z) [that is, the likelihood
## function for a Bernoulli trial between team (x,y,z) and (a,b,c,d).
## Note the appearance of (a) in the denominator: he [surely!] was on
## the losing team.

## The '%<>%' line in the main loop removes the killed player from the
## appropriate team.

## The likelihood function for the full deathorder is then

## (x+y+z)/(a+b+c+d+x+y+z) * (x+y+z)/(b+c+d+x+y+z) * (c+d)/(c+d+x+y+z) * (x+z)/(c+d+x+z)

## Where the four terms correspond to the deaths of a,b,y,c
## respectively.  See how the denominator gets shorter as the teams
## die one by one.

## The function does not have strong error-checking functionality.

## The dataset [here, `zachslist`] is a list whose six elements
## correspond to six rounds of play, which are assumed to be
## statistically independent.  Object H corresponds to an overall
## likelihood function for all the rounds combined: it is created by
## iterating through zachslist and incrementing the likelihood
## function for each round.

## File man/counterstrike.Rd has more details on the data's origin.


## Object Hrand is a randomly generated version of H, created by
## running an in-silico deathmatch on the assumption of equal player
## strengths.

`counterstrike_maker` <- function(team1,team2,deathorder){

  if(identical(sort(c(team1,team2)),sort(deathorder))){
    deathorder <- deathorder[-length(deathorder)]  # last player not killed
  }

  H <- hyper2(pnames=c(team1,team2))
  
  for(killed_player in deathorder){
    if(killed_player %in% team1){
      H[team2] <- H[team2] + 1
      H[c(team1,team2)] <- H[c(team1,team2)] - 1
      team1 %<>% "["(team1 != killed_player)      
    } else {
      H[team1] <- H[team1] + 1
      H[c(team1,team2)] <- H[c(team1,team2)] - 1
      team2 %<>% "["(team2 != killed_player)      
    }
  }
  return(H)
}

`rdeath` <- function(team1,team2){
  out <- NULL
  while(length(team1)>0 & length(team2)>0){
    if(runif(1)>length(team1)/(length(c(team1,team2)))){ ## killed person is on team1
      killed_player <- sample(team1,1)
      team1 %<>% "["(team1 != killed_player)      
    } else { ## killed person is on team2
      killed_player <- sample(team2,1)
      team2 %<>% "["(team2 != killed_player) 
    }
    out <- c(out,killed_player)
  }
  return(out)
}

## all rounds from the match, data kindly supplied by Zachary Hankin:
zachslist <- list(
    round1 =
      c("Skadoodle","olofmeister","tarik",
        "GuardiaN", "RUSH", "rain","Stewie2k",
        "karrigan","autimatic","NiKo"
        ),
    round2 =
      c("karrigan", "NiKo", "Stewie2K",
        "RUSH", "rain", "GuardiaN",
        "autimatic", "olofmeister"
        ),
    round3 =
      c("rain","tarik","autimatic",
        "karrigan","RUSH","GuardiaN",
        "Stewie2K","NiKo","olofmeister"
        ),
    round4 =
      c("rain","GuardiaN","karrigan",
        "NiKo","olofmeister"
        ),
    round5 =
      c("olofmeister","rain","karrigan",
        "tarik","Stewie2K","autimatic"
        ),
    round6 =
      c("GuardiaN","karrigan")
)

H <- hyper2(pnames=c(team1,team2))
for(i in zachslist){
  H <- H + counterstrike_maker(team1,team2, deathorder=i)
}

dotchart(maxp(H),pch=16,main='observed data')

Hrand <- hyper2(pnames=c(team1,team2))
for(i in 1:6){
  pp <- rdeath(team1,team2)
  Hrand <- Hrand + counterstrike_maker(team1,team2, pp)
}
dev.new()
dotchart(maxp(Hrand),pch=16,main='synthetic data')

