## --- Partition additive genetic values by loc ---
res <- AlphaPart(x=AlphaPart.ped, colPath="country", colBV=c("bv1", "bv2"))

## Summarize whole population
ret <- summary(res)

## Summarize population by generation (=trend)
ret <- summary(res, by="gen")

## Summarize population by generation (=trend) but only for domestic location
ret <- summary(res, by="gen", subset=res[[1]]$country == "domestic")

## --- Partition additive genetic values by loc and gender ---

AlphaPart.ped$country.gender <- with(AlphaPart.ped, paste(country, gender, sep="-"))
res <- AlphaPart(x=AlphaPart.ped, colPath="country.gender", colBV=c("bv1", "bv2"))

## Summarize population by generation (=trend)
ret <- summary(res, by="gen")

## Summarize population by generation (=trend) but only for domestic location
ret <- summary(res, by="gen", subset=res[[1]]$country == "domestic")
