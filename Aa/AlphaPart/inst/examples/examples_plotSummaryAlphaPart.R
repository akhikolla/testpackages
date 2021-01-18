\donttest{

## Partition additive genetic values by country
(res <- AlphaPart(x=AlphaPart.ped, colPath="country", colBV=c("bv1", "bv2")))

## Summarize population by generation (=trend)
(ret <- summary(res, by="gen"))

## Plot the partitions
p <- plot(ret, ylab=c("bv for trait 1", "bv for trait 2"), xlab="Generation")
print(p[[1]]$abs)
print(p[[2]]$abs)
print(p)

## Partition additive genetic values by country and sex
AlphaPart.ped$country.gender <- with(AlphaPart.ped, paste(country, gender, sep="-"))
(res <- AlphaPart(x=AlphaPart.ped, colPath="country.gender", colBV=c("bv1", "bv2")))

## Summarize population by generation (=trend)
(ret <- summary(res, by="gen"))

## Plot the partitions
p <- plot(ret, ylab=c("BV for trait 1", "BV for trait 2"), xlab="Generation")
print(p)
p <- plot(ret, ylab=c("BV for trait 1", "BV for trait 2"), xlab="Generation",
        lineTypeList=list("-1"=1, "-2"=2, def=3))
print(p)
p <- plot(ret, ylab=c("BV for trait 1", "BV for trait 2"), xlab="Generation",
        lineTypeList=list("-1"=1, "-2"=2, def=3), useGgplot2=FALSE, useDirectLabels = FALSE)
print(p)

## Plot control (color and type of lines + limits)
p <- plot(ret, ylab=c("BV for trait 1", "BV for trait 2"), xlab="Generation",
        useGgplot2=TRUE, color=c("green", "gray"), lineType=c(2, 3),
        sortValue=FALSE, lineSize=4,
        xlim=c(-1, 7))
print(p)
}
