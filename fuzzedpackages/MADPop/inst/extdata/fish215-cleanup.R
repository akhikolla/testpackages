  # A few formatting steps
fish215 <- read.csv(file = "fish_study_215.csv")
fish215.og <- fish215
alleles <- as.matrix(fish215[,-1])
# give random allele names.
allele.names <- table(c(alleles))
nA <- ceiling(sqrt(length(allele.names)))
allele.names <- sort(sample(letters, 2*nA))
allele.names <- expand.grid(1:nA, allele.names[nA+1:nA])
allele.names <- paste0(allele.names[,2], ".", allele.names[,1])
# allele.names <- sample(allele.names)
allele.names <- allele.names[as.matrix(fish215[,-1])]
allele.names[is.na(allele.names)] <- "" # for data.frame later
fish215 <- cbind(Fish = as.character(fish215$Fish),
                 matrix(allele.names, nrow(fish215)))
# column "Fish" has both lake and fish ids. get rid of fish ids
pop.names <- fish215[,1]
pop.names <- sapply(pop.names,
                    function(str) strsplit(x = str, split = " ")[[1]][1])
fish215[,1] <- pop.names
head(fish215)

# these should be identical except for the name of the alleles
s1 <- HW.suff(X = fish215.og[,-1], popId = fish215[,1])
s2 <- HW.suff(X = fish215[,-1], popId = fish215[,1])
identical(s1[-1], s2[-1])

# mix the NAs and the order of the labels
fish215 <- cbind(Fish = fish215[,1],
                      t(apply(fish215[,-1], 1, function(x) sample(x))))
# convert to a data.frame.  most data will be passed in this fashion.
fish215 <- data.frame(fish215)
# Label the columns
# posterior samples should not depend on these
colnames(fish215) <- c("Lake", "A1", "A2", "A3", "A4")
head(fish215)


save("fish215", file = "data/fish215.RData")
