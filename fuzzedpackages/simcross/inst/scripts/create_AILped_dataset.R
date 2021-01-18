# create a dataset from the AIL pedigree in the QTLRel package

library(QTLRel)
data(miscEx)
AILped <- cbind(id=pedF8$id,
                mom=pedF8$dam,
                dad=pedF8$sire,
                sex=ifelse(pedF8$sex=="M", 1, 0),
                generation=as.numeric(substr(pedF8$generation, 2, 2)))

# renumber individuals
newid <- c(0,order(AILped[,"id"]))
names(newid) <- c("0", as.character(AILped[,"id"]))
for(i in 1:3)
    AILped[,i] <- newid[as.character(AILped[,i])]

# check that the change in sex is correct
newsex <- factor(c("F", "M")[AILped[,"sex"]+1], levels=c("F", "M"))
library(testthat)
expect_equal(pedF8$sex, newsex)

# check that the renumbering is correct
mom <- match(AILped[,2], AILped[,1])
momo <- match(pedF8$dam, pedF8$id)
expect_equal(mom, momo)

dad <- match(AILped[,3], AILped[,1])
dado <- match(pedF8$sire, pedF8$id)
expect_equal(dad, dado)

# generation the same
gen <- factor(paste0("F", AILped[,"generation"]))
expect_equal(gen, pedF8$generation)

# save to file
AILped <- as.data.frame(AILped)
save(AILped, file="../../data/AILped.RData")
