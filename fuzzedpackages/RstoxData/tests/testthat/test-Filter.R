# Satisfy R CMD check
options("mc.cores" = 2)

# NMD Biotic
filenames <- system.file("testresources","biotic_v3_example.xml", package="RstoxData")

inputData <- ReadBiotic(filenames)

# test a wrong expression
filterExpression <- list()
filterExpression$`biotic_v3_example.xml`$individual <- c(
	'weight > 2'
)

context("test-Filter: Invalid column")
# Should be warning but can proceed)
expect_warning(out <-filterData(inputData, filterExpression))

# Therefore should be equal
expect_equal(inputData, out)


# fix and add expressions
filterExpression$`biotic_v3_example.xml`$individual <- c(
	'individualweight > 2'
)
filterExpression$`biotic_v3_example.xml`$agedetermination <- c(
	'age < 10',
	'serialnumber %notin% c(99484)'
)

# Should be OK
out <-filterData(inputData, filterExpression)

context("test-Filter: Filtering Biotic Data")
# Should be TRUE
expect_equal(any(out$`biotic_v3_example.xml`$individual$individualweight <= 2), FALSE)
expect_equal(any(out$`biotic_v3_example.xml`$agedetermination$age >= 10), FALSE)
expect_equal(any(out$`biotic_v3_example.xml`$agedetermination$serialnumber == 99484), FALSE)
expect_equal(all(out$`biotic_v3_example.xml`$individual$individualweight > 2), TRUE)
expect_equal(all(out$`biotic_v3_example.xml`$agedetermination$age < 10), TRUE)
expect_equal(all(out$`biotic_v3_example.xml`$agedetermination$serialnumber != c(99484)), TRUE)


## StoxBiotic
st <- StoxBiotic(inputData)

filterExpressionSt <- list()
filterExpressionSt$Individual <- c(
	"IndividualAge >= 10"
)

# Should be OK
outst <-filterData(st, filterExpressionSt)

context("test-Filter: Filtering StoxBiotic data")
# Should be TRUE
expect_equal(any(outst$Individual$IndividualAge < 10), FALSE)
expect_equal(all(outst$Individual$IndividualAge >= 10), TRUE)

# Filter result consistency between level
inputData1 <- ReadBiotic(filenames)
inputData2 <- StoxBiotic(ReadBiotic(filenames))

filterExpression1 <- list()
filterExpression1$`biotic_v3_example.xml`$fishstation <- c(
	'serialnumber %notin% c(99483)'
)

filterExpression2 <- list()
filterExpression2$Haul <- c(
	"HaulKey %notin% c(99483)"
)

out1 <- StoxBiotic(filterData(inputData1, filterExpression1))
out2 <- filterData(inputData2, filterExpression2)

context("test-Filter: Filter downward propagation")
comparison <- all.equal(out1, out2)

# Should be one difference in the Station table
expect_equal(length(comparison), 1)


# Propagate Upwards
inputData <- ReadBiotic(filenames)

filterExpression <- list()
filterExpression$`biotic_v3_example.xml`$catchsample <- c(
	'commonname %notin% c("torsk", "sei", "hyse", "lange", "lysing", "gråsteinbit", "kveite")'
)

context("test-Filter: Filter upward propagation (BioticData)")
outPrup <- filterData(inputData, filterExpression, propagateUpwards = TRUE)
expect_equal(nrow(outPrup$`biotic_v3_example.xml`$fishstation), 1)

context("test-Filter: Filter downward propagation with blank record tables in between")
expect_equal(nrow(outPrup$`biotic_v3_example.xml`$agedetermination), 0)

# Propagate upwards with StoxBiotic
inputData <- StoxBiotic(ReadBiotic(filenames))

filterExpression <- list()
filterExpression$Sample <- c(
	'SpeciesCategoryKey %like% "breiflabb|lyr"'
)

context("test-Filter: Filter downward + upward propagation (StoxBiotic)")
outPrup <- filterData(inputData, filterExpression, propagateUpwards = TRUE)
expect_equal(nrow(outPrup$SpeciesCategory), 2)
expect_equal(nrow(outPrup$Sample), 2)
expect_equal(nrow(outPrup$Haul), 1)
expect_equal(nrow(outPrup$Station), 1)
expect_equal(nrow(outPrup$Individual), 0)





