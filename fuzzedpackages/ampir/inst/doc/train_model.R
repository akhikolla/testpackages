## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE
)

## -----------------------------------------------------------------------------
library(ampir)

## -----------------------------------------------------------------------------
bat_pos <- read_faa(system.file("extdata/bat_positive.fasta.gz", package = "ampir"))
bat_pos$Label <- "Positive"
bat_pos <- remove_nonstandard_aa(bat_pos)

## -----------------------------------------------------------------------------
bat_neg <- read_faa(system.file("extdata/bat_negative.fasta.gz", package = "ampir"))
bat_neg$Label <- "Negative"
bat_neg <- remove_nonstandard_aa(bat_neg)
bat_neg <- bat_neg[!bat_neg$seq_aa %in% bat_pos$seq_aa,]
bat_neg <- bat_neg[sample(nrow(bat_neg),78),]

## -----------------------------------------------------------------------------
bats <- rbind(bat_pos, bat_neg)

## -----------------------------------------------------------------------------
bats_features <- calculate_features(bats)
bats_features$Label <- as.factor(bats$Label)
rownames(bats_features) <- NULL

## -----------------------------------------------------------------------------
library(caret)

## -----------------------------------------------------------------------------
trainIndex <-createDataPartition(y=bats_features$Label, p=.7, list = FALSE)
bats_featuresTrain <-bats_features[trainIndex,]
bats_featuresTest <-bats_features[-trainIndex,]

## -----------------------------------------------------------------------------
trctrl_prob <- trainControl(method = "repeatedcv", number = 10, repeats = 3,
                            classProbs = TRUE)

## -----------------------------------------------------------------------------
my_bat_svm_model <- train(Label~.,
                       data = bats_featuresTrain[,-1], # excluding seq_name column
                       method="svmRadial",
                       trControl = trctrl_prob,
                       preProcess = c("center", "scale"))

## -----------------------------------------------------------------------------
my_bat_pred <- predict(my_bat_svm_model, bats_featuresTest)
cm <- confusionMatrix(my_bat_pred, bats_featuresTest$Label, positive = "Positive")

## -----------------------------------------------------------------------------
bat_test_set <- bats[bats$seq_name %in% bats_featuresTest$seq_name,][,-3]
rownames(bat_test_set) <- NULL

## -----------------------------------------------------------------------------
my_bat_AMPs <- predict_amps(bat_test_set, min_len = 5, model = my_bat_svm_model)

## ---- echo=FALSE--------------------------------------------------------------
my_bat_AMPs$seq_aa <- paste(substring(my_bat_AMPs$seq_aa,1,10),"...",sep="")
my_bat_AMPs$seq_name <- paste(substring(my_bat_AMPs$seq_name,4,9),"...",sep="")
my_bat_AMPs$prob_AMP <- round(my_bat_AMPs$prob_AMP, digits = 3)
my_bat_AMPs <- my_bat_AMPs[c(1:3,44:46),]

knitr::kable(my_bat_AMPs)

