## ---- eval=FALSE--------------------------------------------------------------
#  options(width=100L)
#  fn <- "dbpedia_csv.tar.gz"
#  
#  if ( !file.exists(fn) ) {
#      download.file("https://github.com/le-scientifique/torchDatasets/raw/master/dbpedia_csv.tar.gz",
#                    fn)
#      untar(fn)
#  }

## ---- eval=FALSE--------------------------------------------------------------
#  library("fastTextR")
#  
#  train <- sample(sprintf("__label__%s", readLines("dbpedia_csv/train.csv")))
#  head(train, 2)

## ---- eval=FALSE--------------------------------------------------------------
#  train <- ft_normalize(train)
#  writeLines(train, con = "dbpedia.train")
#  
#  test <- readLines("dbpedia_csv/test.csv")
#  labels <- trimws(gsub(",.*", "", test))
#  table(labels)

## ---- eval=FALSE--------------------------------------------------------------
#  test <- ft_normalize(test)
#  test <- trimws(sub(".*?,", "", test))
#  head(test, 2)

## ---- eval=FALSE--------------------------------------------------------------
#  cntrl <- ft_control(word_vec_size = 10L, learning_rate = 0.1, max_len_ngram = 2L,
#                      min_count = 1L, nbuckets = 10000000L, epoch = 5L, nthreads = 4L)
#  
#  model <- ft_train(file = "dbpedia.train", method = "supervised", control = cntrl)
#  ft_save(model, "dbpedia.bin")

## ---- eval=FALSE--------------------------------------------------------------
#  model <- ft_load("dbpedia.bin")

## ---- eval=FALSE--------------------------------------------------------------
#  test_pred <- ft_predict(model, newdata=test, k = 1L)
#  str(test_pred)

## ---- eval=FALSE--------------------------------------------------------------
#  confusion_matrix <- table(truth=as.integer(labels),
#                            predicted=as.integer(gsub("\\D", "", test_pred$label)))
#  print(confusion_matrix)

## ---- eval=FALSE--------------------------------------------------------------
#  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
#  print(sprintf("Accuracy: %0.4f", accuracy))

