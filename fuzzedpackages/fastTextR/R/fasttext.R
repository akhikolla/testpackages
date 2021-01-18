# -----------------------------------------------------------
#  ft_control
#  ==========
#' @title Default Control Settings
#' @description A auxiliary function for defining the control variables.
#' @param loss a character string giving the name of the loss function 
#'             allowed values are \code{'softmax'}, \code{'hs'} and \code{'ns'}.
#' @param learning_rate a numeric giving the learning rate, the default value is \code{0.05}.
#' @param learn_update an integer giving after how many tokens the learning rate
#'                     should be updated. The default value is \code{100L}, which
#'                     means the learning rate is updated every 100 tokens. 
#' @param word_vec_size an integer giving the length (size) of the word vectors.
#' @param window_size an integer giving the size of the context window.
#' @param epoch an integer giving the number of epochs.
#' @param min_count an integer giving the minimal number of word occurences.
#' @param min_count_label and integer giving the minimal number of label occurences.
#' @param neg an integer giving how many negatives are sampled (only used if loss is \code{"ns"}).
#' @param max_len_ngram an integer giving the maximum length of ngrams used.
#' @param nbuckets an integer giving the number of buckets.
#' @param min_ngram an integer giving the minimal ngram length.
#' @param max_ngram an integer giving the maximal ngram length.
#' @param nthreads an integer giving the number of threads.
#' @param threshold a numeric giving the sampling threshold.
#' @param label a character string specifying the label prefix (default is \code{'__label__'}).
#' @param verbose an integer giving the verbosity level, the default value
#'                is \code{0L} and shouldn't be changed since Rcpp::Rcout 
#'                cann't handle the traffic.
#' @param pretrained_vectors a character string giving the file path
#'                           to the pretrained word vectors which are used 
#'                           for the supervised learning.
#' @param output a character string giving the output file path.
#' @param save_output a logical (default is \code{FALSE})
#' @param seed an integer 
#' @param qnorm  a logical (default is \code{FALSE})
#' @param retrain a logical (default is \code{FALSE})
#' @param qout a logical (default is \code{FALSE})
#' @param cutoff an integer (default is \code{0L})
#' @param dsub an integer (default is \code{2L})
#' @param autotune_validation_file a character string
#' @param autotune_metric a character string (default is \code{"f1"})
#' @param autotune_predictions an integer (default is \code{1L})
#' @param autotune_duration an integer (default is \code{300L})
#' @param autotune_model_size a character string
#' @return a list with the control variables.
#' @examples
#' ft_control(learning_rate=0.1)
ft_control <- function(loss = c("softmax", "hs", "ns"), 
                       learning_rate=0.05, learn_update=100L, word_vec_size=100L, 
                       window_size=5L, epoch=5L, min_count=5L, min_count_label=0L,
                       neg=5L, max_len_ngram=1L, nbuckets=2000000L, min_ngram=3L,
                       max_ngram=6L, nthreads=1L, threshold=1e-4, label="__label__", 
                       verbose=0, pretrained_vectors="", output="", save_output=FALSE, 
                       seed=0L, qnorm=FALSE, retrain=FALSE, qout=FALSE, cutoff=0L, 
                       dsub=2L, autotune_validation_file="", autotune_metric="f1", 
                       autotune_predictions=1L, autotune_duration=300L, 
                       autotune_model_size="") {
    loss <- match.arg(loss)
    cntrl <- as.list(environment())
    class(cntrl) <- c("ft_control", class(cntrl))
    cntrl
}


check_control_arguments <- function(control) {
    control_arguments <- c("input", "output", "learning_rate", "learn_update", 
        "word_vec_size", "window_size", "epoch", "min_count",  "min_count_label", 
        "neg", "max_len_ngram", "nbuckets", "min_ngram", "max_ngram", "nthreads", 
        "threshold", "label", "verbose", "pretrained_vectors", "save_output", "seed",
        "qnorm", "retrain", "qout", "cutoff", "dsub", "autotune_validation_file", 
        "autotune_metric", "autotune_predictions", "autotune_duration", 
        "autotune_model_size")
    if ( !all(control_arguments %in% names(control)) ) {
        i <- which(!control_arguments %in% names(control))
        if ( length(i) == 1) {
            stop("control argument '", control_arguments[i], "' is missing")
        } else {
            stop("control arguments ", 
                 paste(shQuote(control_arguments[i]), collapse=", "), 
                 " are missing")
        }
    }
}


# -----------------------------------------------------------
#  fasttext
#  ========
#' @title Create a New \code{FastText} Model
#' @description Create a new \code{FastText} model. The available methods
#'  are the same as the package functions but with out the prefix \code{"ft_"}
#'  and without the need to provide the model.
#' @examples
#' ft <- fasttext()
fasttext <- function() {
    model <- new.env(parent = emptyenv())

    model$pointer <- NULL
    model$model_type <- "new_model"
    model$nwords <- 0L
    model$ntoken <- 0L
    model$nlabels <- 0L

    update_model <- function(self, new_model) {
        self$pointer <- new_model$pointer
        self$model_type <- new_model$model_type
        self$nwords <- new_model$nwords
        self$ntoken <- new_model$ntoken
        self$nlabels <- new_model$nlabels
        class(self) <- class(new_model)
        return(invisible(NULL))
    }

    model$load <- function(file) {
        self <- parent.env(environment())$model
        update_model(self, ft_load(file))
        return(invisible(NULL))
    }

    model$save <- function(file, what = c("model", "vectors", "output")) {
        ft_save(parent.env(environment())$model, file, what = what)
    }
    
    model$train <- function(file, method = c("supervised", "cbow", "skipgram"), 
                            control = ft_control(), ...) {
        self <- parent.env(environment())$model
        update_model(self, ft_train(file, method, control, ...))
        return(invisible(NULL))
    }

    model$predict <- function(newdata, k = 1L, threshold = 0, 
                              rval = c("sparse", "dense", "slam"), ...) {
        self <- parent.env(environment())$model
        ft_predict(self, newdata, k, threshold, rval, ...)
    }

    model$test <- function(file, k=1L, threshold=0.0) {
        ft_test(parent.env(environment())$model, file, k, threshold)
    }

    model$words <- function() {
        ft_words(parent.env(environment())$model)
    }

    model$word_vectors <- function(words) {
        ft_word_vectors(parent.env(environment())$model, words)
    }

    model$nearest_neighbors <- function(word, k = 10L) {
        ft_nearest_neighbors(parent.env(environment())$model, word, k)
    }

    model$analogies <- function(word_triplets, k = 10L) {
        ft_analogies(parent.env(environment())$model, word_triplets, k)
    }

    class(model) <- "fasttext"
    
    model
}


.wrap_model <- function(env) {
    env$model_type <- Rft_model_type( env$pointer )
    env$nwords <- Rft_dict_get_nwords( env$pointer )
    env$ntoken <- Rft_dict_get_ntokens( env$pointer )
    env$nlabels <- Rft_dict_get_nlabels( env$pointer )
    class(env) <- c(sprintf("%s_model", env$model_type), "fasttext")
    env
}


# -----------------------------------------------------------
#  ft_train
#  ========
#' @title Train a Model
#' @description Train a new word representation model or supervised 
#'              classification model.
#' @param file a character string giving the location of the input file.
#' @param method a character string giving the method, possible values are 
#'               \code{'supervised'}, \code{'cbow'} and \code{'skipgram'}.
#' @param control a list giving the control variables, for more information
#'                see \code{\link{ft_control}}.
#' @param ... additional control arguments inserted into the control list.
#' @examples
#' \dontrun{
#' cntrl <- ft_control(nthreads = 1L)
#' model <- ft_train("my_data.txt", method="supervised", control = cntrl)
#' }
ft_train <- function(file, method = c("supervised", "cbow", "skipgram"), 
                     control = ft_control(), ...) {
    method <- match.arg(method)
    stopifnot( is.character(file), file.exists(file), inherits(control, "ft_control") )
    control <- modifyList(control, list(...))

    control$input <- file
    control$method <- method

    check_control_arguments(control)
    
    model <- fasttext()
    model$pointer <- Rft_train(control)
    .wrap_model(model)
}


# -----------------------------------------------------------
#  ft_predict
#  ==========
#' @title Predict using a Previously Trained Model
#' @description Predict values based on a previously trained model.
#' @param model an object inheriting from \code{'fasttext'}.
#' @param newdata a character vector giving the new data.
#' @param k an integer giving the number of labels to be returned.
#' @param threshold a double withing \code{[0, 1]} giving lower bound
#'                  on the probabilities. Predictions with probabilities
#'                  below this lower bound are not returned. The default
#'                  is \code{0} which means all predictions are returned.
#' @param rval a character string controlling the return value, allowed 
#'      values are \code{"sparse"}, \code{"dense"} and \code{"slam"}.
#'      The default is sparse, here the values are returned as a \code{data.frame}
#'      in a format similar to a simple triplet matrix (sometimes refereed to as the
#'      coordinate format). If \code{rval} is set to \code{"dense"}, a matrix
#'      of the probabilities is returned. Similarly if \code{rval} is set to 
#'      \code{"slam"}, a matrix in the simple triplet sparse format from the 
#'      \pkg{slam} package is returned.
#' @param ... currently not used. 
#' @return \code{NULL} if a \code{'result_file'} is given otherwise 
#'         if \code{'prob'} is true a \code{data.frame} with the predicted labels 
#'         and the corresponding probabilities, if \code{'prob'} is false a 
#'         character vector with the predicted labels.
#' @examples
#' \dontrun{
#' ft_predict(model, newdata)
#' }
#' @name predict.supervised_model
#' @rdname predict.supervised_model
ft_predict <- function(model, newdata, k = 1L, threshold = 0, 
                       rval = c("sparse", "dense", "slam"), ...) {
    stopifnot( inherits(model, "fasttext"), inherits(newdata, "character") )

    rval <- match.arg(rval)
    pred <- Rft_predict_vec(model$pointer, newdata, as.integer(k), threshold = threshold)

    if (rval == "sparse") return(data.frame(pred, stringsAsFactors = FALSE))

    labels <- unique(pred$label)
    labels <- labels[order(as.integer(gsub("\\D", "", labels)))]
    j <- match(pred$label, labels)
    probs <- simple_triplet_matrix(i = pred$id, j = j, v = pred$prob) 
    colnames(probs) <- labels
    if (rval == "dense") as.matrix(probs) else probs
}


# -----------------------------------------------------------
#  ft_test
#  =======
#' @title Evaluate the Model
#' @description Evaluate the quality of the predictions.
#'              For the model evaluation precision and recall are used.
#' @param model an object inheriting from \code{'fasttext'}.
#' @param file a character string giving the location of the validation file.
#' @param k an integer giving the number of labels to be returned.
#' @param threshold a double giving the threshold.
#' @examples
#' \dontrun{
#' ft_test(model, file)
#' }
ft_test <- function(model, file, k=1L, threshold=0.0) {
    stopifnot(is.character(file), file.exists(file))
    stopifnot( inherits(model, "fasttext") )
    x <- Rft_test(model$pointer, file, as.integer(k), threshold)
    x
}


# -----------------------------------------------------------
#  ft_load_model
#  =============
#' @title Load Model
#' @description Load a previously saved model from file.
#' @param file a character string giving the name of the file
#'             to be read in.
#' @return an object inheriting from \code{"fasttext"}.
#' @examples
#' \dontrun{
#' model <- ft_load("dbpedia.bin")
#' }
ft_load <- function(file) {
    if ( !file.exists(file) ) {
        stop("cannot open file '", file, "': No such file or directory")
    }
    model <- fasttext()
    model$pointer <- try(Rft_load_model(file), silent = TRUE)
    if ( inherits(model, "try-error") ) {
        stop("version doesn't fit. The model was created by a different version than 'fastTextR' uses.")
    }
    .wrap_model(model)
}


# -----------------------------------------------------------
#  ft_save
#  =======
#' @title Write Model
#' @description Write a previously saved model from file.
#' @param model an object inheriting from \code{'fasttext'}.
#' @param file a character string giving the name of the file.
#' @param what a character string giving what should be saved.
#' @examples
#' \dontrun{
#' ft_save(model, "my_model.bin", what = "model")
#' }
ft_save <- function(model, file, what = c("model", "vectors", "output")) {
    what <- match.arg(what)
    if (what == "model") {
        Rft_save_model(model$pointer, file)
    } else if (what == "vectors") {
        Rft_save_vectors(model$pointer, file)
    } else {
        Rft_save_output(model$pointer, file)
    }
    invisible(NULL)
}


# -----------------------------------------------------------
#  ft_words
#  ========
#' @title Get Words
#' @description Obtain all the words from a previously trained model.
#' @param model an object inheriting from \code{"fasttext"}.
#' @return a character vector.
#' @examples
#' \dontrun{
#' ft_words(model)
#' }
ft_words <- function(model) {
    stopifnot( inherits(model, "fasttext") )
    Rft_all_words(model$pointer)
}


# -----------------------------------------------------------
#  get_word_vectors
#  ================
#' @title Get Word Vectors
#' @description Obtain word vectors from a previously trained model.
#' @param model an object inheriting from \code{"fasttext"}.
#' @param words a character vector giving the words.
#' @return a matrix containing the word vectors.
#' @examples
#' \dontrun{
#' ft_word_vectors(model, c("word", "vector"))
#' }
ft_word_vectors <- function(model, words) {
    stopifnot(is.character(words), all(nchar(words) > 0))
    stopifnot( inherits(model, "fasttext") )
    word_vec <- Rft_word_vectors(model$pointer, words)
    word_vec <- do.call(rbind, word_vec)
    rownames(word_vec) <- words
    word_vec
}


# -----------------------------------------------------------
#  ft_nearest_neighbors
#  =================
#' @title Get Nearest Neighbors
#' @description TODO
#' @param model an object inheriting from \code{"fasttext"}.
#' @param word a character string giving the word.
#' @param k an integer giving the number of nearest neighbors to be returned.
#' @return .
#' @examples
#' \dontrun{
#' ft_nearest_neighbors(model, "enviroment", k = 6L)
#' }
ft_nearest_neighbors <- function(model, word, k = 10L) {
    stopifnot(is.character(word), isTRUE(length(word) == 1L))
    stopifnot( inherits(model, "fasttext") )
    Rft_nearest_neighbors(model$pointer, word, as.integer(k))
}

# -----------------------------------------------------------
#  ft_analogies
#  ============
#' @title Get Analogies
#' @description TODO
#' @param model an object inheriting from \code{"fasttext"}.
#' @param word_triplets a character vector of length \code{}string giving the word.
#' @param k an integer giving the number of nearest neighbors to be returned.
#' @return .
#' @examples
#' \dontrun{
#' ft_analogies(model, c("berlin", "germany", "france"), k = 6L)
#' }
ft_analogies <- function(model, word_triplets, k = 10L) {
    stopifnot(is.character(word_triplets), isTRUE(length(word_triplets) == 3L))
    stopifnot( inherits(model, "fasttext") )
    Rft_analogies(model$pointer, word_triplets[1L], word_triplets[2L], 
                  word_triplets[3L], as.integer(k))
}


print_fasttext_methods <- function(x) {
    hide <- c("nlabels",  "ntoken", "nwords", "pointer", "update")
    ft_methods <- setdiff(ls(x), hide)
    for (method in ft_methods) {
        args <- gsub("function ", "", trimws(capture.output(str(get(method, envir = x)))))
        cat("  $", method, args, "\n", sep="")
    }
}


print.fasttext <- function(x, ...) {
    cat("Empty fastText model:\n")
    print_fasttext_methods(x)
}

print.cbow_model <- function(x, ...) {
    cat("fastText", shQuote(x$model_type), "model:", "\n")
    cat(sprintf("    %s tokens, %s words", 
        format(x$ntoken, big.mark=","),
        format(x$nwords, big.mark=",")), "\n")
}

print.skipgram_model <- function(x, ...) {
    cat("fastText", shQuote(x$model_type), "model:", "\n")
    cat(sprintf("    %s tokens, %s words", 
        format(x$ntoken, big.mark=","),
        format(x$nwords, big.mark=",")), "\n")
}

print.supervised_model <- function(x, ...) {
    cat("fastText", shQuote(x$model_type), "model:", "\n")
    cat(sprintf("    %s tokens, %s words, %s labels", 
        format(x$ntoken, big.mark=","),
        format(x$nwords, big.mark=","),
        format(x$nlabels, big.mark=",")), "\n")
}

# -----------------------------------------------------------
#  ft_normalize
#  ============
#' @title Normalize
#' @description Applies normalization to a given text.
#' @param txt a character vector to be normalized.
#' @return a character vector.
#' @examples
#' \dontrun{
#' ft_normalize(some_text)
#' }
ft_normalize <- function(txt) {
    clean_text(txt)
}
