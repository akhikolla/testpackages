#' Step corpus with annotation.
#'
#' Given a \code{tm}-package VCorpus of original text, 
#' returns a VCorpus of stemmed text with '+' appended to all stemmed words. 
#'
#' This is non-optimized code that is expensive to run.
#' First the stemmer chops words.  Then this method passes through and adds a "+"
#' to all chopped words, and builds a list of stems.
#' Finally, the method passes through and adds a "+" to all stems found without a
#' suffix.
#' 
#' So, e.g., goblins and goblin will both be transformed to "goblin+".
#' 
#' Adding the '+' makes stemmed text more readible.
#'
#' Code based on code from Kevin Wu, UC Berkeley Undergrad Thesis 2014.
#'
#' Requires, via the tm package, the SnowballC package.
#'
#' Warning: Do not use this on a \code{\link{textreg.corpus}} object.  Do to text before
#' building the \code{\link{textreg.corpus}} object.
#' 
#' @export
#' @param corpus Original text
#' @param verbose True means print out text progress bar so you can watch progress.
#' @import NLP
#' @import tm
#' @importFrom utils setTxtProgressBar txtProgressBar 
#' @examples
#' \donttest{ 
#' library( tm )
#' texts <- c("texting goblins the dagger", "text these goblins", 
#'             "texting 3 goblins appl daggers goblining gobble")
#' corpus <- VCorpus(VectorSource(texts))
#' stemmed_corpus<-stem.corpus(corpus, verbose=FALSE)
#' inspect( stemmed_corpus[[2]] )
#' }
stem.corpus <- function(corpus, verbose = TRUE) {
    if ( is.character( corpus ) ) {
        warning( "Converting character vector corpus to tm VCorpus (VectorSource)" )
        corpus = VCorpus( VectorSource( corpus ) )
    } else if ( is.textreg.corpus(corpus) ) {
        stop( "Unfortunately, textreg.corpus objects can only be passed to textreg() and not other functions" )
    }
    
	stemmed_corpus <- tm_map(corpus, stemDocument)
	
	new_corpus <- rep(NA, length(corpus))
	stemmed_words <- c()
	if (verbose) {
		pb <- txtProgressBar(min = 0, max = length(corpus), style = 3)
	}

	for (i in 1:length(corpus)) {
		if (verbose) {
			setTxtProgressBar(pb, i/2)
		}

		orig <- scan_tokenizer( content(corpus[[i]]) )
		stem <- scan_tokenizer( content(stemmed_corpus[[i]]) )
		if (length(orig) != length(stem)) {
			if (length(stem) == 0) {
				cat( "missing doc is |", content(corpus[[i]]), "|\n", sep="" )
				stem = c("")
			} else {
			    browser()
				cat( "missing doc is |", content(corpus[[i]]), "|\n", sep="" )
				print(orig)
				cat("\n")
				print(stem)
	
				stop("stem.corpus failure: Corpuses do not match up!")
				#If this is the case, check encodings of corpus and inputs
			}
		}
		dels = orig == stem
		new_corpus[i] = paste0(ifelse(dels, orig, paste0(stem, "+")), collapse = " ")
		stemmed_words = unique( c(stemmed_words, stem[!dels]) )
	}

	K = length(corpus)/2
	for (i in 1:length(new_corpus)) {
		if (verbose) {
			setTxtProgressBar(pb, K + i/2)
		}

		new_words <- strsplit(new_corpus[[i]], " ")[[1]]
		orig <- strsplit(content(corpus[[i]]), " ")[[1]]

		ids = new_words %in% stemmed_words
		new_words[ids] = paste0(new_words[ids], "+")
		new_corpus[i] = paste0(new_words, collapse = " ")
	}

	if (verbose) {
		close(pb)
	}
	rm( stemmed_corpus )
	return(VCorpus(VectorSource(new_corpus)))
}



if (FALSE) {

	texts <- c("texting goblins the dagger", "text the goblin", "texting 3 goblins appl daggers goblining gobble")
	texts = rep(texts, 1000)

	corpus <- VCorpus(VectorSource(texts))
	aa = tm_map(corpus, stemDocument)
	content( aa[[55]] )

	stemmed_corpus <- stem.corpus( corpus )
	stemmed_corpus
	content( stemmed_corpus[[55]] )
	content( stemmed_corpus[[4]] )
	content( stemmed_corpus[[5]] )
}
