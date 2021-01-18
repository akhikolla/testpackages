
// [[Rcpp::plugins("cpp11")]]

#include "interface.h"


namespace fasttext { 

void Args::init_from_list(Rcpp::List control) {
    std::string method = control["method"];
    if ( method == "supervised" ) {
        model = model_name::sup;
        loss = loss_name::softmax;
        minCount = 1;
        minn = 0;
        maxn = 0;
        lr = 0.1;
    } else if ( method == "cbow" ) {
        model = model_name::cbow;
    } else if ( method == "skipgram" ) {
        model = model_name::sg;
    } else {
        Rcpp::stop("unkown method!");
    }
    
    std::string xloss = control["loss"];
    if ( xloss == "softmax" ) {
        loss = loss_name::softmax;
    } else if ( xloss == "hs" ) {
        loss = loss_name::hs;
    } else if ( xloss == "ns" ) {
        loss = loss_name::ns;
    } else {
        Rcpp::stop("unkown loss!");
    }

    input             = Rcpp::as<std::string>( control["input"] );
    // test              = Rcpp::as<std::string>( control["test"] );
    output            = Rcpp::as<std::string>( control["output"] );
    lr                = control["learning_rate"];    
    lrUpdateRate      = control["learn_update"];
    dim               = control["word_vec_size"];
    ws                = control["window_size"];
    epoch             = control["epoch"];
    minCount          = control["min_count"];
    minCountLabel     = control["min_count_label"];
    neg               = control["neg"];
    wordNgrams        = control["max_len_ngram"];        
    bucket            = control["nbuckets"];
    minn              = control["min_ngram"];
    maxn              = control["max_ngram"];
    thread            = control["nthreads"];
    t                 = control["threshold"];
    label             = Rcpp::as<std::string>( control["label"] );
    verbose           = control["verbose"];
    pretrainedVectors = Rcpp::as<std::string>( control["pretrained_vectors"] );

    saveOutput             = control["save_output"]; // bool
    seed                   = control["seed"]; // int
    qnorm                  = control["qnorm"]; // bool
    retrain                = control["retrain"]; // bool
    qout                   = control["qout"]; // bool
    cutoff                 = control["cutoff"]; // size_t
    dsub                   = control["dsub"]; // size_t
    autotuneValidationFile = Rcpp::as<std::string>(control["autotune_validation_file"]); // string
    autotuneMetric         = Rcpp::as<std::string>(control["autotune_metric"]); // string
    autotunePredictions    = control["autotune_predictions"]; // int
    autotuneDuration       = control["autotune_duration"]; // int
    autotuneModelSize      = Rcpp::as<std::string>(control["autotune_model_size"]); // string

}


} // END fasttext NAMESPACE


// [[Rcpp::export]]
std::string Rft_model_type(SEXP ft) {
    Rcpp::XPtr<FastText>fast_text(ft);
    const Args a = fast_text->getArgs();
    if (a.model == fasttext::model_name::cbow) {
        return "cbow";
    } else if (a.model == fasttext::model_name::sg) {
        return "skipgram";
    } else if (a.model == fasttext::model_name::sup) {
        return "supervised";
    }
    return "unkown";
}


// [[Rcpp::export]]
SEXP Rft_new_model() {
    Rcpp::XPtr<FastText> fast_text(new FastText(), true);
    return fast_text;
}


//
// Save/Load Model
// ==========
// [[Rcpp::export]]
SEXP Rft_load_model(std::string file_name) {
    Rcpp::XPtr<FastText> fast_text(new FastText(), true);
    fast_text->loadModel( file_name );
    return fast_text;
}

// [[Rcpp::export]]
SEXP Rft_save_model(SEXP ft, std::string file_name) {
    Rcpp::XPtr<FastText>fast_text(ft);
    fast_text->saveModel( file_name );
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP Rft_save_output(SEXP ft, std::string file_name) {
    Rcpp::XPtr<FastText>fast_text(ft);
    fast_text->saveOutput( file_name );
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP Rft_save_vectors(SEXP ft, std::string file_name) {
    Rcpp::XPtr<FastText>fast_text(ft);
    fast_text->saveVectors( file_name );
    return R_NilValue;
}


//
// Train
// =====
// [[Rcpp::export]]
SEXP Rft_train(SEXP control) {
    Rcpp::XPtr<FastText> fast_text(new FastText(), true);

    Args a = Args();
    a.init_from_list(control);

    fast_text->train(a);
    return fast_text;
}


// [[Rcpp::export]]
int Rft_dict_get_nwords(SEXP ft) {
    Rcpp::XPtr<FastText>fast_text(ft);
    std::shared_ptr<const Dictionary> dict = fast_text->getDictionary();
    int32_t n = dict->nwords();
    return n;
}

// [[Rcpp::export]]
int Rft_dict_get_nlabels(SEXP ft) {
    Rcpp::XPtr<FastText>fast_text(ft);
    std::shared_ptr<const Dictionary> dict = fast_text->getDictionary();
    int32_t n = dict->nlabels();
    return n;
}

// [[Rcpp::export]]
double Rft_dict_get_ntokens(SEXP ft) {
    Rcpp::XPtr<FastText>fast_text(ft);
    std::shared_ptr<const Dictionary> dict = fast_text->getDictionary();
    int64_t n = dict->ntokens();
    return static_cast<double>(n);
}

//
// Predict
// =======
// [[Rcpp::export]]
SEXP Rft_predict_vec(SEXP ft, std::vector<std::string> newdata, int32_t k, float threshold) {
    Rcpp::XPtr<FastText>fast_text(ft);

    std::vector<int32_t> pred_id;
    std::vector<std::string> pred_label;
    std::vector<float> pred_prob;

    for (size_t i=0; i < newdata.size(); i++) {
        std::vector<std::pair<real, std::string>> predictions;
        std::istringstream in_string (newdata[i]);
        bool status = fast_text->predictLine(in_string, predictions, k, threshold);

        if ( predictions.empty() || !status ) continue;

        for (auto it = predictions.cbegin(); it != predictions.cend(); it++) {
            pred_id.push_back(i + 1);
            pred_prob.push_back(it->first);
            pred_label.push_back(it->second);
        }
    }

    return List::create(_["id"] = pred_id, _["label"] = pred_label, _["prob"] = pred_prob);
}


//
// Test
// ====
// [[Rcpp::export]]
SEXP Rft_test(SEXP ft, std::string file_name, int32_t k, float threshold) {
    Rcpp::XPtr<FastText>fast_text(ft);

    std::ifstream infile(file_name);
    fasttext::Meter meter(false);
    fast_text->test(infile, k, threshold, meter);
    
    return List::create(_["nexamples"] =(double)meter.nexamples(), 
                        _["precision"] = meter.precision(), 
                        _["recall"] = meter.recall());
}


// [[Rcpp::export]]
std::vector<std::string> Rft_all_words(SEXP ft) {
    Rcpp::XPtr<FastText>fast_text(ft);
    std::shared_ptr<const Dictionary> dict = fast_text->getDictionary();

    std::vector<std::string> words;
    int32_t nwords = dict->nwords();
    for ( int32_t i=0; i < nwords; i++ ) {
        words.push_back( dict->getWord( i ) );
    }

    return words;
}


// [[Rcpp::export]]
Rcpp::List Rft_word_vectors(SEXP ft, std::vector<std::string> words) {
    Rcpp::XPtr<FastText>fast_text(ft);

    fasttext::Vector vec(fast_text->getDimension());
    Rcpp::List retval(words.size());

    for (int32_t i = 0; i < words.size(); i++) {
        fast_text->getWordVector(vec, words[i]); // vec is set to 0 in getWordVector
        retval[i] = std::vector<float>(vec.data(), vec.data() + vec.size());
    }

    return retval;
}


// [[Rcpp::export]]
Rcpp::NumericVector Rft_nearest_neighbors(SEXP ft, const std::string& word, int32_t k = 10) {
    Rcpp::XPtr<FastText>fast_text(ft);
    Rcpp::NumericVector x(k);
    Rcpp::CharacterVector x_names(k);    

    std::vector<std::pair<real, std::string>> knn = fast_text->getNN(word, k);

    for (int32_t i = 0; i < knn.size(); i++) {
        x[i] = knn[i].first;
        x_names[i] = knn[i].second;
    }
    x.names() = x_names;

    return x;
}


// [[Rcpp::export]]
Rcpp::NumericVector Rft_analogies(SEXP ft, const std::string& wordA, const std::string& wordB, 
                                  const std::string& wordC, int32_t k = 10) {
    Rcpp::XPtr<FastText>fast_text(ft);
    Rcpp::NumericVector x(k);
    Rcpp::CharacterVector x_names(k);
    
    std::vector<std::pair<real, std::string>> analogies;
    analogies = fast_text->getAnalogies(k, wordA, wordB, wordC);

    for (int32_t i = 0; i < analogies.size(); i++) {
        x[i] = analogies[i].first;
        x_names[i] = analogies[i].second;
    }
    x.names() = x_names;

    return x;
}

