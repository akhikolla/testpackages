/*
Copyright (c) 2009 Yahoo! Inc.  All rights reserved.  The copyrights
embodied in the content of this file are licensed under the BSD
(revised) open source license
 */
#include <cstring>		// bzero()
#include <stdio.h>
#include <float.h>

#include "cache.h"
#include "io.h"
#include "parse_regressor.h"
#include "parser.h"
#include "parse_args.h"
#include "sender.h"
#include "network.h"
#include "global_data.h"

#include <Rcpp.h>
#define VWCOUT Rcpp::Rcout

using namespace std;
//
// Does string end with a certain substring?
//
bool ends_with(string const &fullString, string const &ending)
{
    if (fullString.length() > ending.length()) {
        return (fullString.compare(fullString.length() - ending.length(), ending.length(), ending) == 0);
    } else {
        return false;
    }
}

size_t next_pow2(size_t x) {
  int i = 0;
  x = x > 0 ? x - 1 : 0;
  while (x > 0) {
    x >>= 1;
    i++;
  }
  return 1 << i;
}

const float default_decay = 1.;

po::variables_map parse_args(int argc, char *argv[], boost::program_options::options_description& desc,
			     gd_vars& vars,
			     regressor &r, parser* par,
			     string& final_regressor_name)
{
  vars.init();
  global.program_name = argv[0];
  // Declare the supported options.
  desc.add_options()
    ("active_learning", "active learning mode")
    ("active_simulation", "active learning simulation mode")
    ("active_mellowness", po::value<float>(&global.active_c0)->default_value(8.f), "active learning mellowness parameter c_0. Default 8")
    ("adaptive", "use adaptive, individual learning rates.")
    ("exact_adaptive_norm", "use a more expensive exact norm for adaptive learning rates.")
    ("audit,a", "print weights of features")
    ("bit_precision,b", po::value<size_t>(),
     "number of bits in the feature table")
    ("backprop", "turn on delayed backprop")
    ("bfgs", "use bfgs optimization")
    ("cache,c", "Use a cache.  The default is <data>.cache")
    ("cache_file", po::value< vector<string> >(), "The location(s) of cache_file.")
    ("compressed", "use gzip format whenever appropriate. If a cache file is being created, this option creates a compressed cache file. A mixture of raw-text & compressed inputs are supported if this option is on")
    ("conjugate_gradient", "use conjugate gradient based optimization")
    ("regularization", po::value<float>(&global.regularization)->default_value(1.0), "l_2 regularization for bfgs")
    ("corrective", "turn on corrective updates")
    ("data,d", po::value< string >()->default_value(""), "Example Set")
    ("daemon", "persistent daemon mode on port 26542")
    ("num_children", po::value<size_t>(&global.num_children)->default_value(10), "number of children for persistent daemon mode")
    ("pid_file", po::value< string >(), "Write pid file in persistent daemon mode")
    ("decay_learning_rate",    po::value<float>(&global.eta_decay_rate)->default_value(default_decay),
     "Set Decay factor for learning_rate between passes")
    ("input_feature_regularizer", po::value< string >(&global.per_feature_regularizer_input), "Per feature regularization input file")
    ("final_regressor,f", po::value< string >(), "Final regressor")
    ("readable_model", po::value< string >(), "Output human-readable final regressor")
    ("global_multiplier", po::value<float>(&global.global_multiplier)->default_value(1.0), "Global update multiplier")
    ("delayed_global", "Do delayed global updates")
    ("hash", po::value< string > (), "how to hash the features. Available options: strings, all")
    ("help,h","Output Arguments")
    ("hessian_on", "use second derivative in line search")
    ("version","Version information")
    ("ignore", po::value< vector<unsigned char> >(), "ignore namespaces beginning with character <arg>")
    ("initial_weight", po::value<float>(&global.initial_weight)->default_value(0.), "Set all weights to an initial value of 1.")
    ("initial_regressor,i", po::value< vector<string> >(), "Initial regressor(s)")
    ("initial_pass_length", po::value<size_t>(&global.pass_length)->default_value((size_t)-1), "initial number of examples per pass")
    ("initial_t", po::value<float>(&(par->t))->default_value(1.), "initial t value")
    ("l1", po::value<float>(&global.l_1_regularization)->default_value(0.), "l_1 regularization level")
    ("lda", po::value<size_t>(&global.lda), "Run lda with <int> topics")
    ("lda_alpha", po::value<float>(&global.lda_alpha)->default_value(0.1), "Prior on sparsity of per-document topic weights")
    ("lda_rho", po::value<float>(&global.lda_rho)->default_value(0.1), "Prior on sparsity of topic distributions")
    ("lda_D", po::value<float>(&global.lda_D)->default_value(10000.), "Number of documents")
    ("minibatch", po::value<size_t>(&global.minibatch)->default_value(1), "Minibatch size, for LDA")
    ("span_server", po::value<string>(&global.span_server)->default_value(""), "Location of server for setting up spanning tree")
    ("min_prediction", po::value<double>(&global.min_label), "Smallest prediction to output")
    ("max_prediction", po::value<double>(&global.max_label), "Largest prediction to output")
    ("mem", po::value<int>(&global.m)->default_value(15), "memory in bfgs")
    ("multisource", po::value<size_t>(), "multiple sources for daemon input")
    ("noconstant", "Don't add a constant feature")
    ("noop","do no learning")
    ("output_feature_regularizer_binary", po::value< string >(&global.per_feature_regularizer_output), "Per feature regularization output file")
    ("output_feature_regularizer_text", po::value< string >(&global.per_feature_regularizer_text), "Per feature regularization output file, in text")
    ("port", po::value<size_t>(),"port to listen on")
    ("power_t", po::value<float>(&vars.power_t)->default_value(0.5), "t power value")
    ("predictto", po::value< string > (), "host to send predictions to")
    ("learning_rate,l", po::value<float>(&global.eta)->default_value(10),
     "Set Learning Rate")
    ("passes", po::value<size_t>(&global.numpasses)->default_value(1),
     "Number of Training Passes")
    ("predictions,p", po::value< string >(), "File to output predictions to")
    ("quadratic,q", po::value< vector<string> > (),
     "Create and use quadratic features")
    ("quiet", "Don't output diagnostics")
    ("rank", po::value<size_t>(&global.rank)->default_value(0), "rank for matrix factorization.")
    ("random_weights", po::value<bool>(&global.random_weights), "make initial weights random")
    ("raw_predictions,r", po::value< string >(),
     "File to output unnormalized predictions to")
    ("save_per_pass", "Save the model after every pass over data")
    ("sendto", po::value< vector<string> >(), "send example to <hosts>")
    ("testonly,t", "Ignore label information and just test")
    ("thread_bits", po::value<size_t>(&global.thread_bits)->default_value(0), "log_2 threads")
    ("loss_function", po::value<string>()->default_value("squared"), "Specify the loss function to be used, uses squared by default. Currently available ones are squared, classic, hinge, logistic and quantile.")
    ("quantile_tau", po::value<double>()->default_value(0.5), "Parameter \\tau associated with Quantile loss. Defaults to 0.5")

    ("unique_id", po::value<size_t>(&global.unique_id)->default_value(0),"unique id used for cluster parallel")
    ("total", po::value<size_t>(&global.total)->default_value(1),"total number of nodes used in cluster parallel")    
    ("node", po::value<size_t>(&global.node)->default_value(0),"node number used for cluster parallel")    

    ("sort_features", "turn this on to disregard order in which features have been defined. This will lead to smaller cache sizes")
    ("ngram", po::value<size_t>(), "Generate N grams")
    ("skips", po::value<size_t>(), "Generate skips in N grams. This in conjunction with the ngram tag can be used to generate generalized n-skip-k-gram.");


  global.queries = 0;
  global.example_number = 0;
  global.weighted_examples = 0.;
  global.old_weighted_examples = 0.;
  global.backprop = false;
  global.bfgs = false;
  global.corrective = false;
  global.delayed_global = false;
  global.bfgs = false;
  global.hessian_on = false;
  global.stride = 1;
  global.weighted_labels = 0.;
  global.total_features = 0;
  global.sum_loss = 0.0;
  global.sum_loss_since_last_dump = 0.0;
  global.dump_interval = exp(1.);
  global.num_bits = 18;
  global.default_bits = true;
  global.daemon = false;
  global.final_prediction_sink.begin = global.final_prediction_sink.end=global.final_prediction_sink.end_array = NULL;
  global.raw_prediction = -1;
  global.local_prediction = -1;
  global.print = print_result;
  global.min_label = 0.;
  global.max_label = 1.;
  global.update_sum = 0.;
  global.lda = 0;
  global.random_weights = false;
  global.per_feature_regularizer_input = "";
  global.per_feature_regularizer_output = "";
  global.per_feature_regularizer_text = "";
  global.ring_size = 1 << 8;

  global.adaptive = false;
  global.add_constant = true;
  global.exact_adaptive_norm = false;
  global.audit = false;
  global.active = false;
  global.active_simulation =false;
  global.reg = &r;

  global.save_per_pass = false;

  po::positional_options_description p;
  // Be friendly: if -d was left out, treat positional param as data file
  p.add("data", -1);

  po::variables_map vm;

  po::store(po::command_line_parser(argc, argv).
	    style(po::command_line_style::default_style ^ po::command_line_style::allow_guessing).
	    options(desc).positional(p).run(), vm);
  po::notify(vm);

  global.weighted_unlabeled_examples = par->t;
  global.initial_t = par->t;
  global.partition_bits = global.thread_bits;

  if (vm.count("help") || argc == 1) {
    /* upon direct query for help -- spit it out to stdout */
    VWCOUT << "\n" << desc << "\n";
    Rf_error("");
  }

  if (vm.count("quiet"))
    global.quiet = true;
  else
    global.quiet = false;

  if (vm.count("active_simulation"))
      global.active_simulation = true;

  if (vm.count("active_learning") && !global.active_simulation)
    global.active = true;

  if (vm.count("adaptive") || vm.count("exact_adaptive_norm")) {
      global.adaptive = true;
      if (vm.count("exact_adaptive_norm"))
	global.exact_adaptive_norm = true;
      global.stride = 2;
      vars.power_t = 0.0;
      if (global.thread_bits != 0)
	{
	  Rf_error("adaptive code isn't correct with multiple learning cores");
	}
  }

  if (vm.count("backprop")) {
      global.backprop = true;
      VWCOUT << "enabling backprop updates" << endl;
  }

  if (vm.count("corrective")) {
      global.corrective = true;
      VWCOUT << "enabling corrective updates" << endl;
  }

  if (vm.count("delayed_global")) {
      global.delayed_global = true;
      VWCOUT << "enabling delayed_global updates" << endl;
  }
  
  if (vm.count("bfgs") || vm.count("conjugate_gradient")) {
    global.bfgs = true;
    global.stride = 4;
    if (vm.count("hessian_on") || global.m==0) {
      global.hessian_on = true;
    }
    if (!global.quiet) {
      if (global.m>0)
	VWCOUT << "enabling BFGS based optimization ";
      else
	VWCOUT << "enabling conjugate gradient optimization via BFGS ";
      if (global.hessian_on)
	VWCOUT << "with curvature calculation" << endl;
      else
	VWCOUT << "**without** curvature calculation" << endl;
    }
    if (global.numpasses < 2)
      {
	Rf_error("you must make at least 2 passes to use BFGS");
      }
  }


  if (vm.count("version") || argc == 1) {
    /* upon direct query for version -- spit it out to stdout */
    VWCOUT << version << "\n";
    Rf_error("");
  }


  if(vm.count("ngram")){
    global.ngram = vm["ngram"].as<size_t>();
    if(!vm.count("skip_gram")) VWCOUT << "You have chosen to generate " << global.ngram << "-grams" << endl;
    if(vm.count("sort_features"))
      {
	VWCOUT << "ngram is incompatible with sort_features.  " << endl;
	Rf_error("");
      }
  }
  if(vm.count("skips"))
    {
    global.skips = vm["skips"].as<size_t>();
    if(!vm.count("ngram"))
      {
	Rf_error("You can not skip unless ngram is > 1");
      }
    VWCOUT << "You have chosen to generate " << global.skips << "-skip-" << global.ngram << "-grams" << endl;
    if(global.skips > 4)
      {
      VWCOUT << "*********************************" << endl;
      VWCOUT << "Generating these features might take quite some time" << endl;
      VWCOUT << "*********************************" << endl;
      }
    }
  if (vm.count("bit_precision"))
    {
      global.default_bits = false;
      global.num_bits = vm["bit_precision"].as< size_t>();
      if (global.num_bits > 29)
	{
	  Rf_error("Only 29 or fewer bits allowed.  If this is a serious limit, speak up.");
	}
    }
  
  if (vm.count("daemon") || vm.count("pid_file")) {
    global.daemon = true;

    // allow each child to process up to 1e5 connections
    global.numpasses = (size_t) 1e5;
  }

  string data_filename = vm["data"].as<string>();
  if (vm.count("compressed") || ends_with(data_filename, ".gz"))
    set_compressed(par);

  if(vm.count("sort_features"))
    par->sort_features = true;

  if (global.num_bits > 30) {
    Rf_error("The system limits at 30 bits of precision!");
  }

  if (vm.count("quadratic"))
    {
      global.pairs = vm["quadratic"].as< vector<string> >();
      if (!global.quiet)
	{
	  VWCOUT << "creating quadratic features for pairs: ";
	  for (vector<string>::iterator i = global.pairs.begin(); i != global.pairs.end();i++) {
	    VWCOUT << *i << " ";
	    if (i->length() > 2)
	      VWCOUT << endl << "warning, ignoring characters after the 2nd.\n";
	    if (i->length() < 2) {
	      Rf_error("error, quadratic features must involve two sets.");
	    }
	  }
	  VWCOUT << endl;
	}
    }

  for (size_t i = 0; i < 256; i++)
    global.ignore[i] = false;
  global.ignore_some = false;

  if (vm.count("ignore"))
    {
      vector<unsigned char> ignore = vm["ignore"].as< vector<unsigned char> >();
      for (vector<unsigned char>::iterator i = ignore.begin(); i != ignore.end();i++)
	{
	  global.ignore[*i] = true;
	  global.ignore_some = true;
	}
      if (!global.quiet)
	{
	  VWCOUT << "ignoring namespaces beginning with: ";
	  for (vector<unsigned char>::iterator i = ignore.begin(); i != ignore.end();i++)
	    VWCOUT << *i << " ";

	  VWCOUT << endl;
	}
    }

  // matrix factorization enabled
  if (global.rank > 0) {
    // store linear + 2*rank weights per index, round up to power of two
    float temp = ceilf(logf((float)(global.rank*2+1)) / logf (2.f));
    global.stride = 1 << (int) temp;
    global.random_weights = true;
  }

  if (vm.count("noconstant"))
    global.add_constant = false;

  if (vm.count("lda"))
    {
      par->sort_features = true;
      float temp = ceilf(logf((float)(global.lda*2+1)) / logf (2.f));
      global.stride = 1 << (int) temp;
      global.random_weights = true;
      global.add_constant = false;
    }

  if (vm.count("lda") && global.eta > 1.)
    {
      VWCOUT << "your learning rate is too high, setting it to 1" << endl;
      global.eta = min(global.eta,1.f);
    }
  if (!vm.count("lda")) 
    global.eta *= pow(par->t, vars.power_t);

  if (vm.count("minibatch")) {
    size_t minibatch2 = next_pow2(global.minibatch);
    global.ring_size = global.ring_size > minibatch2 ? global.ring_size : minibatch2;
  }
  
  parse_regressor_args(vm, r, final_regressor_name, global.quiet);
  parse_source_args(vm,par,global.quiet,global.numpasses);
  if (vm.count("readable_model"))
    global.text_regressor_name = vm["readable_model"].as<string>();
  
  if (vm.count("active_c0"))
    global.active_c0 = vm["active_c0"].as<float>();
  
  if (vm.count("save_per_pass"))
    global.save_per_pass = true;

  if (vm.count("min_prediction"))
    global.min_label = vm["min_prediction"].as<double>();
  if (vm.count("max_prediction"))
    global.max_label = vm["max_prediction"].as<double>();
  if (vm.count("min_prediction") || vm.count("max_prediction") || vm.count("testonly"))
    set_minmax = noop_mm;

  string loss_function;
  if(vm.count("loss_function"))
	  loss_function = vm["loss_function"].as<string>();
  else
	  loss_function = "squaredloss";
  double loss_parameter = 0.0;
  if(vm.count("quantile_tau"))
    loss_parameter = vm["quantile_tau"].as<double>();

  if (global.rank != 0) {
    loss_function = "classic";
    VWCOUT << "Forcing classic squared loss for matrix factorization" << endl;
  }

  r.loss = getLossFunction(loss_function, loss_parameter);
  global.loss = r.loss;

//   global.eta *= pow(par->t, vars.power_t);

  if (global.eta_decay_rate != default_decay && global.numpasses == 1)
    VWCOUT << "Warning: decay_learning_rate has no effect when there is only one pass" << endl;

  if (pow((double)global.eta_decay_rate, (double)global.numpasses) < 0.0001 )
    VWCOUT << "Warning: the learning rate for the last pass is multiplied by: " << pow((double)global.eta_decay_rate, (double)global.numpasses)
	 << " adjust --decay_learning_rate larger to avoid this." << endl;

  //parse_source_args(vm,par,global.quiet,global.numpasses);


  if (!global.quiet)
    {
      VWCOUT << "Num weight bits = " << global.num_bits << endl;
      VWCOUT << "learning rate = " << global.eta << endl;
      VWCOUT << "initial_t = " << par->t << endl;
      VWCOUT << "power_t = " << vars.power_t << endl;
      if (global.numpasses > 1)
	VWCOUT << "decay_learning_rate = " << global.eta_decay_rate << endl;
      if (global.rank > 0)
	VWCOUT << "rank = " << global.rank << endl;
      if (global.regularization > 0 && global.bfgs)
	VWCOUT << "regularization = " << global.regularization << endl;
    }

  if (vm.count("predictions")) {
    if (!global.quiet)
      VWCOUT << "predictions = " <<  vm["predictions"].as< string >() << endl;
    if (strcmp(vm["predictions"].as< string >().c_str(), "stdout") == 0)
      {
	push(global.final_prediction_sink, (size_t) 1);//stdout
      }
    else
      {
	const char* fstr = (vm["predictions"].as< string >().c_str());
	int f = fileno(fopen(fstr,"w"));
	if (f < 0)
	  VWCOUT << "Error opening the predictions file: " << fstr << endl;
	push(global.final_prediction_sink, (size_t) f);
      }
  }

  if (vm.count("raw_predictions")) {
    if (!global.quiet)
      VWCOUT << "raw predictions = " <<  vm["raw_predictions"].as< string >() << endl;
    if (strcmp(vm["raw_predictions"].as< string >().c_str(), "stdout") == 0)
      global.raw_prediction = 1;//stdout
    else
      global.raw_prediction = fileno(fopen(vm["raw_predictions"].as< string >().c_str(), "w"));
  }

  if (vm.count("audit"))
    global.audit = true;

  parse_send_args(vm, global.pairs);

  if (vm.count("testonly"))
    {
      if (!global.quiet)
	VWCOUT << "only testing" << endl;
      global.training = false;
      if (global.lda > 0)
        global.eta = 0;
    }
  else
    {
      global.training = true;
      if (!global.quiet)
	VWCOUT << "learning_rate set to " << global.eta << endl;
    }

  if (vm.count("predictto"))
    {
      if (!global.quiet)
	VWCOUT << "predictto = " << vm["predictto"].as< string >() << endl;
      global.local_prediction = open_socket(vm["predictto"].as< string > ().c_str());
    }

  return vm;
}

