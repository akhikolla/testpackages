/*
Written by:
Daniel Marcous, Yotam Sandbank
*/

#include "deepboost_C.h"
#include "boost.h"
#include "tree.h"

#include <Rcpp.h>

using namespace Rcpp;


// Train a deepboost model on the given examples, using
// numIter iterations (which not necessarily means numIter trees)
void Train(vector<Example>* train_examples, Model* model, int tree_depth,
 int num_iter, float beta, float lambda, char loss_type, bool verbose) {


	// Train the model
	for (int iter = 1; iter <= num_iter; ++iter) {
		AddTreeToModel(*train_examples, model, loss_type, beta, lambda, tree_depth);
		if (verbose) {
			float error, avg_tree_size;
			int num_trees;
			EvaluateModel(*train_examples, *model, &error, &avg_tree_size,
						  &num_trees);
			Rcpp::Rcout << "Iteration: " << iter
			            << ", error: " << error
			            << ", avg tree size: " << avg_tree_size
			            << ", num trees: " << num_trees
			            << std::endl;
		}
	}
}


// Classify examples using model
vector<Label> Predict(const vector<Example>& examples, const Model& model){
	//TODO::initiate labels
	vector<Label> labels;
    labels.resize(examples.size(), 0);
	for (unsigned i=0; i<examples.size(); i++){
		labels[i] = ClassifyExample(examples[i], model);
    }
	return labels;
}
vector<Probability> PredictProbabilities(const vector<Example>& examples, const Model& model){
  vector<Probability> probabilities;
  probabilities.resize(examples.size(), 0);
  for (unsigned i=0; i<examples.size(); i++){
    probabilities[i] = ComputeExampleClassProbability(examples[i], model);
  }
  return probabilities;
}


// Compute the error of model on examples. Also compute the number of trees in
// model and their average size.
void Evaluate(const vector<Example>& examples, const Model& model,
                   float* error, float* avg_tree_size, int* num_trees){
	EvaluateModel(examples, model, error, avg_tree_size, num_trees);
}
