/*
Copyright 2015 Google Inc. All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include "boost.h"

#include <float.h>
//#include <math.h>
#include <cmath>

#include "tree.h"


float ComputeEta(float wgtd_error, float tree_size, float alpha, float beta, float lambda) {
  wgtd_error = std::fmax(wgtd_error, kTolerance);  // Helps with division by zero.
  const float error_term =
      (1 - wgtd_error) * std::exp(alpha) - wgtd_error * std::exp(-alpha);
  const float complexity_penalty = ComplexityPenalty(tree_size, beta, lambda);
  const float ratio = complexity_penalty / wgtd_error;
  float eta;
  if (std::fabs(error_term) <= 2 * complexity_penalty) {
    eta = -alpha;
  } else if (error_term > 2 * complexity_penalty) {
    eta = std::log(-ratio + std::sqrt(ratio * ratio + (1 - wgtd_error)/wgtd_error));
  } else {
    eta = std::log(ratio + std::sqrt(ratio * ratio + (1 - wgtd_error)/wgtd_error));
  }
  return eta;
}

// TODO(usyed): examples is passed by non-const reference because the example
// weights need to be changed. This is bad style.
void AddTreeToModel(vector<Example>& examples, Model* model, char loss_type, float beta, float lambda, int tree_depth) {
  // Initialize normalizer
  static float normalizer;
  if (model->empty()) {
    if (loss_type == 'e') {
      normalizer = std::exp(1) * static_cast<float>(examples.size());
    } else if (loss_type == 'l') {
      normalizer =
          static_cast<float>(examples.size()) / (std::log(2) * (1 + std::exp(-1)));
    }
  }
  InitializeTreeData(examples, normalizer);
  int best_old_tree_idx = -1;
  float wgtd_error, gradient, best_wgtd_error = 0, best_gradient = 0;

  // Find best old tree
  bool old_tree_is_best = false;
  for (int i = 0; i < model->size(); ++i) {
    const float alpha = (*model)[i].first;
    if (std::fabs(alpha) < kTolerance) continue;  // Skip zeroed-out weights.
    const Tree& old_tree = (*model)[i].second;
    wgtd_error = EvaluateTreeWgtd(examples, old_tree);
    int sign_edge = (wgtd_error >= 0.5) ? 1 : -1;
    gradient = Gradient(wgtd_error, old_tree.size(), alpha, sign_edge, beta, lambda);
    if (std::fabs(gradient) >= std::fabs(best_gradient)) {
      best_gradient = gradient;
      best_wgtd_error = wgtd_error;
      best_old_tree_idx = i;
      old_tree_is_best = true;
    }
  }

  // Find best new tree
  Tree new_tree = TrainTree(examples, beta, lambda, tree_depth);
  wgtd_error = EvaluateTreeWgtd(examples, new_tree);
  gradient = Gradient(wgtd_error, new_tree.size(), 0, -1, beta, lambda);
  if (model->empty() || std::fabs(gradient) > std::fabs(best_gradient)) {
    best_gradient = gradient;
    best_wgtd_error = wgtd_error;
    old_tree_is_best = false;
  }

  // Update model weights
  float alpha = 0;
  const Tree* tree;
  if (old_tree_is_best) {
    alpha = (*model)[best_old_tree_idx].first;
    tree = &((*model)[best_old_tree_idx].second);
  } else {
    alpha = 0;
    tree = &(new_tree);
  }
  const float eta = ComputeEta(best_wgtd_error, tree->size(), alpha, beta, lambda);
  if (old_tree_is_best) {
    (*model)[best_old_tree_idx].first += eta;
  } else {
    model->push_back(make_pair(eta, new_tree));
  }

  // Update examples weights and compute normalizer
  const float old_normalizer = normalizer;
  normalizer = 0;
  for (Example& example : examples) {
    const float u = eta * example.label * ClassifyExample(example, *tree);
    if (loss_type == 'e') {
      example.weight = example.weight * std::exp(-u) * old_normalizer;
    } else if (loss_type == 'l') {
      const float z = (1 - std::log(2) * example.weight * old_normalizer) /
                      (std::log(2) * example.weight * old_normalizer);
      example.weight = 1 / (std::log(2) * (1 + z * std::exp(u)));
    }
    normalizer += example.weight;
  }

  // Renormalize example weights
  // TODO(usyed): Two loops is inefficient.
  for (Example& example : examples) {
    example.weight /= normalizer;
  }
}

Probability ComputeExampleClassProbability(const Example& example, const Model& model) {
  float score = 0;
  float sumOfWeights = 0;
  float probability = 0;
  for (const pair<Weight, Tree>& wgtd_tree : model) {
    score += wgtd_tree.first * ClassifyExample(example, wgtd_tree.second);
    sumOfWeights += wgtd_tree.first;
  }
  probability = ((score/sumOfWeights) + 1) / 2.0;
  return probability;
}

Label ClassifyExample(const Example& example, const Model& model) {
  float score = 0;
  score = ComputeExampleClassProbability(example, model);
  if (score < 0.5) {
    return -1;
  } else {
    return 1;
  }
}

void EvaluateModel(const vector<Example>& examples, const Model& model,
                   float* error, float* avg_tree_size, int* num_trees) {
  float incorrect = 0;
  for (const Example& example : examples) {
    if (example.label != ClassifyExample(example, model)) {
      ++incorrect;
    }
  }
  *num_trees = 0;
  int sum_tree_size = 0;
  for (const pair<Weight, Tree>& wgtd_tree : model) {
    if (std::fabs(wgtd_tree.first) >= kTolerance) {
      ++(*num_trees);
      sum_tree_size += wgtd_tree.second.size();
    }
  }
  *error = (incorrect / examples.size());
  *avg_tree_size = static_cast<float>(sum_tree_size) / *num_trees;
}
