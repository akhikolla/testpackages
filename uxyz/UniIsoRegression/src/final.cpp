#include <Rcpp.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <stack>
#include <cmath>
#include <algorithm>

//using namespace std;
using namespace Rcpp;

void swapDouble(double& a, double& b) {
  double temp = a;
  a = b;
  b = temp;
}
void swapDoubleV(std::vector<double>*& a, std::vector<double>*& b) {
  std::vector<double>* temp = a;
  a = b;
  b = temp;
}

//Fix up for min heap, starting from index k
void fixUpMin(std::vector<double>& heapValue, std::vector<double>& heapWeight, int k) {
  while (k && heapValue[k / 2] > heapValue[k]) {
    swapDouble(heapValue[k / 2], heapValue[k]);
    swapDouble(heapWeight[k / 2], heapWeight[k]);
    k /= 2;
  }
}

//Fix up for max heap, starting from index k
void fixUpMax(std::vector<double>& heapValue, std::vector<double>& heapWeight, int k) {
  while (k && heapValue[k / 2] < heapValue[k]) {
    swapDouble(heapValue[k / 2], heapValue[k]);
    swapDouble(heapWeight[k / 2], heapWeight[k]);
    k /= 2;
  }
}

//Fix down for min heap, starting from index k
void fixDownMin(std::vector<double>& heapValue, std::vector<double>& heapWeight, int k) {
  int heapSize = heapValue.size();
  while (2 * k + 1 < heapSize && (heapValue[2 * k + 1] < heapValue[k] || (2 * k + 2 < heapSize && heapValue[2 * k + 2] < heapValue[k]))) {
    if (2 * k + 2 < heapSize) {
      if (heapValue[2 * k + 1] < heapValue[2 * k + 2]) {
        swapDouble(heapValue[2 * k + 1], heapValue[k]);
        swapDouble(heapWeight[2 * k + 1], heapWeight[k]);
        k = 2 * k + 1;
      }
      else {
        swapDouble(heapValue[2 * k + 2], heapValue[k]);
        swapDouble(heapWeight[2 * k + 2], heapWeight[k]);
        k = 2 * k + 2;
      }
    }
    else {
      swapDouble(heapValue[2 * k + 1], heapValue[k]);
      swapDouble(heapWeight[2 * k + 1], heapWeight[k]);
      k = 2 * k + 1;
    }
  }
}

//Fix down for the max heap, starting from index k
void fixDownMax(std::vector<double>& heapValue, std::vector<double>& heapWeight, int k) {
  int heapSize = heapValue.size();
  while (2 * k + 1 < heapSize && (heapValue[2 * k + 1] > heapValue[k] || (2 * k + 2 < heapSize && heapValue[2 * k + 2] > heapValue[k]))) {
    if (2 * k + 2 < heapSize) {
      if (heapValue[2 * k + 1] > heapValue[2 * k + 2]) {
        swapDouble(heapValue[2 * k + 1], heapValue[k]);
        swapDouble(heapWeight[2 * k + 1], heapWeight[k]);
        k = 2 * k + 1;
      }
      else {
        swapDouble(heapValue[2 * k + 2], heapValue[k]);
        swapDouble(heapWeight[2 * k + 2], heapWeight[k]);
        k = 2 * k + 2;
      }
    }
    else {
      swapDouble(heapValue[2 * k + 1], heapValue[k]);
      swapDouble(heapWeight[2 * k + 1], heapWeight[k]);
      k = 2 * k + 1;
    }
  }
}

void addMin(std::vector<double>& heapValue, std::vector<double>& heapWeight, double toAddVal, double toAddWeight) {
  heapValue.push_back(toAddVal);
  heapWeight.push_back(toAddWeight);
  fixUpMin(heapValue, heapWeight, heapValue.size() - 1);
}

void addMax(std::vector<double>& heapValue, std::vector<double>& heapWeight, double toAddVal, double toAddWeight) {
  heapValue.push_back(toAddVal);
  heapWeight.push_back(toAddWeight);
  fixUpMax(heapValue, heapWeight, heapValue.size() - 1);
}

void popMin(std::vector<double>& heapValue, std::vector<double>& heapWeight) {
  heapValue[0] = heapValue.back();
  heapWeight[0] = heapWeight.back();
  heapValue.pop_back();
  heapWeight.pop_back();
  fixDownMin(heapValue, heapWeight, 0);
}

void popMax(std::vector<double>& heapValue, std::vector<double>& heapWeight) {
  heapValue[0] = heapValue.back();
  heapWeight[0] = heapWeight.back();
  heapValue.pop_back();
  heapWeight.pop_back();
  fixDownMax(heapValue, heapWeight, 0);
}

void addToLevelSet(std::vector<double>& ltHeapValue, std::vector<double>& ltHeapWeight,
                   std::vector<double>& gtHeapValue, std::vector<double>& gtHeapWeight, double& medVal,
                   double& medWeight, double& ltWeight, double& gtWeight, double toAddVal, double toAddWeight) {
  if (toAddVal < medVal) {
    ltWeight += toAddWeight;
    addMax(ltHeapValue, ltHeapWeight, toAddVal, toAddWeight);
  }
  else {
    gtWeight += toAddWeight;
    addMin(gtHeapValue, gtHeapWeight, toAddVal, toAddWeight);
  }
  double totalWeight = medWeight + ltWeight + gtWeight;
  while (ltWeight > totalWeight / 2 || gtWeight > totalWeight / 2) {
    if (ltWeight > gtWeight) {
      addMin(gtHeapValue, gtHeapWeight, medVal, medWeight);
      gtWeight += medWeight;
      medVal = ltHeapValue[0];
      medWeight = ltHeapWeight[0];
      ltWeight -= ltHeapWeight[0];
      popMax(ltHeapValue, ltHeapWeight);
    }
    else {
      addMax(ltHeapValue, ltHeapWeight, medVal, medWeight);
      ltWeight += medWeight;
      medVal = gtHeapValue[0];
      medWeight = gtHeapWeight[0];
      gtWeight -= gtHeapWeight[0];
      popMin(gtHeapValue, gtHeapWeight);
    }
  }
}

void merge(std::vector<double>*& ltHeapValue_a, std::vector<double>*& ltHeapWeight_a, std::vector<double>*& gtHeapValue_a, std::vector<double>*& gtHeapWeight_a,
           double& medVal_a, double& medWeight_a, double& ltWeight_a, double& gtWeight_a, std::vector<double>*& ltHeapValue_b, std::vector<double>*& ltHeapWeight_b,
           std::vector<double>*& gtHeapValue_b, std::vector<double>*& gtHeapWeight_b, double& medVal_b, double& medWeight_b, double& ltWeight_b, double& gtWeight_b) {

  int a_size = ltHeapValue_a->size() + gtHeapValue_a->size() + 1;
  int b_size = ltHeapValue_b->size() + gtHeapValue_b->size() + 1;

  if (b_size > a_size) {
    swapDoubleV(ltHeapValue_a, ltHeapValue_b);
    swapDoubleV(ltHeapWeight_a, ltHeapWeight_b);
    swapDoubleV(gtHeapValue_a, gtHeapValue_b);
    swapDoubleV(gtHeapWeight_a, gtHeapWeight_b);
    swapDouble(medVal_a, medVal_b);
    swapDouble(medWeight_a, medWeight_b);
    swapDouble(ltWeight_a, ltWeight_b);
    swapDouble(gtWeight_a, gtWeight_b);
  }

  for (int i = 0; i < int(ltHeapValue_b->size()); ++i) {
    addToLevelSet(*ltHeapValue_a, *ltHeapWeight_a, *gtHeapValue_a, *gtHeapWeight_a, medVal_a,
                  medWeight_a, ltWeight_a, gtWeight_a, (*ltHeapValue_b)[i], (*ltHeapWeight_b)[i]);
  }
  for (int i = 0; i < int(gtHeapValue_b->size()); ++i) {
    addToLevelSet(*ltHeapValue_a, *ltHeapWeight_a, *gtHeapValue_a, *gtHeapWeight_a, medVal_a,
                  medWeight_a, ltWeight_a, gtWeight_a, (*gtHeapValue_b)[i], (*gtHeapWeight_b)[i]);
  }
  addToLevelSet(*ltHeapValue_a, *ltHeapWeight_a, *gtHeapValue_a, *gtHeapWeight_a, medVal_a,
                medWeight_a, ltWeight_a, gtWeight_a, medVal_b, medWeight_b);
}

// [[Rcpp::export(name = ".reg_1d_l1")]]
std::vector<double> reg_1d_l1(std::vector<double>& ycoords,
                              std::vector<double>& weights, std::vector<double>& error, bool decreasing = false) {
  if (decreasing) {
    std::reverse(ycoords.begin(), ycoords.end());
    std::reverse(weights.begin(), weights.end());
  }

  size_t data_size = ycoords.size();

  std::vector< std::vector<double>* > ltHeapValues;
  std::vector< std::vector<double>* > ltHeapWeights;
  ltHeapValues.reserve(data_size);
  ltHeapWeights.reserve(data_size);

  std::vector< std::vector<double>* > gtHeapValues;
  std::vector< std::vector<double>* > gtHeapWeights;
  gtHeapValues.reserve(data_size);
  gtHeapWeights.reserve(data_size);

  std::vector<double> medValues;
  std::vector<double> medWeights;
  medValues.reserve(data_size);
  medWeights.reserve(data_size);

  std::vector<double> ltWeights;
  std::vector<double> gtWeights;
  ltWeights.reserve(data_size);
  gtWeights.reserve(data_size);

  std::vector<int> lefts;
  lefts.reserve(data_size);

  for (int i = 0; i < int(data_size); ++i) {
    //initialize level set
    ltHeapValues.push_back(0);
    ltHeapWeights.push_back(0);
    ltHeapValues[ltHeapValues.size() - 1] = new std::vector<double>;
    ltHeapWeights[ltHeapWeights.size() - 1] = new std::vector<double>;
    gtHeapValues.push_back(0);
    gtHeapWeights.push_back(0);
    gtHeapValues[gtHeapValues.size() - 1] = new std::vector<double>;
    gtHeapWeights[gtHeapWeights.size() - 1] = new std::vector<double>;
    medValues.push_back(ycoords[i]);
    medWeights.push_back(weights[i]);
    ltWeights.push_back(0);
    gtWeights.push_back(0);
    lefts.push_back(i);

    int back = medValues.size() - 1;
    int stBack = medValues.size() - 2; //second to back
    while (back > 0 &&
           medValues[back] <= medValues[stBack]) {

      //merge the level sets


      merge(ltHeapValues[stBack], ltHeapWeights[stBack], gtHeapValues[stBack],
            gtHeapWeights[stBack], medValues[stBack], medWeights[stBack],
                                                                ltWeights[stBack], gtWeights[stBack], ltHeapValues[back],
                                                                                                                  ltHeapWeights[back], gtHeapValues[back], gtHeapWeights[back],
                                                                                                                                                                        medValues[back], medWeights[back], ltWeights[back], gtWeights[back]);

      delete ltHeapValues[back];
      delete gtHeapValues[back];
      delete ltHeapWeights[back];
      delete gtHeapWeights[back];
      ltHeapValues.pop_back();
      ltHeapWeights.pop_back();
      gtHeapValues.pop_back();
      gtHeapWeights.pop_back();
      medValues.pop_back();
      medWeights.pop_back();
      ltWeights.pop_back();
      gtWeights.pop_back();
      lefts.pop_back();

      --back;
      --stBack;
    }

  }
  for (size_t i = 0; i < ltHeapValues.size(); ++i) {
    delete ltHeapValues[i];
    delete gtHeapValues[i];
    delete ltHeapWeights[i];
    delete gtHeapWeights[i];
  }


  std::vector<double> reg_data(data_size);

  size_t curLSet = 0;
  reg_data[0] = medValues[0];
  error[0] = fabs(ycoords[0] - reg_data[0]);
  for (int i = 1; i < int(data_size); ++i) {
    if (int(curLSet) != int(medValues.size()) - 1 && lefts[curLSet + 1] == i) {
      ++curLSet;
    }
    reg_data[i] = medValues[curLSet];
    error[i] = error[i - 1] + fabs(ycoords[i] - reg_data[i]);
  }

  if (decreasing)
    std::reverse(reg_data.begin(), reg_data.end());

  return reg_data;
}

// [[Rcpp::export(name = ".uni_1d_l1")]]
std::vector<double> uni_1d_l1(std::vector<double>& ycoords, std::vector<double>& weights) {
  size_t size = ycoords.size();

  std::vector<double> incErrors(size, 0);
  std::vector<double> inc_reg = reg_1d_l1(ycoords, weights, incErrors);

  //reverse std::vectors, then call prefix L1 on it again

  std::reverse(ycoords.begin(), ycoords.end());
  std::reverse(weights.begin(), weights.end());

  std::vector<double> decErrors(size, 0);
  std::vector<double> dec_reg = reg_1d_l1(ycoords, weights, decErrors);


  //now find the best point for the mode, and copy the data used in
  //the optimal unimodal regression into a new std::vector

  //initialized to be the min of fully increasing or fully decreasing
  double min_error = std::min(decErrors.back(), incErrors.back());
  //if all are decreasing, set last increasing to be size to denote none
  size_t last_increasing = min_error == incErrors.back() ? size - 1 : size;

  //the prefix regression did the decreasing regression in reverse order,
  //so instead of reversing the resulting std::vector we just go from the back
  for (size_t i = 0; i < size - 2; ++i) {
    double error_at_i = incErrors[i] + decErrors[size - 2 - i];
    if (min_error > error_at_i) {
      min_error = error_at_i;
      last_increasing = i;
    }
  }

  if (last_increasing == size) {
    return dec_reg;
  }
  else if (last_increasing == size - 1) {
    return inc_reg;
  }

  std::vector<double> unimodal(size);

  for (size_t i = 0; i <= last_increasing; ++i) {
    unimodal[i] = inc_reg[i];
  }
  for (size_t i = last_increasing + 1; i < size; ++i) { //double check bounds here***
    unimodal[i] = dec_reg[size - 1 - i];
  }

  return unimodal;

}

class region {
public:
  std::vector <std::pair<int, std::pair<int, int> > > range;
  int m, l, r;
};

// [[Rcpp::export(name = ".reg_1d_l2")]]
std::vector<double> reg_1d_l2(std::vector<double> &y_vec, std::vector<double> &w_vec, bool decreasing = false) {
  if (decreasing) {
    std::reverse(y_vec.begin(), y_vec.end());
    std::reverse(w_vec.begin(), w_vec.end());
  }
  size_t size = y_vec.size();
  std::vector<double> y_new(size, 0);
  std::vector<double> mean_vec(y_vec);
  std::vector<double> sumwy_vec(y_vec);
  std::vector<double> levelerr_vec(size, 0);
  std::vector<double> sumwy2_vec(size, 0);
  std::vector<double> sumw_vec(size, 1);
  std::vector<size_t> left_vec(size, 0);
  std::vector<size_t> right_vec(size, 0);

  for (size_t i = 0; i < size; ++i) {
    sumwy2_vec[i] = y_vec[i] * y_vec[i] * w_vec[i] * w_vec[i];
    sumwy_vec[i] = y_vec[i] * w_vec[i];
    sumw_vec[i] = w_vec[i];
    left_vec[i] = i;
    right_vec[i] = i;
  }

  size_t vec_back = 0;
  for (size_t j = 0; j < size; ++j) {
    bool flag = 0;
    while (mean_vec[j] <= mean_vec[vec_back]) {
      left_vec[j] = left_vec[vec_back];
      sumwy_vec[j] = sumwy_vec[j] + sumwy_vec[vec_back];
      sumwy2_vec[j] = sumwy2_vec[j] + sumwy2_vec[vec_back];
      sumw_vec[j] = sumw_vec[j] + sumw_vec[vec_back];
      mean_vec[j] = sumwy_vec[j] / sumw_vec[j];
      levelerr_vec[j] = sumwy2_vec[j] - sumwy_vec[j] * sumwy_vec[j] / sumw_vec[j];
      if (!vec_back) {
        flag = 1;
        break;
      }
      vec_back = vec_back - 1;
    }
    if (!flag)  vec_back = vec_back + 1;
    left_vec[vec_back] = left_vec[j];
    right_vec[vec_back] = right_vec[j];
    sumwy_vec[vec_back] = sumwy_vec[j];
    sumwy2_vec[vec_back] = sumwy2_vec[j];
    sumw_vec[vec_back] = sumw_vec[j];
    mean_vec[vec_back] = mean_vec[j];
    levelerr_vec[vec_back] = levelerr_vec[j];
  }
  if (levelerr_vec[0] > 0) {}
  for (size_t k = 0; k <= vec_back; ++k) {
    for (size_t l = left_vec[k]; l <= right_vec[k]; ++l) {
      y_new[l] = mean_vec[k];
    }
  }
  if (decreasing) {
    std::reverse(y_new.begin(), y_new.end());
  }
  return y_new;
}




// [[Rcpp::export(name = ".uni_1d_l2")]]
std::vector<double> uni_1d_l2(std::vector<double>& y_vec, std::vector<double>& w_vec) {
  int size = y_vec.size();
  std::vector<double> error((size + 1), 0);
  std::vector<double> err_asc(size, 0);
  std::vector<double> err_dec(size, 0);
  std::vector<double> y_new(size, 0);

  //ascending part
  std::vector<double> mean_vec(y_vec);
  std::vector<double> sumwy_vec(y_vec);
  std::vector<double> levelerr_vec(size, 0);
  std::vector<double> sumwy2_vec(size, 0);
  std::vector<double> sumw_vec(size, 0);

  for (int i = 0; i < size; ++i) {


    sumwy2_vec[i] = y_vec[i] * y_vec[i] * w_vec[i] * w_vec[i];
    sumwy_vec[i] = y_vec[i] * w_vec[i];
    sumw_vec[i] = w_vec[i];
  }
  int vec_back = 0;
  err_asc[0] = 0;
  for (int j = 1; j < size; ++j) {
    err_asc[j] = err_asc[j - 1];
    while (vec_back >= 0 && mean_vec[j] <= mean_vec[vec_back]) {
      sumwy_vec[j] = sumwy_vec[j] + sumwy_vec[vec_back];
      sumwy2_vec[j] = sumwy2_vec[j] + sumwy2_vec[vec_back];
      sumw_vec[j] = sumw_vec[j] + sumw_vec[vec_back];
      mean_vec[j] = sumwy_vec[j] / sumw_vec[j];
      levelerr_vec[j] = sumwy2_vec[j] - sumwy_vec[j] * sumwy_vec[j] / sumw_vec[j];
      err_asc[j] = err_asc[j] - levelerr_vec[vec_back];
      vec_back = vec_back - 1;
    }
    vec_back = vec_back + 1;
    sumwy_vec[vec_back] = sumwy_vec[j];
    sumwy2_vec[vec_back] = sumwy2_vec[j];
    sumw_vec[vec_back] = sumw_vec[j];
    mean_vec[vec_back] = mean_vec[j];
    levelerr_vec[vec_back] = levelerr_vec[j];
    err_asc[j] = err_asc[j] + levelerr_vec[vec_back];
  }

  //decending part
  for (int i = 0; i < size; ++i) {
    mean_vec[i] = y_vec[size - i - 1];
    levelerr_vec[i] = 0;
    sumwy2_vec[i] = y_vec[size - i - 1] * y_vec[size - i - 1] * w_vec[size - i - 1] * w_vec[size - i - 1];
    sumwy_vec[i] = y_vec[size - i - 1] * w_vec[size - i - 1];
    sumw_vec[i] = w_vec[size - i - 1];
  }

  vec_back = 0;
  err_dec[0] = 0;
  for (int j = 1; j < size; ++j) {
    err_dec[j] = err_dec[j - 1];
    while (vec_back >= 0 && mean_vec[j] <= mean_vec[vec_back]) {
      sumwy_vec[j] = sumwy_vec[j] + sumwy_vec[vec_back];
      sumwy2_vec[j] = sumwy2_vec[j] + sumwy2_vec[vec_back];
      sumw_vec[j] = sumw_vec[j] + sumw_vec[vec_back];
      mean_vec[j] = sumwy_vec[j] / sumw_vec[j];
      levelerr_vec[j] = sumwy2_vec[j] - sumwy_vec[j] * sumwy_vec[j] / sumw_vec[j];
      err_dec[j] = err_dec[j] - levelerr_vec[vec_back];
      vec_back = vec_back - 1;
    }
    vec_back = vec_back + 1;
    sumwy_vec[vec_back] = sumwy_vec[j];
    sumwy2_vec[vec_back] = sumwy2_vec[j];
    sumw_vec[vec_back] = sumw_vec[j];
    mean_vec[vec_back] = mean_vec[j];
    levelerr_vec[vec_back] = levelerr_vec[j];
    err_dec[j] = err_dec[j] + levelerr_vec[vec_back];
  }
  //find minimum error
  error[0] = err_dec[size - 1];
  error[size] = err_asc[size - 1];
  for (int i = 1; i < size; ++i) {
    error[i] = err_asc[i - 1] + err_dec[size - i];
  }

  int pos_min = 0;
  for (int i = 0; i < size; ++i) {
    if (error[i] < error[pos_min]) {
      pos_min = i;
    }
  }


  std::vector<double> left_vec(y_vec.size(), 0);
  std::vector<double> right_vec(y_vec.size(), 0);
  for (int i = 0; i < int(y_vec.size()); ++i) {
    left_vec[i] = i;
    right_vec[i] = i;
  }


  //regression on ascending part
  if (pos_min != 0) {

    for (int i = 0; i < pos_min; ++i) {
      mean_vec[i] = y_vec[i];
      sumwy_vec[i] = y_vec[i] * w_vec[i];
      sumw_vec[i] = w_vec[i];
    }

    vec_back = 0;
    for (int j = 1; j < pos_min; ++j) {
      while (vec_back >= 0 && mean_vec[j] <= mean_vec[vec_back]) {
        left_vec[j] = left_vec[vec_back];
        sumwy_vec[j] = sumwy_vec[j] + sumwy_vec[vec_back];
        sumw_vec[j] = sumw_vec[j] + sumw_vec[vec_back];
        mean_vec[j] = sumwy_vec[j] / sumw_vec[j];
        vec_back = vec_back - 1;
      }
      vec_back = vec_back + 1;
      left_vec[vec_back] = left_vec[j];
      right_vec[vec_back] = right_vec[j];
      sumwy_vec[vec_back] = sumwy_vec[j];
      sumw_vec[vec_back] = sumw_vec[j];
      mean_vec[vec_back] = mean_vec[j];
    }

    for (int k = 0; k <= vec_back; ++k) {
      for (int l = (int)left_vec[k]; l <= (int)right_vec[k]; ++l) {
        y_new[l] = mean_vec[k];
      }
    }

  }

  //regression on decending part
  if (pos_min != size) {

    for (int i = 0; i < (size - pos_min); ++i) {
      mean_vec[i] = y_vec[size - i - 1];
      left_vec[i] = i;
      right_vec[i] = i;
      sumwy_vec[i] = y_vec[size - i - 1] * w_vec[size - i - 1];
      sumw_vec[i] = w_vec[size - i - 1];
    }

    vec_back = 0;
    for (int j = 1; j < (size - pos_min); ++j) {
      while (vec_back >= 0 && mean_vec[j] <= mean_vec[vec_back]) {
        left_vec[j] = left_vec[vec_back];
        sumwy_vec[j] = sumwy_vec[j] + sumwy_vec[vec_back];
        sumw_vec[j] = sumw_vec[j] + sumw_vec[vec_back];
        mean_vec[j] = sumwy_vec[j] / sumw_vec[j];
        vec_back = vec_back - 1;
      }
      vec_back = vec_back + 1;
      left_vec[vec_back] = left_vec[j];
      right_vec[vec_back] = right_vec[j];
      sumwy_vec[vec_back] = sumwy_vec[j];
      sumw_vec[vec_back] = sumw_vec[j];
      mean_vec[vec_back] = mean_vec[j];
    }

    for (int k = 0; k <= vec_back; ++k) {
      for (int l = (int)left_vec[k]; l <= (int)right_vec[k]; ++l) {
        y_new[size - l - 1] = mean_vec[k];
      }
    }
  }


  return y_new;

}




// [[Rcpp::export(name = ".reg_1d_linf")]]
std::vector<double> reg_1d_linf(std::vector<double>& y, bool decreasing = false) {

  int size = (int)y.size();



  std::vector<double> error(size + 1, 0), y_new(size, 0), mean(size + 1, 0), maxy(size + 1, 0), miny(size + 1, 0), y_vec(y);
  std::vector<int> left(size + 1, 0);


  if (decreasing) {
    std::reverse(y_vec.begin(), y_vec.end());
  }

  //cout << error.size() << endl;

  mean[0] = std::numeric_limits<double>::min();

  //regression
  for (int i = 1; i <= size; ++i) {
    maxy[i] = y_vec[i - 1];
    miny[i] = y_vec[i - 1];
    mean[i] = y_vec[i - 1];
    left[i] = i;

    while ((left[i] - 1) > 0 && mean[i] <= mean[left[i] - 1]) {

      int j = left[i] - 1;
      miny[i] = miny[i]<miny[j] ? miny[i] : miny[j];
      maxy[i] = maxy[i]>maxy[j] ? maxy[i] : maxy[j];
      mean[i] = (miny[i] + maxy[i]) / double(2);
      left[i] = left[j];
    }
    double levelerror = (maxy[i] - miny[i]) / double(2);
    error[i] = levelerror>error[left[i] - 1] ? levelerror : error[left[i] - 1];
  }

  //recover means
  for (int i = size; i > 0;) {
    if (left[i] != i) {
      for (int j = left[i]; j <= i; ++j) {
        y_new[j - 1] = mean[i];
      }
      i = left[i] - 1;
    }

    else {
      y_new[i - 1] = mean[i];
      --i;
    }
  }

  if (decreasing) {
    std::reverse(y_new.begin(), y_new.end());
  }

  return y_new;
}



// [[Rcpp::export(name = ".uni_1d_linf")]]
std::vector<double> uni_1d_linf(std::vector<double>& y_vec) {

  int size = (int)y_vec.size();
  std::vector<double> error(size + 1, 0), error_asc(size + 1, 0), error_dec(size + 1, 0), y_new(size, 0), mean_asc(size + 1, 0), mean_dec(size + 1, 0), maxy(size + 1, 0), miny(size + 1, 0), y(size + 1, 0), y_rev(size + 1, 0);
  std::vector<int> left_asc(size + 1, 0), left_dec(size + 1, 0);

  for (int i = 1; i <= size; ++i) {
    y[i] = y_vec[i - 1];
  }

  //reversed
  for (int i = 1; i <= size; ++i) {
    y_rev[i] = y_vec[size - i];
  }

  mean_asc[0] = std::numeric_limits<double>::min();
  mean_dec[0] = std::numeric_limits<double>::min();


  //error calculation asc
  for (int i = 1; i <= size; ++i) {

    //setting the default max min mean left as the point itself
    maxy[i] = y[i];
    miny[i] = y[i];
    mean_asc[i] = y[i];
    left_asc[i] = i;

    //while the current point has mean less or equal to prev level, then keep merging
    // until the whole current level set
    while ((left_asc[i] - 1) > 0 && mean_asc[i] <= mean_asc[left_asc[i] - 1]) {

      int j = left_asc[i] - 1;

      //min{prev_level, current_point}, max{prev_level, current_point}
      miny[i] = miny[i]<miny[j] ? miny[i] : miny[j];
      maxy[i] = maxy[i]>maxy[j] ? maxy[i] : maxy[j];

      //the mean only keeps to the current point, even if it's prev level is lowered
      // so the stats in any point before will not change.
      // i.e. 1,3,4 and add 2.
      // 1,3,4 has stats same as themselves. 2 will have mean=3, left=position at 3, max=4, min=2
      // so 2 level sets are 1 and 3,4,2 where information of the later is kept by point at 2
      mean_asc[i] = (miny[i] + maxy[i]) / double(2);
      left_asc[i] = left_asc[j];
    }

    //error for new point, nonzero if it is lower than prev
    double levelerror = (maxy[i] - miny[i]) / double(2);
    error_asc[i] = levelerror>error[left_asc[i] - 1] ? levelerror : error_asc[left_asc[i] - 1];
  }

  //error calculation asc
  for (int i = 1; i <= size; ++i) {
    maxy[i] = y_rev[i];
    miny[i] = y_rev[i];
    mean_dec[i] = y_rev[i];
    left_dec[i] = i;

    while ((left_dec[i] - 1) > 0 && mean_dec[i] <= mean_dec[left_dec[i] - 1]) {

      int j = left_dec[i] - 1;
      miny[i] = miny[i]<miny[j] ? miny[i] : miny[j];
      maxy[i] = maxy[i]>maxy[j] ? maxy[i] : maxy[j];
      mean_dec[i] = (miny[i] + maxy[i]) / double(2);
      left_dec[i] = left_dec[j];
    }
    double levelerror = (maxy[i] - miny[i]) / double(2);
    error_dec[i] = levelerror>error[left_dec[i] - 1] ? levelerror : error_dec[left_dec[i] - 1];
  }

  //diving to 2 halves based on minimum error. points_count+1 possibilities
  // i.e. 1342 into dec1324, inc1+dec342,inc13+dec42, inc134+dec2, asc1342
  error[0] = error_dec[size];
  double min = error[0];
  int pos_min = 0;

  error[size] = error_asc[size];
  for (int i = 1; i < size; ++i) {
    error[i] = error_asc[i] + error_dec[size - i];
  }

  for (int i = 1; i <= size; ++i) {
    if (min > error[i]) {
      pos_min = i;
      min = error[i];
    }
  }

  std::vector<double> y_new_aux(size - pos_min, 0);
  //recover means
  //if the result is not straight descend
  // if it is then goes to the next if condition
  if (pos_min != 0) {
    for (int i = pos_min; i > 0;) {
      if (left_asc[i] != i) {
        for (int j = left_asc[i]; j <= i; ++j) {

          //y_new has the input y_vec size
          y_new[j - 1] = mean_asc[i];
        }
        i = left_asc[i] - 1;
      }

      else {
        y_new[i - 1] = mean_asc[i];
        --i;
      }
    }
  }

  //if the result is not straight ascend
  if (pos_min != size) {
    for (int i = size - pos_min; i > 0;) {
      if (left_dec[i] != i) {
        for (int j = left_dec[i]; j <= i; ++j) {
          y_new_aux[j - 1] = mean_dec[i];
        }
        i = left_dec[i] - 1;
      }

      else {
        y_new_aux[i - 1] = mean_dec[i];
        --i;
      }
    }
  }

  //concatenate
  for (int i = 0; i < size - pos_min; ++i) {
    y_new[pos_min + i] = y_new_aux[size - pos_min - 1 - i];
  }

  return y_new;
}



// [[Rcpp::export(name = ".pre_2d_l1_inc")]]
NumericMatrix pre_2d_l1_inc(NumericMatrix& w, NumericMatrix& data) {

  int row_size = int(data.nrow()), col_size = int(data.ncol());

  //sorted reference
  std::vector <double> sorted_data(row_size*col_size);
  for (int i = 0; i < row_size; ++i) {
    for (int j = 0; j < col_size; ++j) {
      sorted_data[i*col_size + j]=data(i,j);
    }
  }
  sort(sorted_data.begin(), sorted_data.end());

  //output matrix
  NumericMatrix out(row_size,col_size );

  //start point
  region* temp = new region;
  temp->range.reserve(row_size);
  for (int i = 0; i < row_size; ++i) {
    temp->range.push_back(make_pair(i, std::make_pair(0, col_size - 1)));
  }
  temp->l = 0;
  temp->r = row_size*col_size - 1;
  temp->m = ((temp->l + temp->r + 1) / 2);

  //queue for recursion
  std::deque<region*> collection;
  collection.push_back(temp);

  //auxiliary matrices that record actual matrix positions
  std::vector <std::vector<double> > s;
  std::vector <std::vector<int> > t;
  s.resize(row_size);
  t.resize(row_size);

for(int i=0;i<row_size;++i){
  s[i].resize(col_size + 1, 0);
  t[i].resize(col_size + 1, -1);
}
  //m1, m2 is index
  int m1 = 0, m2 = 0;

  while (!collection.empty()) {
    region* cur = collection[0];
    collection.pop_front();

    //indicate if cur needed to be deleted
    bool good_to_delete = false;

    //if only one small grid, or l==r, then fill in output matrix
    if ((cur->range.size() == 1 && cur->range[0].second.second == cur->range[0].second.first)||cur->l>=cur->r) {
      for (int i = 0; i < int(cur->range.size()); ++i) {
        int row = cur->range[i].first;
        int start_col = cur->range[i].second.first, end_col = cur->range[i].second.second;
        for (int col = start_col; col <= end_col; ++col) {

          out(row,col) = sorted_data[cur->m];
        }
      }
      good_to_delete = true;
    }

    else {
      m2 = (cur->l + cur->r + 1) / 2;
      m1 = m2 - 1;

      //calculate error
      for (int pos = 0; pos < int(cur->range.size()); ++pos) {
        int row = cur->range[pos].first, start_col = cur->range[pos].second.first, end_col = cur->range[pos].second.second;
        double rowsum = 0;
        std::vector<double> r(end_col - start_col + 2, 0), l(end_col - start_col + 2, 0);

        //calculate r
        for (int col = start_col, i = 1; col <= end_col; ++col, ++i) {
          rowsum += w(row,col) * fabs(sorted_data[m1] - data(row,col));
          r[i] = rowsum;
        }

        rowsum = 0;

        //calculate l
        for (int col = end_col, i = end_col - start_col; col >= start_col; --col, --i) {
          rowsum += w(row,col) * fabs(sorted_data[m2] - data(row,col));
          l[i] = rowsum;
        }

        //fill in s with r. If not first row, then also add min error from above right
        for (int col = start_col, i = 0; col <= end_col + 1; ++col, ++i) {
          s[row][col] = r[i] + l[i];

          //second row and so on
          if (pos != 0) {
            if (t[cur->range[pos - 1].first][col] == -1) {
              s[row][col] += s[cur->range[pos - 1].first][t[cur->range[pos - 1].first][cur->range[pos - 1].second.first]];
            }

            else s[row][col] += s[cur->range[pos - 1].first][t[cur->range[pos - 1].first][col]];
          }
        }

        //find t matrix for current row
        t[row][end_col + 1] = end_col + 1;
        for (int col = end_col; col >= start_col; --col) {
          if (s[row][col]<s[row][t[row][col + 1]]) t[row][col] = col;
          else t[row][col] = t[row][col + 1];
        }

      }


      //backtrace
      int cur_row = cur->range.back().first, cur_col = t[cur_row][cur->range.back().second.first];
      std::vector<int> trace;
      trace.reserve(cur->range.size());
      trace.push_back(cur_col);

      //they indicate if there are splits
      bool all_m1 = cur_col == cur->range.back().second.second+1, all_m2 = cur_col == cur->range.back().second.first;
      for (int i = int(cur->range.size()) - 2; i >= 0; --i) {
        cur_row = cur->range[i].first;

        //if no data at top, pick the min for whole row
        if (t[cur_row][cur_col] == -1)
          cur_col = t[cur_row][cur->range[i].second.first];
        else cur_col = t[cur_row][cur_col];
        trace.push_back(cur_col);
        if (cur_col != cur->range[i].second.second+1)  all_m1 = false;
        if (cur_col != cur->range[i].second.first)  all_m2 = false;
      }

      //if all are one m and no more m1, m2 retrieving at that direction, then fill in output
      if ((all_m1&&cur->l == m1)||( all_m2&&cur->r == m2)) {
        for (int i = 0; i < int(cur->range.size()); ++i) {
          int row = cur->range[i].first;
          int start_col = cur->range[i].second.first, end_col = cur->range[i].second.second;
          for (int col = start_col; col <= end_col; ++col) {
            out(row,col) = all_m1? sorted_data[m1]: sorted_data[m2];
          }
        }
        good_to_delete = true;
      }

      //else take new m1, m2 and add in this cur back
      else if (all_m1) {
        cur->m= cur->r = m1;
        collection.push_back(cur);
      }
      else if (all_m2) {
        cur->l = m2;
        collection.push_back(cur);
      }

      //split
      else {
        //a is left, b is right
        region* a = new region; a->range.reserve(cur->range.size());
        region* b = new region; b->range.reserve(cur->range.size());
        for (int i = 0; i <= int(cur->range.size()) - 1; ++i) {

          //if all m2
          if (trace[cur->range.size() - i - 1] == cur->range[i].second.first) {
            b->range.push_back(std::make_pair(cur->range[i].first, std::make_pair(cur->range[i].second.first, cur->range[i].second.second)));
          }

          //otherwise
          else {
            a->range.push_back(std::make_pair(cur->range[i].first, std::make_pair(cur->range[i].second.first, trace[cur->range.size() - i - 1] - 1)));
            if (trace[cur->range.size() - i - 1]  <= cur->range[i].second.second) {
              b->range.push_back(std::make_pair(cur->range[i].first, std::make_pair(trace[cur->range.size() - i - 1], cur->range[i].second.second)));
            }
          }
        }

        //the right of left region is m1, vice verse
        a->l = cur->l;
        a->m = a->r = m1;
        b->r = cur->r;
        b->m = b->l = m2;
        collection.push_back(a);
        collection.push_back(b);
        good_to_delete = true;
      }
    }

    //recover s and t
    for (int i = 0; i < int(cur->range.size()); ++i) {
      int row = cur->range[i].first;
      int start_col = cur->range[i].second.first, end_col = cur->range[i].second.second;
      for (int col = start_col; col <= end_col + 1; ++col) {
        s[row][col] = 0;
        t[row][col] = -1;
      }
    }
    if(good_to_delete) delete cur;
  }
  return out;
}



// [[Rcpp::export(name = ".pre_2d_l1_inc")]]
NumericMatrix pre_2d_l2_inc(NumericMatrix & w, NumericMatrix &data) {


  int row_size = int(data.nrow()), col_size = int(data.ncol());

  //output matrix
  NumericMatrix out(row_size, col_size);

  //start point
  region* temp = new region;
  temp->range.reserve(row_size);
  for (int i = 0; i < row_size; ++i) {
    temp->range.push_back(std::make_pair(i, std::make_pair(0, col_size - 1)));
  }

  //queue for recursion
  std::deque<region*> collection;
  collection.push_back(temp);

  //auxiliary matrices that record actual matrix positions
  std::vector <std::vector<double> > s;
  std::vector <std::vector<int> > t;
  s.resize(row_size);
  t.resize(row_size);
  for (int i = 0; i < int(s.size()); ++i) {
    s[i].resize(col_size, 0);
    t[i].resize(col_size, -1);
  }

  while (!collection.empty()) {
    region* cur = collection[0];
    collection.pop_front();
    double x;
    double sumw = 0;
    double sum = 0;

    for (int j = 0; j < int(cur->range.size()); ++j) {
      std::pair<int, std::pair<int, int> > i = cur->range[j];

      int row = i.first;
      int start_col = i.second.first, end_col = i.second.second;
      for (int col = start_col; col <= end_col; ++col) {
        sumw += w(row,col);
        sum += (w(row,col) * data(row,col));

      }
    }
    x = sum / (double)sumw;
    //weighted_avg(*cur, w, data,x);

    //min_index is the index in the region std::vector where smallest happen
    int min_row = 0, min_col = 0, min_index = 0;

    //compute s and t and min_row and min_col
    for (int i = 0; i < int(cur->range.size()); ++i) {

      int row = cur->range[i].first;
      int start_col = cur->range[i].second.first, end_col = cur->range[i].second.second;
      double rowsum = 0;

      //fill in s. If not first row, then also add min error from above right
      for (int col = start_col; col <= end_col; ++col) {

        rowsum += w(row,col) * (data(row,col) - x);
        s[row][col] = rowsum;

        //second row and so on
        if (i != 0) {
          if (t[cur->range[i - 1].first][col] == -1) {
            s[row][col] += s[cur->range[i - 1].first][t[cur->range[i - 1].first][cur->range[i - 1].second.first]];
          }

          else s[row][col] += s[cur->range[i - 1].first][t[cur->range[i - 1].first][col]];
        }
      }

      //find t matrix for current row
      t[row][end_col] = end_col;
      for (int col = end_col - 1; col >= start_col; --col) {
        if (s[row][col]<s[row][t[row][col + 1]]) t[row][col] = col;
        else t[row][col] = t[row][col + 1];
      }

      //find min_col and min_row
      // for first row min_col is t[row][start_col]
      if (i == 0) {
        min_row = row;
        min_col = t[row][start_col];
      }

      //for further rows min
      else {
        if (s[row][t[row][start_col]] <= s[min_row][min_col]) {
          min_row = row;
          min_col = t[row][start_col];
          min_index = i;
        }
      }
    }

    if (cur->range.size() == 1) {
      //first_row_situation(cur, collection, t, out, x);
      int row = cur->range[0].first;
      int start_col = cur->range[0].second.first, end_col = cur->range[0].second.second;

      //check if k==m
      if (t[row][start_col] == end_col) {
        for (int col = start_col; col <= end_col; ++col) {
          out(row,col) = x;
        }
      }

      else {
        region* a = new region;
        region* b = new region;
        a->range.push_back(std::make_pair(row, std::make_pair(start_col, t[row][start_col])));
        b->range.push_back(std::make_pair(row, std::make_pair(t[row][start_col] + 1, end_col)));
        collection.push_back(a);
        collection.push_back(b);
      }

    }

    else {
      //changed
      //if no negative entry, treat the region as a single levelset
      if (s[min_row][min_col] >= 0) {

        for (int j = 0; j < int(cur->range.size()); ++j) {
          std::pair<int, std::pair<int, int> > i = cur->range[j];
          for (int col = i.second.first; col <= i.second.second; ++col) {
            out(i.first,col) = x;
          }
        }
      }

      else {
        int cur_row = min_row, cur_col = min_col;
        std::vector<int> trace;
        trace.reserve(min_index + 1);
        trace.push_back(cur_col);

        //backtrace
        for (int i = min_index - 1; i >= 0; --i) {
          cur_row = cur->range[i].first;
          if (t[cur_row][cur_col] == -1) {
            cur_col = t[cur_row][cur->range[i].second.first];
          }
          else cur_col = t[cur_row][cur_col];
          trace.push_back(cur_col);
        }

        //a is left, b is right, a is never empty
        region* a = new region;
        region* b = new region;
        a->range.resize(0);
        b->range.resize(0);
        a->range.reserve(cur->range.size());
        b->range.reserve(cur->range.size());
        for (int i = 0; i <= min_index; ++i) {
          a->range.push_back(std::make_pair(cur->range[i].first, std::make_pair(cur->range[i].second.first, trace[min_index - i])));
          if (trace[min_index - i] + 1 <= cur->range[i].second.second) {
            b->range.push_back(std::make_pair(cur->range[i].first, std::make_pair(trace[min_index - i] + 1, cur->range[i].second.second)));
          }
        }

        //if the left subregion contain the region's last row
        if (min_index == int(cur->range.size() - 1)) {

          //if left subregion is the entire region
          if (b->range.empty()) {
            delete a;
            delete b;

            for (int j = 0; j < int(cur->range.size()); ++j) {
              std::pair<int, std::pair<int, int> > i = cur->range[j];
              for (int col = i.second.first; col <= i.second.second; ++col) {
                out(i.first,col) = x;
              }
            }
          }
          else {
            collection.push_back(a);
            collection.push_back(b);
          }

        }
        else {
          //still few rows beneath left region
          collection.push_back(a);

          for (int i = min_index + 1; i < int(cur->range.size()); ++i) {
            b->range.push_back(std::make_pair(cur->range[i].first, std::make_pair(cur->range[i].second.first, cur->range[i].second.second)));
          }
          collection.push_back(b);
        }
      }
    }


    //recover s and t
    for (int i = 0; i < int(cur->range.size()); ++i) {
      int row = cur->range[i].first;
      int start_col = cur->range[i].second.first, end_col = cur->range[i].second.second;
      for (int col = start_col; col <= end_col; ++col) {
        s[row][col] = 0;
        t[row][col] = -1;
      }
    }
    delete cur;
  }
  return out;
}

// [[Rcpp::export]]
std::vector<double> reg_1d(std::vector<double>& y_vec, std::vector<double>& w_vec, int metric,
                           bool unimodal = false, bool decreasing = false) {

  if (y_vec.size() == 0) stop("Empty data");
  if (w_vec.empty()) w_vec.resize(y_vec.size(), 1);
  if (w_vec.size() != y_vec.size()) stop("Data and weight have different number of entries");
  for (int i = 0; i < int(w_vec.size()); ++i) {
    if (w_vec[i]<0) {
      stop("Negative weight detected");
    }
  }
  std::vector<double> out;
  if (metric == 1) {
    std::vector<double> error;
    error.resize(y_vec.size(), 0);
    if (!unimodal) {
      if (decreasing) {
        out = reg_1d_l1(y_vec, w_vec, error, true);
      }
      else {
        out = reg_1d_l1(y_vec, w_vec, error);
      }
    }
    else {
      out = uni_1d_l1(y_vec, w_vec);
    }
  }
  else if (metric == 2) {
    if (!unimodal) {
      if (decreasing) {
        out = reg_1d_l2(y_vec, w_vec, true);
      }
      else {
        out = reg_1d_l2(y_vec, w_vec);
      }
    }
    else {
      out = uni_1d_l2(y_vec, w_vec);
    }
  }
  else if (metric == 3) {
    if (!unimodal) {
      if (decreasing) {
        out = reg_1d_linf(y_vec, true);
      }
      else {
        out = reg_1d_linf(y_vec);
      }
    }
    else {
      out = uni_1d_linf(y_vec);
    }
  }
  else {
    stop("Metric does not exist");
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix reg_2d(NumericMatrix& y_vec,
                                NumericMatrix& w_vec, int metric
) {

  //error check
  if (y_vec.nrow() == 0 || y_vec.ncol()== 0) stop("Empty data");
  if (w_vec.nrow() == 0 || w_vec.ncol()== 0) {
    w_vec=NumericMatrix(y_vec.nrow(),y_vec.ncol());
  }

  if (w_vec.nrow() != y_vec.nrow()) stop("Data and weight have different number of rows");
  if (w_vec.ncol() != y_vec.ncol()) stop("Data and weight have different number of columns");

  for (int i = 0; i < w_vec.nrow(); ++i) {
    for (int j = 0; j < w_vec.ncol(); ++j) {
      if (w_vec(i,j) < 0) {
        stop("Negative weight detected");
      }
    }
  }

  NumericMatrix out;
  if (metric == 1) {
    out = pre_2d_l1_inc(w_vec, y_vec);
    return out;
  }
  else if (metric == 2) {
    out = pre_2d_l2_inc(w_vec, y_vec);
    return out;
  }
  else {
    stop("metric does not exist");
  }



}
