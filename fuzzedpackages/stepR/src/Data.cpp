#include "Data.h"

/***************
* class Data
* abstract class maintaining the data
* Florian Pein, 2015
***************/
NumericVector Data::criticalValues_;

void Data::setCriticalValues(const List &input) {
  criticalValues_ = input["q"];
}

void Data::cleanUpLocalVariables() {}

void Data::setLocal(const List &input) {}

double Data::computeSingleStat(unsigned int startIndex, unsigned int leftSegment, unsigned int rightSegment) const {
  return computeSingleStat(0.0);
}

Data::~Data() {}
