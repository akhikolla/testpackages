#include "Node.h"
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

Node::Node()
{

}

Node::Node(vector<int> Rs, int eInt) {
	Rsum = Rs;
	equiInt = eInt;
	nodeTable.resize(199);
}

Node::Node(vector<int> Rs)
{
	Rsum = Rs;
	std::sort(Rs.begin(), Rs.end());

	int eInt = 0;
	for (size_t x = 0; x < Rs.size(); x++) {
		eInt *= 127;
		eInt += Rs[x];
	}
	equiInt = eInt;

}

void Node::show()
{
    /*
	for (size_t i = 0; i < Rsum.size(); i++) {
		cout << Rsum[i] << " ";
	}
	cout << " :: ";
	for (size_t i = 0; i < ChildrenIndex.size(); i++) {
		cout << ChildrenIndex[i] << " ";
	}
	cout << "::" << ub << " :: " << lb << endl;
    */ 
}

void Node::addRsum(int Rs)
{
	Rsum.push_back(Rs);
}

vector<int> Node::getRsum()
{
	return Rsum;
}

int Node::getChildrenIndex(int ind)
{
	return ChildrenIndex[ind];
}

double Node::getLengthToChildren(int ind) {
	return lengthToChildren[ind];
}

double Node::getColChisqToChildren(int ind) {
	return colChisqToChildren[ind];
}

void Node::setUB(double val) {
	ub = val;
}
void Node::setLB(double val) {
	lb = val;
}

double Node::getLB() {
	return lb;
}

double Node::getUB() {
	return ub;
}

void Node::setEquiInt(int x) {
	equiInt = x;
}

int Node::getEquiInt() {
	return equiInt;
}

int Node::getSize() {
	return (int) ChildrenIndex.size();
}

void Node::addPastLen(double len, double chisq) {
	long long int y = (long long int) chisq;

	int i = (int) y % nodeTable.size();
	size_t j = 0;
	while (j < nodeTable[i].size() && nodeTable[i][j].first != y) j++;

	if (j < nodeTable[i].size()) pastLen[nodeTable[i][j].second] += len;

	else {
		pastLen.push_back(len);
		pastChisq.push_back(y);

		//update hash table:
		nodeTable[i].push_back(std::make_pair(y, pastChisq.size()-1));
	}

	/*
	int i = 0;
	while (i < pastChisq.size() && (pastChisq[i] != y)) i++;

	if (i < pastChisq.size()) {
		pastLen[i] += x;
	}
	else {
		pastLen.push_back(x);
		pastChisq.push_back(y);
	}
	*/
}

double Node::getPastLen(int ind) {
	return pastLen[ind];
}

int Node::getPastSize() {
	return (int) pastLen.size();
}

double Node::getPastChisq(int ind) {
	return pastChisq[ind];
}

void Node::setLengthToEnd(double val) {
	lengthToEnd = val;
}

double Node::getLengthToEnd() {
	return lengthToEnd;
}

int Node::getPastChisqSize() {
	return (int) pastChisq.size();
}

void Node::quicksort(int left, int right) {

	double pivot = pastChisq[(left + right) / 2];
	int i, j;
	//double temp;
	i = left;
	j = right;

	while (i <= j) {
		while (pastChisq[i] < pivot) i++;
		while (pastChisq[j] > pivot) j--;

		if (i <= j) {
			//swap arr[i] and arr[j]
			std::swap(pastChisq[i], pastChisq[j]);
			std::swap(pastLen[i], pastLen[j]);
			i++;
			j--;
		}
	}

	if (left < j) quicksort(left, j);
	if (i < right) quicksort(i, right);
}


int Node::bSearch(int chisq) {
	auto it = std::lower_bound(pastChisq.begin(), pastChisq.end(), chisq);
	std::size_t index = std::distance(pastChisq.begin(), it);
	return (int) index;
}


int Node::isChildInList(int x) {
	size_t i = 0;
	while (i < ChildrenIndex.size() && x != ChildrenIndex[i]) i++;
	if (i < ChildrenIndex.size()) 	return (int) i;
	return -1;
}

void Node::addLength(double x, int index) {
	lengthToChildren[index] += x;
}


void Node::addChildLink(int index, double len, double colchisq) {
	ChildrenIndex.push_back(index);
	lengthToChildren.push_back(len);
	colChisqToChildren.push_back(colchisq);
}

void Node::setColChisqToChildren(int ind, double colchisq) {
	colChisqToChildren[ind] = colchisq;
}


void Node::setMinPastChisq(double val) {
  minPastChisq = val;
}

double Node::getMinPastChisq() {
  return minPastChisq;
}

void Node::setMaxPastChisq(double val) {
  maxPastChisq = val;
}

double Node::getMaxPastChisq() {
  return maxPastChisq;
}
