
#ifndef NODE_H
#define NODE_H

#include <vector>

using namespace std;

class Node {

public:
	Node();
	Node(vector<int> Rs);
	Node(vector<int> Rs, int eInt);

	void show();

	void addRsum (int Rs);
	vector<int> getRsum();

	int getChildrenIndex(int ind);
	int isChildInList(int ind);

	void addLength(double val, int ind);
	double getLengthToChildren(int ind);

	void setColChisqToChildren(int ind,  double colchisq);
	double getColChisqToChildren(int ind);

	void addPastLen(double len,  double chisq);
	double getPastLen(int ind);

	double getPastChisq(int ind);
	int bSearch(int chisq);
	int getPastChisqSize();

	int getPastSize();

	void setUB(double val);
	double getUB();

	void setLB(double val);
	double getLB();

	void setEquiInt(int x);
	int getEquiInt();

	int getSize();

	void setLengthToEnd(double val);
	double getLengthToEnd();

	void quicksort(int left, int right);

	void addChildLink(int index, double len, double colchisq);

	void setMinPastChisq(double val);
	double getMinPastChisq();

	void setMaxPastChisq(double val);
	double getMaxPastChisq();

private:
	vector<int> Rsum;	// the remaining rowsums after the previous columns are enumerated
	int equiInt;		// the equivalent integer converted from Rsum to compare the nodes faster
	double lengthToEnd; // the length from this node to the end node, used when the whole branch is counted

	double ub;			// the upper bound for this node
	double lb;			// the lower bound for this node

	vector<int> ChildrenIndex;			// to store the indices of the children node
	vector<double> lengthToChildren;	// to store the lenghths to the children node
	vector<double> colChisqToChildren;	// to store the weights to the children node (weight = partial funchisq for the column enumerated)

	vector<double> pastLen;		// the cummulative lengths from the start node to this node, each entry of this list will be called lengthSoFar
	vector<double> pastChisq;	// the cummulative weights from the start node to this node, each entry of this list will be called chisqSoFar or weightSoFar

	vector<vector<pair <long long int, int>>> nodeTable;
	// first long long int is the chisq, second int is the index of pastLen

	double minPastChisq;
	double maxPastChisq;
};

#endif
