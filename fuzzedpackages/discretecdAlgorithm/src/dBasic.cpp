#include "dBasic.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>
#define MATHLIB_STANDALONE
#define RMATH_DLL
#include <Rmath.h>
#include <chrono>
#include <random>
#include <time.h>
#include <stdlib.h>

#include <string>
#include <cstdlib>

using namespace Eigen;
using namespace std;

//fix seed
// std::default_random_engine generator(14837);
//random seed
// std:: default_random_engine generator((unsigned int) time(nullptr));

// std::default_random_engine generator;
// std:: uniform_real_distribution<double> distribution(0.0,1.0);

// Check whether adding a directed edge from node a to node b induces cycles in G.

const int Cycle(const int& node, const int* G, const int& a, const int& b)
{
	if(a==b) return 1;
    int i,j,lo,nBot=0,nTop=0,SLeng=1,bCycle=0;
    int* color=new int[node]();
    int* S=new int[node]();
    color[a-1]=1;
    S[0]=a;
    while(SLeng>0)
    {
        i=S[nBot];
		SLeng--;
		nBot++;
		for(j=0;j<node;j++)
		{
			lo=(i-1)*node;
			if(G[lo+j]==1)
			{
				if((j+1)==b)
				{
					bCycle=1;
					break;
				}
				else if(color[j]==0)
				{
					nTop++;
					S[nTop]=j+1;
					SLeng++;
					color[j]=1;
				}
			}
		}
		if(bCycle==1) break;
    }
    delete [] color;
    delete [] S;
    return bCycle;
}



// Check whether a graph G has cycles.
bool check(const int& node, const int* G) {
    int* empty = new int[node*node]();
	for(int i = 0; i < node; ++i)
	{
		for(int j = 0; j < node; ++j)
		{
			int pos = i*node + j;
			if(G[pos] == 1)
			{
				if(Cycle(node, empty, j+1, i+1))
				{
					delete [] empty;
					return true;
				}
				empty[pos] = 1;
			}
		}
	}
	delete [] empty;
	return false;
}



// Sorts a into decreasing order, and applies the same permutation to ib.
void revsort(double *a, int *ib, int n)  // a -> prob_array ib -> permu_array
{
    /* Sort a[] into descending order by "heapsort";
     * sort ib[] alongside;
     * if initially, ib[] = 1...n, it will contain the permutation finally
     */

    int l, j, ir, i;
    double ra;
    int ii;

    if (n <= 1) return; // return if there is only one node.

    a--; ib--; // to make index of the ith node "i". first element is 0.

    l = (n >> 1) + 1; // l-1 is the last node which is not a leaf
    ir = n; // record the index to be set.

    for (;;) {
		if (l > 1) { // when the mother node is not the first node.
			l = l - 1; // move forward
			ra = a[l]; // temporary record for the parent node
			ii = ib[l]; // ib goes together with a
		}
		else { // the 1st element in the que is the smallest, move the first element to irth position. sorting part.
			ra = a[ir];
			ii = ib[ir];
			a[ir] = a[1];
			ib[ir] = ib[1];
			if (--ir == 1) {
				a[1] = ra;
				ib[1] = ii;
				return;
			}
		}
		i = l; // parent node
		j = l << 1; // left child: j*2
		while (j <= ir) { // percolate down for j
			if (j < ir && a[j] > a[j + 1]) ++j; //choose the smaller one between the two child (if there exists two child)
			if (ra > a[j]) { //make the parent smaller than child
				a[i] = a[j];
				ib[i] = ib[j];
				j += (i = j); // let i equals to j, and j the left child of ti self
			}
			else
				j = ir + 1; // terminate the loop
		}
		a[i] = ra;
		ib[i] = ii;
    }
}



// Unequal probability sampling without-replacement.
template <class T>
vector<T> ProbSampleNoReplace(const vector<T>& value, int nans, vector<double>& prob)
{
    double rT, mass, totalmass;
    int i, j, k, n1, n = value.size();
	vector<int> perm;
	vector<T> ans;

    /* Record element identities */
    for (i = 0; i < n; i++)
		perm.push_back(i);

    /* Sort probabilities into descending order */
    /* Order element identities in parallel */
	revsort(&prob[0], &perm[0], n);

    /* Compute the sample */
    totalmass = 1;
    for (i = 0, n1 = n-1; i < nans; i++, n1--) {
        	rT = totalmass * unif_rand(); // unif_rand is from r.h
        // rT = totalmass * distribution(generator);
		mass = 0;
		for (j = 0; j < n1; j++) {
			mass += prob[j];
			if (rT <= mass)
                break;
		}
		ans.push_back(value[perm[j]]);
		totalmass -= prob[j];
		for(k = j; k < n1; k++) {
			prob[k] = prob[k + 1];
			perm[k] = perm[k + 1];
		}
    }

	return ans;
}




// Equal probability sampling without-replacement
template <class T>
vector<T> SampleNoReplace(const vector<T>& value, int nans)
{
	int n = value.size();
	vector<T> ans;

	if (nans < 2)
	{
        	ans.push_back(value[static_cast<int>(n * unif_rand())]);
        // ans.push_back(value[static_cast<int>(n * distribution(generator))]);
		return ans;
	}

	int i, j;
	vector<int> x;

    for (i = 0; i < n; ++i)
		x.push_back(i);
    for (i = 0; i < nans; ++i) {
        	j = static_cast<int>(n * unif_rand());
        // j = static_cast<int>(n * distribution(generator));
		ans.push_back(value[x[j]]);
		x[j] = x[--n];
    }
	return ans;
}



vector<int> seq(const int& from, const int& to)
{
	vector<int> ans;
	for (int i = from; i <= to; ++i)
	{
		ans.push_back(i);
	}
	return ans;
}



// Generate the disdtribution of in-degrees (i.e., the number of parents of a node) given the maximal in-degree and the total number of edges in a graph.
vector<int> degreeG(const int& maxdeg, const int& node, const int& nedge)
{
	vector<int> degree_dist(maxdeg + 1, 0);
	int nx = node, nex = nedge, sum = 0, lower, upper;
	if (maxdeg >= 2)
	{
		for (int i = maxdeg; i > 1; --i)
		{
			lower = max(0, static_cast<int>(ceil(nex - (nx - i / 2.0) * (i - 1)))); // lower * i + (nx - lower - (i - 1)) * (i - 1) + (i - 1) * (i - 2) / 2 >= nex
			upper = min(nx - i, static_cast<int>(floor(nex / static_cast<double>(i))));
			if (lower < upper)	degree_dist[i] = SampleNoReplace(seq(lower, upper), 1)[0];
			else if (lower == upper)	degree_dist[i] = lower;
			else continue;
			nx -= degree_dist[i];
            nex -= i * degree_dist[i];
		}
	}
	degree_dist[1] = nex;
	for (int j = 1; j <= maxdeg; ++j)
	{
		sum += degree_dist[j];
	}
    degree_dist[0] = node - sum;
    return degree_dist;
}



// Generate DAGs given the maximum in-degree and the total number of edges.
gStruct GGen(const int& maxdeg, const int& node, const int& nedge)
{
	gStruct ans;
	ans.ordex = MatrixXi::Zero(maxdeg, node); // matrix listing the parents of each node (similar to the adjacency list)
	ans.trueG = MatrixXi::Zero(node, node);
	vector<int> ts = SampleNoReplace(seq(1, node), node); // topological sort
	ans.ts = ts;
	vector<int> degree_dist = degreeG(maxdeg, node, nedge); // true DAG
	vector<int> indegree(node, 0);	// indegree is matched to ts, to record in-degree for each node.
	vector<bool> ind(node, true);  // to record wether the node has been assigned to an in-degree. true for not been assigned to.

	vector<int> candidate, pos, pa;
	int i, j, np, ppa;
	vector<int>::iterator iter;

	for (i = maxdeg; i >= 0; --i) // generate indegree for each node, start from the maximum in-degree, use: "indegree", "ind", "iter", "candidate", "pos".
	{
		if(degree_dist[i] > 0)
		{
			for (j = i; j < node; ++j) // start from i, since parents can only from nodes that are previous than the child in vector ts.
			{
				if (ind[j])	candidate.push_back(j);
			}
			if(candidate.size() == 1)
			{
				indegree[candidate[0]] = i;
				ind[candidate[0]] = false;
			}
			else
			{
				pos = SampleNoReplace(candidate, degree_dist[i]);
				for (iter = pos.begin(); iter < pos.end(); ++iter)
				{
					indegree[*iter] = i;
					ind[*iter] = false;
				}
			}
            //			candidate.~vector();
            while (!candidate.empty()) // instead of using destructor
            {
                candidate.pop_back();
            }

		}
	}

	for (i = 0; i < node; ++i) // generate ordex： i is the index for ts.
	{
		np = indegree[i];
		if(np == 0)	continue;
		else if(i == 1)	ans.ordex(0,ts[i] - 1) = ts[0];
		else
		{
			pa = SampleNoReplace(vector<int>(ts.begin(), ts.begin() + i), np);
			sort(pa.begin(), pa.end());
			for (j = 0; j < np; ++j)
				ans.ordex(j, ts[i] - 1) = pa[j];
		}
	}

	for (i = 0; i < node; ++i) // generate G
	{
		for (j = 0; j < maxdeg; ++j)
		{
			ppa = ans.ordex(j, i);
			if (ppa)	ans.trueG(ppa - 1, i) = 1;
		}
	}

	return ans;
}



// Generate Markov chains.
gStruct MCGen(const int& node)
{
	gStruct ans;
	ans.ordex = MatrixXi::Zero(1, node); // matrix listing the parents of each node (similar to the adjacency list)
	ans.trueG = MatrixXi::Zero(node, node);
	vector<int> ts = SampleNoReplace(seq(1, node), node); // topological sort
	ans.ts = ts;

	int ppa;

	for (int i = 0; i < node; ++i)
	{
		if(i == 0)	continue;
		else ans.ordex(0,ts[i] - 1) = ts[i - 1];
	}

	for (int i = 0; i < node; ++i)
	{
		for (int j = 0; j < 1; ++j)
		{
			ppa = ans.ordex(j, i);
			if (ppa)	ans.trueG(ppa - 1, i) = 1;
		}
	}

	return ans;
}

bool ifIntervene(int node, vector<int> subIvn)
{
  for (int i=0; i<subIvn.size(); i++) {
    if (node==subIvn[i])
      return true;
  }
  return false;
}

// generate mixture of interventionl and observational data.
// ordex, list of parent set
// ts, topological sort
// ivn, list of intervention
// nlevels, number of levels for each node
// data, data set that need to be generated
// coef,
void DatGen(const MatrixXi& ordex, const vector<int>& ts, const vector< vector<int> >& ivn, const vector< vector<int> >& ivn_vals, const bool ivn_rand, const VectorXi& nlevels, MatrixXi& data, const vector<VectorXMXd>& coef)
{

  int node = ordex.cols(), maxdeg = ordex.rows();

  vector< vector<int> > levels; // node by 2 matrix, stores level Index for each node.
  for (int it1 = 0; it1 < node; ++it1)
  {
    levels.push_back(seq(0, nlevels(it1) - 1));
  }

  int dataSize = ivn.size(), cur, pa;// counter -> i, cur -> j
  vector< vector<int> > fX(dataSize);
  if (!ivn_rand) {
    for (int it1=0; it1 < dataSize; ++it1)
    {
      for (int it2=0; it2<ivn_vals[it1].size(); ++it2) {
        fX[it1].push_back(ivn_vals[it1][it2]);
      }
    }
  }
  else {
    for (int it1=0; it1 < dataSize; ++it1)
    {
      for (int it2=0; it2<ivn[it1].size(); ++it2) {
        if (ivn[it1][it2] != -1)
          fX[it1].push_back(SampleNoReplace(levels[ivn[it1][it2]], 1)[0]);
      }
    }
  }

  double rowStat;

  for (int itr = 0; itr < dataSize; ++itr)
  {
    vector<int> iNode = ivn[itr];
    for (int itc = 0; itc < node; ++itc) {
      cur = ts[itc] - 1;
      if (ifIntervene(cur, iNode))
        data(itr, cur) = fX[itr][0];
      else if ((ordex.col(cur).array() == 0).all())
        data(itr, cur) = SampleNoReplace(levels[cur], 1)[0];
      else
      {
        VectorXd logit = VectorXd::Zero(nlevels(cur));
        for (int itp = 0; itp < maxdeg; ++itp)
        {
          pa = ordex(itp, cur);
          if (pa!=0)
          {
            for (int l=0; l<(nlevels[cur]); l++) {
              if (data(itr, pa-1)) {
                logit(l) += coef[cur](itp+1)(l, (data(itr, pa-1)-1)); // for coefficient matrix, row is parent, col is child
              }
              logit(l) += coef[cur](0)(l, 0); // intercept
            }
          }
        }
        rowStat = logit.maxCoeff();
        logit = logit.array() - rowStat;
        logit = logit.array().exp();
        rowStat = logit.sum();
        logit /= rowStat;
        vector<double> temp = vector<double> (logit.data(), logit.data() + nlevels(cur)); // By Jean
        data(itr, cur) = ProbSampleNoReplace(levels[cur], 1, temp)[0]; // By Jean
      }
    }
  }
}

void DatGen_obs(const MatrixXi& ordex, const vector<int>& ts, int nobs, const VectorXi nlevels,
                MatrixXi& data, MatrixXVXi& levelIndex, double coef)
{

  int node = ordex.cols(), maxdeg = ordex.rows();
  //	VectorXi nlevels = VectorXi::Constant(node, 2); // there are only two levels for each node.
  //vector<int> ivn = SampleNoReplace(seq(0, node - 1), noi); //
  //sort(ivn.begin(), ivn.begin() + noi);

  vector< vector<int> > levels; // node by 2 matrix, stores level Index for each node.
  for (int it1 = 0; it1 < node; ++it1)
  {
    levels.push_back(seq(0, nlevels(it1) - 1));
  }

  int dataSize = nobs, Counter = -1, cur, pa;

  double rowStat;

  for (int itr = 0; itr < nobs; ++itr)
  {
    ++Counter; // update i
    for (int itc = 0; itc < node; ++itc)
    {
      cur = ts[itc] - 1; // update j
      if ((ordex.col(cur).array() == 0).all())
        data(Counter, cur) = SampleNoReplace(levels[cur], 1)[0];
      else
      {
        VectorXd logit = VectorXd::Zero(nlevels(cur));
        for (int itp = 0; itp < maxdeg; ++itp)
        {
          pa = ordex(itp, cur);
          if (pa != 0)
          {
            //                        logit(0) += coef - coef * data(Counter, pa - 1);      // coef denotes the magnitude of influence
            //                        logit(1) += coef * data(Counter, pa - 1);
            logit(data(Counter, pa-1)) += coef;
          }
        }
        rowStat = logit.maxCoeff();
        logit = logit.array() - rowStat;
        logit = logit.array().exp();
        rowStat = logit.sum();
        logit /= rowStat;
        vector<double> temp = vector<double> (logit.data(), logit.data() + nlevels(cur)); // By Jean
        data(Counter, cur) = ProbSampleNoReplace(levels[cur], 1, temp)[0]; // By Jean
        //					data(Counter, cur) = ProbSampleNoReplace(levels[cur], 1, vector<double> (logit.data(), logit.data() + nlevels(cur)))[0];
      }
    }
  }
  VectorXi dummy(nlevels.maxCoeff());
  for (int it1 = 0; it1 < nlevels.maxCoeff(); ++it1)
    dummy(it1) = it1;
  for (int it1 = 0; it1 < node; ++it1)
  {
    for (int it2 = 0; it2 < node; ++it2)
    {
      if (it2 != it1)
        levelIndex(it2, it1) = dummy.head(nlevels(it2) - 1);
      else
        levelIndex(it2, it1) = dummy.head(nlevels(it2));
    }
  }
}
