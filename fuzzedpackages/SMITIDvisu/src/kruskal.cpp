/* This file is part of SMITIDvisu package.
* Copyright (C) 2018-2019 Jean-François Rey <jean-francois.rey@inra.fr>

* SMITIDvisu is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SMITIDvisu is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with SMITIDvisu. If not, see <https://www.gnu.org/licenses/>.
*/


#include "kruskal.h"

// compute mst minimum spanning tree of a matrix representing edges between nodes
NumericMatrix mstCompute( const NumericMatrix &mat)
{

  NumericMatrix out(mat.nrow(), mat.ncol());
  colnames(out) = colnames(mat);
  rownames(out) = rownames(mat);

  long long cost = 0;

  vector< pair <long double, pair<int, int> > > edges;

  
  for(int l=0; l<mat.nrow(); l++)
    for(int c=l+1; c<mat.ncol(); c++)
    {
      if(mat(l,c) != 0) edges.push_back(make_pair(mat(l,c), make_pair(l,c)));
    }

  sort(edges.begin(), edges.end(), edgesCompare);
  
  cost = kruskal(edges, out);
   
  return out;
}

// kruskal algorithm
double kruskal(vector< pair <long double, pair<int, int> > > edges, NumericMatrix &out) 
{
  double minCost = 0.0;
  vector< pair <long double, pair<int, int> > >::iterator it_edges;
  
  int * tree = new int[out.ncol()];
  for(int i=0; i<out.ncol(); i++) tree[i] = i;
  
  for( it_edges=edges.begin(); it_edges!=edges.end(); it_edges++ )
  {
    if( rootTree(tree, it_edges->second.first) != rootTree(tree, it_edges->second.second) )
    {
      minCost += it_edges->first;
      addTree(tree, it_edges->second.first, it_edges->second.second);
      out(it_edges->second.first,it_edges->second.second) = out(it_edges->second.second,it_edges->second.first) = 1;
    }
  }
  
  delete[] tree;

  return minCost;
}

// search in disjoint-set data the root node
int rootTree( int * tree, int node )
{
  while(tree[node] != node) node = tree[node];

  return node;
}

// add disjoint-set data a new node
void addTree( int * tree, int node1, int node2)
{
  tree[rootTree(tree, node2)] = tree[rootTree(tree, node1)];
}

// compare two node by their weight
bool edgesCompare(const pair <long double, pair <int, int> > &i, const pair <long double, pair <int, int> > &j)
{
  return i.first < j.first;
}

