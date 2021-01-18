/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
vector<Element>
Matrix::trunk(vector<Element>& edges) {
  mdsize sznan = medusa::snan();
  mdreal rlnan = medusa::rnan();

  /* Determine the number of vertices. */
  mdsize nvert = 0;
  for(mdsize i = 0; i < edges.size(); i++) {
    Element& e = edges[i];
    if(e.row >= nvert) nvert = (e.row + 1);
    if(e.column >= nvert) nvert = (e.column + 1);
  }

  /* Prepare data structures. */
  vector<mdsize> trees(nvert, sznan);
  unordered_map<mdsize, vector<mdsize> > forest;

  /* Kruskal's algorithm. */
  mdsize ntree = (nvert - 1);
  unordered_set<mdsize> mask;
  for(mdsize i = 0; i < edges.size(); i++) {
    Element& e = edges[i];
    if(e.row == e.column) panic("Bad edge.", __FILE__, __LINE__);
    if(e.value == rlnan) panic("Bad value.", __FILE__, __LINE__);
    if(mask.size() >= ntree) break;
  
    /* New tree. */
    mdsize a = e.row;
    mdsize b = e.column; 
    mdsize treeA = trees[a];
    mdsize treeB = trees[b];
    if((treeA == sznan) && (treeB == sznan)) {
      forest[i].push_back(a);
      forest[i].push_back(b);
      trees[a] = i;
      trees[b] = i;
      mask.insert(i);
      continue;
    }

    /* Append to an existing tree. */
    if(treeA == sznan) {
      forest[treeB].push_back(a);
      trees[a] = treeB;
      mask.insert(i);
      continue;
    }
    if(treeB == sznan) {
      forest[treeA].push_back(b);
      trees[b] = treeA;
      mask.insert(i);
      continue;
    }

    /* Check that smaller tree is merged into the larger. */
    if(treeA == treeB) continue;
    if(forest[treeA].size() < forest[treeB].size()) {
      mdsize tmp = treeA;
      treeA = treeB;
      treeB = tmp;
    }

    /* Merge the two trees. */
    vector<mdsize>::iterator it;
    vector<mdsize>& membA = forest[treeA];
    vector<mdsize>& membB = forest[treeB];
    for(it = membB.begin(); it != membB.end(); it++) {
      membA.push_back(*it);
      trees[*it] = treeA;
    }

    /* Remove the smaller tree. */
    forest.erase(treeB);
    mask.insert(i);
  }

  /* Remove the tree from graph. */
  mdsize nremain = 0;
  vector<Element> array;
  for(mdsize i = 0; i < edges.size(); i++) {
    if(mask.count(i) > 0) array.push_back(edges[i]);
    else {edges[nremain] = edges[i]; nremain++;}
  }
  edges.resize(nremain);
  return array;
}
