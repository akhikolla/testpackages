/*******************************************************************************
uspr_neighbors.h

Unrooted SPR neighbor computation and data structures

Copyright 2018 Chris Whidden
cwhidden@fredhutch.org
https://github.com/cwhidden/uspr
May 1, 2018
Version 1.0.1

This file is part of uspr.

uspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

uspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with uspr.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************/


#ifndef INCLUDE_USPR_NEIGHBORS
#define INCLUDE_USPR_NEIGHBORS

// INCLUDES
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <list>
#include <time.h>

#include "utree.h"

using namespace std;

// FUNCTIONS

list<utree> get_neighbors(utree *T, set<string> *known_trees = NULL);
void get_neighbors(utree *T, unode *prev, unode *current, list<utree> &neighbors, set<string> *known_trees = NULL);
void get_neighbors(utree *T, unode *x, unode *y, unode *prev, unode *current, list<utree> &neighbors, set<string> *known_trees = NULL);
void add_neighbor(utree *T, unode *x, unode *y, unode *w, unode *z, list<utree> &neighbors, set<string> *known_trees = NULL);


list<utree> get_neighbors(utree *T, set<string> *known_trees) {
	list<utree> neighbors = list<utree>();
	unode *root = T->get_node(T->get_smallest_leaf());
	get_neighbors(T, NULL, root, neighbors, known_trees);
	return neighbors;
}

// enumerate the source edges
void get_neighbors(utree *T, unode *prev, unode *current, list<utree> &neighbors, set<string> *known_trees) {
	// continue enumerating choices of the first edge
	list<unode *> c_neighbors = current->get_neighbors();
	for (unode *next : c_neighbors) {
		if (next != prev) {
			get_neighbors(T, current, next, neighbors, known_trees);
		}
	}
	// try moving both sides of the edge (prev, current)
	if (prev != NULL) {
		get_neighbors(T, prev, current, prev, current, neighbors, known_trees);
		get_neighbors(T, current, prev, current, prev, neighbors, known_trees);
	}
}

// enumerate the target edges
void get_neighbors(utree *T, unode *x, unode *y, unode *prev, unode *current, list<utree> &neighbors, set<string> *known_trees) {
	// continue enumerating choices of the second edge
	// copy the neighbor list as it may change
	list<unode *> c_neighbors = current->get_neighbors();
	for (unode *next : c_neighbors) {
		if (next != prev) {
			get_neighbors(T, x, y, current, next, neighbors, known_trees);
		}
	}
	// test the spr move (T, x, y, prev, current)
	if (prev != NULL) {
		add_neighbor(T, x, y, prev, current, neighbors, known_trees);
	}
}

void add_neighbor(utree *T, unode *x, unode *y, unode *w, unode *z, list<utree> &neighbors, set<string> *known_trees) {
	// check for duplicate SPR moves
	if (x == y ||
			y == w ||
			y == z ) {
		return;
	}
	// TODO: other duplicates? probably NNIs, like with rooted SPR?

	if (w == y->get_parent()->get_parent() &&
			z == y->get_parent()) {
		return;
	}
	if (z == y->get_parent()->get_parent() &&
			w == y->get_parent()) {
		return;
	}
	if (z->get_parent() == y &&
			w->get_parent() == z) {
		return;
	}
	if (w->get_parent() == y &&
			z->get_parent() == w) {
		return;
	}

	// node info so the uspr can be reversed
	unode *yprime = NULL;
	unode *y1 = NULL;
	unode *y2 = NULL;
	// apply the spr
/*	Rcout << endl;
	Rcout << "T: " << T->str(true) << endl;
	Rcout << "\tx: " << x->get_label() << endl;
	Rcout << "\ty: " << y->get_label() << endl;
	Rcout << "\tw: " << w->get_label() << endl;
	Rcout << "\tz: " << z->get_label() << endl;
*/
	T->uspr(x, y, w, z, &yprime, &y1, &y2);
	// normalize the tree
	distances_from_leaf_decorator(*T, T->get_smallest_leaf());
	T->normalize_order();
	// print the tree
	//	Rcout << "neighbor: " << T->str() << endl;
	string tree_string = T->str();
	bool add_tree = true;
	if (known_trees != NULL) {
		if (known_trees->find(tree_string) == known_trees->end()) {
			known_trees->insert(tree_string);
		}
		else {
			add_tree = false;
		}
	}
	if (add_tree) {
		neighbors.push_back(utree(*T));
	}


	// revert the SPR
	T->uspr(x, yprime, y1, y2);
	distances_from_leaf_decorator(*T, T->get_smallest_leaf());
	T->normalize_order();
//	Rcout << "T: " << T->str() << endl;
//	Rcout << endl;
	return;
}
#endif
