/*******************************************************************************
uspr.h

Unrooted SPR distance computation and data structures

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

#ifndef INCLUDE_USPR
#define INCLUDE_USPR


#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <memory>
#include <ctime>
#include <cstdlib>
#include "utree.h"
#include "unode.h"
#include "uforest.h"
#include "tbr.h"
#include "uspr_neighbors.h"

//#define DEBUG_USPR 1
#ifdef DEBUG_USPR
	#define debug_uspr(x) x
#else
	#define debug_uspr(x)
#endif

// options
bool USE_TBR_APPROX_ESTIMATE = true;
bool USE_TBR_ESTIMATE = true;
bool USE_REPLUG_ESTIMATE = true;

// classes

typedef enum {REPLUG, TBR, TBR_APPROX, BFS} estimator_t;
string estimator_t_name[] = {"REPLUG", "TBR", "TBR_APPROX", "BFS"};

// distance estimate
class tree_distance {
	public:
	// cost to reach this tree
	int cost;
	// estimate of distance to destination
	int estimate;
	// cost + estimate
	int distance;
	// the tree
	string tree;
	// estimator used
	estimator_t estimator;

	tree_distance(int c, int d, string t, estimator_t e) {
		cost = c;
		estimate = d;
		distance = c + d;
		tree = t;
		estimator = e;
	}

	friend bool operator < (tree_distance a, tree_distance b);
	friend bool operator <= (tree_distance a, tree_distance b);

};

// prefer better estimates when equal
bool operator < (tree_distance a, tree_distance b) {
	if (a.distance == b.distance) {
		if (a.estimate == b.estimate) {
			return a.estimator < b.estimator;
		}
		else {
			return a.estimate < b.estimate;
		}
	}
	return a.distance < b.distance; }
bool operator <= (tree_distance a, tree_distance b) {
	if (a.distance == b.distance) {
		if (a.estimate == b.estimate) {
			return a.estimator <= b.estimator;
		}
		else {
			return a.estimate <= b.estimate;
		}
	}
	return a.distance <= b.distance;
}


// function prototypes
int uspr_distance(uforest &T1, uforest &T2);

// functions
int uspr_distance(uforest &T1_original, uforest &T2_original) {

	uforest T1 = uforest(T1_original);
	uforest T2 = uforest(T2_original);

	// normalize tree order
	T1.normalize_order();
	T2.normalize_order();

	// check if trees are equal
	if (utree(T1).str() == utree(T2).str()) {
		return 0;
	}

	// leaf reduction
	map<string, int> label_map = map<string, int>();
	map<int, string> reverse_label_map = map<int, string>();
	leaf_reduction(&T1, &T2, &label_map, &reverse_label_map);
	T1.normalize_order();
	T2.normalize_order();

	debug_uspr(
		Rcout << "T1R: " << T1 << endl;
		Rcout << "T2R: " << T2 << endl;
	)


	// set of visited trees
	set<string> visited_trees = set<string>();

	// target string
	string target = utree(T2).str();


	// priority queue of trees
	multiset<tree_distance> distance_priority_queue = multiset<tree_distance>();


	// start with the first distance
	visited_trees.insert(T1.str());
	distance_priority_queue.insert(tree_distance(0, 1, utree(T1).str(), BFS));



	// final estimator
	estimator_t final_estimator = BFS;
	if (USE_TBR_APPROX_ESTIMATE) {
		final_estimator = TBR_APPROX;
	}
	if (USE_TBR_ESTIMATE) {
		final_estimator = TBR;
	}
	if (USE_REPLUG_ESTIMATE) {
		final_estimator = REPLUG;
	}

	// explore the next tree
	while (!distance_priority_queue.empty()) {
		multiset<tree_distance>::iterator it = distance_priority_queue.begin();

		// debugging
		debug_uspr(
			Rcout << it->distance << ": " << it->cost << " + " << it->estimate << " using " << estimator_t_name[it->estimator] << endl;
			Rcout << "\t" << it->tree << endl;
		)

		// remove the old entry
		int cost = it->cost;
		string tree = it->tree;
		estimator_t prev_estimator = it->estimator;
		distance_priority_queue.erase(it);

		// build the tree
		uforest T = uforest(tree);
		distances_from_leaf_decorator(T, T.get_smallest_leaf());
		T.normalize_order();

		// check if the distance estimate is final
		if (prev_estimator != final_estimator) {
			// if not, compute the next estimate and insert it into the queue
			int distance = 1;
			if (prev_estimator > TBR_APPROX &&
					USE_TBR_APPROX_ESTIMATE) {
				distance = tbr_high_lower_bound(T, T2);
				distance_priority_queue.insert(tree_distance(cost, distance, tree, TBR_APPROX));
			}
			else if (prev_estimator > TBR &&
					USE_TBR_ESTIMATE) {
				distance = tbr_distance(T, T2);
				distance_priority_queue.insert(tree_distance(cost, distance, tree, TBR));
			}
			else if (prev_estimator > REPLUG &&
					USE_REPLUG_ESTIMATE) {
				distance = replug_distance(T, T2);
				distance_priority_queue.insert(tree_distance(cost, distance, tree, REPLUG));
			}
			continue;
		}

		// if final, get the tree neighborhood
		// TODO: enumerate SPRs
		// note: valid SPRs - move one endpoint to anywhere within its subtree
		// i.e. enumerate edges, then enumerate over each subtree
		// what about duplicates? similar to RSPR, NNIs?
		// for each tree
			// check if it has already been seen (visited_trees)
			// if not then insert it into the queue (initial BFS - cost + 1)
			// TODO: stop immediately if we find T2?

		list<utree> neighbors = get_neighbors(&T, &visited_trees);
		debug_uspr(
			Rcout << "examining " << neighbors.size() << " neighbors" << endl;
		)
		for (utree tree : neighbors) {
			string tree_string = tree.str();
//			Rcout << "neighbor: " << tree_string << endl;
//			Rcout << "target: " << target << endl;
				if (tree_string == target) {
//					Rcout << "returning " << cost+1 << endl;
					debug_uspr(
					  Rcout << "examined " << visited_trees.size() << " trees" << endl;
					)
					return cost+1;
				}
//				else {
//					Rcout << "cond: " << (tree_string == target) << endl;
//				}
				distance_priority_queue.insert(tree_distance(cost+1, 1, tree_string, BFS));
		}

	}

	return -1;
}

#endif
