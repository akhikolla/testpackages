/*******************************************************************************
utree.h

Unrooted tree data structure

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

#ifndef INCLUDE_UTREE
#define INCLUDE_UTREE

#include <vector>
#include <iostream>
#include "unode.h"

using namespace std;

class utree;

// options
bool KEEP_LABELS = false;

// prototypes
bool build_utree(utree &t, string &s, map<string, int> *label_map = NULL, map<int, string> *reverse_label_map = NULL);
unsigned int build_utree_helper(utree &t, string &s, int start, unode *parent, bool &valid, map<string, int> *label_map = NULL, map<int, string> *reverse_label_map = NULL);
void find_sibling_pairs_hlpr(utree &t, map<int, int> &sibling_pairs);
map<int, int> distances_from_leaf(utree &T1, int leaf);
void distances_from_leaf_hlpr(utree &T1, map<int, int> &distances, unode *prev, unode *current, int distance);
void distances_from_leaf_decorator(utree &T1, int leaf);
void distances_from_leaf_decorator_hlpr(utree &T1, unode *prev, unode *current, int distance);

class utree {
	protected:
		vector <unode*> internal_nodes;
		vector <unode*> leaves;
		int smallest_leaf;
	public:
		// create the tree
		utree(string &newick, map<string, int> *label_map = NULL, map<int, string> *reverse_label_map = NULL) {
			internal_nodes = vector<unode *>();
			leaves = vector<unode *>();
			build_utree(*this, newick, label_map, reverse_label_map);
		}
		utree(const utree &T) {
			// copy vectors of pointers
			int internal_nodes_size = T.internal_nodes.size();
			int leaves_size = T.leaves.size();
			internal_nodes = vector<unode *>(internal_nodes_size);
			leaves = vector<unode *>(leaves_size);
			smallest_leaf = T.smallest_leaf;
			// create new nodes
			for(int i = 0; i < internal_nodes_size; i++) {
				if (T.internal_nodes[i] != NULL) {
					internal_nodes[i] = new unode(*(T.internal_nodes[i]), false);
				}
			}
			for(int i = 0; i < leaves_size; i++) {
				if (T.leaves[i] != NULL) {
					leaves[i] = new unode(*(T.leaves[i]), false);
				}
			}
			// update neighbor pointers
			for(int i = 0; i < internal_nodes_size; i++) {
				if (internal_nodes[i] != NULL) {
					list<unode *> &old_neighbors = T.internal_nodes[i]->get_neighbors();
					for (unode *u : old_neighbors) {
						internal_nodes[i]->add_neighbor(get_node(u->get_label()));
					}
					list<unode *> &old_contracted_neighbors = T.internal_nodes[i]->get_contracted_neighbors();
					for (unode *u : old_contracted_neighbors) {
						internal_nodes[i]->add_contracted_neighbor(get_node(u->get_label()));
					}
				}
			}
			for(int i = 0; i < leaves_size; i++) {
				if (leaves[i] != NULL) {
					list<unode *> &old_neighbors = T.leaves[i]->get_neighbors();
					for (unode *u : old_neighbors) {
						leaves[i]->add_neighbor(get_node(u->get_label()));
					}
					list<unode *> &old_contracted_neighbors = T.leaves[i]->get_contracted_neighbors();
					for (unode *u : old_contracted_neighbors) {
						leaves[i]->add_contracted_neighbor(get_node(u->get_label()));
					}
				}
			}
		}
		~utree() {
			int end = internal_nodes.size();
			for(int i = 0; i < end; i++) {
				if (internal_nodes[i] != NULL) {
					delete internal_nodes[i];
				}
			}
			end = leaves.size();
			for(int i = 0; i < end; i++) {
				if (leaves[i] != NULL) {
					delete leaves[i];
				}
			}
		}
		utree& operator=(utree T) {
			swap(*this, T);
			return *this;
		}
		friend void swap(utree &first, utree &second) {
			swap(first.internal_nodes, second.internal_nodes);
			swap(first.leaves, second.leaves);
			swap(first.smallest_leaf, second.smallest_leaf);
		}
	friend ostream& operator<<(ostream &os, const utree& t);

	 int add_internal_node() {
		int label = -(internal_nodes.size() + 2);
		internal_nodes.push_back( new unode(label));
		return label;
	}

	 unsigned int add_leaf(unsigned int label) {
	  unsigned int start = leaves.size();
		if (leaves.size() <= label) {
			leaves.resize(label + 1);
		}
		for(unsigned int i = start; i < label; i++) {
			leaves[i] = NULL;
		}
		leaves[label] = new unode(label);
		return label;
	}

	int add_phi_node() {
		int new_node_label = add_leaf(leaves.size());
		unode *new_node = get_node(new_node_label);
		new_node->set_phi(true);
		return new_node_label;
	}

	 unode *get_internal_node(int label) const {
		 return internal_nodes[-(label) - 2];
	 }

	 vector<unode *> &get_internal_nodes() {
		 return internal_nodes;
	 }

	 vector<unode *> &get_leaves() {
		 return leaves;
	 }

	 list<unode *> get_node_list() {
		 list<unode *> L = list<unode *>();
		 for (unode *i: leaves) {
			 if (i != NULL) {
				 L.push_back(i);
			 }
		 }
		 for (unode *i: internal_nodes) {
			 if (i != NULL) {
				 L.push_back(i);
			 }
		 }
			 return L;
	 }

	 unode *get_leaf(int label) const {
		 return leaves[label];
	 }

	 unode *get_node(int label) const {
		 if (label < 0) {
			 return get_internal_node(label);
		 }
		 else {
			 return get_leaf (label);
		 }
	 }

	int num_leaves() const {
		return leaves.size();
	}

	string str(bool print_internal = false, map<int, string> *reverse_label_map = NULL) const{
		stringstream s;
		int start = smallest_leaf;
		if (start == -1) {
			return "empty tree";
		}
		unode *root = leaves[start]->get_neighbors().front();
		str_subtree(s, root, root, print_internal, reverse_label_map);
		return s.str();
	}
	string str(int start, string contracted_sep = ",", bool print_internal = false, map<int, string> *reverse_label_map = NULL) const{
		stringstream s;
		if (start == -1) {
			return "empty tree";
		}
		unode *root = get_node(start);
		str_subtree(s, root, root, contracted_sep, print_internal, reverse_label_map);
		return s.str();
	}

	list<int> find_leaves() {
		list<int> leaf_list = list<int>();
		for (unode *i : leaves) {
			if (i != NULL) {
				leaf_list.push_back(i->get_label());
			}
		}
		return leaf_list;
	}

	map<int, int> find_sibling_pairs() {
		map<int, int> sibling_pairs = map<int, int>();
		find_sibling_pairs_hlpr(*this, sibling_pairs);
		return sibling_pairs;
	}

	void set_smallest_leaf(int l) {
		smallest_leaf = l;
	}

	int get_smallest_leaf() {
		return smallest_leaf;
	}

	void root() {
		root(smallest_leaf);
	}

	void root(int l) {
		unode *n = get_node(l);
		if (n != NULL) {
			n->root(n->get_label());
		}
	}

	void root(unode *n) {
		if (n != NULL) {
			n->root(n->get_label());
		}
	}

	string str_subtree(unode *n, bool print_internal_labels = false, map<int, string> *reverse_label_map = NULL) {
		stringstream ss;
		str_subtree(ss, n, n->get_parent(), print_internal_labels, reverse_label_map);
		return ss.str();
	}

	string str_subtree(unode *n, unode *p, bool print_internal_labels = false, map<int, string> *reverse_label_map = NULL) {
		stringstream ss;
		str_subtree(ss, n, p, print_internal_labels, reverse_label_map);
		return ss.str();
	}

	string str_subtree(unode *n, unode *p, string contracted_sep, bool print_internal_labels = false, map<int, string> *reverse_label_map = NULL) {
		stringstream ss;
		str_subtree(ss, n, p, contracted_sep, print_internal_labels, reverse_label_map);
		return ss.str();
	}

	void str_subtree(stringstream &s, unode *n, unode *prev, bool print_internal_labels = false, map<int, string> *reverse_label_map = NULL) const {
		// only leaf labels
		if (print_internal_labels || n->get_label() >= 0) {
			s << n->str(reverse_label_map);
		}
		list<unode *>::const_iterator i;
		int count = 0;
		bool has_contracted = false;
		for(unode *i : n->const_neighbors()) {
			if (prev == NULL || (*i).get_label() != prev->get_label()) {
				if (count == 0) {
					s << "(";
				}
				else {
					s << ",";
				}
				count++;
				str_subtree(s, i, n, print_internal_labels, reverse_label_map);
			}
		}
		for(unode *i : n->const_contracted_neighbors()) {
			if (prev == NULL || (*i).get_label() != prev->get_label()) {
				if (count == 0) {
					s << "<";
				}
				else {
					s << ",";
				}
				count++;
				has_contracted = true;
				str_subtree(s, i, n, print_internal_labels, reverse_label_map);
			}
		}
		if (has_contracted) {
			s << ">";
		}
		else if (count > 0) {
			s << ")";
		}
	}

	void str_subtree(stringstream &s, unode *n, unode *prev, string contracted_sep, bool print_internal_labels = false, map<int, string> *reverse_label_map = NULL) const {
		// only leaf labels
		if (print_internal_labels || n->get_label() >= 0) {
			s << n->str(reverse_label_map);
		}
		list<unode *>::const_iterator i;

		int count = 0;
		bool has_contracted = false;
		for(unode *i : n->const_neighbors()) {
			if (prev == NULL || (*i).get_label() != prev->get_label()) {
				if (count == 0) {
					s << "(";
				}
				else {
					s << ",";
				}
				count++;
				str_subtree(s, i, n, contracted_sep, print_internal_labels, reverse_label_map);
			}
		}
		for(unode *i : n->const_contracted_neighbors()) {
			if (prev == NULL || (*i).get_label() != prev->get_label()) {
				if (count == 0) {
					s << "<";
				}
				else {
					s << contracted_sep;
				}
				count++;
				has_contracted = true;
				str_subtree(s, i, n, contracted_sep, print_internal_labels, reverse_label_map);
			}
		}
		if (has_contracted) {
			s << ">";
		}
		else if (count > 0) {
			s << ")";
		}
	}

	void str_subtree_with_depths(stringstream &s, unode *n, unode *prev, bool print_internal_labels = false) const {
		// only leaf labels
		if (print_internal_labels || n->get_label() >= 0) {
			s << n->str();
		}
		list<unode *>::const_iterator i;

		int count = 0;
		bool has_contracted = false;
		for(unode *i : n->const_neighbors()) {
			if (prev == NULL || (*i).get_label() != prev->get_label()) {
				if (count == 0) {
					s << "(";
				}
				else {
					s << ",";
				}
				count++;
				str_subtree_with_depths(s, i, n, print_internal_labels);
			}
		}
		for(unode *i : n->const_contracted_neighbors()) {
			if (prev == NULL || (*i).get_label() != prev->get_label()) {
				if (count == 0) {
					s << "<";
				}
				else {
					s << ",";
				}
				count++;
				has_contracted = true;
				str_subtree_with_depths(s, i, n, print_internal_labels);
			}
		}
		if (has_contracted) {
			s << ">";
		}
		else if (count > 0) {
			s << ")";
		}
		s << ":" << n->get_distance();
	}

	void normalize_order() {
		get_node(get_smallest_leaf())->get_parent()->normalize_order();
	}
	void normalize_order(int n) {
		get_node(n)->normalize_order();
	}

	// apply a USPR operation moving (x,y) to (x,yprime) where yprime is adjacent to x, w, and z
	bool uspr(unode *x, unode *y, unode *w, unode *z, unode **yprime = NULL, unode **y1 = NULL, unode **y2 = NULL) {
		// y must have 3 neighbors
		if (y->get_num_neighbors() != 3) {
			return false;
		}
		// remove (x,y)
		x->remove_neighbor(y);
		y->remove_neighbor(x);

		// remove y's first other neighbor
		unode *y1_real = y->get_neighbors().front();
		y->remove_neighbor(y1_real);
		y1_real->remove_neighbor(y);
		if (y1 != NULL) {
			*y1 = y1_real;
		}
		unode *y2_real = y->get_neighbors().front();
		y->remove_neighbor(y2_real);
		y2_real->remove_neighbor(y);
		if (y2 != NULL) {
			*y2 = y2_real;
		}

		// connect y's previous neighbors
		y1_real->add_neighbor(y2_real);
		y2_real->add_neighbor(y1_real);

		//remove (w,z)
		w->remove_neighbor(z);
		z->remove_neighbor(w);

		// add the nodes adjacent to y
		y->add_neighbor(x);
		x->add_neighbor(y);
		y->add_neighbor(w);
		w->add_neighbor(y);
		y->add_neighbor(z);
		z->add_neighbor(y);

		if (yprime != NULL) {
			*yprime = y;
		}

		// cleanup the tree
		// TODO: optional? this is probably slow
//		unode *root = get_node(get_smallest_leaf());
//		distances_from_leaf_decorator(this, root);
//		normalize_order();

		return true;
	}

};

ostream& operator<<(ostream &os, const utree& t) {
	os << t.str() << ";";
	return os;
}

bool build_utree(utree &t, string &s, map<string, int> *label_map, map<int, string> *reverse_label_map) {
	bool valid = true;
	unode dummy = unode(-1);
	build_utree_helper(t, s, 0, &dummy, valid, label_map, reverse_label_map);
	unode *root = dummy.get_parent();
	root->remove_neighbor(&dummy);
	root->contract();

	int end = t.num_leaves();
	int start = -1;
	for(int i = 0; i < end; i++) {
		if (t.get_leaf(i) != NULL) {
			start = i;
			break;
		}
	}
	t.set_smallest_leaf(start);
	return valid;
}

unsigned int build_utree_helper(utree &t, string &s, int start, unode *parent, bool &valid, map<string, int> *label_map, map<int, string> *reverse_label_map) {
	// next special char
	unsigned int loc = s.find_first_of("(,)", start);
	/*if (loc == string::npos) { // MS note: Compiler reports that this is always false
		return loc;
	}*/
	if (s[loc] != '(') {
		// leaf
		while(s[start] == ' ' || s[start] == '\t')
			start++;
		int end = loc;
		while(s[end] == ' ' || s[end] == '\t')
			end--;
		string name = s.substr(start, end - start);
		int label = -1;
		if (label_map != NULL) {
			map<string, int>::iterator m = label_map->find(name);
			if (m != label_map->end()) {
				label = m->second;
			}
			else {
				label = label_map->size();
				if (KEEP_LABELS) {
					label = atoi(name.c_str());
				}
				label_map->insert(make_pair(name, label));
				reverse_label_map->insert(make_pair(label, name));
			}
		}
		else {
			// auto keep labels when no label map is given
			label = atoi(name.c_str());
		}
		unode *new_node = t.get_leaf(t.add_leaf(label));
		parent->add_neighbor(new_node);
		new_node->add_neighbor(parent);
	}
	else {
		// internal node
	int l = t.add_internal_node();
		unode *new_node = t.get_internal_node(l);
		loc = build_utree_helper(t, s, loc + 1, new_node, valid, label_map, reverse_label_map);
		while(s[loc] == ',') {
			loc = build_utree_helper(t, s, loc + 1, new_node, valid, label_map, reverse_label_map);
		}
		if (s[loc] != ')') {
			valid = false;
			return s.size()-1;
		}
	//	if (parent->get_label() != -1) {
			new_node->add_neighbor(parent);
			parent->add_neighbor(new_node);
	//	}
		loc++;
	}
	return loc;
}


void find_sibling_pairs_hlpr(utree &t, map<int, int> &sibling_pairs) {
	for(int l : t.find_leaves()) {
		unode *n = t.get_leaf(l);
		unode *p = n->get_neighbors().front();
		for (unode *u : p->get_neighbors()) {
			int ul = u->get_label();
			if (u->is_leaf() && ul > l) {
				sibling_pairs.insert(make_pair(l,ul));
				sibling_pairs.insert(make_pair(ul,l));
			}
		}
	}
	return;
}

map<int, int> distances_from_leaf(utree &T1, int leaf) {
	map<int, int> distances = map<int, int>();
	unode *node = T1.get_leaf(leaf);
	distances_from_leaf_hlpr(T1, distances, node, node, 0);
	return distances;
}

void distances_from_leaf_hlpr(utree &T1, map<int, int> &distances, unode *prev, unode *current, int distance) {
	distances.insert(make_pair(current->get_label(), distance));
	for(unode *n : current->get_neighbors()) {
		if (n != prev) {
			distances_from_leaf_hlpr(T1, distances, current, n, distance+1);
		}
	}
}

void distances_from_leaf_decorator(utree &T1, int leaf) {
	unode *node = T1.get_leaf(leaf);
	distances_from_leaf_decorator_hlpr(T1, node, node, 0);
}

void distances_from_leaf_decorator_hlpr(utree &T1, unode *prev, unode *current, int distance) {
	current->set_distance(distance);
	for(unode *n : current->get_neighbors()) {
		if (n != prev) {
			distances_from_leaf_decorator_hlpr(T1, current, n, distance+1);
		}
	}
}

#endif
