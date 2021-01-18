/*******************************************************************************
unode.h

Unrooted tree node data structure and functions

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

#ifndef INCLUDE_UNODE
#define INCLUDE_UNODE

//#define DEBUG 1
#ifdef DEBUG
	#define debug(x) x
#else
	#define debug(x)
#endif

#include <list>
#include <sstream>
#include <cstdio>
#include <climits>

using namespace std;

class unode;

class unode {
	private:
	int label;
	list<unode *> neighbors;
	list<unode *> contracted_neighbors;
	int num_neighbors;
	int component;
	bool terminal;
	int distance;
	bool b_protected;
	bool phi;

	public:
	unode() {
		label = -1;
		neighbors = list<unode *>();
		contracted_neighbors = list<unode *>();
		num_neighbors = 0;
		component = -1;
		terminal = false;
		distance = -1;
		b_protected = false;
		phi = false;
	}
	unode(int l) {
		label = l;
		neighbors = list<unode *>();
		contracted_neighbors = list<unode *>();
		num_neighbors = 0;
		component = -1;
		terminal = false;
		distance = -1;
		b_protected = false;
		phi = false;
	}
	unode(const unode &n, bool include_neighbors = true) {
		label = n.label;
		// don't include neighbors when copying as they will be updated later
		if (include_neighbors) {
			neighbors = list<unode *>(n.neighbors);
			contracted_neighbors = list<unode *>(n.contracted_neighbors);
			num_neighbors = n.num_neighbors;
		}
		else {
			neighbors = list<unode *>();
			contracted_neighbors = list<unode *>();
			num_neighbors = 0;
		}
		component = n.component;
		terminal = n.terminal;
		distance = n.distance;
		b_protected = n.b_protected;
		phi = n.phi;
	}
	~unode() {
		neighbors.clear();
		contracted_neighbors.clear();
	}

	void add_neighbor(unode *n) {
		if (num_neighbors > 0 && neighbors.front()->get_distance() > n->get_distance()) {
			neighbors.push_front(n);
		}
		else {
			neighbors.push_back(n);
		}
		num_neighbors++;
	}

	void add_contracted_neighbor(unode *n) {
		contracted_neighbors.push_back(n);
	}

	void add_parent(unode *n) {
		neighbors.push_front(n);
		num_neighbors++;
	}

	bool remove_neighbor(unode *n) {
		list<unode *>::iterator i;
		bool result = false;
		for(i = neighbors.begin(); i != neighbors.end(); i++) {
			if ((*i) == n) {
				result = true;
				break;
			}
		}
		if (result) {
			neighbors.remove(*i);
			num_neighbors--;
		}
		return result;
	}

	bool remove_contracted_neighbor(unode *n) {
		list<unode *>::iterator i;
		bool result = false;
		for(i = contracted_neighbors.begin(); i != contracted_neighbors.end(); i++) {
			if ((*i) == n) {
				result = true;
				break;
			}
		}
		if (result) {
			contracted_neighbors.remove(*i);
		}
		return result;
	}

	bool contract_neighbor(unode *n) {
		bool ret = remove_neighbor(n);
		if (ret) {
			contracted_neighbors.push_back(n);
		}
		return ret;
	}

	string str(map<int, string> *reverse_label_map = NULL) const {
		stringstream ss;
		if (phi) {
			ss << "*";
		}
		else {
			if (reverse_label_map != NULL &&
					reverse_label_map->find(label) != reverse_label_map->end()) {
				ss << (*reverse_label_map)[label];
			}
			else {
				ss << label;
			}
		}
		return ss.str();
	}

	bool operator ==(const unode n) const {
		return this->label == n.label;
	}
	bool operator !=(const unode n) const {
		return this->label != n.label;
	}


	int get_label() const {
		return label;
	}

	const list<unode *> &const_neighbors() const {
		return neighbors;
	}

	const list<unode *> &const_contracted_neighbors() const {
		return contracted_neighbors;
	}

	list<unode *> &get_neighbors() {
		return neighbors;
	}

	list<unode *> &get_contracted_neighbors() {
		return contracted_neighbors;
	}

	int get_num_neighbors() {
		return num_neighbors;
	}

	int get_num_all_neighbors() {
		return num_neighbors + contracted_neighbors.size();
	}

	bool is_leaf() {
		return (num_neighbors == 1);
	}

	void set_component(int c) {
		component = c;
	}

	int get_component() {
		return component;
	}

	void set_phi(bool b) {
		phi = b;
	}

	bool is_phi() {
		return phi;
	}

	void root(int l) {
		unode *p = NULL;
		for (unode *n : get_neighbors()) {
			if (n->get_label() == l) {
				p = n;
			}
			else {
				n->root(get_label());
			}
		}
		if (p != NULL) {
			neighbors.remove(p);
			neighbors.push_front(p);
		}
	}

	void rotate(int l) {
		unode *p = NULL;
		for (unode *n : get_neighbors()) {
			if (n->get_label() == l) {
				p = n;
			}
		}
		if (p != NULL) {
			neighbors.remove(p);
			neighbors.push_front(p);
		}
	}

	unode *get_parent() {
		if (neighbors.empty()) {
			return NULL;
		}
		return neighbors.front();
	}

	bool is_adjacent(unode *a) {
//		Rcout << label << "->is_adjacent(" << a->get_label() << ")" << endl;
		for (unode *n : neighbors) {
			if (n == a) {
//				Rcout << "true" << endl;
				return true;
			}
		}
		for (unode *n : contracted_neighbors) {
			if (n == a) {
//				Rcout << "true" << endl;
				return true;
			}
		}
//		Rcout << "false" << endl;
		return false;
	}

	unode *get_neighbor_not(unode *a) {
		return get_neighbor_not(a, a);
	}

	unode *get_neighbor_not(unode *a, unode *b) {
		list<unode *>::reverse_iterator x;
		for (x = neighbors.rbegin(); x != neighbors.rend(); x++) {
			if (*x != a && *x != b) {
				return *x;
			}
		}
		return NULL;
	}

	void set_terminal(bool t) {
		terminal = t;
	}

	bool get_terminal() {
		return terminal;
	}

	unode *get_sibling() {
		unode *parent = get_parent();
		int count = 0;
		for (unode *x : parent->get_neighbors()) {
			if (count > 0 && x != this) {
				return x;
			}
			count++;
		}
		return parent; // Never reached, but appeases compiler
	}

	void clear_neighbors() {
		neighbors.clear();
		num_neighbors = 0;
	}


	void clear_contracted_neighbors() {
		contracted_neighbors.clear();
	}

	void uncontract_neighbors() {
		for (unode *x: contracted_neighbors) {
			add_neighbor(x);
		}
		clear_contracted_neighbors();
	}

	void uncontract_subtree(unode *last = NULL) {
		for (unode *n : neighbors) {
			if (last == NULL || n != last) {
				n->uncontract_subtree(this);
			}
		}
		for (unode *n : contracted_neighbors) {
			if (last == NULL || n != last) {
				n->uncontract_subtree(this);
			}
		}
		uncontract_neighbors();
	}

	unode *contract_degree_two_subtree(unode *last = NULL) {
		debug(
			Rcout << label << ".contract_degree_two_subtree()" << endl;
		)
		list<unode *> neighbor_copy = list<unode *>(get_neighbors());
		for (unode *n : neighbor_copy) {
			if (last == NULL || n != last) {
				n->contract_degree_two_subtree(this);
			}
		}
		return contract();
	}

	unode *contract() {
		debug(
			Rcout << "n: " << num_neighbors << endl;
			Rcout << "c_n: " << contracted_neighbors.size() << endl;
		)
		if (num_neighbors == 1 && contracted_neighbors.empty()) {
			unode *p = neighbors.front();
			if (p->is_leaf() && this->get_label() < -1) {
				p->remove_neighbor(this);
				this->remove_neighbor(p);
				if (component > -1) {
					p->set_component(component);
				}
				if (is_protected()) {
					p->set_protected(true);
				}
				return p;
			}
		}
//		/*
		else if (num_neighbors == 0 && contracted_neighbors.size() == 2) {
//			uncontract_neighbors();
			unode *p = contracted_neighbors.front();
			unode *c = *(next(contracted_neighbors.begin(), 1));
			debug(
				Rcout << "contracting:" << endl;
				Rcout << p << "\t" << p->get_num_all_neighbors() << endl;
				Rcout << c << "\t" << c->get_num_all_neighbors() << endl;
			)
			if (p->get_num_all_neighbors() < c->get_num_all_neighbors()) {
				unode *temp = p;
				p = c;
				c = temp;
			}
			if (p->get_num_all_neighbors() > 1) {
				clear_contracted_neighbors();
				p->remove_neighbor(this);
				p->remove_contracted_neighbor(this);
				c->remove_neighbor(this);
				c->remove_contracted_neighbor(this);
				c->add_parent(p);
				p->add_contracted_neighbor(c);
				if (p->get_distance() > distance &&
						c->get_distance() > distance) {
					p->set_distance(distance-1);
					c->set_distance(distance);
				}
				else {
					c->set_distance(p->get_distance()+1);
				}
				if (!get_terminal()) {
					p->set_terminal(false);
				}
				else {
					p->set_terminal(true);
				}
				if (component > -1) {
					p->set_component(component);
				}
				if (is_protected()) {
					c->set_protected(true);
				}
				return p;
			}
		}
//		*/
		else if (num_neighbors == 2 && contracted_neighbors.empty()) {
			unode *p = neighbors.front();
			unode *c = *(next(neighbors.begin(), 1));
			debug(
				Rcout << "contracting:" << endl;
				Rcout << p << endl;
				Rcout << c << endl;
			)
			if (!p->is_leaf() ||
						!(p->get_contracted_neighbors().empty()) ||
						!c->is_leaf()) {
				clear_neighbors();
				p->remove_neighbor(this);
				c->remove_neighbor(this);
				c->add_parent(p);
				p->add_neighbor(c);
				if (p->get_distance() > distance &&
						c->get_distance() > distance) {
					p->set_distance(distance-1);
					c->set_distance(distance);
				}
				else {
					c->set_distance(p->get_distance()+1);
				}
				if (!get_terminal()) {
					p->set_terminal(false);
				}
				if (component > -1) {
					p->set_component(component);
				}
				if (is_protected()) {
					c->set_protected(true);
				}
				return p;
			}
		}
		return this;
	}

	void set_distance(int d) {
		distance = d;
	}

	int get_distance() {
		return distance;
	}
	bool is_singleton() {
		if (num_neighbors == 0) {
			return true;
		}
		return false;
	}

	bool is_protected() {
		return b_protected;
	}

	void set_protected(bool b) {
		b_protected = b;
	}

	void get_connected_nodes(list<unode *> &connected_nodes, unode *last = NULL) {
		for (unode *n : neighbors) {
			if (last == NULL || n != last) {
				n->get_connected_nodes(connected_nodes,this);
			}
		}
		for (unode *n : contracted_neighbors) {
			if (last == NULL || n != last) {
				n->get_connected_nodes(connected_nodes,this);
			}
		}
		connected_nodes.push_back(this);
	}

	int normalize_order_hlpr(unode *prev = NULL) {
		// return leaf label
		if (label >= 0 && prev != NULL) {
			return label;
		}
		map<int, unode *> ordered_children = map<int, unode *>();
		int min_descendant = INT_MAX;
		unode *parent = NULL;
		for (unode *n : neighbors) {
			if (n != prev) {
				int n_min_descendant = n->normalize_order_hlpr(this);
				ordered_children.insert(make_pair(n_min_descendant, n));
				if (n_min_descendant < min_descendant) {
					min_descendant = n_min_descendant;
				}
			}
			else {
				parent = n;
			}
		}
		// remove all neighbors
		clear_neighbors();

		// re-add in correct order
		if (parent != NULL) {
			add_neighbor(parent);
		}
		while(!ordered_children.empty()) {
			map<int, unode *>::iterator next = ordered_children.begin();
			add_neighbor(next->second);
			ordered_children.erase(next);
		}

		// contracted neighbors
		ordered_children.clear();
		for (unode *n : contracted_neighbors) {
			int n_min_descendant = n->normalize_order_hlpr(this);
			ordered_children.insert(make_pair(n_min_descendant, n));
			if (n_min_descendant < min_descendant) {
				min_descendant = n_min_descendant;
			}
		}
		// remove all neighbors
		clear_contracted_neighbors();

		// re-add in correct order
		while(!ordered_children.empty()) {
			map<int, unode *>::iterator next = ordered_children.begin();
			add_contracted_neighbor(next->second);
			ordered_children.erase(next);
		}
		return min_descendant;
	}



	// normalize branching order by smallest subtree leaf
	// guaranteed unique if started at the smallest leaf
	void normalize_order() {
		// normalize order
		normalize_order_hlpr();
	}

	unode *find_uncontracted_node() {
		unode *ret = this;
		unode *prev = this;
		while (ret->is_leaf()) {
			unode *next = ret->get_parent();
			if (prev == next) {
				return ret;
			}
			prev = ret;
			ret = next;
		}
		return ret;
	}

};

#endif
