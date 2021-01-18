/*******************************************************************************
uforest.h

Forest data structure and functions

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

#ifndef INCLUDE_UFOREST
#define INCLUDE_UFOREST

#include <iostream>
#include "unode.h"
#include "utree.h"

class uforest: public utree {
	private:
		vector<unode *> components;

	public:
		uforest(string &newick, map<string, int> *label_map, map<int, string> *reverse_label_map) : utree(newick, label_map, reverse_label_map) {
			components = vector<unode *>();
			if (leaves.size() > 0) {
				components.push_back(leaves[smallest_leaf]);
				leaves[smallest_leaf]->set_component(0);
			}
		}
		uforest(string &newick) : utree(newick) {
			components = vector<unode *>();
			if (leaves.size() > 0) {
				components.push_back(leaves[smallest_leaf]);
				leaves[smallest_leaf]->set_component(0);
			}
		}
		uforest(const uforest &F) : utree(F) {
			// copy vector of pointers
			int components_size = F.components.size();
			components = vector<unode *>(components_size);
			// update with the new nodes
			for(int i = 0; i < components_size; i++) {
//				Rcout << "components[" << i << "]" << endl;
				components[i] = get_node(F.components[i]->get_label());
//				Rcout << "\t" << F.components[i] << ", " << components[i] << endl;
			}
		}
		uforest& operator=(uforest F) {
			swap(*this, F);
			return *this;
		}
		friend void swap(uforest &first, uforest &second) {
			swap(static_cast<utree&>(first), static_cast<utree&>(second));
			swap(first.components, second.components);
		}
		string str(bool print_internal = false, map<int, string> *reverse_label_map = NULL) const {

			stringstream ss;
			for(int i = 0; i < (int) components.size(); i++) {
				if (i > 0) {
					ss << " ";
				}
				unode *root = components[i];
				if (root->get_component() != i) {
					ss << "@";
				}
				if (root->get_label() > -1) {
					if (root->is_leaf()) {
						root = root->get_neighbors().front();
					}
					else if (!(root->get_contracted_neighbors().empty())){
						root = root->get_contracted_neighbors().front();
					}
				}
				str_subtree(ss, root, root, print_internal, reverse_label_map);
				ss << ";";
			}
			return ss.str();
		}
		string str_with_depths(bool print_internal = false) const {
			stringstream ss;
			for(unsigned int i = 0; i < components.size(); i++) {
				if (i > 0) {
					ss << " ";
				}
				unode *root = components[i];
				if (root->get_label() > -1) {
					if (root->is_leaf()) {
						root = root->get_neighbors().front();
					}
					else if (!(root->get_contracted_neighbors().empty())){
						root = root->get_contracted_neighbors().front();
					}
				}
				str_subtree_with_depths(ss, root, root, print_internal);
				ss << ";";
			}
			return ss.str();
		}
		friend ostream& operator<<(ostream &os, const uforest& f);

		pair<int,int> cut_edge(int x, int y) {
			unode *X, *Y;
			X = get_node(x);
			Y = get_node(y);
			//Rcout << "d_X: " << X->get_distance() << endl;
			//Rcout << "d_Y: " << Y->get_distance() << endl;
			bool swapped = false;
			if (Y->get_distance() > X->get_distance()) {
				//Rcout << "AHH!" << endl;
				X = get_node(y);
				Y = get_node(x);
				swapped = true;
			}
			bool cut_x = X->remove_neighbor(Y);
			bool cut_y = Y->remove_neighbor(X);

			if (!cut_x || !cut_y) {
				return make_pair(-1,-1);
			}

			unode *Xprime = X->contract();
			unode *Yprime = Y->contract();


			if (Xprime->get_component() > -1) {
				//Rcout << "boo" << endl;
				add_component(Yprime);
				update_component(Xprime->get_component(), Xprime);
				if (Yprime->get_component() > -1) {
					update_component(Yprime->get_component(), Yprime);
				}
			}
			else {
				//Rcout << "urns" << endl;
				add_component(Xprime);
				if (Yprime->get_component() > -1) {
					update_component(Yprime->get_component(), Yprime);
				}
			}
			if (swapped) {
				return make_pair(Yprime->get_label(),Xprime->get_label());
			}
			else {
				return make_pair(Xprime->get_label(),Yprime->get_label());
			}
		}

		void update_component(int c, int l) {
			components[c] = get_node(l);
		}

		void update_component(int c, unode *n) {
			components[c] = n;
		}

		void add_component(unode *C) {
			C->set_component(components.size());
			components.push_back(C);
		}

		int num_components() {
			return components.size();
		}

		void uncontract() {
			for (unode *c : components) {
				unode *root = c;
				if (root->get_label() > -1) {
					if (root->is_leaf()) {
						root = root->get_neighbors().front();
					}
					else if (!(root->get_contracted_neighbors().empty())){
						root = root->get_contracted_neighbors().front();
					}
				}
				root->uncontract_subtree();
			}
		}
		void contract_degree_two() {
			for(unsigned int i = 0; i < components.size(); i++) {
				unode *c = components[i];
				unode *result = c->contract_degree_two_subtree();
				if (result != c) {
					components[i] = result;
				}
			}
		}

	list<unode *> get_alive_nodes() {
		list<unode *> alive = list<unode *>();
		for (unode *c : components) {
			c->get_connected_nodes(alive);
		}
		return alive;
	}

	vector<unode *> get_components() {
		return components;
	}

	void normalize_order() {
		for (unode *c : components) {
			c->normalize_order();
		}
	}
};

ostream& operator<<(ostream &os, const uforest& f) {
	os << f.str();
	return os;
}

#endif
