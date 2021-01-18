#include "jti_types.h"
#include "set_ops.h"

struct hypergraph {

  // Member vars
  std::vector<int> nodes;                                // clique nodes
  std::vector<edge_pair> edges;                          // weighted edges (w, {u,v}) of the clique nodes
  std::vector<bool> visited;                             // for dfs_*
  std::unordered_map<int, int>  parent;                  // for dfs_*
  std::unordered_map<int, std::vector<int>> clique_tree; // product of kruskal
  bool has_cycle;                                        // does clique_tree contain cycles
  arma::Mat<int> rooted_junction_tree;

  // Constructor
  hypergraph(Rcpp::List cliques) {

    // Init nodes
    std::vector<int> ns(cliques.size());
    std::iota(ns.begin(), ns.end(), 0);
    nodes = ns;
    
    has_cycle = false;
    rooted_junction_tree = arma::zeros<arma::Mat<int>>(nodes.size(), nodes.size());

    // Init visited
    visited.resize(nodes.size());
    std::fill(visited.begin(), visited.end(), false);

    // Sort the elements so int_set_intersect dosent need to do it every single time
    int nc = cliques.size();
    for (int i = 0; i < nc; i++) {
      std::vector<int> cs = cliques[i];
      std::sort(cs.begin(), cs.end());
      cliques[i] = cs;
    }

    // Init weighted edges as vectors of pair<w, pair<v,u>>
    for (int i = 0; i < nc; i++) {
      for (int j = i; j < nc; j++) {
	if (j == i) continue;
	std::vector<int> ci  = cliques[i];
	std::vector<int> cj  = cliques[j];
	std::vector<int> sij = int_set_intersect(ci, cj);
	if (!sij.empty()) {
	  edge_pair new_edge = std::make_pair(-sij.size(), std::make_pair(i,j));
	  edges.push_back(new_edge);
	}
      }
    }

    // Sort the edges in decending order to prepare for kruskal
    std::sort(edges.begin(), edges.end()); 
  }

  // Member funcs
  void reset_dfs_variable();
  void dfs_detect_cycle(int src);
  void dfs_root_clique_tree(int src); 
  void kruskal();

  // For debugging
  // void show_clique_tree();
  // void show_edges();
  
};

// void hypergraph::show_edges() {
//   for (auto & e : edges) {
//     auto p = e.second;
//     auto w = e.first;
//     std::cout << "(" << p.first << ", " << p.second << ") : " << w << "\n";
//   }
// }

// void hypergraph::show_clique_tree() {
//   for (auto & e : clique_tree) {
//     std::cout << e.first << " : ";
//     for (auto & u : e.second) {
//       std::cout << u << ", ";
//     }
//     std::cout << "\n";
//   }
// }

void hypergraph::reset_dfs_variable() {
  parent.clear();
  std::fill(visited.begin(), visited.end(), false);
  has_cycle = false;
}

void hypergraph::dfs_detect_cycle(int src) {
  visited[src] = true;
  std::vector<int> src_neighbors = clique_tree[src];
  for (auto & clique_tree_node : src_neighbors) {
    if (!visited[clique_tree_node]) {
      parent[clique_tree_node] = src;
      dfs_detect_cycle(clique_tree_node);
    } else if (parent[src] != clique_tree_node) {
      has_cycle = true;
      return;
    }
  }
}


void hypergraph::dfs_root_clique_tree(int src) {
  visited[src] = true;
  std::vector<int> src_neighbors = clique_tree[src];
  for (auto & clique_tree_node : src_neighbors) {
    if (!visited[clique_tree_node]) {
      rooted_junction_tree(clique_tree_node, src) = 1;
      dfs_root_clique_tree(clique_tree_node);
    }
  }
}


void add_edge(std::unordered_map<int, std::vector<int>>& clique_tree, edge_pair edge) {
  clique_tree[edge.second.first].push_back(edge.second.second);
  clique_tree[edge.second.second].push_back(edge.second.first);
}

void hypergraph::kruskal() {
  int n_edges_added = 0;
  int n_edges = edges.size();
  for (int i = 0; i < n_edges; i++) {
    auto edge = edges[i];
    auto node = edge.second.first;
    auto tmp_clique_tree = clique_tree;
    add_edge(clique_tree, edge);
    dfs_detect_cycle(node);
    if (has_cycle) {
      clique_tree = tmp_clique_tree;
    } else { 
     n_edges_added += 1;
    }
    reset_dfs_variable();
    if (n_edges_added == nodes.size() - 1) return;
  }
}


void root_as_largest_clique(Rcpp::List & cliques, int & root) {
  int biggest = 0;
  for (int i = 0; i < cliques.size(); i++) {
    std::vector<int> ci = cliques[i];
    if (ci.size() > biggest) {
      biggest = ci.size();
      root = i + 1;
    }
  }  
}
  


// [[Rcpp::export]]
Rcpp::List rooted_junction_tree(Rcpp::List cliques, int root = 0) {

  // select the biggest clique if no root is given
  if (!root) root_as_largest_clique(cliques, root);

  if (root > cliques.size() || root < 1) Rcpp::stop("root must be in {1, 2, ..., #cliques}");
  
  hypergraph g(cliques); // - 1 to compensate for the R side
  g.kruskal();
  g.dfs_root_clique_tree(root - 1);
  return Rcpp::List::create(_["collect"]     = g.rooted_junction_tree,
			    _["distribute"]  = g.rooted_junction_tree.t(),
			    _["clique_root"] = root);
}
