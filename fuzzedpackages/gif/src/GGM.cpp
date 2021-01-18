// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <RcppEigen.h>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <tuple>
#include <exception>
#include <iostream>
#include <list>

using namespace Rcpp;
using namespace std;

typedef Eigen::Triplet<double> triplet;

int sgn(double x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

double matrix_relative_difference(const Eigen::MatrixXd &S1, const Eigen::MatrixXd &S2) {
  double mean_value = (S1.norm() + S2.norm()) / 2;
  Eigen::MatrixXd mat = S1 - S2;
  return mat.norm() / mean_value;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> bcd(Eigen::MatrixXd &S, Eigen::MatrixXi &act_set,
                                int iter_max = 10, double eps = 1e-3) {
  int p = (int) S.rows();
  int iter_time = 0;
  Eigen::MatrixXd W = S, last_W = S;
  Eigen::VectorXd updated_w12;
  Eigen::VectorXi dense_number_vec = Eigen::VectorXi::Zero(p);
  int s = (int) act_set.rows();
  int nzs = 2 * s;
  Eigen::VectorXd beta_hat(p - 1);
  double temp;
  Eigen::MatrixXi non_zero_index(nzs, 2);
  Eigen::MatrixXi diag_index(p, 2);
  std::vector<triplet> trp;
  std::vector<Eigen::VectorXi> vec_non_zero_index_list(p);

  Eigen::MatrixXi act_set2 = act_set;
  act_set2.col(0).swap(act_set2.col(1));
  non_zero_index << act_set,
                    act_set2;

  for(int i = 0; i < nzs; i++) {
    dense_number_vec(non_zero_index(i, 1))++;
  }

  for(int i = 0; i < p; i++) {
    if(dense_number_vec(i) > 0) {
      int dense_number = dense_number_vec(i);
      Eigen::VectorXi vec_non_zero_index((int) dense_number);
      int t = 0;
      for(int j = 0; j < 2 * s; j++) {
        if(non_zero_index(j, 1) == i) {
          if(non_zero_index(j, 0) == (p - 1)) {
            vec_non_zero_index(t++) = i;
          } else {
            vec_non_zero_index(t++) = non_zero_index(j, 0);
          }
        }
      }
      sort(vec_non_zero_index.data(), vec_non_zero_index.data() + vec_non_zero_index.size());
      vec_non_zero_index_list[i] = vec_non_zero_index;
    }
  }

  while (iter_time < iter_max) {
    for (int i = 0; i < p; ++i) {
      // When there is any element in the active set, the update would work
      // else, we should set the corresponding diagonal element in y_inv 1 / covariance(i,i)
      // auto start = high_resolution_clock::now();

      if (dense_number_vec(i) > 0) {
        int dense_number = dense_number_vec(i);

        Eigen::VectorXi vec_non_zero_index = vec_non_zero_index_list[i];

        Eigen::MatrixXd W11 = Eigen::MatrixXd::Zero(dense_number, dense_number);
        for (int j = 0; j < dense_number; j++) {
          int k = vec_non_zero_index[j];
          if (k != i) {
            W11(j, j) = W(k, k);
          } else {
            W11(j, j) = W(p - 1, p - 1);
          }
        }
        for (int j = 0; j < dense_number; j++) {
          int index_j = vec_non_zero_index[j];
          for(int k = j + 1; k < dense_number; k++) {
            int index_k = vec_non_zero_index[k];
            if (index_j != i && index_k != i) {
              W11(j, k) = W(index_j, index_k);
              W11(k, j) = W11(j, k);
            } else if (index_j == i && index_k != i) {
              W11(j, k) = W(p - 1, index_k);
              W11(k, j) = W11(j, k);
            } else {
              W11(j, k) = W(index_j, p - 1);
              W11(k, j) = W11(j, k);
            }
          }
        }

        Eigen::VectorXd s12 = Eigen::VectorXd::Zero(dense_number);
        for(int j = 0; j < dense_number; j++) {
          int k = vec_non_zero_index[j];
          if(k != i) {
            s12(j) = S(k, i);
          } else {
            s12(j) = S(p - 1, i);
          }
        }

        Eigen::LLT<Eigen::MatrixXd> lltOfW(W11); // compute the Cholesky decomposition of W11
        Eigen::VectorXd beta_star;
        if(lltOfW.info() == Eigen::NumericalIssue)
        {
          throw 1;
        } else {
          // Eigen::VectorXd beta_star = W11.colPivHouseholderQr().solve(s12);
          beta_star = W11.llt().solve(s12);
        }

        updated_w12 = Eigen::VectorXd::Zero(p - 1);
        for(int k = 0; k < dense_number; k++) {
          int index_k = vec_non_zero_index[k];
          for(int j = 0; j < p - 1; j++) {
            if(j != i) {
              if(index_k != i) {
                updated_w12(j) += beta_star(k) * W(j, index_k);
              } else {
                updated_w12(j) += beta_star(k) * W(j, p - 1);
              }
            } else {
              if(index_k != i) {
                updated_w12(i) += beta_star(k) * W(p - 1, index_k);
              } else {
                updated_w12(i) += beta_star(k) * W(p - 1, p - 1);
              }
            }
          }
        }

        for(int j = 0; j < p - 1; j++) {
          if(j != i) {
            W(i, j) = updated_w12(j);
            W(j, i) = W(i, j);
          } else {
            W(p - 1, j) = updated_w12(j);
            W(j, p - 1) = W(p - 1, j);
          }
        }

        if(iter_time == (iter_max - 1)) {
          beta_hat = Eigen::VectorXd::Zero(p - 1);
          for(int j = 0; j < dense_number; j++) {
            beta_hat(vec_non_zero_index[j]) = beta_star(j);
          }
          temp = updated_w12.dot(beta_hat);

          double temp_diag = 1 / (S(i, i) - temp);
          trp.push_back(triplet(i, i, temp_diag));
          for(int j = 0; j < dense_number; j++) {
            int k = vec_non_zero_index[j];
            double temp_entry = - 0.5 * beta_hat(k) * temp_diag;
            if(k != i) {
              trp.push_back(triplet(k, i, temp_entry));
              trp.push_back(triplet(i, k, temp_entry));
            } else {
              trp.push_back(triplet(p - 1, k, temp_entry));
              trp.push_back(triplet(k, p - 1, temp_entry));
            }
          }

          // Here for push_back, it would count twice, so we should divide the corresponding entris by 2,
          // but for y_inv, it would not be the case, the corresponding entries would be updated twice
          // and only the last time would be the result
        }
      } else {
        if(iter_time == (iter_max - 1)) {
          trp.push_back(triplet(i, i, 1 / S(i, i)));
        }
      }

    }

    if (iter_time == (iter_max - 1)){
      break;
    } else if (matrix_relative_difference(W, last_W) < eps) {
      iter_time = iter_max - 1;
    } else {
      last_W = W;
      iter_time++;
    }
  }

  for(int i = 0; i < nzs + p; i++) {
    if(!std::isfinite(trp[i].value())) throw 1;
  }

  Eigen::SparseMatrix<double> Omega(p, p);
  Omega.setFromTriplets(trp.begin(), trp.end());

  return Omega;
}

// Class for an undirected graph
class Graph {
  int V;    // No. of vertices
  list<int> *adj;    // Pointer to an array containing adjacency lists
  bool isCyclicUtil(int v, bool visited[], int parent);
public:
  Graph(int V);   // Constructor
  void addEdge(int v, int w);   // to add an edge to graph
  bool isCyclic();   // returns true if there is a cycle
};

Graph::Graph(int V) {
  this->V = V;
  adj = new list<int>[V];
}

void Graph::addEdge(int v, int w) {
  adj[v].push_back(w); // Add w to v’s list.
  adj[w].push_back(v); // Add v to w’s list.
}

// A recursive function that uses visited[] and parent to detect
// cycle in subgraph reachable from vertex v.
bool Graph::isCyclicUtil(int v, bool visited[], int parent) {
  // Mark the current node as visited
  visited[v] = true;

  // Recur for all the vertices adjacent to this vertex
  list<int>::iterator i;
  for (i = adj[v].begin(); i != adj[v].end(); ++i)
  {
    // If an adjacent is not visited, then recur for that adjacent
    if (!visited[*i])
    {
      if (isCyclicUtil(*i, visited, v))
        return true;
    }

    // If an adjacent is visited and not parent of current vertex,
    // then there is a cycle.
    else if (*i != parent)
      return true;
  }
  return false;
}

// Returns true if the graph contains a cycle, else false.
bool Graph::isCyclic() {
  // Mark all the vertices as not visited and not part of recursion
  // stack
  bool *visited = new bool[V];
  for (int i = 0; i < V; i++)
    visited[i] = false;

  // Call the recursive helper function to detect cycle in different
  // DFS trees
  for (int u = 0; u < V; u++)
    if (!visited[u]) // Don't recur for u if it is already visited
      if (isCyclicUtil(u, visited, -1))
        return true;

      return false;
}

bool is_acyclic(Eigen::MatrixXd &S, Eigen::MatrixXi &act_set) {
  int p = S.cols();
  int s = act_set.rows();
  Graph g(p);
  for(int i = 0; i < s; i++) {
    g.addEdge(act_set(i, 0), act_set(i, 1));
  }
  if(g.isCyclic()) return false; else return true;
}

// [[Rcpp::export]]
List soft_GT(Eigen::MatrixXd &S, double lambda, Eigen::MatrixXi &act_set) {
  // Variable declaration
  int p = S.cols();
  Eigen::MatrixXd sigma_res = Eigen::MatrixXd::Zero(p, p);
  std::vector<triplet> trp;
  Eigen::SparseMatrix<double> S_opt(p, p);
  bool acyclic = is_acyclic(S, act_set);

  // Compute sigma_res
  for(int i = 0; i < p; i++) {
    for(int j = i + 1; j < p; j++) {
      if(fabs(S(i, j)) > lambda) {
        sigma_res(i, j) = S(i, j) - lambda * sgn(S(i, j));
        sigma_res(j, i) = sigma_res(i, j);
      }
    }
  }

  // Set the diagonal entries for S_opt
  for(int i = 0; i < p; i++) {
    double tmp = 0;
    for(int m = 0; m < p; m++) {
      if(m != i && sigma_res(i, m) != 0) {
        tmp += (sigma_res(i, m) * sigma_res(i, m)) / (S(i, i) * S(m, m) - sigma_res(i, m) * sigma_res(i, m));
      }
    }
    trp.push_back(triplet(i, i, (1 + tmp) / S(i, i)));
  }

  // Set the off-diagonal entries for S_opt
  for(int i = 0; i < p; i++) {
    for(int j = i + 1; j < p; j++) {
      if(sigma_res(i ,j) != 0) {
        double tmp = - sigma_res(i, j) / (S(i, i) * S(j, j) - sigma_res(i, j) * sigma_res(i, j));
        trp.push_back(triplet(i, j, tmp));
        trp.push_back(triplet(j, i, tmp));
      }
    }
  }

  S_opt.setFromTriplets(trp.begin(),trp.end());

  return List::create(Named("Omega") = S_opt, Named("is.acyclic") = acyclic);
}
