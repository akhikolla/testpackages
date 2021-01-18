#ifndef SAMPLER_HPP_
#define SAMPLER_HPP_

unsigned randmult(const arma::rowvec & vecproba);
void rand_label_edges(PixelGraph & G, VertexMap & label_vertices,
                      EdgeMap * les_poids);

struct SW_filter {
  
  SW_filter() : maxi(0.) {};
  SW_filter(EdgeMap_SW & eg_, EdgeMap_dir & eg_dir_,
            arma::rowvec & maxi_): eg(&eg_), eg_dir(&eg_dir_), maxi(maxi_) {};
  
  virtual ~SW_filter() {};
  
  template <typename Edge>
  bool operator()(const Edge& e) const {
    return ((*eg)[e] >= 0. and (*eg)[e] < maxi[(*eg_dir)[e]]);
  };
  
  EdgeMap_SW * eg;
  EdgeMap_dir * eg_dir;
  arma::rowvec maxi;
};



class Visitor_cc : public boost::default_bfs_visitor {
public:
  Visitor_cc(VertexMap_pot & mylabel_pot, arma::rowvec * myproba) {
    label_pot = &mylabel_pot;
    proba = myproba;
  };
  
  // Si je rencontre un nouveau sommet, je lance ca
  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex v, const Graph &) {
    for (size_t k = 0; k < proba->size(); ++k){
      (*proba)[k] *= exp((*label_pot)[v][k]);
    }
  };
private:
  arma::rowvec * proba;
  VertexMap_pot * label_pot;
};

class Visitor_cc_cond : public boost::default_bfs_visitor {
public:
  Visitor_cc_cond(arma::rowvec & maxi_,
                  VertexMap_pot & mylabel_pot,
                  VertexMap & mylabel_vert,
                   PixelGraph & myGborder,
                  VertexMap & my_vertices_border,
                  EdgeMap & my_edges_border,
                  EdgeMap_dir & my_dir,
                  arma::rowvec * myproba) {
    label_pot = &mylabel_pot;
    label_vertices = &mylabel_vert;
    G_border = &myGborder;
    vertices_border = &my_vertices_border;
    edges_border = &my_edges_border;
    label_dir = &my_dir;
    proba = myproba;
    maxi = maxi_;
  };
  
  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex v, const Graph &) {
    for (size_t k = 0; k < proba->size(); ++k){
      (*proba)[k] *= exp((*label_pot)[v][k]);
    }
    PixelGraph::out_edge_iterator out, out_end;
    for (tie(out, out_end) = out_edges(v, *G_border); out != out_end; ++out){
      if ((*label_vertices)[v] == (*vertices_border)[target(*out, *G_border)]
            and R::runif(0,1) < maxi[(*label_dir)[*out]]){
        (*proba)[(*vertices_border)[target(*out, *G_border)]] *= exp((*edges_border)[*out]);
      }
    }
  };
private:
  arma::rowvec * proba;
  VertexMap_pot * label_pot;
  VertexMap * label_vertices;
  PixelGraph * G_border;
  VertexMap * vertices_border;
  EdgeMap * edges_border;
  EdgeMap_dir * label_dir;
  arma::rowvec maxi;
};


class Visitor_color : public boost::default_bfs_visitor {
public:
  Visitor_color(unsigned new_color_,
                std::set<unsigned> * myset_p,
                VertexMap * img_p_): new_color(new_color_),
                img_p(img_p_), pointeur_set(myset_p){};
  
  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex v, const Graph &) {
    put(*img_p, v, new_color);
    pointeur_set->erase(v);
  };
private:
  unsigned new_color;
  VertexMap * img_p;
  std::set<unsigned> * pointeur_set;
};


#endif /* SAMPLER_HPP_ */
