#include <GiRaF.hpp>

using namespace boost;
using namespace arma;
using namespace std;
using namespace Rcpp;

#include <RcppArmadilloExtensions/sample.h>


unsigned randmult(const rowvec & vecproba) {
  
  unsigned int couleur = 0;
  double sumproba = vecproba[0];
  double p = R::runif(0,1);
  
  while (p > sumproba) {
    ++couleur;
    sumproba += vecproba[couleur];
  }
  return couleur;
}

void Block::exact_sample (const long double & Z) {
  
  VertexMap label_vertices = get(vertex_color, G);
  const VertexMap_pot label_pot = get(vertex_potential, G);
  
  rowvec vecproba_init = pow(g, double(num_edges(G))) / Z * (z_rec).row((z_rec).n_rows - 1);
  
  unsigned col_opt = randmult(vecproba_init);
  vector<unsigned> col_config = config_base_K(col_opt, h, K);
  
  unsigned lc = nb_pixel-h;
  for (size_t i = 0 ; i < h ; ++i){
    put(label_vertices, lc +i, col_config[i]);
  }
  
  
  if (lc > 0){
    unsigned f_pixel = lc - 1;
    unsigned given_col = col_opt;
    rowvec q_cond = conv_to< rowvec >::from(factor_ll);
    q_cond = factor_lc[given_col]*q_cond;
    
    // Probability to swap color
    rowvec vec_proba = zeros<rowvec>(K);
    for (size_t k = 0 ; k < K ; ++k){
      vec_proba[k] = z_rec(f_pixel-1,
                           k + K * ref(given_col, 1))/z_rec(f_pixel, given_col);
      vec_proba[k] *= q_cond[k + ref(given_col, 0)] * exp(label_pot[f_pixel][k]);
    }
    
    put(label_vertices, f_pixel, randmult(vec_proba));
    
    for (size_t pixel = f_pixel - 1 ; pixel > 0 ; --pixel){
      
      given_col = label_vertices[pixel+1] + K*ref(given_col, 1); // Update given_col
      
      unsigned row = pixel%h;
      
      if (row == 0 and h != 1){
        q_cond = conv_to< rowvec >::from(factor_fl);
      } else if (row == h -1){
        q_cond = conv_to< rowvec >::from(factor_ll);
      } else {
        q_cond = conv_to< rowvec >::from(factor);
      }
      
      
      for (size_t k = 0 ; k < K ; ++k){
        vec_proba[k] = z_rec(pixel-1,
                             k + K * ref(given_col, 1))/z_rec(pixel, given_col);
        vec_proba[k] *= q_cond[k + ref(given_col, 0)] * exp(label_pot[pixel][k]);
      }
      put(label_vertices, pixel, randmult(vec_proba));
    }
    
    // Update first pixel
    given_col = label_vertices[1] + K*ref(given_col, 1);
    
    for (unsigned k = 0 ; k < K ; ++k){
      vec_proba[k] = factor_fl[k+ref(given_col, 0)]/z_rec(0, given_col) * exp(label_pot[0][k]);
    }
    put(label_vertices, 0, randmult(vec_proba));
  }
}


void Block::exact_sample_cond (const long double & Z, Border & border) {
  
  VertexMap label_vertices = get(vertex_color, G);
  const VertexMap_pot label_pot = get(vertex_potential, G);
  const EdgeMap label_edges_border = get(edge_weight, border.G_border);
  const VertexMap & label_vertices_border = get(vertex_color, border.G_border);
  
  rowvec vecproba_init = (z_rec).row((z_rec).n_rows - 1) *pow(g, double(num_edges(G))) / Z;
  
  unsigned col_opt = randmult(vecproba_init);
  vector<unsigned> col_config = config_base_K(col_opt, h, K);
  
  unsigned lc = nb_pixel-h;
  for (size_t i = 0 ; i < h ; ++i){
    put(label_vertices, lc +i, col_config[i]);
  }
  
  if (lc > 0){
    
    unsigned f_pixel = lc - 1;
    
    unsigned given_col = col_opt;
    
    rowvec q_cond = conv_to< rowvec >::from(factor_ll);
    q_cond = factor_lc_cor[given_col]*q_cond;
    
    
    // Initialisation du vecteur de proba
    rowvec vec_proba = zeros<rowvec>(K);
    for (size_t k = 0 ; k < K ; ++k){
      vec_proba[k] = z_rec(f_pixel-1,
                           k + K * ref(given_col, 1))/z_rec(f_pixel, given_col);
      vec_proba[k] *= q_cond[k + ref(given_col, 0)] * exp(label_pot[f_pixel][k]);
    }
    
    PixelGraph::out_edge_iterator out, out_end;
    for (tie(out, out_end) = out_edges(f_pixel, border.G_border); out != out_end; ++out){
      vec_proba[label_vertices_border[target(*out, border.G_border)]] *= exp(label_edges_border[*out]);
    }
    
    vec_proba = vec_proba/sum(vec_proba);
    
    put(label_vertices, f_pixel, randmult(vec_proba));
    
    for (size_t pixel = f_pixel - 1 ; pixel > 0 ; --pixel){
      
      given_col = label_vertices[pixel+1] + K*ref(given_col, 1);
      
      unsigned row = pixel%h;
      
      if (row == 0 and h != 1){
        q_cond = conv_to< rowvec >::from(factor_fl);
      } else if (row == h -1){
        q_cond = conv_to< rowvec >::from(factor_ll);
      } else {
        q_cond = conv_to< rowvec >::from(factor);
      }
      
      for (size_t k = 0 ; k < K ; ++k){
        vec_proba[k] = z_rec(pixel-1,
                             k + K * ref(given_col, 1))/z_rec(pixel, given_col);
        vec_proba[k] *= q_cond[k + ref(given_col, 0)] * exp(label_pot[pixel][k]);
      }
      
      for (tie(out, out_end) = out_edges(pixel, border.G_border); out != out_end; ++out){
        vec_proba[label_vertices_border[target(*out, border.G_border)]] *= exp(label_edges_border[*out]);
      }
      
      vec_proba = vec_proba/sum(vec_proba);
      
      put(label_vertices, pixel, randmult(vec_proba));
    }
    
    given_col = label_vertices[1] + K*ref(given_col, 1);
    
    
    for (unsigned k = 0 ; k < K ; ++k){
      vec_proba[k] = factor_fl[k+ref(given_col,0)]/z_rec(0, given_col) * exp(label_pot[0][k]);
    }
    
    for (tie(out, out_end) = out_edges(0, border.G_border); out != out_end; ++out){
      vec_proba[label_vertices_border[target(*out, border.G_border)]] *= exp(label_edges_border[*out]);
    }
    
    vec_proba = vec_proba/sum(vec_proba);
    
    put(label_vertices, 0, randmult(vec_proba));
  }
}


void Lattice::GibbsSampler (unsigned niter, bool random, bool initialize_z) {
  
  VertexMap label_vertices = get(vertex_color, G);
  const VertexMap_pot label_pot = get(vertex_potential, G);
  const EdgeMap label_edges = get(edge_weight, G);
  
  if (initialize_z) {
    PixelGraph::vertex_iterator it, it_end;
    for (tie(it, it_end) = vertices(G) ; it != it_end ; ++it){
      put(label_vertices, *it, floor(K*R::runif(0,1)));
    }
  }
  
  unsigned nb_pixel = num_vertices(G);
  Row<unsigned> seq = linspace< Row<unsigned> >(0, nb_pixel-1, nb_pixel);
  rowvec vecproba = zeros<rowvec> (K);
  
  for (unsigned iter = 0; iter < niter; ++iter) {
    
    if (random == true){
      seq = Rcpp::RcppArmadillo::sample(seq, nb_pixel, false);
    }
    
    for (unsigned i = 0 ; i < nb_pixel ; ++i) {
      
      
      for (size_t k = 0 ; k < vecproba.size() ; ++k) {
        vecproba[k] = label_pot[seq[i]][k];
      }
      
      PixelGraph::out_edge_iterator out, out_end;
      for (tie(out, out_end) = out_edges(seq[i], G); out != out_end; ++out){
        vecproba[label_vertices[target(*out, G)]] += label_edges[*out];
      }
      
      for (size_t k = 0 ; k < vecproba.size() ; ++k) {
        vecproba[k] = exp(vecproba[k]);
      }
      
      vecproba = vecproba/sum(vecproba);
      put(label_vertices, seq[i], randmult(vecproba));
    }
  }
}




void Lattice::GibbsSamplerCond (unsigned niter, Border & border,
                                bool random, bool initialize_z) {
  
  
  const VertexMap label_vertices_border = get(vertex_color, border.G_border);
  const EdgeMap label_edges_border = get(edge_weight, border.G_border);
  
  VertexMap label_vertices = get(vertex_color, G);
  const VertexMap_pot label_pot = get(vertex_potential, G);
  const EdgeMap label_edges = get(edge_weight, G);
  
  
  if (initialize_z) {
    PixelGraph::vertex_iterator it, it_end;
    for (tie(it, it_end) = vertices(G) ; it != it_end ; ++it){
      put(label_vertices, *it, floor(K*R::runif(0,1)));
    }
  }
  
  unsigned nb_pixel = num_vertices(G);
  Row<unsigned> seq = linspace< Row<unsigned> >(0, nb_pixel-1, nb_pixel);
  
  rowvec vecproba = zeros<rowvec>(K);
  
  for (unsigned iter = 0; iter < niter; ++iter) {
    
    if (random == true){
      seq = Rcpp::RcppArmadillo::sample(seq, nb_pixel, false);
    }
    
    for (size_t i = 0 ; i < nb_pixel ; ++i) {
      
      for (size_t k = 0 ; k < vecproba.size() ; ++k) {
        vecproba[k] = label_pot[seq[i]][k];
      }
      
      PixelGraph::out_edge_iterator out, out_end;
      for (tie(out, out_end) = out_edges(seq[i], G); out != out_end; ++out){
        vecproba[label_vertices[target(*out, G)]] += label_edges[*out];
      }
      
      for (tie(out, out_end) = out_edges(seq[i], border.G_border); out != out_end; ++out){
        vecproba[label_vertices_border[target(*out, border.G_border)]] += label_edges_border[*out];
      }
      
      for (size_t k = 0 ; k < vecproba.size() ; ++k) {
        vecproba[k] = exp(vecproba[k]);
      }
      
      vecproba = vecproba/sum(vecproba);
      put(label_vertices, seq[i], randmult(vecproba));
    }
  }
}


void rand_label_edges(PixelGraph & G, VertexMap & label_vertices, EdgeMap_SW * les_poids){
  
  PixelGraph::edge_iterator it, it_end;
  for (tie(it, it_end) = edges(G); it != it_end; ++it) {
    if (label_vertices[source(*it, G)] == label_vertices[target(*it, G)]){
      put(*les_poids, *it, R::runif(0,1));
    } else {
      put(*les_poids, *it, -1.);
    }
  }
}

void Lattice::SWSampler(unsigned niter, bool initialize_z) {
  
  VertexMap label_vertices = get(vertex_color, G);
  EdgeMap_SW les_poids = get(edge_update, G);
  EdgeMap_dir label_dir = get(edge_weight2, G);
  VertexMap_pot label_pot = get(vertex_potential, G);
  
  if (initialize_z) {
    PixelGraph::vertex_iterator it, it_end;
    for (tie(it, it_end) = vertices(G) ; it != it_end ; ++it){
      put(label_vertices, *it, floor(K*R::runif(0,1)));
    }
  }
  
  rand_label_edges(G, label_vertices, &les_poids);
  
  rowvec proba = zeros<rowvec>(param_inter.size());
  for (size_t i = 0 ; i < param_inter.size() ; ++i){
    proba[i] = 1. - exp(-param_inter[i]);
  }
  
  /*Filtered graph to access the connected component with all edges
  * label lower than proba*/
  
  SW_filter my_filter_sw(les_poids, label_dir, proba);
  
  for (unsigned iter = 0; iter < niter; ++iter) {
    
    filtered_graph< PixelGraph, SW_filter> g_SW (G, my_filter_sw);
    
    set<unsigned> index_pixel;
    
    PixelGraph::vertex_iterator it, it_end;
    for (tie(it, it_end) = vertices(G) ; it != it_end ; ++it){
      index_pixel.insert(*it);
    }
    
    set<unsigned>::iterator it_s;
    while (!index_pixel.empty()){
      it_s = index_pixel.begin();
      rowvec vec_proba = ones<rowvec>(K);
      Visitor_cc vis_cc(label_pot, &vec_proba);
      breadth_first_search(g_SW, *it_s, visitor(vis_cc));
      
      vec_proba = vec_proba/sum(vec_proba);
      unsigned new_color = randmult(vec_proba);
      
      Visitor_color vis_color(new_color, &index_pixel, &label_vertices);
      breadth_first_search(g_SW, *it_s, visitor(vis_color));
    }
    
    rand_label_edges(G, label_vertices, &les_poids);
  }
}


void Lattice::SWSamplerCond(unsigned niter, Border & border, bool initialize_z) {
  
  VertexMap label_vertices = get(vertex_color, G);
  EdgeMap_SW les_poids = get(edge_update, G);
  EdgeMap_dir label_dir = get(edge_weight2, G);
  VertexMap_pot label_pot = get(vertex_potential, G);
  
  if (initialize_z) {
    PixelGraph::vertex_iterator it, it_end;
    for (tie(it, it_end) = vertices(G) ; it != it_end ; ++it){
      put(label_vertices, *it, floor(K*R::runif(0,1)));
    }
  }
  
  rand_label_edges(G, label_vertices, &les_poids);
  
  rowvec proba = zeros<rowvec>(param_inter.size());
  for (size_t i = 0 ; i < param_inter.size() ; ++i){
    proba[i] = 1. - exp(-param_inter[i]);
  }
  
  SW_filter my_filter_sw(les_poids, label_dir, proba);
  
  VertexMap label_vertices_border = get(vertex_color, border.G_border);
  EdgeMap label_edges_border = get(edge_weight, border.G_border);
  EdgeMap_dir label_dir_border = get(edge_weight2, border.G_border);
  
  
  for (unsigned iter = 0; iter < niter; ++iter) {
    
    filtered_graph< PixelGraph, SW_filter> g_SW (G, my_filter_sw);
    
    set<unsigned> index_pixel;
    
    PixelGraph::vertex_iterator it, it_end;
    for (tie(it, it_end) = vertices(G) ; it != it_end ; ++it){
      index_pixel.insert(*it);
    }
    
    set<unsigned>::iterator it_s;
    while (!index_pixel.empty()){
      it_s = index_pixel.begin();
      rowvec vec_proba = ones<rowvec>(K);
      Visitor_cc_cond vis_cc(proba, label_pot, label_vertices,
                             border.G_border, label_vertices_border,
                             label_edges_border, label_dir_border,
                             &vec_proba);
      breadth_first_search(g_SW, *it_s, visitor(vis_cc));
      
      vec_proba = vec_proba/sum(vec_proba);
      unsigned new_color = randmult(vec_proba);
      
      Visitor_color vis_color(new_color, &index_pixel, &label_vertices);
      breadth_first_search(g_SW, *it_s, visitor(vis_color));
    }
    
    rand_label_edges(G, label_vertices, &les_poids);
  }
}
