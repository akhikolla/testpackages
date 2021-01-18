#include <GiRaF.hpp>

using namespace boost;
using namespace arma;
using namespace std;
using namespace Rcpp;


Lattice::Lattice(unsigned height_img, unsigned width_img,
                 unsigned nb_color, unsigned nb_neigh,
                 arma::vec beta):
  h(height_img), w(width_img),
  K(nb_color), nb_nei(nb_neigh), nb_pixel(height_img*width_img),
  param_inter(beta), param_pot(arma::zeros< arma::vec >(nb_color)), tropism(4,0.), g(1.){
  if (nb_nei == 4) {
    init_graph_4(h, w, param_inter, &G);
    
  } else if (nb_nei == 8) {
    init_graph_8(h, w, param_inter, &G);
  } else {
    throw(std::runtime_error("Wrong number of neighbors"));
  }
  set_potential(param_pot, &G);
}


Lattice::Lattice(unsigned height_img, unsigned width_img,
                 unsigned nb_color, unsigned nb_neigh,
                 arma::vec beta, arma::vec gamma): h(height_img), w(width_img),
                 K(nb_color), nb_nei(nb_neigh),
                 nb_pixel(height_img*width_img),
                 param_inter(beta), param_pot(gamma), tropism(4,0.), g(1.){
  if (nb_nei == 4) {
    init_graph_4(h, w, param_inter, &G);
  } else if (nb_nei == 8) {
    init_graph_8(h, w, param_inter, &G);
  } else {
    throw(std::runtime_error("Wrong number of neighbors"));
  }
  set_potential(param_pot, &G);
}

Lattice::Lattice(unsigned height_img, unsigned width_img,
                 unsigned nb_color, unsigned nb_neigh,
                 arma::vec beta, arma::vec gamma, std::vector<unsigned> tropisme_dir):
  h(height_img), w(width_img),
  K(nb_color), nb_nei(nb_neigh),
  nb_pixel(height_img*width_img),
  param_inter(beta), param_pot(gamma), tropism(tropisme_dir), g(1.){
  if (nb_nei == 4) {
    init_graph_4(h, w, param_inter, &G);
  } else if (nb_nei == 8) {
    init_graph_8(h, w, param_inter, &G);
  } else {
    throw(std::runtime_error("Wrong number of neighbors"));
  }
  set_potential(param_pot, &G);
}

Row<unsigned> Lattice::view_vertices(){
  
  Row<unsigned> ans(nb_pixel);
  
  VertexMap label_vertices = get(boost::vertex_color, G);
  
  PixelGraph::vertex_iterator it, it_end;
  for (tie(it, it_end) = vertices(G) ; it != it_end ; ++it){
    ans[*it] = label_vertices[*it];
  }
  
  return ans;
}


Border::Border(unsigned height_block, unsigned width_block, unsigned nb_neigh,
               arma::vec beta): h(height_block), w(width_block), nb_nei(nb_neigh), param_inter(beta){
  if (nb_nei == 4) {
    init_graph_border_4 (h, w, param_inter, &G_border);
  } else if (nb_nei == 8) {
    init_graph_border_8 (h, w, param_inter, &G_border);
  } else {
    throw(std::runtime_error("Wrong number of neighbors"));
  }
}


void Border::set_borders (const vector<unsigned> & top, const vector<unsigned> & left,
                          const vector<unsigned> & bottom, const vector<unsigned> & right,
                          const vector<unsigned> & corner) {
  
  /* Initialize vertices labels on the borders */
  
  VertexMap label_vertex = get(vertex_color, G_border);
  
  for (unsigned i = 0 ; i < w ; ++i){
    put(label_vertex, h*w + i, top[w-1-i]);
    put(label_vertex, h*w + w + h +i, bottom[i]);
  }
  
  for (unsigned i = 0 ; i < h ; ++i){
    put(label_vertex, h*w + w + i, left[i]);
    put(label_vertex, h*w + 2*w + h +i, right[h-1-i]);
  }
  
  for (unsigned i = 0 ; i < 4 ; ++i){
    put(label_vertex, h*w + 2*w + 2*h +i, corner[i]);
  }
}

void Block::initFactor () {
  
  Mat<unsigned> dico_factor;
  dictionnary_factor(K, nb_nei, &dico_factor);
  
  switch(nb_nei){
  case 4:
    dictionnary(h, K, nb_nei, &ref);
    break;
  case 8:
    dictionnary(h+1, K, nb_nei, &ref);
    factor.resize(K*K*K*K*K, 1.);
    factor_fl.resize(K*K*K*K*K, 1.);
    factor_ll.resize(K*K*K*K*K, 1.);
    break;
  }
  
  PixelGraph G_, G_fl, G_ll, G_lc;
  init_graph_factor(h, nb_nei, param_inter, &G_, &G_fl, &G_ll, &G_lc);
  
  Model_Factor(dico_factor, G_, &factor, g);
  Model_Factor(dico_factor, G_fl, &factor_fl, g);
  Model_Factor(dico_factor, G_ll, &factor_ll, g);
  
  VertexMap_pot pot_on_singletons = get(vertex_potential, G);
  Model_Factor_lc(h, w, K, g, pot_on_singletons, G_lc, &factor_lc);
  
  if (h == 1){
    factor_fl = factor_ll;
  }
}

void Block::correctFactor(Border & border) {
  
  // Same than Factor_Cor but for the last column
  
  VertexMap label_vertex_border = get(vertex_color, border.G_border);
  EdgeMap label_edge_border = get(edge_weight, border.G_border);
  
  factor_lc_cor = factor_lc;
  
  for (size_t i = 0 ; i < factor_lc.n_cols ; ++i) {
    
    vector<unsigned> data = config_base_K(i, h, K);
    
    for (size_t j = 0 ; j < h ; ++j) {
      
      double modif = 0.;
      
      PixelGraph::vertex_descriptor Vpixel= vertex(h * (w-1) + j, border.G_border);
      PixelGraph::out_edge_iterator out, out_end;
      
      for (tie(out, out_end) = out_edges(Vpixel, border.G_border); out != out_end; ++out){
        modif += label_edge_border[*out]
        * Model_Pot(data[j], label_vertex_border, out, border.G_border);
      }
      factor_lc_cor[i] *= exp(modif);
    }
  }
}


void init_graph_4 (unsigned h, unsigned w, const vec & param_inter, PixelGraph *G) {
  
  add_vertex(0, *G); // For lattices of size 1
  
  EdgeMap label_edge = get(edge_weight, *G);
  EdgeMap_dir dir = get(edge_weight2, *G);
  
  unsigned pixel = 0 ;
  for(unsigned col = 0; col < w -1; ++col, ++pixel) {
    for (unsigned row = 0; row < h -1; ++row, ++pixel) {
      add_edge(pixel, pixel+1, *G);
      put(dir, edge(pixel, pixel+1, *G).first, 0);
      put(label_edge, edge(pixel, pixel+1, *G).first, param_inter[0]);
      
      add_edge(pixel, pixel+h, *G);
      put(dir, edge(pixel, pixel+h, *G).first, 1);
      put(label_edge, edge(pixel, pixel+h, *G).first, param_inter[1]);
    }
    // Last line
    add_edge(pixel, pixel + h, *G);
    put(dir, edge(pixel, pixel+h, *G).first, 1);
    put(label_edge, edge(pixel, pixel+h, *G).first, param_inter[1]);
  }
  // Last column
  for (unsigned row = 0; row < h -1; ++row, ++pixel){
    add_edge(pixel, pixel + 1, *G);
    put(dir, edge(pixel, pixel+1, *G).first, 0);
    put(label_edge, edge(pixel, pixel+1, *G).first, param_inter[0]);
  }
}

void init_graph_8 (unsigned h, unsigned w, const vec & param_inter, PixelGraph *G) {
  
  init_graph_4 (h, w, param_inter, G);
  
  EdgeMap label_edge = get(edge_weight, *G);
  EdgeMap_dir dir = get(edge_weight2, *G);
  
  unsigned pixel = 0 ;
  
  for (unsigned col = 0 ; col < w -1 ; ++col, ++pixel) {
    ++pixel;
    
    for (unsigned row = 1 ; row < h -1 ; ++row, ++pixel) {
      
      add_edge(pixel, pixel+h -1, *G);
      put(dir, edge(pixel, pixel+h-1, *G).first, 2);
      put(label_edge, edge(pixel, pixel+h-1, *G).first, param_inter[2]);
      
      add_edge(pixel, pixel+h +1, *G);
      put(dir, edge(pixel, pixel+h +1, *G).first, 3);
      put(label_edge, edge(pixel, pixel+h +1, *G).first, param_inter[3]);
    }
    
    // Last line
    add_edge(pixel, pixel+h -1,*G);
    put(dir, edge(pixel, pixel+h -1, *G).first, 2);
    put(label_edge, edge(pixel, pixel+h -1, *G).first, param_inter[2]);
    
    
    // First line
    add_edge(pixel-(h -1), pixel+2, *G);
    put(dir, edge(pixel-(h -1), pixel+2, *G).first, 3);
    put(label_edge, edge(pixel-(h -1), pixel+2, *G).first, param_inter[3]);
  }
}

void init_graph_border_4 (unsigned h, unsigned w, const vec & param_inter, PixelGraph *G) {
  
  // Graph of the borders (does not include the edges within the block)
  
  EdgeMap label_edge = get(edge_weight, *G);
  EdgeMap_dir label_dir = get(edge_weight2, *G);
  
  for (unsigned i = 0 ; i < w ; ++i){
    unsigned cible_top = h*w + w -1 -i;
    unsigned cible_bottom = h*w + w + h +i;
    add_edge(i*h, cible_top, *G);
    add_edge((i+1)*h-1, cible_bottom, *G);
    
    put(label_dir, edge(i*h, cible_top, *G).first, 0);
    put(label_dir, edge((i+1)*h-1, cible_bottom, *G).first, 0);
    put(label_edge, edge(i*h, cible_top, *G).first, param_inter[0]);
    put(label_edge, edge((i+1)*h-1, cible_bottom, *G).first, param_inter[0]);
  }
  
  for (unsigned i = 0 ; i < h ; ++i){
    unsigned cible_left = h*w + w + i;
    unsigned cible_right = h*w + 2*w + 2*h-1 -i;
    add_edge(i, cible_left, *G);
    add_edge(h*(w-1)+i, cible_right, *G);
    
    put(label_dir, edge(i, cible_left, *G).first, 1);
    put(label_dir, edge(h*(w-1)+i, cible_right, *G).first, 1);
    put(label_edge, edge(i, cible_left, *G).first, param_inter[1]);
    put(label_edge, edge(h*(w-1)+i, cible_right, *G).first, param_inter[1]);
  }
  
  // Ajout des corners
  add_vertex(h*w+2*w+2*h, *G);
  add_vertex(h*w+2*w+2*h+1, *G);
  add_vertex(h*w+2*w+2*h+2, *G);
  add_vertex(h*w+2*w+2*h+3, *G);
}


void init_graph_border_8 (unsigned h, unsigned w, const vec & param_inter, PixelGraph *G) {
  
  init_graph_border_4 (h, w, param_inter, G);
  
  EdgeMap label_edge = get(edge_weight, *G);
  EdgeMap_dir label_dir = get(edge_weight2, *G);
  
  for (unsigned i = 0 ; i < w-1 ; ++i){
    unsigned cible_top = h*w + w -1 -i;
    unsigned cible_bottom = h*w + w + h +i;
    add_edge(i*h, cible_top-1, *G);
    add_edge((i+1)*h-1, cible_bottom+1, *G);
    
    put(label_dir, edge(i*h, cible_top-1, *G).first, 2);
    put(label_dir, edge((i+1)*h-1, cible_bottom+1, *G).first, 3);
    
    put(label_edge, edge(i*h, cible_top-1, *G).first, param_inter[2]);
    put(label_edge, edge((i+1)*h-1, cible_bottom+1, *G).first, param_inter[3]);
  }
  
  for (unsigned i = 1 ; i < w ; ++i){
    unsigned cible_top = h*w + w -i;
    unsigned cible_bottom = h*w + w + h-1 +i;
    add_edge(i*h, cible_top, *G);
    add_edge((i+1)*h-1, cible_bottom, *G);
    
    put(label_dir, edge(i*h, cible_top, *G).first, 3);
    put(label_dir, edge((i+1)*h-1, cible_bottom, *G).first, 2);
    
    put(label_edge, edge(i*h, cible_top, *G).first, param_inter[3]);
    put(label_edge, edge((i+1)*h-1, cible_bottom, *G).first, param_inter[2]);
  }
  
  for (unsigned i = 0 ; i < h-1 ; ++i){
    unsigned cible_left = h*w + w + i;
    unsigned cible_right = h*w + 2*w + 2*h-1 -i;
    add_edge(i, cible_left, *G);
    add_edge(i, cible_left+1, *G);
    add_edge(h*(w-1)+i, cible_right, *G);
    add_edge(h*(w-1)+i, cible_right-1, *G);
    
    put(label_dir, edge(i, cible_left, *G).first, 1);
    put(label_dir, edge(i, cible_left+1, *G).first, 2);
    put(label_dir, edge(h*(w-1)+i, cible_right, *G).first, 1);
    put(label_dir, edge(h*(w-1)+i, cible_right-1, *G).first, 3);
    
    put(label_edge, edge(i, cible_left, *G).first, param_inter[1]);
    put(label_edge, edge(i, cible_left+1, *G).first, param_inter[2]);
    put(label_edge, edge(h*(w-1)+i, cible_right, *G).first, param_inter[1]);
    put(label_edge, edge(h*(w-1)+i, cible_right-1, *G).first, param_inter[3]);
  }
  
  for (unsigned i = 1 ; i < h ; ++i){
    unsigned cible_left = h*w + w-1 + i;
    unsigned cible_right = h*w + 2*w + 2*h -i;
    add_edge(i, cible_left, *G);
    add_edge(h*(w-1)+i, cible_right, *G);
    
    put(label_dir, edge(i, cible_left, *G).first, 3);
    put(label_dir, edge(h*(w-1)+i, cible_right, *G).first, 2);
    
    put(label_edge, edge(i, cible_left, *G).first, param_inter[3]);
    put(label_edge, edge(h*(w-1)+i, cible_right, *G).first, param_inter[2]);
  }
  
  // Corners
  add_edge(0, h*w+2*w+2*h, *G);
  add_edge(h-1, h*w+2*w+2*h+1, *G);
  add_edge(h*w-1, h*w+2*w+2*h+2, *G);
  add_edge(h*(w-1), h*w+2*w+2*h+3, *G);
  
  put(label_dir, edge(0, h*w+2*w+2*h, *G).first, 3);
  put(label_dir, edge(h-1, h*w+2*w+2*h+1, *G).first, 2);
  put(label_dir, edge(h*w-1, h*w+2*w+2*h+2, *G).first, 3);
  put(label_dir, edge(h*(w-1), h*w+2*w+2*h+3, *G).first, 2);
  
  put(label_edge, edge(0, h*w+2*w+2*h, *G).first, param_inter[3]);
  put(label_edge, edge(h-1, h*w+2*w+2*h+1, *G).first, param_inter[2]);
  put(label_edge, edge(h*w-1, h*w+2*w+2*h+2, *G).first, param_inter[3]);
  put(label_edge, edge(h*(w-1), h*w+2*w+2*h+3, *G).first, param_inter[2]);
}

void set_potential(const vec & potential, PixelGraph *G) {
  
  VertexMap_pot label_pot = get(boost::vertex_potential, *G);
  
  PixelGraph::vertex_iterator it, it_end;
  for (tie(it, it_end) = vertices(*G) ; it != it_end ; ++it){
    put(label_pot, *it , potential);
  }
}

void rand_label_edges(PixelGraph *G){
  
  EdgeMap les_poids = get(edge_weight, *G);
  const VertexMap label_vertices = get(boost::vertex_color, *G);
  
  PixelGraph::edge_iterator it, it_end;
  for (tie(it, it_end) = edges(*G); it != it_end; ++it) {
    if (label_vertices[source(*it, *G)] == label_vertices[target(*it, *G)]){
      put(les_poids, *it, R::runif(0,1));
    } else {
      les_poids[*it] = -1. ;
    }
  }
}

void dictionnary (unsigned length, unsigned K, unsigned nb_neigh, Mat<unsigned> *ref) {
  
  /* Store quantities to link a word (column of the image) with the corresponding factor */
  
  unsigned int opt = pow(double(K), int(length));
  ref->set_size(opt, 2);
  
  unsigned int ante_p = pow(double(K), int(length-3));
  unsigned int av_der = pow(double(K), int(length-2));
  unsigned int fin = pow(double(K), int(length-1));
  
  switch(nb_neigh) {
  case 4:
    for (size_t i = 0 ; i < opt ; ++i) {
      (*ref)(i,0) = (i%K)*K + (i/fin)*K*K; // Links for factors
      (*ref)(i,1) = i%fin ; // Word minus its last letter
    }
    break;
  case 8:
    for (size_t i = 0; i < opt ; ++i) {
      (*ref)(i,0) = (i%K)*K + ((i%av_der)/ante_p)*K*K
      + ((i%fin)/av_der)*K*K*K + (i/fin)*K*K*K*K ;
      (*ref)(i,1) = i%fin ;
    }
    break;
  }
}


void dictionnary_factor (unsigned K, unsigned nb_neigh, Mat<unsigned> *dico_factor) {
  
  /* Save useful configurations for factors computation */
  
  unsigned sec = K*K;
  unsigned ter = K*K*K;
  unsigned fin = sec;
  
  switch(nb_neigh){
  case 4:
    dico_factor->set_size(fin*K, 3);
    for (size_t i = 0 ; i < dico_factor->n_rows ; ++i) {
      (*dico_factor)(i,0)= i%K ; // First letter
      (*dico_factor)(i,1)= (i%fin)/K ; // 2nd letter
      (*dico_factor)(i,2)= i/fin ; // 3rd letter
    }
    break;
  case 8:
    fin = ter*K;
    dico_factor->set_size(fin*K, 5);
    for (size_t i = 0 ; i < dico_factor->n_rows ; ++i) {
      (*dico_factor)(i,0) = i%K ;
      (*dico_factor)(i,1) = (i%sec)/K ;
      (*dico_factor)(i,2) = (i%ter)/sec ;
      (*dico_factor)(i,3) = (i%fin)/ter ;
      (*dico_factor)(i,4) = i/fin;
    }
    break;
  }
}

vector<unsigned> config_base_K (unsigned z, unsigned length, unsigned K) {
  
  /* Number z written in base K */
  
  vector<unsigned> ans(length, 0);
  unsigned q = z;
  for(size_t i = 0 ; i < length ; ++i){
    ans[i] = q%K ;
    q = q/K;
  }
  return ans ;
}

void Model_Factor(const Mat<unsigned> & dico_factor,
                  PixelGraph & G, vector<double> *factor,
                  double g) {
  
  /* Computation of all the possible factors */
  
  const EdgeMap label_edge = get(boost::edge_weight, G);
  
  // On calcule les facteurs correspondants a un mot
  for (size_t i_opt = 0 ; i_opt < factor->size() ; ++i_opt){
    
    vector<unsigned> word(dico_factor.n_cols, 0);
    for (size_t i = 0 ; i < word.size() ; ++i){
      word[i] = dico_factor(i_opt,i);
    }
    
    PixelGraph::edge_iterator it, it_end;
    for (tie(it, it_end) = edges(G); it != it_end; ++it) {
      (*factor)[i_opt] *= exp(-log(g) + label_edge[*it] * Model_Pot(word, it, G));
    }
  }
}

void Model_Factor_lc(unsigned h, unsigned w, unsigned K, double g,
                     const VertexMap_pot & pot_on_singletons,
                     PixelGraph & G, rowvec *factor_lc) {
  
  /* Computation of all the possible factors for the last column */
  
  const EdgeMap label_edge = get(boost::edge_weight, G);
  
  for (size_t i = 0 ; i < factor_lc->size() ; ++i){
    vector<unsigned> data = config_base_K(i, h, K);
    
    PixelGraph::edge_iterator it, it_end;
    for (tie(it, it_end) = edges(G); it != it_end; ++it) {
      (*factor_lc)[i] *= exp(-log(g)
                               + label_edge[*it] * Model_Pot(data, it, G));
    }
    
    PixelGraph::vertex_iterator vit, vit_end;
    for (tie(vit, vit_end) = vertices(G); vit != vit_end; ++vit) {
      (*factor_lc)[i] *= exp(pot_on_singletons[*vit+h*(w-1)][data[*vit]] );
    }
  }
}

void init_graph_factor_4 (const vec & param, bool last_line, PixelGraph *G) {
  
  /* 1st order dependency structure of a node for factors computation */
  
  EdgeMap label_edges = get(edge_weight, *G);
  EdgeMap_dir dir = get(edge_weight2, *G);
  
  if (!last_line){
    add_edge(0, 1, *G);
    put(dir, edge(0, 1, *G).first, 0);
    put(label_edges, edge(0, 1, *G).first, param[0]);
  }
  add_edge(0, 2, *G);
  put(dir, edge(0, 2, *G).first, 1);
  put(label_edges, edge(0, 2, *G).first, param[1]);
}

void init_graph_factor_8 (const vec & param, bool first_line, bool last_line, PixelGraph *G) {
  
  /* 2nd order dependency structure of a node for factors computation */
  
  EdgeMap label_edges = get(edge_weight, *G);
  EdgeMap_dir dir = get(edge_weight2, *G);
  
  if (!first_line){
    add_edge(0, 2, *G);
    put(dir, edge(0, 2, *G).first, 2);
    put(label_edges, edge(0, 2, *G).first, param[2]);
  }
  if (!last_line){
    add_edge(0, 1, *G);
    put(dir, edge(0, 1, *G).first, 0);
    put(label_edges, edge(0, 1, *G).first, param[0]);
    
    add_edge(0, 4, *G);
    put(dir, edge(0, 4, *G).first, 3);
    put(label_edges, edge(0, 4, *G).first, param[3]);
  }
  add_edge(0, 3, *G);
  put(dir, edge(0, 3, *G).first, 1);
  put(label_edges, edge(0, 3, *G).first, param[1]);
}


void init_graph_lc (const vec & param, unsigned h, PixelGraph *G) {
  
  /* Dependency structure for the last column */
  
  EdgeMap label_edges = get(edge_weight, *G);
  EdgeMap_dir dir = get(edge_weight2, *G);
  
  add_vertex(0, *G);
  for (unsigned i = 0 ; i < h-1 ; ++i){
    add_edge(i, i+1,  *G);
    put(dir, edge(i, i+1,  *G).first, 0);
    put(label_edges, edge(i, i+1,  *G).first, param[0]);
  }
}


void init_graph_factor (unsigned h, unsigned nb_nei, const vec & param_inter,
                        PixelGraph *G, PixelGraph *G_fl,
                        PixelGraph *G_ll, PixelGraph *G_lc) {
  
  if (nb_nei == 4){
    init_graph_factor_4(param_inter, false, G);
    init_graph_factor_4(param_inter, false, G_fl);
    init_graph_factor_4(param_inter, true, G_ll);
  } else {
    init_graph_factor_8(param_inter, false, false, G);
    init_graph_factor_8(param_inter, true, false, G_fl);
    init_graph_factor_8(param_inter, false, true, G_ll);
  }
  
  init_graph_lc(param_inter, h, G_lc);
}
