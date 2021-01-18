#include <GiRaF.hpp>

using namespace std;
using namespace boost;
using namespace arma;


long double Block::recursion() {
  
  unsigned opt = ref.n_rows;
  rowvec z_cond = zeros<rowvec> (opt);
  
  VertexMap_pot pot_on_singletons = get(vertex_potential, G);
  
  if (nb_pixel != 1){
    for (size_t i = 0 ; i < opt ; ++i){
      for (size_t k = 0 ; k < K ; ++k){
        z_cond[i] += factor_fl[k+ref(i,0)] * exp(pot_on_singletons[0][k]);
      }
    }
  } else {
    z_cond.ones(opt);
  }
  
  for (size_t pixel = 1 ; pixel < nb_pixel-h ; ++pixel) {
    
    unsigned row = pixel%(h);
    
    rowvec z_old = z_cond;
    z_cond.zeros(opt);
    
    for (size_t i = 0 ; i < opt ; ++i) {
      
      if (row == 0 and h != 1) { // First Line
        for (size_t k = 0 ; k < K ; ++k){
          
          z_cond[i] += factor_fl[k+ref(i,0)]
          * exp(pot_on_singletons[pixel][k])
          * z_old[k+K*ref(i, 1)];
        }
      } else if (row == h-1) { // Last Line
        for (unsigned k = 0 ; k < K ; ++k){
          z_cond[i] += factor_ll[k+ref(i,0)]
          * exp(pot_on_singletons[pixel][k])
          * z_old[k+K*ref(i, 1)];
        }
      } else {
        for (size_t k = 0 ; k < K ; ++k){
          z_cond[i] += factor[k+ref(i,0)]
          * exp(pot_on_singletons[pixel][k])
          * z_old[k+K*ref(i, 1)];
        }
      }
    }
  }
  
  return sum(z_cond.subvec(0,factor_lc.n_cols-1)%factor_lc) * pow(g, double(num_edges(G)));
}


long double Block::recursion_mem () {
  
  unsigned opt = ref.n_rows;
  rowvec z_cond = zeros<rowvec>(opt);
  z_rec.zeros(max(1, int(nb_pixel-h)), opt);
  
  VertexMap_pot pot_on_singletons = get(vertex_potential, G);
  
  if (nb_pixel != 1){
    for (size_t i = 0 ; i < opt ; ++i){
      for (size_t k = 0 ; k < K ; ++k){
        z_cond[i] += factor_fl[k+ref(i,0)] * exp(pot_on_singletons[0][k]);
      }
    }
    z_rec.row(0) = z_cond;
  } else {
    z_rec.ones();
  }
  
  for (size_t pixel = 1 ; pixel < nb_pixel-h ; ++pixel) {
    
    unsigned row = pixel%(h);
    
    rowvec z_old = z_cond;
    z_cond.zeros();
    
    for (size_t i = 0 ; i < opt ; ++i) {
      
      if (row == 0 and h != 1) { // First Line
        for (size_t k = 0 ; k < K ; ++k){
          
          z_cond[i] += factor_fl[k+ref(i,0)]
          * exp(pot_on_singletons[pixel][k])
          * z_old[k+K*ref(i, 1)];
        }
      } else if (row == h-1) { // Last Line
        for (unsigned k = 0 ; k < K ; ++k){
          z_cond[i] += factor_ll[k+ref(i,0)]
          * exp(pot_on_singletons[pixel][k])
          * z_old[k+K*ref(i, 1)];
        }
      } else {
        for (size_t k = 0 ; k < K ; ++k){
          z_cond[i] += factor[k+ref(i,0)]
          * exp(pot_on_singletons[pixel][k])
          * z_old[k+K*ref(i, 1)];
        }
      }
    }
    
    z_rec.row(pixel) = z_cond;
  }
  
  long double Z = 0.;
  unsigned end_z = z_rec.n_rows -1;
  for (size_t i = 0 ; i < factor_lc.size() ; ++i) {
    z_rec(end_z,i) *= factor_lc[i];
    Z += z_rec(end_z,i);
  }
  
  Z *= pow(g, double(num_edges(G)));
  return Z ;
}


double Factor_Cor(unsigned pixel, Border & border, rowvec *factor) {
  
  // Correct factors when there are borders
  
  VertexMap label_vertex = get(vertex_color, border.G_border);
  EdgeMap label_edge_border = get(edge_weight, border.G_border);
  PixelGraph::vertex_descriptor Vpixel = vertex(pixel, border.G_border);
  
  PixelGraph::out_edge_iterator out, out_end;
  for (tie(out, out_end) = out_edges(Vpixel, border.G_border); out != out_end; ++out){
    
    (*factor)[label_vertex[target(*out, border.G_border)]] *= exp(label_edge_border[*out]);
    
  }
  
  return sum(*factor);
}

long double Block::recursion_cond(Border & border) {
  
  unsigned opt = ref.n_rows;
  rowvec z_cond = zeros<rowvec> (opt);
  rowvec factor_temp = zeros<rowvec>(K);
  
  VertexMap_pot pot_on_singletons = get(vertex_potential, G);
  
  if(nb_pixel != 1){
    for (size_t i = 0 ; i < opt ; ++i) {
      
      for (unsigned k = 0 ; k < K ; ++k){
        factor_temp[k] = factor_fl[k+ref(i,0)]* exp(pot_on_singletons[0][k]);
      }
      
      z_cond[i] = Factor_Cor (0, border, &factor_temp);
    }
  } else {
    // When there is only a single pixel, it is enough to correct the last column
    z_cond.ones(opt);
  }
  
  for (size_t pixel = 1 ; pixel < nb_pixel-h ; ++pixel) {
    
    unsigned row = pixel%h;
    rowvec z_old = z_cond;
    
    for (size_t i = 0 ; i < opt ; ++i) {
      
      if (row == 0 and h != 1) { // First line
        for (size_t k = 0 ; k < K ; ++k){
          factor_temp[k] = factor_fl[k+ref(i,0)]
          * z_old[k+K*ref(i, 1)]
          * exp(pot_on_singletons[pixel][k]);
        }
      }
      else if (row == h-1) { // Last line
        for (size_t k = 0 ; k < K ; ++k){
          factor_temp[k] = factor_ll[k+ref(i,0)]
          * z_old[k+K*ref(i, 1)]
          * exp(pot_on_singletons[pixel][k]);
        }
      }
      else {
        for (size_t k = 0 ; k < K ; ++k){
          factor_temp[k] = factor[k+ref(i,0)]
          * z_old[k+K*ref(i, 1)]
          * exp(pot_on_singletons[pixel][k]);
        }
      }
      z_cond[i] = Factor_Cor (pixel, border, &factor_temp);
    }
  }
  
  return sum(z_cond.subvec(0,factor_lc_cor.n_cols-1)%factor_lc_cor) * pow(g, double(num_edges(G)));
  
}

long double Block::recursion_cond_mem (Border & border) {
  
  unsigned opt = ref.n_rows;
  rowvec z_cond = zeros<rowvec> (opt);
  rowvec factor_temp = zeros<rowvec> (K);
  z_rec.zeros(max(1, int(nb_pixel-h)), opt);
  
  VertexMap_pot pot_on_singletons = get(vertex_potential, G);
  
  if(nb_pixel != 1){
    for (size_t i = 0 ; i < opt ; ++i) {
      for (unsigned k = 0 ; k < K ; ++k){
        factor_temp[k] = factor_fl[k+ref(i,0)]* exp(pot_on_singletons[0][k]);
      }
      z_cond[i] = Factor_Cor (0, border, &factor_temp);
    }
    z_rec.row(0) = z_cond;
  } else {
    // When there is only a single pixel, it is enough to correct the last column
    z_rec.ones();
  }
  
  for (size_t pixel = 1 ; pixel < nb_pixel-h ; ++pixel) {
    
    unsigned row = pixel%h;
    rowvec z_old = z_cond;
    
    for (size_t i = 0 ; i < opt ; ++i) {
      
      if (row == 0 and h != 1) { // First line
        for (size_t k = 0 ; k < K ; ++k){
          factor_temp[k] = factor_fl[k+ref(i,0)]
          * exp(pot_on_singletons[pixel][k])
          * z_old[k+K*ref(i, 1)];
        }
      }
      else if (row == h-1) { // Last line
        for (size_t k = 0 ; k < K ; ++k){
          factor_temp[k] = factor_ll[k+ref(i,0)]
          * exp(pot_on_singletons[pixel][k])
          * z_old[k+K*ref(i, 1)];
        }
      }
      else {
        for (size_t k = 0 ; k < K ; ++k){
          factor_temp[k] = factor[k+ref(i,0)]
          * exp(pot_on_singletons[pixel][k])
          * z_old[k+K*ref(i, 1)];
        }
      }
      z_cond[i] = Factor_Cor (pixel, border, &factor_temp);
    }
    z_rec.row(pixel) = z_cond;
  }
  
  unsigned end_z = z_rec.n_rows -1;
  z_rec.row(end_z).subvec(0,factor_lc_cor.n_cols-1) = z_rec.row(end_z).subvec(0,factor_lc_cor.n_cols-1)%(factor_lc_cor);
  
  return sum(z_rec.row(end_z).subvec(0,factor_lc_cor.n_cols-1))*pow(g, double(num_edges(G)));
  
}



