#ifndef LATTICE_HPP_
#define LATTICE_HPP_

class Border {
public:
  Border(unsigned height_block, unsigned width_block, unsigned nb_neigh,
         arma::vec beta);
  
  virtual ~Border() {};
  
  void set_borders(const std::vector<unsigned> & top, const std::vector<unsigned> & left,
                   const std::vector<unsigned> & bottom, const std::vector<unsigned> & right,
                   const std::vector<unsigned> & corner);
  
public:
  unsigned h, w, nb_nei;
  arma::vec param_inter;
  PixelGraph G_border;
};

class Lattice {
public:
  Lattice(unsigned height_img, unsigned width_img,
          unsigned nb_color, unsigned nb_neigh,
          arma::vec beta);
  
  Lattice(unsigned height_img, unsigned width_img,
          unsigned nb_color, unsigned nb_neigh,
          arma::vec beta, arma::vec gamma);
  
  Lattice(unsigned height_img, unsigned width_img,
          unsigned nb_color, unsigned nb_neigh,
          arma::vec beta, arma::vec gamma, std::vector<unsigned> tropisme_dir);
  
  virtual ~Lattice() {};
  
  arma::Row<unsigned> view_vertices();
  
  void GibbsSampler(unsigned niter, bool random, bool initialize_z = true);
  void GibbsSamplerCond(unsigned niter, Border & border,
                        bool random, bool initialize_z = true);
  void SWSampler(unsigned niter, bool initialize_z = true);
  void SWSamplerCond(unsigned niter, Border & border, bool initialize_z = true);
  
public:
  unsigned h, w, K, nb_nei, nb_pixel;
  arma::vec param_inter, param_pot;
  std::vector<unsigned> tropism;
  double g;
  PixelGraph G;
};


class Block : public Lattice {
public:
  Block(unsigned height_block, unsigned width_block,
        unsigned nb_color, unsigned nb_neigh,
        arma::vec beta): Lattice(height_block, width_block,
        nb_color, nb_neigh, beta),
        nb_letters(nb_neigh/2 + 1),
        factor(K*K*K, 1.),
        factor_fl(K*K*K, 1.),
        factor_ll(K*K*K, 1.),
        factor_lc(arma::ones< arma::rowvec >(unsigned(pow(double(K), height_block)))),
        factor_lc_cor(factor_lc),
        ref(arma::zeros<arma::Mat<unsigned> >(0)),
        z_rec(arma::zeros<arma::mat>(0)){ };
  
  Block(unsigned height_block, unsigned width_block,
        unsigned nb_color, unsigned nb_neigh,
        arma::vec beta, arma::vec gamma): Lattice(height_block, width_block,
        nb_color, nb_neigh, beta, gamma),
        nb_letters(nb_neigh/2 + 1),
        factor(K*K*K, 1.),
        factor_fl(K*K*K, 1.),
        factor_ll(K*K*K, 1.),
        factor_lc(arma::ones< arma::rowvec >(unsigned(pow(double(K), height_block)))),
        factor_lc_cor(factor_lc),
        ref(arma::zeros<arma::Mat<unsigned> >(0)),
        z_rec(arma::zeros<arma::mat>(0)){};
  
  Block(unsigned height_block, unsigned width_block,
        unsigned nb_color, unsigned nb_neigh,
        arma::vec beta, arma::vec gamma, std::vector<unsigned> tropisme_dir):
    Lattice(height_block, width_block, nb_color, nb_neigh,
            beta, gamma, tropisme_dir),
            nb_letters(nb_neigh/2 + 1),
            factor(K*K*K, 1.),
            factor_fl(K*K*K, 1.),
            factor_ll(K*K*K, 1.),
            factor_lc(arma::ones< arma::rowvec >(unsigned(pow(double(K), height_block)))),
            factor_lc_cor(factor_lc),
            ref(arma::zeros<arma::Mat<unsigned> >(0)),
            z_rec(arma::zeros<arma::mat>(0)){};
  
  virtual ~Block() {};
  
  
  void initFactor();
  void correctFactor(Border & border);
  
  long double recursion();
  long double recursion_mem();
  
  long double recursion_cond(Border & border);
  long double recursion_cond_mem(Border & border);
  
  
  void exact_sample (const long double & Z);
  void exact_sample_cond (const long double & Z, Border & border);
  
public:
  unsigned nb_letters;
  std::vector<double> factor, factor_fl, factor_ll;
  arma::rowvec factor_lc, factor_lc_cor;
  arma::Mat<unsigned> ref;
  arma::mat z_rec;
};
/***************************************/





/***************************************/
void init_graph_4 (unsigned h, unsigned w, const arma::vec & param_inter, PixelGraph *G);
void init_graph_8 (unsigned h, unsigned w, const arma::vec & param_inter, PixelGraph *G);
void init_graph_border_4 (unsigned h, unsigned w, const arma::vec & param_inter, PixelGraph *G);
void init_graph_border_8 (unsigned h, unsigned w, const arma::vec & param_inter, PixelGraph *G);

//void set_vertices(const std::vector<unsigned> & word, PixelGraph *G);
void set_potential(const arma::vec & potential, PixelGraph *G);
//void set_label_edges(const arma::vec & param, PixelGraph *G);


inline double Model_Pot(const std::vector<unsigned> & word,
                        const PixelGraph::edge_iterator & it,
                        const PixelGraph & G){
  return word[source(*it, G)] == word[target(*it, G)];
};

inline double Model_Pot(unsigned k,
                        const VertexMap & label_vertex,
                        const PixelGraph::out_edge_iterator & out,
                        const PixelGraph & g_border){
  return k == label_vertex[target(*out, g_border)];
}

void dictionnary (unsigned length, unsigned K, unsigned nb_neigh, arma::Mat<unsigned> *ref);
void dictionnary_factor (unsigned K, unsigned nb_neigh, arma::Mat<unsigned> *dico_factor);
std::vector<unsigned> config_base_K(unsigned z, unsigned length, unsigned K);

void Model_Factor(const arma::Mat<unsigned> & dico_factor,
                  PixelGraph & G,
                  std::vector<double> *factor,
                  double g);
void Model_Factor_lc(unsigned h_block, unsigned w_block, unsigned K, double g,
                     const VertexMap_pot & pot_on_singletons,
                     PixelGraph & G, arma::rowvec *factor_lc);

void init_graph_factor_4 (const arma::vec & param, bool last_line,
                          PixelGraph *G);
void init_graph_factor_8 (const arma::vec & param, bool first_line, bool last_line,
                          PixelGraph *G);
void init_graph_lc (const arma::vec & param, unsigned h_block,
                    PixelGraph *G);
void init_graph_factor (unsigned h_block, unsigned nb_nei, const arma::vec & param_inter,
                        PixelGraph *G, PixelGraph *G_fl,
                        PixelGraph *G_ll, PixelGraph *G_lc);
/***************************************/



#endif /* LATTICE_HPP_ */
