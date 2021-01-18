#include "fwpopsim.h"
#include "fwsim_types.h"
#include <iostream>
#include <fstream>
#include <string>

//#include "gperftools/profiler.h"


using namespace Rcpp;


std::ostream& operator<<(std::ostream &strm, const SimulatedGenealogy& simres) {
  strm << "Size = " << simres.m_pop_size << std::endl;
  return strm;
}

/*
RcppExport SEXP start_profiler(SEXP str) {
  ProfilerStart(as<const char*>(str));
  return R_NilValue;
}

RcppExport SEXP stop_profiler() {
  ProfilerStop();
  return R_NilValue;
}
*/



  
/*
http://eli.thegreenplace.net/2009/11/23/visualizing-binary-trees-with-graphviz
*/

/*
void genealogy_to_dot_draw_node(Individual* node, std::ostringstream& stream)
{
    if (!node) {
      return;
    } 
    
    static int nullcount = 0;
    
    stream << "  " << node->get_id() << " [label=\"" << node->get_label() << "\"];" << std::endl;

    std::vector<Individual*> children = node->get_children();
    
    for (auto &child: children) { 
      stream << "    " << node->get_id() << " -> " << child->get_id() << ";" << std::endl;
      genealogy_to_dot_draw_node(child, stream);
    }
}

void genealogy_to_dot(std::vector<Individual*> population, std::ostringstream& stream)
{
    stream << "digraph BST {" << std::endl;
    stream << "    node [fontname=\"Arial\"];" << std::endl;
    
    for (auto &child: population) { 
      //Rcout << child->get_label() << " has " << child->get_children().size() << " children" << std::endl;
      
      if (!child || child->get_children().size() == 0) {
        continue;
      }
      
      genealogy_to_dot_draw_node(child, stream);
    }
    
    stream << "}" << std::endl;
}
*/


void genealogy_to_dot_draw_node(Individual* node, std::ostringstream& stream, std::vector<int>& mark_ids) {
  if (!node) {
    return;
  } 
  
  static int nullcount = 0;
  int node_id = node->get_id();
  bool mark = false;
  
  for (auto &mark_id: mark_ids) { 
    if (node_id == mark_id) {
      mark = true;
      break;
    }
  }
  
  //stream << "  " << node_id << " [label=\"" << node->get_label() << "\"];" << std::endl;
  stream << "  " << node_id << " [label=\"" << node->get_label() << "\"" << (mark ? " fillcolor=yellow style=filled" : "") << "];" << std::endl;

  std::vector<Individual*> children = node->get_children();
  
  for (auto &child: children) { 
    stream << "    " << node_id << " -> " << child->get_id() << ";" << std::endl;
    genealogy_to_dot_draw_node(child, stream, mark_ids);
  }
}

void genealogy_to_dot(std::vector<Individual*> population, std::ostringstream& stream, bool skip_empty, std::vector<int>& mark_ids) {
  stream << "digraph BST {" << std::endl;
  stream << "    node [fontname=\"Arial\"];" << std::endl;
  
  for (auto &child: population) { 
    //Rcout << child->get_label() << " has " << child->get_children().size() << " children" << std::endl;
    
    if (skip_empty) {
      if (!child || child->get_children().size() == 0) {
        continue;
      }
    }
    
    genealogy_to_dot_draw_node(child, stream, mark_ids);
  }
  
  stream << "}" << std::endl;
}

// Overload
void genealogy_to_dot(std::vector<Individual*> population, std::ostringstream& stream, bool skip_empty) {
  std::vector<int> mark_ids;
  genealogy_to_dot(population, stream, skip_empty, mark_ids);
}

Individual* find_MRCA_with_lineage(Individual* i1, Individual* i2, std::vector<Individual*>& lineage_ids) {
  // FIXME: Not strictly technically required, but population theoretically it is easier to handle
  if (i1->get_generation() != i2->get_generation()) {
    throw std::invalid_argument("i1 and i2 must be individuals from same generation");
  }

  Individual* p1 = i1->get_parent();
  Individual* p2 = i2->get_parent();
  
  if (!p1 || !p2) {
    throw std::invalid_argument("went back to founders, no MRCA found");
  }
  
  if (p1->get_id() == p2->get_id()) {
    lineage_ids.push_back(p1);
    return p1;
  }
  
  // p1 != p2
  lineage_ids.push_back(p1);
  lineage_ids.push_back(p2);
  
  return find_MRCA_with_lineage(p1, p2, lineage_ids);
}

Individual* find_MRCA(Individual* i1, Individual* i2) {
  // FIXME: Not strictly technically required, but population theoretically it is easier to handle
  if (i1->get_generation() != i2->get_generation()) {
    throw std::invalid_argument("i1 and i2 must be individuals from same generation");
  }

  Individual* p1 = i1->get_parent();
  Individual* p2 = i2->get_parent();
  
  if (!p1 || !p2) {
    throw std::invalid_argument("went back to founders, no MRCA found; consider more generations");
  }
  
  if (p1->get_id() == p2->get_id()) {
    return p1;
  }
  
  return find_MRCA(p1, p2);
}

std::vector<int> all_pairwise_MRCA(std::vector<Individual*> population) {
  std::vector<int> gs;
  
  int pop_size = population.size();
  
  if (pop_size <= 1) {
    throw std::invalid_argument("expected pop_size of at least 2");
  }
  
  Rcout << "Considers " << pop_size*(pop_size-1)/2 << " pairs of individuals" << std::endl;
  
  for (int idx1 = 0; idx1 < (pop_size - 1); ++idx1) {
    Individual* i1 = population[idx1];
    
    for (int idx2 = (idx1 + 1); idx2 < pop_size; ++idx2) {
      Individual* i2 = population[idx2];
      
      // i1 and i2 can originate from different founders, warn about more generations required
      try {
        Individual* mrca = find_MRCA(i1, i2);
        int g = i1->get_generation() - mrca->get_generation();
        gs.push_back(g);
      } catch( const std::invalid_argument& e ) {
        std::ostringstream warningstream;
        warning(e.what());
      }
    }
  }
  
  Rcout << "Got " << gs.size() << " actual pairs of individuals with common founder" << std::endl;
  
  return gs;
}



std::vector<int> sample_pairwise_MRCA(std::vector<Individual*> population, int n) {
  std::vector<int> gs;

  if (n <= 0) {
    throw std::invalid_argument("expected n of at least 1 random pair");
  }
  
  int pop_size = population.size();
    
  if (pop_size <= 1) {
    throw std::invalid_argument("expected pop_size of at least 2");
  }
  
  Rcout << "Considers " << n << " random pairs of individuals" << std::endl;
  
  for (int i = 0; i < n; ++i) {
    int idx1 = (int)(R::runif(0, 1)*pop_size);
    int idx2 = (int)(R::runif(0, 1)*pop_size);
    
    while (idx2 == idx1) {
      idx2 = (int)(R::runif(0, 1)*pop_size); // sure to exit as pop_size >= 2
    }
    
    //Rcout << "(" << idx1 << ", " << idx2 << ")" << std::endl;
    
    Individual* i1 = population[idx1];
    Individual* i2 = population[idx2];

    try {
      Individual* mrca = find_MRCA(i1, i2);
      int g = i1->get_generation() - mrca->get_generation();
      gs.push_back(g);
    } catch( const std::invalid_argument& e ) {
      std::ostringstream warningstream;
      warning(e.what());
    }
  }
  
  Rcout << "Got " << gs.size() << " actual pairs of individuals with common founder" << std::endl;
  
  return gs;
}


// [[Rcpp::export]]
Rcpp::XPtr<SimulatedGenealogy> fwpopsim_fixed_genealogy(int G, IntegerVector H0, int pop_size, 
  List mutmodel, bool progress, bool trace, bool cleanup_haplotypes = true, bool cleanup_lineages = true, 
  bool plot = false, bool all_pairs = false, int random_pairs = 0, bool continue_to_one_founder = false) {
  
  Function Rprint("print");

  int mutation_model_type = as<int>(mutmodel["modeltype"]);
  NumericMatrix mutpars = mutmodel["mutpars"];
  
  int loci = H0.length();
   
  /****************************************************************************/
  /* PRINT PARAMETERS */
  /****************************************************************************/
  if (trace) {
    Rcout << "#--- G ------------------------#" << std::endl;
    Rcout << G << std::endl;
    Rcout << std::endl;

    Rcout << "#--- r -----------------------#" << std::endl;
    Rcout << loci << std::endl;
    Rcout << std::endl;
    
    Rcout << "#--- H0 -----------------------#" << std::endl;
    Rprint(H0);
    Rcout << std::endl;
    
    Rcout << "#--- pop size------------------#" << std::endl;
    Rprint(pop_size);
    Rcout << std::endl;

    Rcout << "#--- mutmodel -----------------#" << std::endl;
    Rprint(mutmodel);
  }
  
  /****************************************************************************/
  /* DETERMINE MUTATION MODEL */
  /****************************************************************************/
  MutationModel* mutation_model = NULL;
  
  SMM mutation_smm;
  LMM mutation_lmm;
  EMM mutation_emm;
  
  if (mutation_model_type == 1) {
    mutation_smm = SMM(mutpars);
    mutation_model = &mutation_smm;
  } else if (mutation_model_type == 2) {
    mutation_lmm = LMM(mutpars);
    mutation_model = &mutation_lmm;
  } else if (mutation_model_type == 3) {
    mutation_emm = EMM(mutpars);
    mutation_model = &mutation_emm;
  } else {
    throw std::invalid_argument("The mutation model was not recognized!");
  }
  
  /****************************************************************************/
  /* INITIAL POPULATION */
  /****************************************************************************/

  std::vector<Individual*> init_population(pop_size);
  std::vector<Individual*> population(pop_size);
  
  int id = 0;
  
  std::vector<int> H0vec(loci);
  for (int locus = 0; locus < loci; locus++) {
    H0vec[locus] = H0[locus];
  }
  
  for (int individual = 0; individual < pop_size; ++individual) {
    Individual* ind = new Individual(++id, 0, individual, NULL, H0vec);
    population[individual] = ind;
    init_population[individual] = ind;
  }
  
  int founders_with_descendants = 0;
  
  //for (int generation = 1; generation <= G; generation++) {
  //for (int generation = 1; generation <= G && (!continue_to_one_founder || (continue_to_one_founder && founders_with_descendants == 1)); generation++) {
  //for (int generation = 1; generation <= G && (continue_to_one_founder && founders_with_descendants == 1); generation++) {
  for (int generation = 1; generation <= G || (continue_to_one_founder && founders_with_descendants != 1); generation++) {
    if (trace) {
      Rcout << "===============================================" << std::endl;
      Rcout << "Generation " << generation << " (out of " << G << ")" << std::endl;
      Rcout << "===============================================" << std::endl;
    }
    
    std::vector<Individual*> new_population(pop_size);
    
    // FIXME: Smarter data structure! (Hash)Map?
    std::vector<int> founders(pop_size);
    
    for (int individual = 0; individual < pop_size; individual++) {
      int parent_index = (int)(pop_size*Rf_runif(0, 1));
      Individual* parent = population[parent_index];
      
      std::vector<int> haplotype = parent->get_haplotype();
      
      if (trace) {
        Rcout << "Before mutation: " << sprint_vector(haplotype) << std::endl;
      }

      for (int locus = 0; locus < loci; locus++) {
        double locus_mut_prob[3];
        mutation_model->mutation_table(haplotype[locus], locus, locus_mut_prob);
        
        double mut_down = locus_mut_prob[0];
        double mut_up_cum = locus_mut_prob[0] + locus_mut_prob[1];
        
        double u = Rf_runif(0, 1);

        if (u <= mut_down) {
          haplotype[locus] -= 1;
        } else if (u <= mut_up_cum) {
          haplotype[locus] += 1;
        } 
      }
      
      if (trace) {
        Rcout << "After mutation:  " << sprint_vector(haplotype) << std::endl;
        Rcout << std::endl;
      }
      
      Individual* ind = new Individual(++id, generation, individual, parent, haplotype);
      new_population[individual] = ind;      
      parent->add_child(ind);
      founders[ind->get_founder_id()] += 1;
    } /* for (individual in pop_tree) */
  
    if (cleanup_haplotypes) {
      for (int individual = 0; individual < pop_size; individual++) {
        population[individual]->cleanup_haplotype();
      }
    }
    
    if (cleanup_lineages) {
      for (int individual = 0; individual < pop_size; individual++) {
        Individual::cleanup_lineage(population[individual]);
      }
    }    
    
    population = new_population;

    founders_with_descendants = 0;
    
    for (int individual = 0; individual < pop_size; individual++) {
      if (founders[individual] > 0) {
        founders_with_descendants++;
      }
    }
        
    if (progress && !trace) {
      Rcout << "Generation " << generation << " / " << G << " done (" << founders_with_descendants << " founders left)\r" << 
        std::flush;
    }
  } /* for (generation in 1:G) */

  if (progress && !trace) {
    Rcout << std::endl;
  }  
  
  SimulatedGenealogy* simres = new SimulatedGenealogy(population, init_population);
  
  if (all_pairs) {
    std::vector<int> gs_all = all_pairwise_MRCA(population);
    double mean_all = std::accumulate(gs_all.begin(), gs_all.end(), 0.0) / gs_all.size();
    Rcout << "MRCA mean is " << mean_all << " generations back based on all pairs of individuals" << std::endl;
  }
  
  if (random_pairs > 0) {
    std::vector<int> gs_sample = sample_pairwise_MRCA(population, random_pairs);
    double mean_sample = std::accumulate(gs_sample.begin(), gs_sample.end(), 0.0) / gs_sample.size();
    Rcout << "MRCA mean is " << mean_sample << " generations back based on sample of " << gs_sample.size() << " random pairs" << std::endl;
  }
  
  
  if (plot) {      
    std::ofstream outfile;
    
    
    std::ostringstream dotstream;
    genealogy_to_dot(init_population, dotstream, cleanup_lineages);
    

    outfile.open("tmp-proto.dot", std::ios::out | std::ios::trunc );    
    outfile << dotstream.str();
    outfile.close();
    
    
    ////////////

    int idx1 = 1;
    int idx2 = 8;
    Individual* i1 = population[idx1];
    Individual* i2 = population[idx2];
    std::vector<Individual*> lineage_ids;
    lineage_ids.push_back(i1);
    lineage_ids.push_back(i2);
    Individual* mrca = find_MRCA_with_lineage(i1, i2, lineage_ids);
    int g = i1->get_generation() - mrca->get_generation();
    
    Rcout << "MRCA is " << g << " generations back" << std::endl;
    //Rcout << "MRCA is " << (i1->get_generation() - find_MRCA(i1, i2)->get_generation()) << " generations back" << std::endl;


    std::vector<int> mark_ids;
    for (auto &node: lineage_ids) { 
      mark_ids.push_back(node->get_id());
    }
    
    //mark_ids.push_back(3);
    //mark_ids.push_back(50);
    
    std::ostringstream dotstream_marked;
    genealogy_to_dot(init_population, dotstream_marked, cleanup_lineages, mark_ids);
    
    outfile.open("tmp-proto-marked.dot", std::ios::out | std::ios::trunc );    
    outfile << dotstream_marked.str();
    outfile.close();

  }
  
  //http://www.r-bloggers.com/external-pointers-with-rcpp/
  //http://r.789695.n4.nabble.com/reinterpreting-externalptr-in-R-td4653908.html
  Rcpp::XPtr<SimulatedGenealogy> res(simres);
  
  /*
  CharacterVector classes = res.attr("class") ;
  classes.push_back("myclass") ;
  res.attr("class") = classes;
  */
  res.attr("class") = CharacterVector::create("fwsim_fixed_sim_genealogy", "externalptr");

  return(res);
}

// [[Rcpp::export]]
void print_simulation_info(Rcpp::XPtr<SimulatedGenealogy> object) {  
  // treat object as SimulatedGenealogy*
  Rcout << *object;
}


