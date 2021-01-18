#include "fwpopsim.h"
//#include "gperftools/profiler.h"
//#include <iostream>
//#include <fstream>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

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
  _G:        Guarantied to be an IntegerVector with length 1
  _H0:       Guarantied to be an n x r IntegerMatrix
  _N0:       Guarantied to be an IntegerVector with length n
  _mutmodel: Guarantied to have class "haptools_mutmodel", thus a list with names
               modeltype
               modeldescription
               mutdw.description
               mutup.description
               mutpars
               mutlogic.mut.dw
               mutlogic.mut.up
               mutlogic.not.mut
               prob.mut.per.locus
  
  _save_gs:   Guarantied to be an IntegerVector with length G of 0's and 1's
  _progress:  Guarantied to be an LogicalVector with length 1
*/
// [[Rcpp::export]]
List Cpp_fwpopsim_fixed(int G, IntegerMatrix H0, IntegerVector N0, 
  List mutmodel, bool SNP, IntegerVector save_gs, bool progress, bool trace) {
  
  Function Rprint("print");

  int mutation_model_type = as<int>(mutmodel["modeltype"]);
  NumericMatrix mutpars = mutmodel["mutpars"];
  
  int n = H0.nrow();
  int loci = H0.ncol();
  
  IntegerVector popsizes(G);
  
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
    
    Rcout << "#--- N0 -----------------------#" << std::endl;
    Rprint(N0);
    Rcout << std::endl;

    Rcout << "#--- mutmodel -----------------#" << std::endl;
    Rprint(mutmodel);
    Rcout << "#--- save.gs ------------------#" << std::endl;
    print_save_gs(save_gs, G);
    Rcout << std::endl;
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
  int pop_size = 0;
  
  for (int i = 0; i < n; i++) {
    pop_size += N0(i);
  }
  
  IntegerMatrix pop_tree(pop_size, loci);
  IntegerVector parent_vec(pop_size);
  IntegerVector founder_vec(pop_size);

  for (int individual = 0; individual < pop_size; individual++) {
    parent_vec[individual] = R_NaInt;
    founder_vec[individual] = individual;
  }
      
  int pop_tree_index = 0;
  
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < N0(i); k++) {
      for (int locus = 0; locus < loci; locus++) {
        pop_tree(pop_tree_index, locus) = H0(i, locus);
      }
      pop_tree_index += 1;
    }
  }
  
  /*
  List saved_populations(G - 1);  
  List saved_parents(G - 1); 
  List saved_founders(G - 1);
  */

  List saved_populations(G-1);
  IntegerMatrix saved_parents(G-1, pop_size); 
  IntegerMatrix saved_founders(G-1, pop_size);
  
  IntegerVector NAVEC(pop_size);
  for (int individual = 0; individual < pop_size; individual++) {
    NAVEC[individual] = R_NaInt;
  }


  for (int generation = 0; generation < G; generation++) {
    if (trace) {
      Rcout << "===============================================" << std::endl;
      Rcout << "Generation " << generation + 1 << " (out of " << G << ")" << std::endl;
      Rcout << "===============================================" << std::endl;
    }

    IntegerMatrix new_pop_tree(pop_size, loci);
    IntegerVector new_parent_vec(pop_size);
    IntegerVector new_founder_vec(pop_size);

    for (int individual = 0; individual < pop_size; individual++) {
      int parent = (int)(pop_size*Rf_runif(0, 1));
      IntegerVector person = pop_tree(parent, Rcpp::_);
      new_pop_tree(individual, Rcpp::_) = person;
      new_parent_vec[individual] = parent;
      new_founder_vec[individual] = founder_vec[parent];

      for (int locus = 0; locus < loci; locus++) {
        double locus_mut_prob[3];
        mutation_model->mutation_table(person[locus], locus, locus_mut_prob);
        
        double mut_down = locus_mut_prob[0];
        double mut_up_cum = locus_mut_prob[0] + locus_mut_prob[1];
        
        double u = Rf_runif(0, 1);
        
        /*
        Rprint(mut_down);
        Rprint(mut_up_cum);
        Rprint(u);
        */
        
        if (u <= mut_down) {
          new_pop_tree(individual, locus) = new_pop_tree(individual, locus) - 1;
        } else if (u <= mut_up_cum) {
          new_pop_tree(individual, locus) = new_pop_tree(individual, locus) + 1;
        } 
      }
    } /* for (individual in pop_tree) */
    
    popsizes[generation] = pop_size;
    
    //if (generation < (G - 1)) {
    if (generation < (G - 1)) {
      if (save_gs[generation] == 1) {
        saved_populations[generation] = pop_tree;
        saved_parents(generation, Rcpp::_) = parent_vec;
        saved_founders(generation, Rcpp::_) = founder_vec;
      } else {
        saved_populations[generation] = R_NaInt;
        saved_parents(generation, Rcpp::_) = NAVEC;
        saved_founders(generation, Rcpp::_) = NAVEC;
      }
    }
    
    pop_tree = new_pop_tree;
    parent_vec = new_parent_vec;
    founder_vec = new_founder_vec;
    
    if (trace) {
      Rprint(pop_tree);
    }
    
    if (progress && !trace) {
      Rcout << "Generation " << generation + 1 << " / " << G << " done\r" << 
        std::flush;
    }
  } /* for (generation in 1:G) */

  if (progress && !trace) {
    Rcout << std::endl;
  }
  
  //hashes_file.close();
    
  /****************************************************************************/
  /* RETURN VALUE */
  /****************************************************************************/
  List pars;
  pars["G"] = G;
  pars["H0"] = H0;
  pars["N0"] = N0;
  pars["alpha"] = NA_REAL;
  pars["mutmodel"] = mutmodel;
  pars["progress"] = progress;
  pars["trace"] = progress;

  List extra;
  extra["parents"] = parent_vec; 
  extra["founders"] = founder_vec; 
  extra["saved_parents"] = saved_parents; 
  extra["saved_founders"] = saved_founders; 
  
  List ret;
  ret["pars"] = pars;
  ret["saved_populations"] = saved_populations;
  ret["population"] = pop_tree;
  ret["pop_sizes"] = popsizes; 
  ret["extra"] = extra; 

  return(ret);
}

