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

IntegerMatrix unordered_set_to_matrix(RCPP_UNORDERED_SET<haplotype, 
  haplotype_hash> pop_tree, const int loci) {
  
  int n = pop_tree.size();
  Rcpp::IntegerMatrix mat(n, loci + 1);
  int row = 0;
  
  for (RCPP_UNORDERED_SET<haplotype, haplotype_hash>::iterator itr = 
    pop_tree.begin(); itr != pop_tree.end(); ++itr) {
    
    for (int k = 0; k < loci; k++) {
      mat(row, k) = itr->profile[k];
    }
    
    mat(row, loci) = itr->count;
    row += 1;
  }
  
  return mat;
}

/*
  _G:        Guarantied to be an IntegerVector with length 1
  _H0:       Guarantied to be an n x r IntegerMatrix
  _N0:       Guarantied to be an IntegerVector with length n
  _alpha:    Guarantied to be an NumericVector with length G
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
List Cpp_fwpopsim(int G, IntegerMatrix H0, IntegerVector N0, NumericVector alpha, 
  List mutmodel, bool SNP, IntegerVector save_gs, bool progress, bool trace, bool ensure_children = false) {
  
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

    Rcout << "#--- alpha --------------------#" << std::endl;
    print_alpha(alpha, G);
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
  RCPP_UNORDERED_SET<haplotype, haplotype_hash> pop_tree;

  //std::ofstream hashes_file;
  //hashes_file.open("hashes.txt");
  
  for (int i = 0; i < n; i++) {
    std::vector<int> profile(loci);
    
    for (int locus = 0; locus < loci; locus++) {
      profile[locus] = H0(i, locus);
    }
    
    const haplotype h = { profile, N0(i) };
    pop_tree.insert(h);
    
    //hashes_file << sprint_vector(profile) << "\t" << hfunc(h) << std::endl;
  }

  IntegerMatrix pop_mat = unordered_set_to_matrix(pop_tree, loci);

  if (trace) {
    Rcout << "===============================================" << std::endl;
    Rcout << "Initial population (generation 0)" << std::endl;
    Rcout << "===============================================" << std::endl;
    Rprint(pop_mat);
  }
  
  /****************************************************************************/
  /* EVOLUTION */
  /****************************************************************************/
  
  List saved_populations(G - 1);
  std::vector< RCPP_UNORDERED_SET<haplotype, haplotype_hash> > pop_trees_vec(G);

  for (int generation = 0; generation < G; generation++) {
    if (trace) {
      Rcout << "===============================================" << std::endl;
      Rcout << "Generation " << generation + 1 << " (out of " << G << 
        ") with growth " << alpha[generation] << std::endl;
      Rcout << "===============================================" << std::endl;
    }
    
    /*
    if (generation >= 1) {
      stop("DONE");
    }
    */
    
    /*
    The previous generation is no longer needed, only if it should be saved.
    Else, free it.
    
    Note, that we will never catch pop_trees_vec[G-1], but
    as it is the end population (that is returned), this is 
    as it should be.
    */    
    if (generation > 0 && save_gs[generation - 1] == 0) {
      pop_trees_vec[generation- 1].clear();
    }

    RCPP_UNORDERED_SET<haplotype, haplotype_hash> new_pop_tree;

    for (int individual = 0; individual < pop_mat.nrow(); individual++) {
      IntegerVector person = pop_mat(individual, Rcpp::_);
      std::vector<int> profile(loci);
      int count = person[loci];
      
      for (int locus = 0; locus < loci; locus++) {
        profile[locus] = person(locus);
      }
      
      std::vector<haplotype> children;
      
      //int children_count = (int)Rf_rpois(alpha[generation]*(double)count);
      // FIXME: 2016-05-20
      int children_count = (ensure_children) ? (int)Rf_rpois(alpha[generation]*(double)count) + 1 : (int)Rf_rpois(alpha[generation]*(double)count);
      
      if (trace) {
        Rcout << "-----------------------------------------------" 
          << std::endl;
        print_vector(profile);
        Rcout << " x " << count << " gets " << children_count << 
          " children" << std::endl;
        Rcout << "-----------------------------------------------" 
          << std::endl;
      }

      if (children_count == 0) {
        continue;
      }

      const haplotype h = { profile, children_count };
      children.push_back(h);
                        
      for (int locus = 0; locus < loci; locus++) {
        if (trace) {
          Rcout << "Locus " << locus + 1 << std::endl;
        }
        
        std::vector<haplotype> new_children;
        
        /*
        We will be adding to the children list if a mutation happens.
        Thus, we must now only iterate 
        */
        for (std::vector<haplotype>::iterator it = children.begin(); 
          it != children.end(); ++it) {
          
          double locus_mut_prob[3];
          mutation_model->mutation_table(it->profile[locus], locus, locus_mut_prob);
          
          int child_count = it->count;
          int result[3]; /* 3: down, up, not mutate */
          R::rmultinom(child_count, locus_mut_prob, 3, result);

          if (trace) {
            Rcout << "  Children ";
            print_vector(it->profile);
            Rcout << " x " << child_count << " distributed as follows:" << std::endl;

            Rcout << "    #dw mutation = " << result[0] << 
              ";\tP(dw mutation) = " << locus_mut_prob[0] << std::endl;
              
            Rcout << "    #up mutation = " << result[1] << 
              ";\tP(up mutation) = " << locus_mut_prob[1] << std::endl;
              
            Rcout << "    #no mutation = " << result[2] << 
              ";\tP(no mutation) = " << locus_mut_prob[2] << std::endl;
          }
          
          /* dw mutation */
          if (result[0] > 0) {
            std::vector<int> mutated_profile = it->profile;

            if (SNP) {
              /* +1 on purpose to avoid negative values */
              mutated_profile[locus] = (mutated_profile[locus] + 1) % 2; 
            } else {
              mutated_profile[locus] -= 1;
            }
            
            const haplotype mut_h = { mutated_profile, result[0] };
            new_children.push_back(mut_h);
          }
          
          /* up mutation */
          if (result[1] > 0) {
            std::vector<int> mutated_profile = it->profile;
            
            if (SNP) {
              mutated_profile[locus] = (mutated_profile[locus] + 1) % 2;
            } else {
              mutated_profile[locus] += 1;
            }
            
            const haplotype mut_h = { mutated_profile, result[1] };
            new_children.push_back(mut_h);
          }
          
          /* No mutation is updated */
          it->count = result[2];
        }
        
        /* Mutated children during children loop is added to children for next locus */
        for (std::vector<haplotype>::iterator it = 
          new_children.begin(); it != new_children.end(); ++it) {
          
          children.push_back(*it);
        }
      }
      
      for (std::vector<haplotype>::iterator it = children.begin(); 
        it != children.end(); ++it) {
        
        if (it->count > 0) {
          RCPP_UNORDERED_SET<haplotype, haplotype_hash>::iterator got = 
            new_pop_tree.find(*it);

          //hashes_file << sprint_vector(it->profile) << "\t" << hfunc(*it) << std::endl;
          
          if (got == new_pop_tree.end()) {
            new_pop_tree.insert(*it);
          } else {
            const haplotype updated_h = { got->profile, got->count + it->count };
            new_pop_tree.erase(got);
            new_pop_tree.insert(updated_h);
          }
        }
      }

      if (trace) {
        Rcout << "The resulting " << children_count << " children that ";
        print_vector(profile);
        Rcout << " x " << count << " gets are:" << std::endl;
        
        for (std::vector<haplotype>::iterator it = children.begin(); 
          it != children.end(); ++it) {
          
          Rcout << "  " << *it << std::endl;
        }
      }      
    } /* for (individual in pop_mat) */
    
    pop_trees_vec[generation] = new_pop_tree;
    pop_mat = unordered_set_to_matrix(new_pop_tree, loci);
    
    int size = 0;
    for (int individual = 0; individual < pop_mat.nrow(); individual++) {
      size += pop_mat(individual, loci);
    }
    
    popsizes[generation] = size;
    
    if (generation < (G - 1)) {
      if (save_gs[generation] == 1) {
        saved_populations[generation] = pop_mat;
      } else {
        saved_populations[generation] = R_NaInt;
      }
    }
    
    if (trace) {
      Rcout << "###############################################" << std::endl;
      Rcout << "Population after evolution of this generation" << std::endl;
      Rcout << "###############################################" << std::endl;
      Rprint(pop_mat);
      Rcout << "###############################################" << std::endl;
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
  pars["alpha"] = alpha;
  pars["mutmodel"] = mutmodel;
  pars["progress"] = progress;
  pars["trace"] = progress;

  List ret;
  ret["pars"] = pars;
  ret["saved_populations"] = saved_populations;
  ret["population"] = pop_mat;
  ret["pop_sizes"] = popsizes; 

  return(ret);
}

