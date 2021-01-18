// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]


#include <unordered_set>
#include <unordered_map>
#include <Rcpp.h>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>


using namespace Rcpp;

struct idm_int
{
  int id;
  int set_size;
};

typedef boost::multi_index::multi_index_container<
  idm_int,
  boost::multi_index::indexed_by<
    boost::multi_index::hashed_unique<
      boost::multi_index::member<
        idm_int,int,&idm_int::id
      >
    >,
    boost::multi_index::ordered_non_unique<
      boost::multi_index::member<
        idm_int,int, &idm_int::set_size
      >,
      std::greater<int>
    >
  >
> idm_int_multi;


// [[Rcpp::export]]
IntegerMatrix greedy_set_cover2(
    const IntegerVector &i0,
    const IntegerVector &i1,
    const IntegerVector &group_sizes_i0,
    const IntegerVector &group_sizes_i1
) {
  
  
  std::vector<std::unordered_set<int> > i0_to_i1;
  std::unordered_map<int,std::unordered_set<int> > i1_to_i0;

  // Reserve space in i0_to_i1 and i1_to_i0 from groupsizes.
  std::unordered_set<int> tmp;
  int N_i0=group_sizes_i0.size();
  
  
  for(int i=0;i<N_i0;i++) {
    i0_to_i1.push_back(tmp);
    i0_to_i1[i].reserve(group_sizes_i0[i]);
  }
  
  for(int i=0;i<group_sizes_i1.size();i++) {
    i1_to_i0[i]=tmp;
    i1_to_i0[i].reserve(group_sizes_i1[i]);
  }
  
  //std::cout << "Constructing containers...";
  
  int _i0,_i1;
  
  for (int i=0;i<i1.size();i++) {
    
    _i0=i0[i];
    _i1=i1[i];
    
    i0_to_i1[_i0].insert(_i1);
    i1_to_i0[_i1].insert(_i0);
    
    
  }
  
  //std::cout << " done." << '\n';
  
  // Construct container for set_sizes:
  idm_int_multi idm_set_sizes;

  for(int i=0;i<N_i0;i++){
    int set_size=i0_to_i1[i].size();
    idm_set_sizes.insert({i,set_size});
  }
  
  // Nec. iterator objects:
  auto &set_size_index = idm_set_sizes.get<1>();
  auto &idm_id_index = idm_set_sizes.get<0>();
  
  idm_int idm_int_with_max_set_size;

  
  int co=0;
  int co_sets=0;
  int row_out=0;
  double N_covered=0;
  double coverage=0.0;
  double N_to_cover=double (i1_to_i0.size());
  IntegerMatrix Out(group_sizes_i1.size(),2);

  while(!i1_to_i0.empty() && !i0_to_i1.empty()) {
    
    // Find largest set
    // Set size sorts by size, largest first.
    idm_int_with_max_set_size=*set_size_index.begin(); 
    ++co_sets;
    
    int out0_id=idm_int_with_max_set_size.id;
    std::unordered_set<int> i1_to_delete=i0_to_i1[out0_id];
    N_covered += double (i1_to_delete.size());
    coverage = 100 * (N_covered / N_to_cover); 
    
    for(std::unordered_set<int>::iterator i1_to_delete_it=i1_to_delete.begin();
        i1_to_delete_it !=i1_to_delete.end();
        i1_to_delete_it++){

      int i1_curr=*i1_to_delete_it;
      
      // Adjust Out
      Out(row_out,0)=out0_id;
      Out(row_out,1)=i1_curr;
      ++row_out;
      
      std::unordered_set<int>& i0_which_cover_i1_curr=i1_to_i0[i1_curr];
      int i0_curr;
      
      std::unordered_map<int,int> adjustments_set_sizes;
      
      // Delete current_idm1 from all idm0_sets which cover it
      for(std::unordered_set<int>::iterator i0_which_cover_i1_curr_it=i0_which_cover_i1_curr.begin();
          i0_which_cover_i1_curr_it!=i0_which_cover_i1_curr.end();
          i0_which_cover_i1_curr_it++) {
        
        i0_curr=*i0_which_cover_i1_curr_it;
        
        // Erase current_id from set:
        i0_to_i1[i0_curr].erase(i1_curr);

        // Save how many times set sizes are reduced.
        adjustments_set_sizes[i0_curr]=adjustments_set_sizes[i0_curr]+1;
        
      }
      
      // Adjust set sizes:
      for(std::unordered_map<int,int>::iterator it=adjustments_set_sizes.begin();
          it!=adjustments_set_sizes.end();
          it++) {
        
        auto it_2 = idm_id_index.find(it->first);
        int to_substract=it->second; // Number of elements which were taken from the set.
        idm_id_index.modify(it_2, [&to_substract](idm_int &a){ a.set_size = a.set_size-to_substract;});
        
      }
      
      // Erase idm1 from idm1_map.
      i1_to_i0.erase(i1_curr);
      
    }
    
    if(co_sets % 1000 == 0 || i1_to_i0.empty() || i0_to_i1.empty()){
      
      Rcpp::Rcout << coverage << "% covered by " << 
        round( co_sets * 1000.0 ) / 1000.0 <<  " sets.\n";
  
    }
    
    ++co;
    
  }
  
  return Out;
  
}

