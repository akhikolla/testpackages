extern "C" {
	#include "sailr.h"
	#include "sailr_ext.h"
}
#include "datasailr_ext_func.hpp"


#include <fstream>
#include <iostream>
#include <string>
#include <iterator>
#include <tuple>
#include <cstring>

#include "stdlib.h"
#include "stdio.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

/* --------------------------------------------------------- */
// Macro variables
typedef unsigned int SXPTYPE;
#define INTSXP 13
#define REALSXP 14
#define STRSXP 16
#define NILSXP 0

#define INTNUM 0
#define DBLNUM 1

#define ORIGINAL 0
#define UPDATED 1

/* --------------------------------------------------------- */
// Conditions
#define SAILR_NULL_UPDATED

// Macro
#ifdef DEBUG
#define IF_DEBUG(x) do{ x } while(false) 
#define IF_DEBUG_DECL(x) x
#else
#define IF_DEBUG(x) do{ } while( false )
#define IF_DEBUG_DECL(x) do{ } while( false )
#endif

// nullptr support

#ifdef NEEDS_NULLPTR_SUPPORT

const                        // this is a const object...
class {
public:
  template<class T>          // convertible to any type
    operator T*() const      // of null non-member
    { return 0; }            // pointer...
  template<class C, class T> // or any type of null
    operator T C::*() const  // member pointer...
    { return 0; }
private:
  void operator&() const;    // whose address can't be taken
} nullptr = {};              // and whose name is nullptr

#endif

/* --------------------------------------------------------- */

// Each element of tuple corresponds to 
// 0. column name, 1. original vector, 2. SXPTYPE, 3. size, 4. extra vector (doubles for ints. ints for doubles, new STRs for STRSXP.).
// 5. type information(Previously std::vector<int>* ) 6. std::string* class_name 7. std::vector<std::string>* levels 
typedef std::tuple< char* , void* , SXPTYPE , int , void*, void* , std::string*, std::vector<std::string>* > VEC_ELEM;
typedef std::vector< VEC_ELEM > VEC_LIST; 

void vec_list_add_int_vec( VEC_LIST* vec_list, char* var_name, IntegerVector* r_vec , int size);
void vec_list_add_double_vec( VEC_LIST* vec_list, char* var_name, NumericVector* r_vec , int size);
void vec_list_add_string_vec( VEC_LIST* vec_list, char* var_name, StringVector* r_vec , int size);
VEC_LIST* ConvertDataFrame( DataFrame df , char** var_names , int num_of_vars, char** lhs_var_names, int num_of_lhs_vars, char** srv_var_names, int num_of_srv_vars, int* conversion_error);


int vec_list_nrows(VEC_LIST* vec_list);
VEC_ELEM* vec_elem_find(VEC_LIST* vl, const char* var_name);
bool   cstring_exists_in_charactervector(char* var_name, CharacterVector var_vector);
SXPTYPE  vec_elem_type_of(VEC_ELEM* vec_elem);
char* vec_elem_name_of(VEC_ELEM* vec_elem);

void
vec_list_add_int_vec( VEC_LIST* vec_list, char* var_name, IntegerVector* r_vec, int size )
{
  std::vector<int>* cpp_i_vec;
  std::vector<double>* cpp_d_vec;
  std::vector<int>* cpp_type_vec;
  if(r_vec != NULL){
    cpp_i_vec = new std::vector<int> (size);
    cpp_d_vec = new std::vector<double> (size);
    cpp_type_vec = new std::vector<int> (size, INTNUM );
    auto cpp_iter = cpp_i_vec->begin();
    auto cpp_d_iter = cpp_d_vec->begin();
    auto type_iter = cpp_type_vec->begin();
    for( IntegerVector::iterator r_iter = r_vec->begin(); r_iter != r_vec->end() ; ++r_iter){
      if( (! Rcpp::IntegerVector::is_na(*r_iter)) && (! Rcpp::traits::is_nan<INTSXP>(*r_iter) ) ){ 
      // Not NA or NaN in IntegerVector
        *cpp_iter = *r_iter;
        *type_iter = INTNUM;
      }else{
      // NA or NaN in IntegerVector
        *cpp_d_iter = NA_REAL; // -nan in C++
        *type_iter = DBLNUM;
      }
      ++cpp_iter;
      ++cpp_d_iter;
      ++type_iter;
    }
  } else {
    cpp_i_vec = new std::vector<int> (size);
    cpp_d_vec = new std::vector<double> (size, NA_REAL );
    cpp_type_vec = new std::vector<int> (size, DBLNUM );
  }
    
  VEC_ELEM new_vec_elem = VEC_ELEM { strdup(var_name), (void*) cpp_i_vec, INTSXP , size, (void*) cpp_d_vec , (void*) cpp_type_vec , nullptr, nullptr};
  vec_list->push_back( new_vec_elem );
}

void
vec_list_add_double_vec( VEC_LIST* vec_list, char* var_name, NumericVector* r_vec , int size)
{
  std::vector<double>* cpp_d_vec;

  if(r_vec != NULL){
    cpp_d_vec = new std::vector<double> (size, NA_REAL);
    auto cpp_iter = cpp_d_vec->begin();
    for( auto r_iter = r_vec->begin(); r_iter != r_vec->end() ; ++r_iter){
      if( (! Rcpp::NumericVector::is_na(*r_iter)) && (! Rcpp::traits::is_nan<REALSXP>(*r_iter) ) ){ 
      // Not NA_REAL or R_NaN in NumericVector
        *cpp_iter = *r_iter;
		++cpp_iter;
      }else{
      // NA_REAL or R_NaN
        *cpp_iter = NA_REAL; // nan in C++
		++cpp_iter;
      }
    }
  }else{
    cpp_d_vec = new std::vector<double> (size, NA_REAL);
  }
  std::vector<int>* cpp_i_vec = new std::vector<int> (size);
  std::vector<int>* cpp_type_vec = new std::vector<int> (size, DBLNUM );
  VEC_ELEM new_vec_elem = VEC_ELEM { strdup(var_name), (void*) cpp_d_vec, REALSXP , size, (void*) cpp_i_vec, (void*) cpp_type_vec, nullptr, nullptr};
  vec_list->push_back( new_vec_elem );
}

void
vec_list_add_string_vec( VEC_LIST* vec_list, char* var_name, StringVector* r_vec , int size )
{
  std::vector<std::string*>* pvec_pstr ; 
  if(r_vec != NULL){
    pvec_pstr = new std::vector<std::string*>(size);
    int idx;
    for(idx = 0; idx < size ; ++idx ){
      if( ! Rcpp::StringVector::is_na(r_vec->operator[](idx))){
        pvec_pstr->operator[](idx) = new std::string(r_vec->operator[](idx));
      }else{
        pvec_pstr->operator[](idx) = NULL;
      }
    }
  }else{
    pvec_pstr = new std::vector<std::string*>(size, NULL);
  }
  std::vector<std::string* > *new_pvec_pstr = new std::vector<std::string* >(size, NULL);
  std::vector<int>* cpp_updated_vec = new std::vector<int> ( size, ORIGINAL );
  VEC_ELEM new_vec_elem = VEC_ELEM { strdup(var_name), (void*) pvec_pstr , STRSXP, size , (void*) new_pvec_pstr, (void*) cpp_updated_vec, nullptr, nullptr};
  vec_list->push_back( new_vec_elem );
}

void
vec_list_add_factor_vec( VEC_LIST* vec_list, char* var_name, IntegerVector* r_vec , int size )
{
  std::vector<std::string*>* pvec_pstr ; 

  std::string* class_name;
  StringVector r_attr_class = r_vec->attr("class");
  class_name = new std::string(r_attr_class[0]);

  std::vector<std::string>* factor_levels;
  StringVector r_factor_levels = r_vec->attr("levels");
  int level_size = r_factor_levels.size();
  factor_levels= new std::vector<std::string>(level_size);
  int level_idx;
  for(level_idx = 0 ; level_idx < level_size; ++level_idx){
    factor_levels->operator[](level_idx) = r_factor_levels[level_idx];
  }

  if(r_vec != NULL){
    pvec_pstr = new std::vector<std::string*>(size);
    int idx;
    for(idx = 0; idx < size ; ++idx ){
      if( ! Rcpp::IntegerVector::is_na(r_vec->operator[](idx))){
        pvec_pstr->operator[](idx) = new std::string( factor_levels->operator[](r_vec->operator[](idx) - 1));
      }else{
        pvec_pstr->operator[](idx) = NULL;
      }
    }
  }else{
    // This branch should never be run
    pvec_pstr = new std::vector<std::string*>(size, NULL);
  }
  std::vector<std::string* > *new_pvec_pstr = new std::vector<std::string* >(size, NULL);
  std::vector<int>* cpp_updated_vec = new std::vector<int> ( size, ORIGINAL );
  VEC_ELEM new_vec_elem = VEC_ELEM {strdup(var_name), (void*) pvec_pstr , STRSXP, size , (void*) new_pvec_pstr, (void*) cpp_updated_vec, class_name, factor_levels};
  vec_list->push_back( new_vec_elem );
}

void  
vec_list_add_null_vec ( VEC_LIST* vec_list, char* var_name , int size)
{
  std::vector<int>* cpp_vec = new std::vector<int>(size);
  VEC_ELEM new_vec_elem = VEC_ELEM { strdup(var_name), (void*) cpp_vec, NILSXP , size, nullptr, nullptr, nullptr, nullptr};
  vec_list->push_back( new_vec_elem );
}


void
vec_list_free( VEC_LIST* vl){
	IF_DEBUG_DECL( char* var_name; );
	void* column_vec1;
	void* column_vec2;
	void* column_vec3;
	std::vector<int>* int_vec;
	std::vector<double>* double_vec;
	std::vector<int>* type_vec;
	std::vector<std::string* >* strvec_ori;
	std::vector<std::string* >* strvec;
	std::vector<int>* updated_vec;
	
	
	for( auto it = vl->begin(); it!=vl->end() ; ++it){
		IF_DEBUG( var_name = std::get<0>(*it); );
		switch( std::get<2>(*it)){
		case INTSXP:
			IF_DEBUG( Rcpp::Rcout << "Free INTSXP element (" << var_name << ")"  << std::endl; );
			column_vec1 = std::get<1>(*it);
			int_vec = (std::vector<int>*) column_vec1;
			delete int_vec;	
			column_vec2 = std::get<4>(*it);
			double_vec = (std::vector<double>*) column_vec2;
			delete double_vec;	
			column_vec3 = std::get<5>(*it);
			type_vec = (std::vector<int>*) column_vec3;
			delete type_vec;
			break;
		case REALSXP:
			IF_DEBUG( Rcpp::Rcout << "Free REALSXP element (" << var_name << ")"  << std::endl; );
			column_vec1 = std::get<1>(*it);
			double_vec = (std::vector<double>*) column_vec1;
			delete double_vec;	
			column_vec2 = std::get<4>(*it);
			int_vec = (std::vector<int>*) column_vec2;
			delete int_vec;	
			column_vec3 = std::get<5>(*it);
			type_vec = (std::vector<int>*) column_vec3;
			delete type_vec;
			break;
		case STRSXP:
			IF_DEBUG( Rcpp::Rcout << "Free STRSXP element (" << var_name << ")"  << std::endl; );
			column_vec1 = std::get<1>(*it);
			strvec_ori = (std::vector<std::string* >*) column_vec1;
			for( auto str_iter = strvec_ori->begin(); str_iter != strvec_ori->end(); ++str_iter){
				delete (*str_iter) ;
			}
			delete strvec_ori;

			column_vec2 = std::get<4>(*it);
			strvec = (std::vector<std::string* >*)column_vec2;
			for( auto str_iter = strvec->begin(); str_iter != strvec->end(); ++str_iter){
				delete (*str_iter) ;
			}
			delete strvec;

			column_vec3 = std::get<5>(*it);
			updated_vec = (std::vector<int>*)column_vec3;
			delete updated_vec;
			break;
		case NILSXP:
			IF_DEBUG( Rcpp::Rcout << "Free NILSXP element (" << var_name << ")"  << std::endl; );
			column_vec1 = std::get<1>(*it);
			int_vec = (std::vector<int>*) column_vec1;
			delete int_vec;
			break;
		default:
			IF_DEBUG( Rcpp::Rcout << "Unintended type of element is found (" << var_name << ")"  << std::endl; );
			break;
		}
        free(std::get<0>(*it) );
        delete std::get<6>(*it);
        delete std::get<7>(*it);
	}

	delete vl;
}

bool
cstring_exists_in_pchar_array(const char* cstring, char** array, int from_idx, int to_idx)
{
  int idx;
  char* current_cstring;
  for( idx = from_idx; idx <= to_idx; ++idx){
    current_cstring = array[idx];
    if( 0 == std::strcmp(cstring, current_cstring) ){
      return true;
    }
  }
  return false;
}

VEC_LIST*
ConvertDataFrame( DataFrame df , char** var_names , int num_of_vars, char** lhs_var_names, int num_of_lhs_vars, char** srv_var_names, int num_of_srv_vars, int* conversion_error)
{
  // all_vars store var_names + srv_var_names
  IF_DEBUG(Rcpp::Rcout << "all_vars = all the variable names that appear in souce + spceial variables. " << std::endl; );
  int num_of_all_vars = num_of_vars + num_of_srv_vars;
  char** all_vars = (char**) malloc(sizeof(char*) * num_of_all_vars);
  int all_var_idx = 0;
  while(all_var_idx < num_of_vars ){
    all_vars[all_var_idx] = var_names[all_var_idx];
    ++all_var_idx;
  }
  int srv_var_idx = 0;
  while(srv_var_idx < num_of_srv_vars){
    if(! cstring_exists_in_pchar_array( srv_var_names[srv_var_idx] , all_vars, 0, num_of_vars - 1)){
      all_vars[all_var_idx] =  srv_var_names[srv_var_idx];
    }else{
      all_vars[all_var_idx] =  NULL;
    }
    ++all_var_idx;
    ++srv_var_idx;
  }

  // Convert DataFrame columns into VEC_LIST vector elements.
  char* var_name;
  int var_idx;
  VEC_LIST* vec_list = new VEC_LIST;
  CharacterVector df_var_vector = df.attr("names") ;
  CharacterVector lhs_var_vector;
  CharacterVector srv_var_vector;
  int lhs_var_idx;
  for( lhs_var_idx = 0; lhs_var_idx < num_of_lhs_vars; ++lhs_var_idx){
    lhs_var_vector.push_back( lhs_var_names[lhs_var_idx]); 
  }

  // Dataframe column variables
  IF_DEBUG( Rcpp::Rcout << "Dataframe variable names: " << std::endl ; );
  IF_DEBUG( Rcpp::Rcout << df_var_vector << std::endl ; );

  // All the variable names in source (= LHS U RHS) + service variables (which should be always included from dataframe)
  IF_DEBUG( Rcpp::Rcout << "All the variable names in source + service variable names (which should be always included from dataframe): " << std::endl ; );
  // IF_DEBUG( show each char* in var_names );

  // foreach var_name (= All the variable names appear in source code.)
  for( var_idx = 0 ; var_idx < num_of_all_vars; ++var_idx ){
      var_name = all_vars[var_idx] ;
    if( var_name == NULL ){
      IF_DEBUG( Rcpp::Rcout << "SKIP (service name already exists in source, and corresponding vector is included in VEC_LIST.)." << std::endl; ); 
      goto SKIP_ADD_NEW_VEC;
    };

    // Even if the DataFrame does not have the variable name
    // , and the variable exist on LHS
    // then the variable can be allowed.(The variable can be initialized dynamically)
    if(! cstring_exists_in_charactervector(var_name, df_var_vector)){
    	if(cstring_exists_in_charactervector(var_name, lhs_var_vector)){
   			IF_DEBUG( Rcpp::Rcout << "\"" << var_name << "\"" <<  " does not exist in dataframe column name but appears as LHS. " ; );
			IF_DEBUG( Rcpp::Rcout << "Added to VEC_LIST as null vector" << std::endl; );
			vec_list_add_null_vec ( vec_list, var_name, df.nrows() );
			continue;
		}else{
			Rcout << "Error: " << "\"" << var_name << "\"" << " appears in source code, but deoes not exist in dataframe or appear as LHS. Never to be defined." << std::endl;
			*conversion_error = 1;
		}
    }else{
		// This branch means the var_name exists in dataframe.
		IF_DEBUG( Rcpp::Rcout << "\"" << var_name << "\"" << " exists in dataframe column name (=type known) and also appears as variables (LHS or RHS) or as service variables." ; );

    	// Rcpp:GenericVector col_vec = df[var_name];    
		IntegerVector int_vec;
		NumericVector num_vec;
		StringVector  str_vec;

		StringVector r_attr_class;
		std::string r_class_name;

		IF_DEBUG( Rcpp::Rcout << "Type: " << TYPEOF(df[var_name]) << std::endl ;  );
		switch( TYPEOF(df[var_name]) ){
		case INTSXP:
			int_vec = df[var_name];
			if(int_vec.hasAttribute("class")){
				r_attr_class = int_vec.attr("class");
				r_class_name = (std::string) r_attr_class[0];
			}else{
				r_class_name = std::string("");
			}
			IF_DEBUG(Rcpp::Rcout << "class attribute is " << r_class_name << std::endl; );
			if( r_class_name.compare( "factor") != 0 ){  // Normal IntegerVector
				vec_list_add_int_vec ( vec_list, var_name, &int_vec , int_vec.size());
				IF_DEBUG( Rcpp::Rcout << "Added to VEC_LIST as int vector" << std::endl; );
			}else{  // Factor IntegerVector
				vec_list_add_factor_vec ( vec_list, var_name, &int_vec , int_vec.size());
				IF_DEBUG( Rcpp::Rcout << "Added to VEC_LIST as string vector (<= from IntegerVector with 'factor' class attribute)" << std::endl; );
			}
			break;
		case REALSXP:
			num_vec = df[var_name];
			vec_list_add_double_vec ( vec_list, var_name, &num_vec, num_vec.size());
			IF_DEBUG( Rcpp::Rcout << "Added to VEC_LIST as double vector" << std::endl; );
			break;	
		case STRSXP:
			str_vec = df[var_name];
			vec_list_add_string_vec ( vec_list, var_name, &str_vec, str_vec.size());
			IF_DEBUG( Rcpp::Rcout << "Added to VEC_LIST as string vector" << std::endl; );
			break;
		default:
			Rcpp::Rcout << "Error: " << var_name << "'s type of this dataframe is not supported. (TYPEOF value is " <<  TYPEOF(df[var_name]) << ")" << std::endl;
			break;
		}
	}
    SKIP_ADD_NEW_VEC: ;
  }
  IF_DEBUG( Rcpp::Rcout << "DataFrame is converted." << std::endl; );
  free(all_vars);
  return vec_list;
}	

IntegerVector
reorder_intvec(IntegerVector vec, IntegerVector order)
{
    return vec[order];
}
NumericVector
reorder_dblvec(NumericVector vec, IntegerVector order)
{
    return vec[order];
}
StringVector
reorder_strvec(StringVector vec, IntegerVector order)
{
    return vec[order];
}
#define REORDER_INTVEC( ori_vec, order_vec) do{ if( ! IntegerVector::is_na(order_vec[0]) ){ ori_vec = reorder_intvec(ori_vec, order_vec); } }while(0)
#define REORDER_DBLVEC( ori_vec, order_vec) do{ if( ! IntegerVector::is_na(order_vec[0]) ){ ori_vec = reorder_dblvec(ori_vec, order_vec); } }while(0)
#define REORDER_STRVEC( ori_vec, order_vec) do{ if( ! IntegerVector::is_na(order_vec[0]) ){ ori_vec = reorder_strvec(ori_vec, order_vec); } }while(0)

// #define REORDER_INTVEC( ori_vec, order_vec) do{if(order_vec == R_NilValue){Rcpp::Rcout << "REORDER_INTVEC" << std::endl;}}while(0)
// #define REORDER_DBLVEC( ori_vec, order_vec) do{if(order_vec == R_NilValue){Rcpp::Rcout << "REORDER_DBLVEC" << std::endl;}}while(0)
// #define REORDER_STRVEC( ori_vec, order_vec) do{if(order_vec == R_NilValue){Rcpp::Rcout << "REORDER_STRVEC" << std::endl;}}while(0)

DataFrame
ConvertVecList(VEC_LIST* vl, std::vector<std::string> lvars, std::vector<int>* order, std::vector<std::string> vars_excluded)
{
  DataFrame new_df;
  CharacterVector name_vec;
  std::vector<std::string> var_name_list;
  char* var_name;
  std::string cpp_var_name;
  
  void* column_vec1;
  void* column_vec2;
  void* column_vec3;

  IF_DEBUG(Rcpp::Rcout << "Start Convert VEC_LIST. " << std::endl;);

  IntegerVector rcpp_order;
  if(order != NULL){
      rcpp_order = Rcpp::wrap(*order);
  }else{
      rcpp_order = IntegerVector::create(NA_INTEGER);
  }

  IF_DEBUG(Rcpp::Rcout << "R_NilValue is set. " << std::endl;);

  IntegerVector intvec;
  NumericVector dblvec;
  IntegerVector typevec;
  NumericVector numvec;
  IntegerVector new_intvec;
  
  std::vector<std::string* > strvec;
  std::vector<std::string* > strvec_ori;
  std::vector<int> updated_vec;
  StringVector rstrvec;
  LogicalVector dbl_pos;
  LogicalVector nan_pos;

  int index;
  unsigned int idx;
  
  for(auto it = vl->begin(); it != vl->end(); ++it ){
    var_name = std::get<0>(*it);
    cpp_var_name = std::string(var_name);

    if (std::find(lvars.begin(), lvars.end(), cpp_var_name) != lvars.end()){ // Varible name is included in lvars
      if(std::find(vars_excluded.begin(), vars_excluded.end(), cpp_var_name) == vars_excluded.end()){ // Variable name is not included in vars_excluded
	  IF_DEBUG( Rcpp::Rcout << "VEC_LIST column " << cpp_var_name << " : going to be converted" << std::endl;  );
        goto MATCH;
      }
    }
    IF_DEBUG( Rcpp::Rcout << "VEC_LIST column " << cpp_var_name << " : conversion is skipped" << std::endl;  );
 
    continue;  // if the flow is not MATCH, go to next loop.
    
    MATCH: switch( std::get<2>(*it)){
    case INTSXP:
      IF_DEBUG( Rcpp::Rcout << "Convert integer vector (" << var_name << ")"  << std::endl; );
      name_vec.push_back(var_name);
      
      column_vec1 = std::get<1>(*it);
      intvec = Rcpp::wrap( *((std::vector<int>*)column_vec1));
      column_vec2 = std::get<4>(*it);
      dblvec = Rcpp::wrap( *((std::vector<double>*)column_vec2));  
      column_vec3 = std::get<5>(*it);
      typevec = Rcpp::wrap( *((std::vector<int>*)column_vec3));  
      dbl_pos = (typevec == DBLNUM);
      if(is_false(any(dbl_pos))){
        REORDER_INTVEC( intvec, rcpp_order);
        new_df.push_back(intvec, var_name);
        IF_DEBUG( Rcpp::Rcout << "integer vector (" << var_name << ")"  << " is added to R Dataframe" << std::endl; );
      }else if(all(is_na(ifelse(dbl_pos, dblvec, NA_REAL)))){ // All the DBLNUM positions have na/nan/-nan. Return intvec.
        new_intvec = ifelse(!dbl_pos, intvec, NA_INTEGER);
        REORDER_INTVEC(new_intvec, rcpp_order);
        new_df.push_back(new_intvec, var_name);
        IF_DEBUG( Rcpp::Rcout << "integer vector (" << var_name << ")"  << " is added to R Dataframe" << std::endl; );
      }else{
        for( index = 0 ; index < typevec.size() ; ++index ){
          if(typevec(index) == INTNUM){
            if( intvec(index) == NA_INTEGER)
              dblvec(index) = NA_REAL;
            else
              dblvec(index) = (double) intvec(index);
          }else if(typevec(index) == DBLNUM){
              // dblvec is used.
          }else{
            Rcpp::Rcout << "ERROR: type_vec should have INTNUM or DBLNUM" << std::endl;
          }
        }
        REORDER_DBLVEC(dblvec, rcpp_order);
        new_df.push_back(dblvec, var_name);
        IF_DEBUG( Rcpp::Rcout << "numeric(=double) vector (" << var_name << ")"  << " is added to R Dataframe" << std::endl; );
      }
      break;
    case REALSXP:
      IF_DEBUG( Rcpp::Rcout << "Convert real(=double) vector (" << var_name << ")"  << std::endl; );
      name_vec.push_back(var_name);
      
      column_vec1 = std::get<1>(*it);
      dblvec = Rcpp::wrap( *((std::vector<double>*)column_vec1));
      column_vec2 = std::get<4>(*it);
      intvec = Rcpp::wrap( *((std::vector<int>*)column_vec2));  
      column_vec3 = std::get<5>(*it);
      typevec = Rcpp::wrap( *((std::vector<int>*)column_vec3));
      dbl_pos = (typevec == DBLNUM);
      if(is_false(any(dbl_pos))){
        REORDER_INTVEC(intvec, rcpp_order);
        new_df.push_back(intvec, var_name);
        IF_DEBUG( Rcpp::Rcout << "integer vector (" << var_name << ")"  << " is added to R Dataframe" << std::endl; );
      }else if(all(is_na(ifelse(dbl_pos, dblvec, NA_REAL)))){ // All the DBLNUM positions have na/nan/-nan. Return intvec.
        new_intvec = ifelse(!dbl_pos, intvec, NA_INTEGER);
        REORDER_INTVEC(new_intvec, rcpp_order);
        new_df.push_back(new_intvec, var_name);
        IF_DEBUG( Rcpp::Rcout << "integer vector (" << var_name << ")"  << " is added to R Dataframe" << std::endl; );
      }else{
        for( index = 0 ; index < typevec.size() ; ++index ){
          if(typevec(index) == INTNUM){
            if( intvec(index) == NA_INTEGER)
              dblvec(index) = intvec(index);
            else
              dblvec(index) = (double) intvec(index);
          }else if(typevec(index) == DBLNUM){
            // dblvec is used.
          }else{
            Rcpp::Rcout << "ERROR: type_vec should have INTNUM or DBLNUM" << std::endl;
          }
        }
        REORDER_DBLVEC(dblvec, rcpp_order);
        new_df.push_back(dblvec, var_name);
        IF_DEBUG( Rcpp::Rcout << "numeric(=double) vector (" << var_name << ")"  << " is added to R Dataframe" << std::endl; );
      }
      break;
    case STRSXP:
      IF_DEBUG( Rcpp::Rcout << "Convert string vector (" << var_name << ")"  << std::endl; );
      name_vec.push_back(var_name);
      
      column_vec1 = std::get<1>(*it);
      strvec_ori = *((std::vector<std::string* >*)column_vec1);
      column_vec2 = std::get<4>(*it);
      strvec = *((std::vector<std::string* >*)column_vec2);
      column_vec3 = std::get<5>(*it);
      updated_vec = *((std::vector<int>*)column_vec3);
      rstrvec = StringVector(strvec.size());
      
      for(idx = 0; idx < updated_vec.size(); ++idx ){
        if(updated_vec[idx] == UPDATED){
          IF_DEBUG( printf("%d : (UPDATED) Address: %p , Value: ", idx, strvec[idx]); Rcpp::Rcout << *(strvec[idx]) << std::endl; );
//          rstrvec.push_back(*(strvec[idx])) ; // Too inefficient => Deleted.
	      rstrvec[idx] = *(strvec[idx]);
        }else if(updated_vec[idx] == ORIGINAL){
//          rstrvec.push_back(*(strvec_ori[idx])) ; // Too inefficient => Deleted.
		  if(strvec_ori[idx] == NULL){
            IF_DEBUG( printf("%d : (ORIGINAL) NULL", idx); );
		    // Nothing to be assigned. The default element value of StringVector is ""
		  }else{
            IF_DEBUG( printf("%d : (ORIGINAL) Address: %p , Value: ", idx, strvec_ori[idx]); Rcpp::Rcout << *(strvec_ori[idx]) << std::endl; );
			rstrvec[idx] = *(strvec_ori[idx]);
		  }
        }else{
            IF_DEBUG( Rcpp::Rcout << "ERROR: type_vec should have UPDATED or ORIGINAL. TYPE ID: " << updated_vec[idx] << std::endl; );
          }
      }
      REORDER_STRVEC(rstrvec, rcpp_order);
      new_df.push_back(rstrvec, var_name);
      IF_DEBUG( Rcpp::Rcout << "character(=string) vector (" << var_name << ")"  << " is added to R Dataframe" << std::endl; );
      break;
    case NILSXP:
      break;
    default:
      break;
    }
  }
//  new_df.attr("names") = name_vec; // Previsouly used. Now, column name is added with push_back method.

  IF_DEBUG( Rcpp::Rcout << "Column names just after converting VEC_LIST to Rcpp vectors :" ; );
  IF_DEBUG( Rcpp::CharacterVector colname_vec = new_df.names(); for(const auto& iter : colname_vec){ Rcpp::Rcout << iter << " "; }; Rcpp::Rcout << std::endl; );
  
  Rcpp::DataFrame new_df_out = Rcpp::DataFrame::create(new_df, _["stringsAsFactors"] = false);  // This step is required to output proper data.frame. Also, stringsAsFactors attribute is required to be dataframe.
  if(order != NULL){ // This DataFrame is ordered.
      std::vector<int>* rownum_vec = (std::vector<int>*) std::get<1>(*(vec_elem_find(vl, "_n_")));
      IntegerVector rownum_rcpp_vec = Rcpp::wrap(*rownum_vec);
      REORDER_INTVEC(rownum_rcpp_vec, rcpp_order);
      new_df_out.attr("DataSailr_NewOrder") = Rcpp::LogicalVector::create(1); // TRUE in R
      new_df_out.attr("DataSailr_NewOrderVector") = rownum_rcpp_vec ; 
  }
  IF_DEBUG( Rcpp::Rcout << "Column names after converting to DataFrame: " ; );
  IF_DEBUG( Rcpp::CharacterVector colname_vec2 = new_df_out.names(); for(const auto& iter : colname_vec2){ Rcpp::Rcout << iter << " "; }; Rcpp::Rcout << std::endl; );

  return new_df_out;
}

int
vec_list_push_cloned_row(VEC_LIST* vl, int from_idx)
{
  // To access Tuple elemnet
  void* column_vec1;
  void* column_vec2;
  void* column_vec3;

  // For integer & double columns
  std::vector<int> *intvec;
  std::vector<double> *dblvec;
  std::vector<int> *typevec;

  // For string column
  std::vector<std::string* > *strvec;
  std::vector<std::string* > *strvec_ori;
  std::vector<int> *updated_vec;
  std::string* temp_str;

  // For null column
  std::vector<int>* nullvec;
  
  for(auto it = vl->begin(); it != vl->end(); ++it ){
    switch( std::get<2>(*it)){
    case INTSXP:
      IF_DEBUG( Rcpp::Rcout << "Exteniding integer vector (" << ((char*) std::get<0>(*it)) << ")"  << std::endl; );
      column_vec1 = std::get<1>(*it);
      intvec = ((std::vector<int>*)column_vec1);
      column_vec2 = std::get<4>(*it);
      dblvec = ((std::vector<double>*)column_vec2);
      column_vec3 = std::get<5>(*it);
      typevec = ((std::vector<int>*)column_vec3);

      intvec->push_back(intvec->operator[](from_idx));
      dblvec->push_back(dblvec->operator[](from_idx));
      typevec->push_back(typevec->operator[](from_idx));

      std::get<3>(*it) = std::get<3>(*it) + 1; // increment size.
      break;
    case REALSXP:
      IF_DEBUG( Rcpp::Rcout << "Exteniding doubler vector (" << ((char*) std::get<0>(*it)) << ")"  << std::endl; );
      column_vec1 = std::get<1>(*it);
      dblvec = ((std::vector<double>*)column_vec1);
      column_vec2 = std::get<4>(*it);
      intvec = ((std::vector<int>*)column_vec2);
      column_vec3 = std::get<5>(*it);
      typevec = ((std::vector<int>*)column_vec3);

      dblvec->push_back(dblvec->operator[](from_idx));
      intvec->push_back(intvec->operator[](from_idx));
      typevec->push_back(typevec->operator[](from_idx));

      std::get<3>(*it) = std::get<3>(*it) + 1; // increment size.
      break;
    case STRSXP:
      IF_DEBUG( Rcpp::Rcout << "Exteniding string vector (" << ((char*) std::get<0>(*it)) << ")"  << std::endl; );
      
      column_vec1 = std::get<1>(*it);
      strvec_ori = (std::vector<std::string* >*)column_vec1;
      column_vec2 = std::get<4>(*it);
      strvec = (std::vector<std::string* >*)column_vec2;
      column_vec3 = std::get<5>(*it);
      updated_vec = (std::vector<int>*)column_vec3;

      if(updated_vec->operator[](from_idx) == UPDATED){
        temp_str = (std::string*) strvec->operator[](from_idx);
        strvec_ori->push_back(new std::string(*temp_str));
      }else{
        if( strvec->operator[](from_idx) != NULL){
          temp_str = (std::string*) strvec->operator[](from_idx);
          strvec_ori->push_back(new std::string(*temp_str));
        }else{
          strvec_ori->push_back(NULL);
        }
      }
      strvec->push_back(NULL);
      updated_vec->push_back(ORIGINAL);

      std::get<3>(*it) = std::get<3>(*it) + 1; // increment size.
      break;
    case NILSXP:
      nullvec->push_back(0);

      std::get<3>(*it) = std::get<3>(*it) + 1; // increment size.
      break;
    default:
      break;
    }
  }
  return (vec_list_nrows(vl) - 1); // Index for the last row (= size - 1)
}

void
vec_list_show_summary(VEC_LIST* vlist )
{
	char* colname;
	int sexptype;
	for(auto iter = vlist->begin(); iter != vlist->end(); ++iter){
		colname = std::get<0>(*iter);
		sexptype = std::get<2>(*iter);
		switch( sexptype ){
		case INTSXP:
			Rcpp::Rcout << "\"" << colname << "\":integer " ; 
			break;
		case REALSXP:
			Rcpp::Rcout << "\"" << colname << "\":double " ; 
			break;
		case STRSXP:
			Rcpp::Rcout << "\"" << colname << "\":string " ;
			break;
		case NILSXP:
			Rcpp::Rcout << "\"" << colname << "\":null " ; 
			break;
		default:
			Rcpp::Rcout << "\"" << colname << "\":unknown " ; 
			break;
		}
	}
	Rcpp::Rcout << std::endl;
}


bool
cstring_exists_in_charactervector(char* var_name, CharacterVector var_vector){
  std::string string_cpp; 
  bool result = false;
  for(auto iter = var_vector.begin(); iter != var_vector.end() ; ++iter){
    string_cpp = *iter ;
    if (string_cpp.compare(std::string(var_name)) == 0){
      result = true;
    }
  }
  return result;
}

int
vec_list_nrows(VEC_LIST* vec_list)
{
  int size = std::get<3>( vec_list->operator[](0) );
  return size;
}

std::vector<std::string>
vec_list_extract_nil_vars( VEC_LIST* vec_list )
{
	std::vector<std::string> nil_vars ;
	VEC_ELEM vec_elem;
	char* nil_var_name;
	for( auto iter = vec_list->begin(); iter != vec_list->end(); ++iter ){
		vec_elem = *iter;
		if (vec_elem_type_of(&vec_elem) == NILSXP){
			nil_var_name = vec_elem_name_of(&vec_elem);
			nil_vars.push_back(std::string(nil_var_name));
		}
	}
	return nil_vars;
}

int
vec_elem_remove_nil(VEC_LIST* vl, char* nil_var_name)
{
  char* vec_elem_name;
  std::vector<int>* nilvec;
  for( auto it = vl->begin(); it != vl->end(); ++it){
    vec_elem_name = std::get<0>(*it);
    if( strcmp(vec_elem_name, nil_var_name) == 0){
      nilvec = (std::vector<int>*)std::get<1>(*it);
      free(vec_elem_name);
      delete(nilvec);
      vl->erase(it);
      return 0;
    }
  } 
  return 1;
}

std::vector<void*>
update_vec_elem_with_new_type(VEC_LIST* vec_list, char* nil_var_name, char new_type)
{
  int size = ((std::vector<int>*)std::get<1>(*(vec_elem_find(vec_list, nil_var_name))))->size();
  
  void* vec_ori; 
  void* vec_new;
  void* vec_type;
  VEC_ELEM* added_elem ;
  std::vector<void*> new_vec_info = std::vector<void*>(3);
  
  if(new_type == 'i'){
    vec_elem_remove_nil(vec_list, nil_var_name );
    vec_list_add_int_vec( vec_list, nil_var_name, NULL , size );
    added_elem = vec_elem_find(vec_list, nil_var_name);
    vec_ori = std::get<1>(*added_elem);
    vec_new = std::get<4>(*added_elem);
    vec_type = std::get<5>(*added_elem);
    new_vec_info[0] = ((void*) vec_ori);
    new_vec_info[1] = ((void*) vec_new);
    new_vec_info[2] = ((void*) vec_type);
    return new_vec_info ;
  } else if(new_type == 'd'){
    vec_elem_remove_nil(vec_list, nil_var_name );
    vec_list_add_double_vec( vec_list, nil_var_name, NULL , size );
    added_elem = vec_elem_find(vec_list, nil_var_name);
    vec_ori = std::get<1>(*added_elem);
    vec_new = std::get<4>(*added_elem);
    vec_type = std::get<5>(*added_elem);
    new_vec_info[0] = ((void*) vec_ori);
    new_vec_info[1] = ((void*) vec_new);
    new_vec_info[2] = ((void*) vec_type);
    return new_vec_info ;
  } else if(new_type == 's'){
    vec_elem_remove_nil(vec_list, nil_var_name );
    vec_list_add_string_vec(vec_list, nil_var_name, NULL, size);
    added_elem = vec_elem_find(vec_list, nil_var_name);
    vec_ori = std::get<1>(*added_elem);
    vec_new = std::get<4>(*added_elem);
    vec_type = std::get<5>(*added_elem);
    new_vec_info[0] = ((void*) vec_ori);
    new_vec_info[1] = ((void*) vec_new);
    new_vec_info[2] = ((void*) vec_type);
    return new_vec_info ;
  } else if(new_type == 'r'){
    new_vec_info[0] = NULL;
    new_vec_info[1] = NULL;
    new_vec_info[2] = NULL;
    return new_vec_info ;
  } else{
  	if(new_type == 'n'){ 
        IF_DEBUG( printf("WARNING: The variable, %s , on ptr_table does not seem to be updated.\n", nil_var_name ); );
    } else if(new_type == 'x'){ 
        IF_DEBUG( printf("ERROR: The variable, %s , on ptr_table is updated to unknown type. TYPE ID: %d \n", nil_var_name, new_type ); );
    } else { 
        IF_DEBUG( printf("ERROR: The variable, %s , on ptr_table is updated to unknown type.  TYPE ID: %d \n", nil_var_name, new_type ); );
	}
    new_vec_info[0] = NULL;
    new_vec_info[1] = NULL;
    new_vec_info[2] = NULL;
    return new_vec_info ;    	
  }
}


VEC_ELEM*
vec_elem_find(VEC_LIST* vl, const char* var_name)
{
  char* vec_elem_name;
  for( auto it = vl->begin(); it != vl->end(); ++it){
    vec_elem_name = std::get<0>(*it);
    if( strcmp(vec_elem_name, var_name) == 0){
      return (&(*it));
    }
  } 
  return NULL;
}

SXPTYPE
vec_elem_type_of(VEC_ELEM* vec_elem)
{
  SXPTYPE sxptype = std::get<2>(*vec_elem);
  return sxptype;
}

char*
vec_elem_name_of(VEC_ELEM* vec_elem)
{
    char* var_name = std::get<0>(*vec_elem);
    return var_name;
  }


/* --------------------------------------------------------- */


int
update_sailr_ptr_table ( ptr_table_object* table, char** vars, int var_num, VEC_LIST* vl, int row_idx )
{
	int var_idx ;
	char* var_name;
	std::vector<int>* intvec;
	std::vector<double>* dblvec;
	std::vector<int>* typevec;
	std::vector<std::string* >* strvec;

	int* int_currpos;
	double* dbl_currpos;
	int row_type;
	std::string* pstr;

	for( var_idx = 0; var_idx < var_num; ++var_idx){
		var_name = vars[var_idx];
		if(var_name == NULL){
		    IF_DEBUG( Rcpp::Rcout << "NULL is ignored while updating ptr_table." << std::endl;);
		    continue;
		}else{
		    IF_DEBUG(Rcpp::Rcout << "Update var(" << var_name  << ") " ; );
		}
		VEC_ELEM* vec_elem = vec_elem_find (vl, var_name);
	  
		switch(vec_elem_type_of(vec_elem)){
		case INTSXP :
			intvec = (std::vector<int>*) std::get<1>(*vec_elem);
			int_currpos = (intvec->data()) + row_idx;
			dblvec = (std::vector<double>*) std::get<4>(*vec_elem);
			dbl_currpos = (dblvec->data()) + row_idx;
			typevec = (std::vector<int>*) std::get<5>(*vec_elem);
			row_type = typevec->operator[](row_idx);
			if(row_type == INTNUM){
				sailr_ptr_table_create_int_from_ptr(&table, var_name , &int_currpos, &dbl_currpos);
			}else if(row_type == DBLNUM){
				sailr_ptr_table_create_double_from_ptr(&table, var_name , &dbl_currpos, &int_currpos);
			}
		break;
		case REALSXP:
			dblvec = (std::vector<double>*) std::get<1>(*vec_elem);
			dbl_currpos = (dblvec->data()) + row_idx;
			intvec = (std::vector<int>*) std::get<4>(*vec_elem);
			int_currpos = (intvec->data()) + row_idx;
			typevec = (std::vector<int>*) std::get<5>(*vec_elem);
			row_type = typevec->operator[](row_idx);
			if(row_type == INTNUM){
				sailr_ptr_table_create_int_from_ptr(&table, var_name , &int_currpos, &dbl_currpos);
			}else if(row_type == DBLNUM){
				sailr_ptr_table_create_double_from_ptr(&table, var_name , &dbl_currpos, &int_currpos);
			}
		break;
		case STRSXP:
			strvec = (std::vector<std::string* >*) std::get<1>(*vec_elem);
			pstr = strvec->operator[](row_idx);
			if( pstr != NULL){
				sailr_ptr_table_create_string_from_cstring(&table, var_name , pstr->c_str());
			}else{
				sailr_ptr_table_create_string_from_cstring(&table, var_name , "");
			}
		break;
		case NILSXP:
			sailr_ptr_table_create_null(&table, var_name );
		break;
		default:
			Rcpp::Rcout << "ERROR: This type of column is not supported. " << std::endl;
		break; 
		}
	}
	return 0;
}



int
update_sailr_vec_list(VEC_LIST*  vl, std::vector<std::string> vars, ptr_table_object* table, int row_idx ) 
{
  char* var_name;
  std::vector<std::string* >* new_pvec_pstr;
  const char* cstr_on_table;
  unsigned int var_idx;
  
  std::vector<int>* type_vec;
  std::vector<int>* updated_vec;

  
  for( var_idx = 0; var_idx < vars.size() ; ++var_idx){
    var_name = (char*) vars[var_idx].c_str();
    VEC_ELEM* vec_elem = vec_elem_find (vl, var_name);

    switch(vec_elem_type_of(vec_elem)){
    case INTSXP :
      // if ptr_table has PTR_DBL for the corresponding variable, (sailr_ptr_table_get_type(table, nil_var_name);)
      // type should be set to DBLNUM
      if(sailr_ptr_table_get_type(&table, var_name) == 'd'){
        IF_DEBUG( Rcpp::Rcout << "type mismatch (INTSXP on vec_elem; PTR_DBL on ptr_table)" << std::endl; );
        type_vec = (std::vector<int>*)std::get<5>(*vec_elem);
        type_vec->operator[](row_idx) = DBLNUM;
      } else if(sailr_ptr_table_get_type(&table, var_name) == 'i'){
        type_vec = (std::vector<int>*)std::get<5>(*vec_elem);
        type_vec->operator[](row_idx) = INTNUM;     
      }
      break;
    case REALSXP:
      // if ptr_table has PTR_INT for the corresponding variable, (sailr_ptr_record_get_type(table, nil_var_name);)
      // type should be set to DBLNUM
      if(sailr_ptr_table_get_type(&table, var_name) == 'i'){
        IF_DEBUG( Rcpp::Rcout << "type mismatch (DBLSXP on vec_elem; PTR_INT on ptr_table)" << std::endl; );
        type_vec = (std::vector<int>*)std::get<5>(*vec_elem);
        type_vec->operator[](row_idx) = INTNUM;
      } else if(sailr_ptr_table_get_type(&table, var_name) == 'd'){
        type_vec = (std::vector<int>*)std::get<5>(*vec_elem);
        type_vec->operator[](row_idx) = DBLNUM;     
      }
      break;
    case STRSXP:
      new_pvec_pstr = (std::vector<std::string*>*) std::get<4>(*vec_elem);
      cstr_on_table =  sailr_ptr_table_read_string(&table, (char*)var_name);
      new_pvec_pstr->operator[](row_idx) = new std::string(cstr_on_table);
      updated_vec = (std::vector<int>*)std::get<5>(*vec_elem);
      updated_vec->operator[](row_idx) = UPDATED;
      break;
    default:
      IF_DEBUG( Rcpp::Rcout << var_name << ": VEC LIST needs not be updated for this type" << std::endl; );
      break; 
    }
  }
  return 0;
}

void
show_sailr_vec_list_nth(VEC_LIST*  vl, int nth )
{
  VEC_ELEM* vec_elem;
  unsigned int elem_idx;
  char* elem_name;
  std::vector<int>* int_vec;
  std::vector<double>* double_vec;
  std::vector<int>* type_vec;
  std::vector<std::string* >* ori_str_vec;
  std::vector<std::string* >* new_str_vec;
  std::string ori_str;
  std::string new_str;
  for(elem_idx = 0; elem_idx < vl->size() ; ++elem_idx){
    vec_elem = &( vl->operator[](elem_idx));
    switch(vec_elem_type_of(vec_elem)){
    case INTSXP :
      elem_name = (char*)std::get<0>(*vec_elem);
      int_vec = (std::vector<int>*)std::get<1>(*vec_elem);
      double_vec = (std::vector<double>*)std::get<4>(*vec_elem);
      type_vec = (std::vector<int>*)std::get<5>(*vec_elem);
      Rcpp::Rcout << elem_name << ":INTSXP" << " " ;
      Rcpp::Rcout << int_vec->operator[](nth) << " |" ;
      Rcpp::Rcout << double_vec->operator[](nth) << " |" ;
      Rcpp::Rcout << type_vec->operator[](nth) << std::endl;
      break;
    case REALSXP:
      elem_name = (char*)std::get<0>(*vec_elem);
      double_vec = (std::vector<double>*)std::get<1>(*vec_elem);
      int_vec = (std::vector<int>*)std::get<4>(*vec_elem);
      type_vec = (std::vector<int>*)std::get<5>(*vec_elem);
      Rcpp::Rcout << elem_name << ":REALSXP"  << " " ;
      Rcpp::Rcout << double_vec->operator[](nth) << " |" ;
      Rcpp::Rcout << int_vec->operator[](nth) << " |" ;
      Rcpp::Rcout <<  type_vec->operator[](nth) << std::endl;
      break;
    case STRSXP:
      elem_name = (char*)std::get<0>(*vec_elem);
      ori_str_vec = (std::vector<std::string*>*)std::get<1>(*vec_elem);
      new_str_vec = (std::vector<std::string*>*)std::get<4>(*vec_elem);
      if(ori_str_vec->operator[](nth) == NULL){
        ori_str = std::string("");
      }else{
        ori_str = *(ori_str_vec->operator[](nth));
      }
      if(new_str_vec->operator[](nth) == NULL){
        new_str = std::string("");
      }else{
        new_str = *(new_str_vec->operator[](nth));
      }
      type_vec = (std::vector<int>*)std::get<5>(*vec_elem);
      Rcpp::Rcout << elem_name << ":STRSXP" << " " ;
      Rcpp::Rcout << ori_str << "(:ori) " ;
      Rcpp::Rcout << new_str << "(:new) " ;
      Rcpp::Rcout << type_vec->operator[](nth) << std::endl;
      break;
    case NILSXP:
      elem_name = (char*)std::get<0>(*vec_elem);
      Rcpp::Rcout << elem_name << ":NILSXP" << std::endl;
      break; 
    default:
      elem_name = (char*)std::get<0>(*vec_elem);
      Rcpp::Rcout << elem_name << ":OTHER TYPES" << std::endl;
    break; 
    }
  }
}

void
ShowVecList(VEC_LIST* vl, unsigned int to_nth_row )
{
	unsigned int idx;
	unsigned int max_row = vec_list_nrows(vl);
	unsigned int row_to_show;

	if( to_nth_row < max_row ){
		row_to_show = to_nth_row;
	}else{
		row_to_show = max_row;
	}

	for( idx = 0 ; idx < row_to_show; ++idx){
		Rcpp::Rcout << "Row " << (idx + 1 ) << std::endl;
		show_sailr_vec_list_nth(vl, idx);
	}
}

  

/* --------------------------------------------------------- */

// [[Rcpp::export(.data_sailr_cpp_execute)]]
Rcpp::DataFrame
data_sailr_cpp_execute( Rcpp::CharacterVector rchars, Rcpp::DataFrame df)
{
	std::string code = Rcpp::as<std::string> (rchars);
  IF_DEBUG( Rcpp::Rcout << code << std::endl; );

	// Initializing pointer table. Pointer table is implemented using UThash in C.
	ptr_table_object* table = sailr_ptr_table_init() ;
  IF_DEBUG( printf("ptr_table is initialized (table's pointer is %p ) \n", table); );

	// Initializing parser_state, which stores TreeNode*.
	IF_DEBUG( Rcpp::Rcout << "Constructing parse tree."  << std::endl ; );
	int parse_result;
	parser_state_object* ps = sailr_new_parser_state ((char*)"Code from R", table);
	parse_result = sailr_run_parser( code.c_str(), ps );  // Now ps now holds AST tree and ptr_table!!

	// Check whether parsing succeeded or not. When fails, free memory and stop this function.
	if(parse_result == 0 ){
	  IF_DEBUG( Rcpp::Rcout << "Success: sailr script is successfully parsed" << std::endl; );
	}else {
	  /* Free memory */
	  IF_DEBUG( Rcpp::Rcout << "Free parser tree" << std::endl; );
	  sailr_tree_free(ps);
	  IF_DEBUG( Rcpp::Rcout << "Free pointer table" << std::endl; );
	  sailr_ptr_table_del_all(&table);
	  IF_DEBUG( Rcpp::Rcout << "Free parser state object" << std::endl; );
	  sailr_parser_state_free(ps);

	  if(parse_result == 1){
	    Rcpp::stop( "sailr script syntax error (code: %d)\n", parse_result );
	  }else if(parse_result == 2){
	    Rcpp::stop( "memory exhausted during parsing (code: %d)\n", parse_result );
	  }else{
	    Rcpp::stop( "yyparse returned unknown code (code: %d)\n" , parse_result );
	  }
	}

	// Show parse tree
	IF_DEBUG( Rcpp::Rcout << "Show parse tree"  << std::endl ; );
	IF_DEBUG( sailr_tree_dump( ps ); std::cout << std::flush; );

	// Show pointer table after constructing parse tree
	IF_DEBUG( Rcpp::Rcout << "Show pointer table. (At this point, anonymous STRING and REGEXP's should be already added.)" << std::endl; );
	IF_DEBUG( sailr_ptr_table_show_all(&table);  std::cout << std::flush; );

	// Show variables on source code
	IF_DEBUG( Rcpp::Rcout << "Variable names on RHS and LHS were collected during constructing parse tree." << std::endl; );
	IF_DEBUG( Rcpp::Rcout << "Variable names that appear in source code ( = LHS + RHS):" << std::endl; );
	std::vector<std::string> vars;
	char** var_array = sailr_varnames(ps);
	int var_num = sailr_varnames_num(ps);
	int var_idx;
	for(var_idx=0; var_idx < var_num; var_idx++){
	  IF_DEBUG( Rcpp::Rcout << " " << var_array[var_idx] ; );
	  vars.push_back(std::string(var_array[var_idx]));
	}
	IF_DEBUG( Rcpp::Rcout << std::endl; );

	IF_DEBUG( Rcpp::Rcout << "All the LHS varialbes need ptr_record. "; );
	IF_DEBUG( Rcpp::Rcout << "LHS varialbes that exist on dataframe: (type) known varaibles(=non-nil on VEC_LIST). "; );
	IF_DEBUG( Rcpp::Rcout << "LHS varialbes that do not exist on dataframe: (type) unknown varaibles (=nil on VEC_LIST), which types may or may not be defined during execution. " << std::endl; );
	IF_DEBUG( Rcpp::Rcout << "Variable names on LHS:"  << std::endl ; );
	std::vector<std::string> lhs_vars;
	char** lhs_var_array = sailr_lhs_varnames(ps);
	int lhs_var_num = sailr_lhs_varnames_num(ps);
	if(lhs_var_num > 0){
	int lhs_var_idx;
	  for(lhs_var_idx=0; lhs_var_idx < lhs_var_num; lhs_var_idx++){
	    IF_DEBUG( Rcpp::Rcout << " \"" << lhs_var_array[lhs_var_idx] << "\""; );
	    lhs_vars.push_back(std::string(lhs_var_array[lhs_var_idx]));
	  }
	}
	IF_DEBUG( Rcpp::Rcout << std::endl; );

	IF_DEBUG( Rcpp::Rcout << "All the RHS variables should appear on LHS or should exist in dataframe column names" << std::endl; );
	IF_DEBUG( Rcpp::Rcout << "Variable names on RHS:"  << std::endl ; );
	std::vector<std::string> rhs_vars;
	char** rhs_var_array = sailr_rhs_varnames(ps);
	int rhs_var_num = sailr_rhs_varnames_num(ps);
	if(rhs_var_num > 0){
	  int rhs_var_idx;
	  for(rhs_var_idx=0; rhs_var_idx < rhs_var_num; rhs_var_idx++){
	    IF_DEBUG( Rcpp::Rcout << " \"" << rhs_var_array[rhs_var_idx] << "\""; );
	    rhs_vars.push_back(std::string(rhs_var_array[rhs_var_idx]));
	  }
	}
	IF_DEBUG( Rcpp::Rcout << std::endl; );

	// Specify variable names for special purposes
	#define SERVICE_VAR_NUM 2
	char* service_var_array[ SERVICE_VAR_NUM ];
	std::string service_var_one( "_n_" );
	std::string service_var_two( "_discard_" );
	service_var_array[0] = &service_var_one[0];
	service_var_array[1] = &service_var_two[0];

	std::vector<std::string> service_vars;
	if(SERVICE_VAR_NUM > 0){
	  int srv_var_idx;
	  for(srv_var_idx=0; srv_var_idx < SERVICE_VAR_NUM; srv_var_idx++){
	    IF_DEBUG( Rcpp::Rcout << " \"" << service_var_array[srv_var_idx] << "\""; );
	    service_vars.push_back(std::string(service_var_array[srv_var_idx]));
	  }
	}
	
	// Convert R dataframe into C++ Vec_LIST
	IF_DEBUG( Rcpp::Rcout << "Convert Rcpp DataFame to C++ VEC_LIST" << std::endl; );

	int conversion_error = 0;
	VEC_LIST* vec_list = ConvertDataFrame(df, var_array, var_num, lhs_var_array, lhs_var_num , service_var_array, SERVICE_VAR_NUM, &conversion_error );
	if( conversion_error != 0 ){
		// Free variable names
		sailr_varnames_free(var_array, var_num);
		sailr_varnames_free(lhs_var_array, lhs_var_num);
		sailr_varnames_free(rhs_var_array, rhs_var_num);

		// Free parser tree
		sailr_tree_free(ps);

		// Free pointer table
		sailr_ptr_table_del_all(&table);

		// Free parser_state object
		sailr_parser_state_free(ps);

		// Free VEC_LIST
		vec_list_free(vec_list);
		Rcpp::stop("Stopped: check variable names again");
	}

	// Extract nil variables from vec_list
	IF_DEBUG( Rcpp::Rcout << "nil vars are (type) unknown (LHS) variables that do not exist in dataframe. " << std::endl; );
	IF_DEBUG( Rcpp::Rcout << "Nil Variable names on VEC_LIST:" ; );

	std::vector<std::string> nil_vars;
	nil_vars = vec_list_extract_nil_vars( vec_list );
	for(auto iter = nil_vars.begin(); iter != nil_vars.end(); ++iter){
		std::string nil_var_name = *iter;
		IF_DEBUG( Rcpp::Rcout << " \"" << nil_var_name << "\""; );
	}
	IF_DEBUG( Rcpp::Rcout << std::endl; );

	// Variables to update = var_array + service_var_array
	// This is used to update ptr_table
	char** vars_to_update;
	int vars_to_update_size = var_num + SERVICE_VAR_NUM;
	IF_DEBUG(Rcpp::Rcout << "Variables to update = var_array + service_var_array " << std::endl; );
	vars_to_update = (char**) malloc( sizeof(char*) * vars_to_update_size );
	int var_update_idx = 0;
	while(var_update_idx < var_num){
	    vars_to_update[var_update_idx] = var_array[var_update_idx];
	    ++var_update_idx;
	}
	IF_DEBUG( Rcpp::Rcout << "Showing vars_to_update that is copied from var_array: " << std::endl; );
	IF_DEBUG( int show_var_idx = 0; while(show_var_idx < var_num){ Rcpp::Rcout << vars_to_update[show_var_idx] << " " ; ++show_var_idx;}; Rcpp::Rcout << std::endl; );
	int srv_var_idx = 0;
	while(srv_var_idx < SERVICE_VAR_NUM){
	    if(! cstring_exists_in_pchar_array( service_var_array[srv_var_idx] , vars_to_update, 0, var_num - 1)){
		vars_to_update[var_update_idx] =  service_var_array[srv_var_idx];
		IF_DEBUG(Rcpp::Rcout << service_var_array[srv_var_idx] << ": service var does not exist in source var names. added to vars_to_update." << std::endl; );
	    }else{
		vars_to_update[var_update_idx] =  NULL;
		IF_DEBUG(Rcpp::Rcout << service_var_array[srv_var_idx] << ": service var already exists in source var names. vars_to_update is updated with NULL." << std::endl; );
	    }
	    ++var_update_idx;
	    ++srv_var_idx;
	}
	IF_DEBUG( Rcpp::Rcout << "vars_to_update initialized." << std::endl; );

	// For int and double, pass the address of each element of VEC_LIST to ptr_table.
	// The values VEC_LIST points to are updated automatically after each calculation. 
	//
	// String variables are dealt differently from int and double vectors
	// because the data size of the result is different from their inputs. We may need more space to store result.
	// In this implementation, string variables in VEC_LIST store results in a new vector.
	//

	// Initialize ptr_table using vec_list's vectors' first element.
	update_sailr_ptr_table( table, vars_to_update, vars_to_update_size, vec_list, 0 );

	// Collect regular expression ptr_record's & string ptr_record's
	// Regular expression literals on ptr_table need to be reset for every row, but need not be deleted.
	// Regular expressions assigned to some variable need to be deleted. (Next row, it may have different regular expression)
	// String object literals can stay thorough executions.
	// String objects assigned to some varialbe need to be freed (delted) for every row.
	std::vector<ptr_record_object*> rexp_records_for_literals;
	std::vector<ptr_record_object*> rexp_records_for_vars;
	std::vector<ptr_record_object*> str_records_for_vars;
	ptr_record_object* curr_pr;
	curr_pr = sailr_ptr_table_first_record( &table );
	while( curr_pr != NULL ){
		if( sailr_ptr_record_get_type(curr_pr) == 's'){
			if( sailr_ptr_record_is_anonym( curr_pr ) != 1 ){  // not anonym (=not literal)
				str_records_for_vars.push_back(curr_pr);
			}
		}
		if( sailr_ptr_record_get_type(curr_pr) == 'r'){
			if( sailr_ptr_record_is_anonym( curr_pr ) != 1 ){  // not anonym (=not literal)
				rexp_records_for_vars.push_back(curr_pr);
			}else{  // anonym (=literal)
				rexp_records_for_literals.push_back(curr_pr);
			}
		}
		curr_pr = sailr_ptr_record_next(curr_pr);
	}

	// Free string objects on ptr_table that are used just for ptr_table initialization.
	for(auto str_record_iter = str_records_for_vars.begin(); str_record_iter != str_records_for_vars.end(); ++str_record_iter){
		sailr_ptr_record_free_objects( *str_record_iter );
	}
  
	IF_DEBUG( Rcpp::Rcout << "Show ptr table before calculation." << std::endl; );
	IF_DEBUG( sailr_ptr_table_show_all(&table); );

	// ptr_table should have variable names and their types here. 
	// Inside gen_code(). Inside gen_code_ident(), variable types are looked up using ptr_table. 
	vm_inst_object* inst_list = sailr_gen_code( ps, table); 
	IF_DEBUG( Rcout << "VM code list is generated" << std::endl; );
	IF_DEBUG( sailr_vm_inst_list_show_all ( inst_list ); ); 

	// VM instruction list => VM instruction array.
	vm_inst_object* vmcode = sailr_vm_inst_list_to_code(inst_list); 
	int vmcode_size = sailr_vm_inst_list_size( inst_list);
	
	// Variable for virtual machine
	vm_stack_object* vmstack;

	// External Functions
	ext_func_hash_object* extfunchash; 
	extfunchash = sailr_ext_func_hash_init();
	sailr_ext_func_hash_add( &extfunchash, "push!", 0, &sailr_external_push_row );
	sailr_ext_func_hash_add( &extfunchash, "discard!", 0, &sailr_external_discard_row );
	sailr_ext_func_hash_add( &extfunchash, "MYPRINTLN", 1, &sailr_external_println );
	
	// Execute code 
	IF_DEBUG( Rcpp::Rcout << "VM code is to be executed on Virtual machine" << std::endl; );
	
	int row_idx ;
	int ori_row_idx;
	int pushed_row_idx;
	int num_of_rows = vec_list_nrows(vec_list) ;
	char new_type;
	std::vector<void*> new_vec_info;
	void** new_ptr;
	int vm_exec_code_result = 1; // success
	bool vm_exec_success = true;
	bool sort_required = false;
	bool discard_flag = false;


	IF_DEBUG( Rcpp::Rcout << "Calculation started. (" << "Num of rows to be processed: " << num_of_rows << ")" << std::endl; );

	for( ori_row_idx = 0 ; ori_row_idx < num_of_rows ; ++ori_row_idx){  // Conduct for each row of dataframe.

		loop_when_suspended: ;

		if( vm_exec_code_result ==  1){ // if previously success
			row_idx = ori_row_idx;
			vmstack = sailr_vm_stack_init(); // In this case, vm stack must have been freed, and should be created newly.
		} else if( vm_exec_code_result == 2){ // if previously suspended
			if( extfunchash != NULL && ( strcmp( sailr_ext_func_hash_get_last_executed(&extfunchash), "push!") == 0 )){
				row_idx = pushed_row_idx;  // Result needs to be pushed to a new row.
			}
		}
		IF_DEBUG( if(row_idx == 0){ Rcpp::Rcout << "Virtual machine is generated for processing the first row.\n" << std::flush;} );

		if( extfunchash != NULL){
			sailr_ext_func_hash_reset_last_executed( &extfunchash );
		}

		update_sailr_ptr_table( table, vars_to_update, vars_to_update_size, vec_list, row_idx );
		IF_DEBUG( if(row_idx == 0){  Rcpp::Rcout << "ptr_table is updated.\n" << std::flush;}) ;

		if( vm_exec_code_result ==  1){ // if previously success
			vm_exec_code_result = sailr_vm_exec_code(vmcode, vmcode_size , table , vmstack, extfunchash);
		} else if( vm_exec_code_result == 2 ){ // if previously suspended
			vm_exec_code_result = sailr_vm_resume_code(vmcode, vmcode_size , sailr_vm_stack_get_code_position(vmstack), table , vmstack, extfunchash);
		}
		IF_DEBUG( Rcpp::Rcout << "VM code execution is finished (vm_exec_code_result: " << vm_exec_code_result << ")\n" << std::flush; );

		if( vm_exec_code_result == 0 ){
			vm_exec_success = false;  // Runtime error
			goto finalize;
		}else if(vm_exec_code_result == 2){
			// discard!() ends the current row processing, but to obtain the last executed function name (which is held by vm), the vm should not be destroyed.
			// For this purpose, the vm_exec_code_result is set 2 (suspend).
			// When discard!() is called, the following code is executed, and vm_exec_code_result needs to be set 1 at the end to finish the current row.
			if(extfunchash != NULL && sailr_ext_func_hash_get_last_executed(&extfunchash) != NULL && ( strcmp( sailr_ext_func_hash_get_last_executed(&extfunchash), "discard!") == 0 )){
				IF_DEBUG(Rcpp::Rcout << "discard!() sets _discard_ to be 1"  << std::endl;);
				sailr_ptr_table_update_int(&table, "_discard_", 1 ); // 0: not discard, 1: discard
				discard_flag = true;
				sailr_ext_vm_stack_end(vmstack);
				vm_exec_code_result = 1;
			}
		}
		
		IF_DEBUG( Rcpp::Rcout << "VM execution finished." << std::endl; );
		IF_DEBUG( if(row_idx == 0){ Rcpp::Rcout << "Showing ptr_table just after the first row calculation." << std::flush;} );
		IF_DEBUG( if(row_idx == 0){ sailr_ptr_table_show_all(&table); std::cout << std::flush;} );

#ifdef SAILR_NULL_UPDATED
		int null_updated_val = sailr_ptr_table_info_get_null_updated(&table);
		int extract_str_dbl_int_bits = (1<<0)|(1<<1)|(1<<2); // 0b00000111 can also be equivalent after C++14
		if( null_updated_val & extract_str_dbl_int_bits){ // if ( null is updated to int, doulbe or string )
		IF_DEBUG( Rcpp::Rcout << "Null updated" << std::endl; );
#endif
		for(auto iter_nil = nil_vars.begin(); iter_nil != nil_vars.end() ;  )
		{
		  std::string nil_cpp_name = *iter_nil;
		  char* nil_var_name = (char*) nil_cpp_name.c_str();
		  if(sailr_ptr_record_is_ptr_null(&table, nil_var_name) == 0 ){
		    IF_DEBUG( Rcpp::Rcout << "\"" << nil_var_name << "\"" << " becomes non null on ptr_table (= type defined). " << std::endl; );
		    // 1. obtain TYPE  **OK**
		    new_type = sailr_ptr_table_get_type(&table, nil_var_name);  // OK
		    // 2. update null with new VEC_ELEM corresponding to the TYPE.  
        	new_vec_info = update_vec_elem_with_new_type(vec_list, nil_var_name , new_type);  // OK
		    // 3. obtain pointer
		    new_ptr = sailr_ptr_table_get_pptr(&table, nil_var_name);  // OK
		    // 4. copy pointer value to VEC_ELEM
			if(new_type == 'i'){
		    	IF_DEBUG( Rcpp::Rcout << "Its int value on ptr_table is copied into new vector in VEC_LIST." << std::endl;);
				((std::vector<int>*) new_vec_info[0])->operator[](row_idx) = *((int*) *new_ptr) ;
				((std::vector<int>*) new_vec_info[2])->operator[](row_idx) = INTNUM;
				sailr_ptr_table_free_objects(&table, nil_var_name);
		    }else if (new_type == 'd'){
				IF_DEBUG( Rcpp::Rcout << "Its double value on ptr_table is copied into new vector in VEC_LIST." << std::endl;);
				((std::vector<double>*) new_vec_info[0])->operator[](row_idx) = *((double*) *new_ptr);
				((std::vector<int>*) new_vec_info[2])->operator[](row_idx) = DBLNUM;
				sailr_ptr_table_free_objects(&table, nil_var_name);
		    }else if (new_type == 's'){
				IF_DEBUG( Rcpp::Rcout << "New std::string same as the one on ptr_table is created. (Values are not assigned yet)" << std::endl;);
				str_records_for_vars.push_back( sailr_ptr_table_find( &table, nil_var_name ) );
		    }else if (new_type == 'r'){
				IF_DEBUG( Rcpp::Rcout << "NIl var" << nil_var_name << " is updated to regular expression on ptr_table. Nothing is done for VEC_LIST." << std::endl;);
				rexp_records_for_vars.push_back( sailr_ptr_table_find( &table, nil_var_name ) );
			}else{
				IF_DEBUG( Rcpp::Rcout << "WARNING: NIl var" << nil_var_name << " is updated to unknown type on ptr_table. Nothing is done for VEC_LIST." << std::endl;);
			}
		    // 5. From nil_vars, current element should be removed.
		    iter_nil = nil_vars.erase( iter_nil );
		  }else{
		    ++iter_nil;
		  }
		}
#ifdef SAILR_NULL_UPDATED
		}
		sailr_ptr_table_info_reset_null_updated(&table);
#endif

		// VEC_LIST's STRVECTORs should be updated using result of ptr_table.
		update_sailr_vec_list(vec_list, lhs_vars, table, row_idx ) ; 
		IF_DEBUG( show_sailr_vec_list_nth(vec_list, row_idx ); );
		IF_DEBUG( Rcpp::Rcout << "VEC_LIST is updated" << std::endl; );

		// Free string objects on ptr_table
		IF_DEBUG( Rcpp::Rcout << "Num of str records to be cleaned : " << str_records_for_vars.size() << std::endl; );
		for(auto str_record_iter = str_records_for_vars.begin(); str_record_iter != str_records_for_vars.end(); ++str_record_iter){
			sailr_ptr_record_free_objects( *str_record_iter );
		}

		if( vm_exec_code_result == 2 ){ // if suspended
			if( extfunchash != NULL && ( strcmp( sailr_ext_func_hash_get_last_executed(&extfunchash), "push!") == 0 )){
				// Push new row to DF & set the row number to pushed_row_idx
				pushed_row_idx = vec_list_push_cloned_row(vec_list, row_idx);
				IF_DEBUG( Rcpp::Rcout << "Copying row from " << row_idx << " to " << pushed_row_idx << std::endl; );
				// When push!() is called, sort is required
				sort_required = true;
				goto loop_when_suspended;
			}
		}
			     
		// Reset literal regular expressions
		IF_DEBUG( Rcpp::Rcout << "Num of rexp literal records to be reset : " << rexp_records_for_literals.size() << std::endl; );
		for(auto rexp_iter = rexp_records_for_literals.begin(); rexp_iter != rexp_records_for_literals.end(); ++rexp_iter){
			sailr_ptr_record_reset_rexp( *rexp_iter );
		}

		// Free non-literal regular expressions (assigned to some vars)
		IF_DEBUG( Rcpp::Rcout << "Num of rexp literal records to be cleaned : " << rexp_records_for_vars.size() << std::endl; );
		for(auto rexp_iter = rexp_records_for_vars.begin(); rexp_iter != rexp_records_for_vars.end(); ++rexp_iter){
			sailr_ptr_record_free_objects( *rexp_iter );
		}
		IF_DEBUG( Rcpp::Rcout << "Ptr_record's with type of PTR_REXP are reset." << std::endl; );

	}

	IF_DEBUG( Rcout << "Calculation finished successfully" << std::endl; );
	IF_DEBUG( Rcout << "Show ptr_table" << std::endl; );
	IF_DEBUG( sailr_ptr_table_show_all(&table); std::cout << std::flush; );
	
	/* Show Vec list summary */
	IF_DEBUG( Rcpp::Rcout << "Show the summary of vec list." << std::endl; );
	IF_DEBUG( vec_list_show_summary(vec_list); ); 
	IF_DEBUG( Rcpp::Rcout << "Show contents of vec list." << std::endl; );
	IF_DEBUG( ShowVecList(vec_list, 10); );
	
	finalize: ;

	/* Convert VEC_LIST* to DataFrame. */
	IF_DEBUG( Rcpp::Rcout << "Convert vec_list to Rcpp::DataFrame" << std::endl; );
	IF_DEBUG( Rcpp::Rcout << "vm_exec_success code: " << vm_exec_success << std::endl; );
	DataFrame new_df;
	if(vm_exec_success){
		// For sorting
		std::vector<std::tuple<int, int, int > > tup_vec;
		unsigned int pos;
		unsigned int target_pos;
		unsigned int current_nrow;
		std::vector<int> order_vec;
		std::vector<int> final_order_vec;
		if(sort_required || discard_flag){
			IF_DEBUG( Rcpp::Rcout << "Find _n_ column in VEC_LIST" << std::endl;);
			VEC_ELEM* rownum_vec_elem = vec_elem_find( vec_list, "_n_");
			void* rownum_intvec_void = std::get<1>(*rownum_vec_elem);
			std::vector<int>* rownum_intvec = (std::vector<int>*) rownum_intvec_void;
			current_nrow = rownum_intvec->size();

			IF_DEBUG( Rcpp::Rcout << "Find _discard_ column in VEC_LIST" << std::endl;);
			VEC_ELEM* discard_vec_elem = vec_elem_find( vec_list, "_discard_");
			void* discard_intvec_void = std::get<1>(*discard_vec_elem);
			std::vector<int>* discard_intvec = (std::vector<int>*) discard_intvec_void;

			IF_DEBUG( Rcpp::Rcout << "Prepare tuple vec for sorting" << std::endl; );
			for( pos = 0 ; pos < current_nrow ; ++pos ){
				tup_vec.emplace_back( rownum_intvec->operator[](pos), pos, discard_intvec->operator[](pos) );
			}
			// sort by first element
			IF_DEBUG( Rcpp::Rcout << "Sorting" << std::endl; );
			std::stable_sort(tup_vec.begin(), tup_vec.end());
			IF_DEBUG( Rcpp::Rcout << "Sorting vector is obtained" << std::endl; );

			// extract second elements as vector 
			order_vec = std::vector<int>(current_nrow);

			pos = 0;
			target_pos = 0;
			while( pos < current_nrow){
				if(std::get<2>(tup_vec[pos]) != 1){ // discard is not 1 (= keep the row)
					order_vec[target_pos] = std::get<1>(tup_vec[pos]); // std::get<1>() second element
					++target_pos;
				}
				++pos;
			}
			IF_DEBUG( Rcpp::Rcout << "Size of the final dataframe is to be " << target_pos << std::endl;);
			final_order_vec = std::vector<int>( order_vec.begin(), order_vec.begin() + target_pos);
			IF_DEBUG( Rcpp::Rcout << "Size of the final dataframe is " << final_order_vec.size() << std::endl;);
			
			IF_DEBUG( for(auto iter_final: final_order_vec){ Rcpp::Rcout << iter_final << " ";}; Rcpp::Rcout << std::endl; );
			
			new_df = ConvertVecList(vec_list, lhs_vars, &final_order_vec , service_vars);
		}else{
			new_df = ConvertVecList(vec_list, lhs_vars, NULL , service_vars);
		}
	}
	
	// Free temporarily used cstring variable names
	// In the future, it is better to provide iteration for var name access, then no need to free.
	sailr_varnames_free(var_array, var_num);
	sailr_varnames_free(lhs_var_array, lhs_var_num);
	sailr_varnames_free(rhs_var_array, rhs_var_num);

	// Free vars_to_update which is temporarily used to store variable names to update ptr_table
	free(vars_to_update);

	/* Free memory of instructions */
	IF_DEBUG( Rcpp::Rcout << "Free VM instruction list" << std::endl; );
	sailr_vm_inst_list_free( inst_list );
	IF_DEBUG( Rcpp::Rcout << "Free VM instruction code" << std::endl; );
	sailr_vm_inst_code_free( vmcode );

	/* Free memory */
	IF_DEBUG( Rcpp::Rcout << "Free parser tree" << std::endl; );
	sailr_tree_free(ps);
	
	IF_DEBUG( Rcpp::Rcout << "Free pointer table" << std::endl; );
	sailr_ptr_table_del_all(&table);

	IF_DEBUG( Rcpp::Rcout << "Free parser state object" << std::endl; );
	sailr_parser_state_free(ps);

	/* Free external function hash */
	sailr_ext_func_hash_free(&extfunchash);

	/* Free VEC_LIST */
	vec_list_free(vec_list);

	if(vm_exec_success){
		return new_df;  // should return dataframe.
	}else{
	    Rcpp::stop( "sailr virtual machine runtime error. stopped. \n" );
	}
}



