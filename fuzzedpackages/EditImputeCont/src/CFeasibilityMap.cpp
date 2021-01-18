#include "CHeader.h" 
#include "lp_lib.h"

#include "CData.h"
#include "CFeasibilityMap.h"
#include "R.h"

MapManager::MapManager(CData &Data):data(Data) {  
		
	iter = 0;	
	last_trim.resize(data.n_faulty + 1);
	for (int i = 0; i < data.n_faulty + 1; i++) {
	  last_trim[i] = 1;
	}
	
	map_feasible_tau_t.resize(data.n_faulty + 1);
	
	last_trim_fn2.resize(data.n_faulty + 1);
	for (int i = 0; i < data.n_faulty + 1; i++) {
	  last_trim_fn2[i] = 1;
	}
		
	map_count_x_out.resize(data.n_faulty + 1);

}

//Constructor
CFeasibilityMap::CFeasibilityMap(){
  Debug = false;
 	useMap = true;
  // useMap = false; // memory error if this is false when using large dataset
	pmm = NULL;
}

//Destructor
CFeasibilityMap::~CFeasibilityMap(){}

void CFeasibilityMap::Build(CData &Data) {
	
  if (useMap) {
		
    if (pmm!=NULL) {
      delete pmm; pmm=NULL;
    }
    pmm = new MapManager(Data);
		
  } else {
    
		feasibleMap = Matrix(Data.n_tau,Data.n_faulty); 
    for (int i_faulty=1; i_faulty<=Data.n_faulty; i_faulty++){
						
      int i_original = Data.Faulty2Original[i_faulty-1];
			
      ColumnVector x_tilde_i = (Data.D_Observed.row(i_original)).t(); 
			
      for (int i_tau=1; i_tau<=Data.n_tau; i_tau++){			
       ColumnVector s_i = tau_to_s_fn(i_tau,Data.n_var);
			 feasibleMap(i_tau,i_faulty) = feasible_test_fn(Data, x_tilde_i, s_i, i_original, Data.epsilon) ;
	    } 
			 
      if ( ((1.0*i_faulty/100)==(floor(1.0*i_faulty/100))) ){ 
		Rprintf( "Mat_feasible_tau_t for i_sample= %d\n",i_faulty);
  	  }	
    
		} // for (int i_faulty=1; i_faulty<=Data.n_faulty; i_faulty++)  
  
	} // if (useMap) {...} else {
	
} // void CFeasibilityMap::Build

int CFeasibilityMap::feasible_test_fn(CData &Data, ColumnVector &x_tilde_i, ColumnVector &s_i, int i_original, float epsilon) {
    ColumnVector Dummy;
    return feasible_test_fn(Data, x_tilde_i, s_i, i_original, false, epsilon, Dummy);
}

// Commented on 05/21/2015
// void CFeasibilityMap::test(CData &Data) {
//   int i_original = 519;
//   ColumnVector s_i(Data.n_var); s_i = 0;
//   s_i(4) = 1;s_i(5) = 1; s_i(6) = 1;
//   s_i(7) = 1;s_i(8) = 1; s_i(9) = 1;
//   s_i(10) = 1;s_i(11) = 1;
//   s_i(15) = 1;
//   s_i(16) = 1;s_i(18) = 1;
//   s_i(19) = 1;
//   s_i(22) = 1;s_i(23) = 1; s_i(24) = 1;
//   s_i(26) = 1;s_i(27) = 1;
//   cout << " start testing..." << endl;
//   try {
//     ColumnVector x_tilde_i = (pmm->data.D_Observed.row(i_original)).t() ;    // defined in mm
//     cout << x_tilde_i;
//     cout << s_i;
//     test_feasible_test_fn(Data, x_tilde_i, s_i);
//   } catch (...) {
//     cout << "oops" << endl;
//   }
//
//   cout << "Done test" << endl;
// }

int CFeasibilityMap::isCandidateFeasible(ColumnVector &cand_s, int i_faulty) {
  // This function gives 0 for infeasible solution s_i
  //    								 1 for feasible s_i
  if (!useMap) {
    return (int)feasibleMap(s_to_tau_fn(cand_s),i_faulty);
  } else {
  	static int MAX_MAP_SIZE = 30000; 		// WAS 30000;		//the maximum map size allowed
  	static int TRIM_INCREAMENT = 2;   	// 100
  																			// calculate appropriately to remove 30%
                                        
    int i_tau = s_to_tau_fn(cand_s);
  	map<int,int>& current = pmm->map_feasible_tau_t[i_faulty];	// 'current' is the vector for i_faulty
  	// if ( i_tau is not in the vector ) then store it to 'current'
  	// if 'current' is full, eliminate one of elements
  	if(current.find(i_tau) == current.end()) { 
  		// ask result of i_tau is NOT stored in the map
  		// .find returns an iterator if found, o.w. map::end

  		//trim the map
  		if (current.size() > MAX_MAP_SIZE) {
  			int last_trim = pmm->last_trim[i_faulty];
  			while (current.size() > ((int)(double)MAX_MAP_SIZE * 0.7 )) { //reduce map to 70% full
  				map<int,int>::iterator it = current.begin();
  				while (it != current.end()) {
  					int temp = it->second ;
  					if (abs(temp) <= last_trim) {
  						current.erase(it++); //note the post increment, do not use pre-increment
  					} else {
  						++it;	//can use pre-increment in this case
  					}
  				}
  				last_trim += TRIM_INCREAMENT;
  			}
  			if (last_trim > pmm->iter) { last_trim = pmm->iter;}
  			pmm->last_trim[i_faulty] = last_trim;
  		} // if (current.size() > MAX_MAP_SIZE)
      int i_original = pmm->data.Faulty2Original[i_faulty-1];
      
  		ColumnVector x_tilde_i = (pmm->data.D_Observed.row(i_original)).t() ;  	// defined in mm
       
      int map_v = pmm->iter * feasible_test_fn(pmm->data, x_tilde_i, cand_s, i_original,pmm->data.epsilon);
  		if (map_v==0) map_v = -pmm->iter ;
  		current[i_tau] = map_v; 

  	} else {
  		// update iteration number
  		if ( current[i_tau]>0 ){
  			current[i_tau] = pmm->iter;
  		} else {
  			current[i_tau] = -pmm->iter ;
  		} 
  	} // if(current.find(i_tau) == current.end())
  	return current[i_tau]>0?1:0;
  }
}


int CFeasibilityMap::test_feasible_test_fn(CData &Data, ColumnVector &x_tilde_i, ColumnVector &s_i){
  int is_feasible = 0 ;
  
  int n_var = x_tilde_i.nrows() ; 
  int n_EditVec = Data.EditVec.nrows() ;
  int sum_s_1 = s_i.sum() ; int sum_s_0 = n_var - sum_s_1 ;   
  ColumnVector x_0(sum_s_0); x_0 = 0;
  Matrix A_0(n_EditVec,sum_s_0) ; A_0 = 0;
	Matrix A_1(n_EditVec,sum_s_1) ; A_1 = 0;
	
  int count0 = 0;
  int count1 = 0;
	for (int j_var=1; j_var<=n_var; j_var++){
		if ((int)s_i(j_var)==0){
			count0++; 
			x_0(count0) = x_tilde_i(j_var) ; 
			A_0.column(count0) = Data.EditMat.column(j_var) ; 
		} else {
      count1++;
			A_1.column(count1) = Data.EditMat.column(j_var) ; 
		}
	} 
  
	ColumnVector which_rows(n_EditVec); which_rows = 0 ; 
	for (int i_row=1; i_row<=n_EditVec; i_row++){
		int n_zero = 0 ; 
		for (int j=1; j<=sum_s_1; j++){
			if (A_1(i_row,j)==0) n_zero++; 
		}
		if (n_zero<sum_s_1) which_rows(i_row) = 1 ; 
	}

	Matrix A_1_red(which_rows.sum(),sum_s_1); 
	ColumnVector b_1_red(which_rows.sum()); 
	Matrix A_0_red(which_rows.sum(),sum_s_0); 
	for (int i_row=1, count_row = 0; i_row<=n_EditVec; i_row++){
		if (which_rows(i_row)==1){
			count_row++; 
			A_1_red.row(count_row) = A_1.row(i_row) ;  	// Constraints
			b_1_red.row(count_row) = Data.EditVec.row(i_row) ; 
			A_0_red.row(count_row) = A_0.row(i_row) ; 
		} 
	}
  ColumnVector EditVec_shr = b_1_red - A_0_red * x_0 ;
  is_feasible = SolveLP(A_1_red, EditVec_shr);
	return is_feasible?1:0; 
}


int CFeasibilityMap::feasible_test_fn(CData &Data, ColumnVector &x_tilde_i, 
    ColumnVector &s_i, int i_original, bool initD_S, float epsilon, ColumnVector &x){
  int is_feasible = 0 ;
  if (initD_S) {
    if (!Data.PassStep0_for_init(x_tilde_i,s_i,epsilon)) { return is_feasible;}
  } else {
    // Step 0 - check balance edits hold, but only one s_ij has the value 1. 
    if (!Data.PassStep0(s_i,i_original)) { return is_feasible;}
  }
  
  // Step 1
  int n_var = x_tilde_i.nrows() ; 
  int n_EditVec = Data.EditVec.nrows() ;
  int sum_s_1 = s_i.sum() ; int sum_s_0 = n_var - sum_s_1 ;   
	ColumnVector x_0(sum_s_0); x_0 = 0;
  Matrix A_0(n_EditVec,sum_s_0) ; A_0 = 0;
	Matrix A_1(n_EditVec,sum_s_1) ; A_1 = 0;
	
  int count0 = 0;
  int count1 = 0;
	for (int j_var=1; j_var<=n_var; j_var++){
		if ((int)s_i(j_var)==0){
			count0++; 
			x_0(count0) = x_tilde_i(j_var) ; 
			A_0.column(count0) = Data.EditMat.column(j_var) ; 
		} else {
      count1++;
			A_1.column(count1) = Data.EditMat.column(j_var) ; 
		}
	} 
 
  Data.Debug = Debug;
  if (!Data.PassStep1(s_i, A_0, x_0)) { return is_feasible;}
	
  // Step 2: lpSolve to check whether s_i has a feasible solution of x_i
	ColumnVector which_rows(n_EditVec); which_rows = 0 ; 
	for (int i_row=1; i_row<=n_EditVec; i_row++){
		int n_zero = 0 ; 
		for (int j=1; j<=sum_s_1; j++){
			if (A_1(i_row,j)==0) n_zero++; 
		}
		if (n_zero<sum_s_1) which_rows(i_row) = 1 ; 
	}

	Matrix A_1_red(which_rows.sum(),sum_s_1); 
	ColumnVector b_1_red(which_rows.sum()); 
	Matrix A_0_red(which_rows.sum(),sum_s_0); 
	for (int i_row=1, count_row = 0; i_row<=n_EditVec; i_row++){
		if (which_rows(i_row)==1){
			count_row++; 
			A_1_red.row(count_row) = A_1.row(i_row) ;  	// Constraints
			b_1_red.row(count_row) = Data.EditVec.row(i_row) ; 
			A_0_red.row(count_row) = A_0.row(i_row) ; 
		} 
	}
  ColumnVector EditVec_shr = b_1_red - A_0_red * x_0 ;
  
  if (initD_S) {
    is_feasible = SolveLP(A_1_red, EditVec_shr,x);    
  } else {
    is_feasible = SolveLP(A_1_red, EditVec_shr);
  }
	return is_feasible?1:0; 
}


void CFeasibilityMap::initilize_D_and_S(CData &Data) {
	
  ColumnVector order_to_test = get_order_to_test(Data.n_tau, Data.n_var);
  ColumnVector list_feasible_type2 = get_feasible_tau(Data);
	
	for (int i_faulty = 1; i_faulty <= Data.n_faulty; i_faulty++) {  
		
    bool is_pass = false;
    int i_original = Data.Faulty2Original[i_faulty-1];
    ColumnVector x_tilde_i = (Data.D_Observed.row(i_original)).t(); 
    for (int i_order = 1; i_order <= order_to_test.nrows() && !is_pass; i_order++) {
      int i_tau = order_to_test(i_order);
      ColumnVector s_i = tau_to_s_fn(i_tau,Data.n_var);
    
      bool skip_for_type2 = false; 
      if (Data.is_case(i_original,2) && list_feasible_type2(i_tau) == 0) { skip_for_type2 = true;}
      if (!skip_for_type2){
        ColumnVector x_mean;
        int is_feasible = feasible_test_fn(Data, x_tilde_i, s_i, i_original, true, Data.epsilon, x_mean);
        
        if (is_feasible > 0) {
          //copy solution to a temp vector
          ColumnVector temp = x_tilde_i;
          for (int index = 1, count =0; index <= Data.n_var; index++){
            if (s_i(index) == 1) {
              count++;
              temp(index) = x_mean(count);
            }
          }
          Data.Debug = Debug;
          if (Data.PassEdits(temp)) { //then check if it satifies edits          
            is_pass = true;
            Data.initial_S_Mat.row(i_faulty) = s_i.t();   
            Data.D_initial.row(i_original) = temp.t();
          }
        }//is_feasible
        Debug = false;
      }
    } 
  
	}  
	
}

bool CFeasibilityMap::SolveLP(Matrix &A, ColumnVector &b, ColumnVector &x) {
  lprec *lp ;
  int n_row = A.nrows(); int n_col = A.ncols();
  x = ColumnVector(n_col); x = 0;
  lp = make_lp(0,n_col) ; 
  
  double *input_row = new double[1+n_col];
  for (int i_row=1; i_row<=n_row; i_row++){
      input_row[0] = 0 ; // The first zero is for matrix form
      for (int j=1; j<=n_col; j++){
          input_row[j] = A(i_row,j) ;
      }
      add_constraint(lp, input_row, LE, b(i_row)) ;
  }
  delete [] input_row;
  
  double *input_obj = new double[1+n_col];    // The first zero is for matrix form
  input_obj[0] = 0 ;
  for (int j=1; j<=n_col; j++){
      input_obj[j] = 1 ;
  }
  set_obj_fn(lp, input_obj) ;
  delete [] input_obj;
  set_verbose(lp, IMPORTANT); // NEUTRAL (0), IMPORTANT (3), NORMAL (4), FULL (6)
  bool is_feasible = (solve(lp)==0); // 0: feasible solution found,  2: not found
                                     // solution for minimizing objective function                               
  double* x_min = new double[n_col];
  double* x_max = new double[n_col];                      
  if (is_feasible) {
    get_variables(lp, x_min);
    set_maxim(lp);
    is_feasible = (solve(lp)==0); // 0: feasible solution found,  2: not found
    if (is_feasible) {
      get_variables(lp, x_max);
      for (int i = 0; i < n_col; i++) {
        x(i+1) = (x_min[i] + x_max[i]) / 2.0;
      }
    }
  }
  
  delete [] x_min;
  delete [] x_max;
                                     
  delete_lp(lp);
  return is_feasible;
}

//generate s_n in reverse order
ReturnMatrix CFeasibilityMap::tau_to_s_fn2( double tau_i, int n_var ){
    ColumnVector s_i(n_var) ; 
    unsigned int tau = (unsigned int)tau_i;
    for (int i_var = 1; i_var <=n_var; i_var++) {
      s_i(n_var - i_var + 1) = tau & 1; tau >>= 1; 
    }
    s_i.release(); return s_i;
}

ReturnMatrix CFeasibilityMap::tau_to_s_fn( double tau_i, int n_var ){
    ColumnVector s_i(n_var) ; 
    unsigned int tau = (unsigned int)tau_i;
    for (int i_var = 1; i_var <=n_var; i_var++) {
      s_i(i_var) = tau & 1; tau >>= 1; 
    }
    s_i.release(); return s_i;
}

int CFeasibilityMap::s_to_tau_fn( ColumnVector &s_i ){
  int tau_i=0; 
  for (int i_var = 0; i_var<s_i.nrows(); i_var++){
    if (s_i(i_var+1)==1) tau_i += 1 << i_var;
	}
  return tau_i;
}


ColumnVector CFeasibilityMap::get_feasible_tau(CData &Data) {
  ColumnVector list_feasible_type2(Data.n_tau); 
  for (int i_tau = 1; i_tau <=Data.n_tau; i_tau++){
		if (i_tau%10000000==0) {
			// double Prog = 100.00 * i_tau / Data.n_tau ;
			// int Prog_int = Prog ; // Changed 
			Rprintf( "progress = %d percent\n");
			// Changed
		}
    ColumnVector s_i = tau_to_s_fn(i_tau, Data.n_var);
    list_feasible_type2(i_tau) = Data.get_feasible_tau(s_i); 
  } 
  return(list_feasible_type2);
}

ColumnVector CFeasibilityMap::get_order_to_test(int n_tau, int n_var){
  ColumnVector order_to_test(n_tau);
  int count = 0;  
  for (int i = 1; i < n_var; i++) { //1 to n_var-1
      if ( (n_var<=13)||(i>=(n_var-5)) ){ // Added 2015/3/10
				for (int i_tau = n_tau; i_tau >0; i_tau--)  {
	  			ColumnVector s_i = tau_to_s_fn2(i_tau,n_var);
	  			if (s_i.sum() == i) {order_to_test(++count) = s_to_tau_fn(s_i);}
				}
      } // Added 2015/3/10
  }
  return order_to_test; 
}

//return tau_q
int CFeasibilityMap::EvaluateMove(int i_original, CData &Data, ColumnVector &s_i, int iter,
  int &what_type_move, double &g_option_q, double &g_mode_q, bool SampleMove) {
  int i_faulty = Data.Original2Faulty[i_original-1];
  ColumnVector free_s_i = Data.get_free_s_i(i_original, s_i);
  int sum_free_s_i_ones = 0; 
  int sum_free_s_i_zeros = 0;
  for (int j_var=1; j_var<=Data.n_var; j_var++){ // free_s_i can have values 0, 1 and 9.
    if (free_s_i(j_var)==1) sum_free_s_i_ones++; 
		if (free_s_i(j_var)==0) sum_free_s_i_zeros++; 
	}
  
	// propose legitimate s_q
	ColumnVector isMoveOption_q(3) ; isMoveOption_q = 1 ; // for 1:birth, 2:move, 3:death
  if ( sum_free_s_i_ones<=1 ) isMoveOption_q(3) = 0 ; 		// no death option
  if ( sum_free_s_i_zeros<=1 ) isMoveOption_q(1) = 0 ; 		// no birth option
  
  double g_option_birth_q=0, g_option_move_q=0, g_option_death_q=0; 
	int tau_birth_q=0, tau_move_q=0, tau_death_q=0 ;
  
  
	if ( isMoveOption_q(1)==1 ){
		// birth	ex. s_i (10010) -> s_q (11010) 
		ColumnVector var_zeros = whichZeros(free_s_i) ;	// ex. (2 3 5)
		Matrix cand_S_Mat(sum_free_s_i_zeros,Data.n_var); 
		ColumnVector is_cand_tau_feasible(sum_free_s_i_zeros) ; is_cand_tau_feasible = 0 ; 
		for (int j_temp=1; j_temp<=sum_free_s_i_zeros; j_temp++){
			ColumnVector cand_s = s_i ;
			int loc_zero = var_zeros(j_temp) ; // ex. 2, 3, or 5
			cand_s(loc_zero) = 1; 
			
      Data.AddjustCandidate_s_i_ForSatisfiedBalanceEdit(cand_s, i_original);
      if (SampleMove) {
			  cand_S_Mat.row(j_temp) = cand_s.t();
      }
      is_cand_tau_feasible(j_temp) = isCandidateFeasible(cand_s,i_faulty);
		} 

		if ( is_cand_tau_feasible.sum()>0 ){
      if (SampleMove) {
  			// Changed
  			// ex. loc_zero = 2, 3, or 5, if putting 1 for s_i2 and s_i5 is feasible 
  			// 		 sum_s_i_zeros = 3, is_cand_tau_feasible = (1 0 1)
  			ColumnVector cand_tau_index = whichOnes(is_cand_tau_feasible) ; // ex. (1 3)
  			int random_int = runifdiscrete_fn(is_cand_tau_feasible.sum()) ; // ex. 1 or 2 since is_is_cand_tau_feasible_ble.sum()=2
  			int selected_order = cand_tau_index(random_int) ; // ex. 1 or 3
  			ColumnVector selected_s_i = (cand_S_Mat.row(selected_order)).t() ; // ex. s_i corresponding to putting 1 for s_i2 or s_i5
  			tau_birth_q = CFeasibilityMap::s_to_tau_fn(selected_s_i) ;
      }
			g_option_birth_q = 1.0/is_cand_tau_feasible.sum(); 
		} else {
			isMoveOption_q(1) = 0 ; 
		} 
	} 
	if ( isMoveOption_q(2)==1 ){
		// move		(10010) -> (11000) 
		ColumnVector var_ones = whichOnes(free_s_i) ; // ex. (1 4)
		ColumnVector var_zeros = whichZeros(free_s_i) ; // ex. (2 3 5)
		Matrix cand_S_Mat(sum_free_s_i_ones*sum_free_s_i_zeros,Data.n_var) ; 
		ColumnVector is_cand_tau_feasible(sum_free_s_i_ones*sum_free_s_i_zeros) ; is_cand_tau_feasible = 0 ; 
    
		int count_row = 0 ; 
		for (int j_temp=1; j_temp<=sum_free_s_i_ones; j_temp++){
			for (int k_temp=1; k_temp<=sum_free_s_i_zeros; k_temp++){
				count_row = count_row + 1 ; 
				ColumnVector cand_s = s_i ;
				int loc_one = var_ones(j_temp) ; // ex. 1 or 4
				cand_s(loc_one) = 0 ; 
				int loc_zero = var_zeros(k_temp) ; // ex. 2, 3, or 5
				cand_s(loc_zero) = 1 ; 
        
        Data.AddjustCandidate_s_i_ForSatisfiedBalanceEdit(cand_s, i_original);
        if (SampleMove) {
				  cand_S_Mat.row(count_row) = cand_s.t() ;
        }
        is_cand_tau_feasible(count_row) = isCandidateFeasible(cand_s,i_faulty);
			} 
		} 
    
		if ( is_cand_tau_feasible.sum()>0 ){
      if (SampleMove) {
  			// ex. var_ones = 2, 3, and var_ones = 1, 4 if putting changing (1,2) or (3,4) is feasible 
  			// 		 sum_free_s_i_zeros = sum_free_s_i_zeros = 2, count_row:1~4, is_cand_tau_feasible = (1 0 0 1)
  			ColumnVector cand_tau_index = whichOnes(is_cand_tau_feasible) ; // ex. (1, 4)
  			int random_int = runifdiscrete_fn(is_cand_tau_feasible.sum()) ; // ex. 1~2 since is_cand_tau_feasible.sum=2
  			int selected_order = cand_tau_index(random_int) ; // ex. 1 or 4
  			ColumnVector selected_s_i = (cand_S_Mat.row(selected_order)).t() ; // ex. s_i for changing (1,2) or (3,4)
  			tau_move_q = CFeasibilityMap::s_to_tau_fn(selected_s_i) ; 
      }
			g_option_move_q = 1.0/is_cand_tau_feasible.sum() ; 
		} else {
			isMoveOption_q(2) = 0 ;  
		} 
	} 

	if ( isMoveOption_q(3)==1 ) {
		// death ex. s_i (11110) -> s_q (11010) 
		ColumnVector var_ones = whichOnes(free_s_i) ; // ex. (1 2 3 4)
		Matrix cand_S_Mat(sum_free_s_i_ones,Data.n_var) ; 
		// ColumnVector cand_tau(sum_free_s_i_ones) ; 
		ColumnVector is_cand_tau_feasible(sum_free_s_i_ones) ; is_cand_tau_feasible = 0 ; 
		for (int j_temp=1; j_temp<=sum_free_s_i_ones; j_temp++){
			ColumnVector cand_s = s_i ;
			int loc_one = var_ones(j_temp) ; // ex. 1, 2, 3, or 4
			cand_s(loc_one) = 0 ; 

			Data.AddjustCandidate_s_i_ForSatisfiedBalanceEdit(cand_s, i_original);
      if (SampleMove) {
			    cand_S_Mat.row(j_temp) = cand_s.t() ; 
      }
      is_cand_tau_feasible(j_temp) = isCandidateFeasible(cand_s,i_faulty);
		} 

		if ( is_cand_tau_feasible.sum()>0 ){
      if (SampleMove) {
  			ColumnVector cand_tau_index = whichOnes(is_cand_tau_feasible) ; // ex. (2, 3, 5)
  			int random_int = runifdiscrete_fn(is_cand_tau_feasible.sum()) ; // ex. 1~3
  			int selected_order = cand_tau_index(random_int) ; // ex. 2, 3, or 5
  			ColumnVector selected_s_i = (cand_S_Mat.row(selected_order)).t() ; 
  			tau_death_q = CFeasibilityMap::s_to_tau_fn(selected_s_i) ;
      }
			g_option_death_q = 1.0/is_cand_tau_feasible.sum() ; 	
		} else {
			isMoveOption_q(3) = 0 ; 
		} 
	} 
	
  
	if ( isMoveOption_q.sum() == 0 ){
		Rprintf( "Check: case_i(i) > 0 but isMoveOption_q.sum() == 0 \n");
		Rprintf( "iter= %d, i= %d, s_j= \n",iter,i_original);
		for (int j_temp=1; j_temp<=Data.n_var; j_temp++){
			Rprintf( " %d",(int)s_i(j_temp));
		}
		Rprintf( "\n");
		Rprintf( "tau_i =%d\n",s_to_tau_fn(s_i));
	}  

	g_mode_q = 1.0/(isMoveOption_q.sum()) ;
  if (SampleMove) {
  	ColumnVector temp_index = whichOnes(isMoveOption_q) ;	// ex. (1 3)
  	int temp_int = runifdiscrete_fn(isMoveOption_q.sum()) ; 
  	what_type_move = temp_index(temp_int) ;	// ex. 1 or 3
  }
	int tau_q = 0; 
  if (SampleMove) {
  	if (what_type_move==1){ tau_q = tau_birth_q ; g_option_q = g_option_birth_q ; }
  	if (what_type_move==2){ tau_q = tau_move_q ; g_option_q = g_option_move_q ; }
  	if (what_type_move==3){ tau_q = tau_death_q ; g_option_q = g_option_death_q ; }
  } else {
    if (what_type_move==1){ g_option_q = g_option_death_q ; }
  	if (what_type_move==2){ g_option_q = g_option_move_q ; }
  	if (what_type_move==3){ g_option_q = g_option_birth_q ; }
  }
  
  return tau_q;  
}

bool CFeasibilityMap::SolveLP(Matrix & A, ColumnVector &b, bool domax) {
  lprec *lp ;
  int n_row = A.nrows(); int n_col = A.ncols();
  lp = make_lp(0,n_col) ; 

  double *input_row = new double[1+n_col];
  for (int i_row=1; i_row<=n_row; i_row++){
      input_row[0] = 0 ;
      for (int j=1; j<=n_col; j++){
          input_row[j] = A(i_row,j) ;
      }
      add_constraint(lp, input_row, LE, b(i_row)) ;
  }
  delete [] input_row;
  // we don't even need an objective function !!
  double *input_obj = new double[1+n_col];    // The first zero is for matrix form
  input_obj[0] = 0 ;
  for (int j=1; j<=n_col; j++){
      input_obj[j] = 1 ;
  }
  set_obj_fn(lp, input_obj) ;
  if (domax) {
    set_maxim(lp);
  }
  delete [] input_obj;
  set_verbose(lp, IMPORTANT); // NEUTRAL (0), IMPORTANT (3), NORMAL (4), FULL (6)
  bool is_feasible = (solve(lp)==0); // 0: feasible solution found,  2: not found
  	                                 // solution for minimizing objective function
	delete_lp(lp);
  return is_feasible;
}

bool CFeasibilityMap::SolveLP(Matrix &A, ColumnVector &b) {
  lprec *lp ;
  int n_row = A.nrows(); int n_col = A.ncols();
  lp = make_lp(0,n_col) ; 

  double *input_row = new double[1+n_col];
  for (int i_row=1; i_row<=n_row; i_row++){
      input_row[0] = 0 ;
      for (int j=1; j<=n_col; j++){
          input_row[j] = A(i_row,j) ;
      }
      add_constraint(lp, input_row, LE, b(i_row)) ;
  }
  delete [] input_row;
  
  double *input_obj = new double[1+n_col];    // The first zero is for matrix form
  input_obj[0] = 0 ;
  for (int j=1; j<=n_col; j++){
      input_obj[j] = 1 ;
  }
  set_obj_fn(lp, input_obj) ;
  delete [] input_obj;
  set_verbose(lp, IMPORTANT); // NEUTRAL (0), IMPORTANT (3), NORMAL (4), FULL (6)
  bool is_feasible = (solve(lp)==0); // 0: feasible solution found,  2: not found
		                                 // solution for minimizing objective function
	delete_lp(lp);

  return is_feasible;
}

int CFeasibilityMap::count_x_out_fn(CData &Data,int i_tau, int i_original,int n_simul, Uniform &randUnif) {
	
  // double case2_count_out = 0;
	int case2_count_out = 0; // Changed by Hang on 5/16/2015
	
  ColumnVector s_i = tau_to_s_fn( i_tau, Data.n_var );   
  ColumnVector item_by_joint = Data.copy_non_balance_edit(s_i);
  ColumnVector tilde_y_i = Data.log_D_Observed.row(i_original).t();
  
	for (int i_simul=1; i_simul<=n_simul; i_simul++){
			//Generate from uniform distribution
			ColumnVector y_q = tilde_y_i;
			for ( int temp_j=1; temp_j<=Data.n_var; temp_j++ ){
				if ( item_by_joint(temp_j)==1 ){
					y_q(temp_j) = Data.logB_L(temp_j)+Data.logB_U_L(temp_j)*randUnif.Next(); 
				} 
			} 
	
			ColumnVector x_q = exp_ColumnVector(y_q) ;
      Data.update_full_x_for_balance_edit(x_q);
			// if (!Data.PassEdits(x_q)) { case2_count_out += 1.0;}
      if (!Data.PassEdits(x_q)) { case2_count_out += 1;}  // Changed by Hang on 5/16/2015
	} 
  if (case2_count_out ==0) {
    case2_count_out = 1;
  }
	return case2_count_out; // ADDED by Hang on 5/16/2015
}

int CFeasibilityMap::Map_count_x_out_fn(int i_tau, int i_original,int n_simul, Uniform &randUnif) {
  if (!useMap) {
    //should not be called
	Rprintf( "Map_count_x_out_fn should not be called in non-map mode\n");
    return 0;//which will cause an error (segment fault)
  } else {
  	static int MAX_MAP_SIZE_FN2 = 20000;   	// WAS 20000;		//the maximum map size allowed
	  static int TRIM_INCREAMENT_FN2 = 2;   // 100  
	  static int MODE = 100000;
    int i_faulty = pmm->data.Original2Faulty[i_original-1];
    map<int,int>& current_fn2 = pmm->map_count_x_out[i_faulty];  // 'current_fn2' is the vector for i_faulty
  	
  	// if ( i_tau is not in the vector ) then store it to 'current'
  	// if 'current' is full, eliminate one of elements
  	if(current_fn2.find(i_tau) == current_fn2.end()) { 	
  		//trim the map
  		if (current_fn2.size() > MAX_MAP_SIZE_FN2) {
        int last_trim_fn2 = pmm->last_trim_fn2[i_faulty];
  			while (current_fn2.size() > ((int)(double)MAX_MAP_SIZE_FN2 * 0.7 )) { //reduce map to 70% full
  				map<int,int>::iterator it = current_fn2.begin();					
  				while (it != current_fn2.end()) {
  					int temp = it->second;		
            temp = temp % MODE; // get the iter 
  					if (temp <= last_trim_fn2) {
  						current_fn2.erase(it++); //note the post increment, do not use pre-increment
  					} else {
  						++it;	//can use pre-increment in this case
  					}
  				}
  				last_trim_fn2 += TRIM_INCREAMENT_FN2;
  			}
  			if (last_trim_fn2 > pmm->iter) { last_trim_fn2 = pmm->iter;}
  			pmm->last_trim_fn2[i_faulty] = last_trim_fn2;
  		} // if (current.size() > MAX_MAP_SIZE)
  		
      int i_original = pmm->data.Faulty2Original[i_faulty-1];
      int count_out = count_x_out_fn(pmm->data, i_tau,i_original, n_simul, randUnif);
      current_fn2[i_tau] = (count_out-1)*MODE + pmm->iter;
  	} else {
      int temp = current_fn2[i_tau];
      int count_out = (temp-temp%MODE)/MODE+1 ; 
  	  current_fn2[i_tau] = (count_out-1)*MODE + pmm->iter; 
  	}
    int temp = current_fn2[i_tau];
    int print_count_out = (temp-temp%MODE)/MODE+1 ; 
	  return print_count_out;
  }
}

double CFeasibilityMap::Simulate_logUnif_case2(int i_tau, int i_original,int n_simul,Uniform &randUnif) {
  if (!useMap) {return 0;} //should not happen!
  double case2_count_out = Map_count_x_out_fn(i_tau, i_original,n_simul,randUnif);
  ColumnVector s_i = tau_to_s_fn( i_tau, pmm->data.n_var );
  ColumnVector item_by_joint = pmm->data.copy_non_balance_edit(s_i);
	double Area = 1.0; 
	for ( int temp_j=1; temp_j<=pmm->data.n_var; temp_j++ ){
		if ( item_by_joint(temp_j)==1 ){
			Area = Area * pmm->data.logB_U_L(temp_j) ; 
		} 
	} 
  Area = Area * case2_count_out / n_simul ;
  return -log(Area);        
}

void CFeasibilityMap::Simulate_logUnif_case2(int n_simul, Uniform &randUnif, CData &Data) {
  if (useMap) {return;}
  if (feasibleMap.nrows() == 0 || feasibleMap.maximum()==0) {
	  Rprintf( "Feasibility Map need to be set or computed first\n");

    return;
  }
  Data.logUnif_case2 =Matrix(Data.n_faulty,Data.n_tau); Data.logUnif_case2 = 0;
  for (int i_original=1; i_original<=Data.n_sample; i_original++){
    if (Data.is_case(i_original,2)) {
      int i_faulty = Data.Original2Faulty[i_original-1];
      ColumnVector tilde_y_i = Data.log_D_Observed.row(i_original).t();
      for (int i_tau=1; i_tau<=Data.n_tau; i_tau++){
        if ( feasibleMap(i_tau,i_faulty)==1 ){
          double case2_count_out = 0;
          ColumnVector s_i = tau_to_s_fn( i_tau, Data.n_var );
          ColumnVector item_by_joint = Data.copy_non_balance_edit(s_i);
          
    			for (int i_simul=1; i_simul<=n_simul; i_simul++){
    					//Generate from uniform distribution
    					ColumnVector y_q = tilde_y_i;
    					for ( int temp_j=1; temp_j<=Data.n_var; temp_j++ ){
    						if ( item_by_joint(temp_j)==1 ){
    							y_q(temp_j) = Data.logB_L(temp_j)+Data.logB_U_L(temp_j)*randUnif.Next(); 
    						} 
    					} 
    			
    					ColumnVector x_q = exp_ColumnVector(y_q) ;
              Data.update_full_x_for_balance_edit(x_q);
              if (!Data.PassEdits(x_q)) { case2_count_out += 1.0;} 
    			} 
    		
    			double Area = 1.0; 
    			for ( int temp_j=1; temp_j<=Data.n_var; temp_j++ ){
    				if ( item_by_joint(temp_j)==1 ){
    					Area = Area * Data.logB_U_L(temp_j) ; 
    				} 
    			} 
          Area = Area * case2_count_out / n_simul ;
          Data.set_logUnif_case2(i_original, i_tau, -log(Area));
     		} 
    	} 
    	
      if ( ((1.0*i_original/100)==(floor(1.0*i_original/100))) ){ 
		   Rprintf( "logUnif_y_tilde for i_sample= %d\n",i_original);
    	}
    }
  }
}
