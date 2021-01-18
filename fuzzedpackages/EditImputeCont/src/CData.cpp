#include "CHeader.h"  
#include "R.h"
#include "CData.h"

CData::CData() {
  has_initialValue = false;
  has_trueS = false;
  n_sample = 0; n_var = 0; n_EditVec = 0;n_balance_edit = 0;
  seed_no = 0.101;
  epsilon = 0.1;
  data_folder = "";
  Debug = false;
  
}

// Commented on 05/21/2015
// CData::CData(string folder){		
//   has_initialValue = true;
//   //read in data from files
//   data_folder = folder;
//   ReadData();
//   init();
//   //read in optional data if any
//   ReadOptionalData();
// }

// Commented on 05/21/2015
// void CData::init(string folder){
//   has_initialValue = true;
//   //read in data from files
//   data_folder = folder;
//   ReadData();
//   init();
//   //read in optional data if any
//   ReadOptionalData();
// }


//Destructor
CData::~CData(){
}

void CData::init() {
  //pre-compute log_D_Observed
  log_D_Observed = Matrix(n_sample,n_var);
  double *d = D_Observed.data();
  double *log_d = log_D_Observed.data(); 
  for (int i = 0; i < n_sample * n_var; i++) {
    *log_d++ = log(*d++);
  }
  //figure out banalce edit and record violation types
  initilize_obs_fail_matrix();
  initilize_balance_edits();
  //logUnif_case2 = Matrix(n_faulty,n_tau); logUnif_case2 = 0;
  D_initial = D_Observed;
  initial_S_Mat = Matrix(n_faulty,n_var);initial_S_Mat = 0;
}

void CData::SetData(Matrix &X_, Matrix &Edit_, Matrix &LogB_, int n_blanceedit) {
  n_var = X_.ncols();
  n_sample= X_.nrows();
  n_EditVec = Edit_.nrows();
  n_balance_edit = n_blanceedit;
  
  EditMat = Edit_.columns(1,n_var);
  EditVec = Edit_.column(n_var + 1);
  
  logB_L = LogB_.column(1);
  logB_U = LogB_.column(2);
  logB_U_L = logB_U - logB_L;
  log_logB_U_L = log_ColumnVector(logB_U_L);
  
  D_Observed = X_;
  
  n_tau = ipow(2,n_var)-2;
}

// Commented on 05/21/2015
// bool CData::ReadData() {
//   ReadDimInfo();
//   Matrix Edit = ReadMatrix("Edit.txt",n_EditVec,n_var + 1);
//   EditMat = Edit.columns(1,n_var);
//   EditVec = Edit.column(n_var + 1);
//
//   Matrix logB = ReadMatrix("logB.txt",n_var,2);
//   logB_L = logB.column(1);
//   logB_U = logB.column(2);
//   logB_U_L = logB_U - logB_L;
//   log_logB_U_L = log_ColumnVector(logB_U_L);
//
//   D_Observed = ReadMatrix("D_obs.txt",n_sample,n_var);
//   return true;
// }

// Commented on 05/21/2015
// bool CData::ReadOptionalData() {
//   has_trueS = true;
//   try {
//     True_S_Mat = ReadMatrix("True_S_Mat.dat",n_faulty,n_var).t();//optional?!
//   } catch (...) {
//     has_trueS = false;
//   }
//   if (has_initialValue) {
//     //initilize intial imputed data for records
//     D_initial = D_Observed;
//     Matrix temp = ReadMatrix("D_initial.dat",n_faulty,n_var);
//     for (int i = 1; i <= n_faulty; i++ ) {
//       D_initial.row(Faulty2Original[i-1]) = temp.row(i);
//     }
//     // for faulty values, import its initial values from R
//     initial_S_Mat = ReadMatrix("S_initial.dat",n_faulty,n_var);
//   }
//   return true;
// }

int CData::get_feasible_tau(ColumnVector &s) {
  int skip_for_infeasible_type2 = 0;
  for (int be = 1; be <= n_balance_edit; be++) {
    double sum = 0.0;
    double sum1 = 0.0;
    for (int j = 1; j <= BalanceEdits[be-1].nrows(); j++) {
      sum += s(BalanceEdits[be-1](j));
      if (j > 1) {sum1 += s(BalanceEdits[be-1](j));}
    }
    if (sum == 1) { skip_for_infeasible_type2 = 1 ;} 
    if (s(BalanceEdits[be-1](1))==0 && sum1 > 0) {skip_for_infeasible_type2 = 1;}
  }
  ColumnVector s_temp = subvector_by_index(s, Index_for_BE_sum_and_not_in_BE);
  if (s_temp.sum() == 0 ) {skip_for_infeasible_type2 = 1;}
  return (1 - skip_for_infeasible_type2);
}

void CData::initilize_obs_fail_matrix() {
  int p = n_balance_edit+1;
  int n_non_balance_edit = n_EditVec - 2 * n_balance_edit;
  obs_edit_fail = Matrix(n_sample,p+1); obs_edit_fail = 0;
  Matrix mat_fail_edits = EditMat * D_Observed.t();
  n_faulty = 0; 
  for (int i = 1; i <= n_sample; i++){	
    ColumnVector fail = mat_fail_edits.column(i) - EditVec;
    bool balance_edit_fail = false;
    for (int be = 1; be <= n_balance_edit; be++) {
      int base = n_non_balance_edit + be * 2 - 1;
      if (fail(base) > 0 || fail(base+1) > 0) {
	      obs_edit_fail(i,be) = 1;
	      balance_edit_fail = true;
      }
    }
    for (int v = 1; obs_edit_fail(i,p) == 0 && v <= n_non_balance_edit; v++) {
      if (fail(v) > 0) obs_edit_fail(i,p) = 1;
    }
    if (balance_edit_fail) {
      obs_edit_fail(i,p+1) = 1 ;
    } else {
      if (obs_edit_fail(i,p) == 1) {
	      obs_edit_fail(i,p+1) = 2;
      } else {
	      obs_edit_fail(i,p+1) = 0; 
      } 
    }
    if (obs_edit_fail(i,p+1) >0) {
      n_faulty++;
    }
  }    
  //build orignal to faulty map
  Faulty2Original = vector<int>(n_faulty);
  Original2Faulty = vector<int>(n_sample);
  for (int i = 1, count = 0; i <= n_sample; i++){
    if (obs_edit_fail(i,p+1) > 0) {
      Faulty2Original[count++] = i;
      Original2Faulty[i-1] = count;
    } else {
      Original2Faulty[i-1] = 0;
    }
  }
}

void CData::initilize_balance_edits() {
  int NonBalancedEdit = n_EditVec - 2 * n_balance_edit;
  Matrix BalanceEditMat = EditMat.rows(NonBalancedEdit+1,n_EditVec);
    
  BalanceEdits.reserve(n_balance_edit);
  BalanceEdits_coeff.reserve(n_balance_edit); // ADDED 02/02/2015
    
  for (int be = 1; be <= n_balance_edit; be++) {
      
    ColumnVector balance_edit = BalanceEditMat.row((be - 1)*2 + 1).t();
    int n_positive = 0, n_negative =0;
    for (int i = 1; i <= n_var; i++) {
      if (balance_edit(i) > 0) n_positive++;
      if (balance_edit(i) < 0) n_negative++;
    } // for (int i<=n_var)
      
    if ((n_positive ==  1 && n_negative > 1) || (n_positive >  1 && n_negative == 1)) { 
      if (n_negative == 1) { //flip
	for (int i = 1; i <= n_var; i++) {
	  balance_edit(i) = -balance_edit(i);
	}
	n_negative = n_positive; n_positive = 1;
      }
      
      ColumnVector id_bal(n_positive+n_negative);
      ColumnVector coef_bal(n_positive+n_negative); // ADDED 02/02/2015
      for (int i = 1; i <= n_var && n_positive >0; i++) {
	if (balance_edit(i) > 0) {
	  id_bal(1) = i;
	  coef_bal(1) = balance_edit(i); // ADDED 02/02/2015 
	  n_positive--;
	}
      }
      for (int i = n_var; i >=1 && n_negative >0; i--) {
	if (balance_edit(i) < 0) {
	  id_bal(1+n_negative) = i;
	  coef_bal(1+n_negative) = balance_edit(i); // ADDED 02/02/2015
	  n_negative--;
	}
      }
      BalanceEdits.push_back(id_bal);
      BalanceEdits_coeff.push_back(coef_bal); // ADDED 02/02/2015
	
    }  else {
      //should be checked long before here, but just in case
	  Rprintf( "Balance Edit %d  not recognized\n", be);
      //cout << "Balance Edit " << be << " not recognized" << endl; 
    }
      
  }
  build_index_for_balance_edit();
}

double CData::get_logUnif_y_tilde(ColumnVector &s, int i_tau, int i_original) {
  double logUnif_y_tilde = 0;
  if (is_case(i_original,1)) {
    ColumnVector item_by_runif = get_var_by_runif(i_original,s);
    for (int j=1; j<=n_var; j++){
      if (item_by_runif(j)==1){
	      logUnif_y_tilde -= log_logB_U_L(j);
      }
    } 
  } else {
    if (is_case(i_original,2)) {
      logUnif_y_tilde = logUnif_case2(Original2Faulty[i_original-1],i_tau);  
    }
  }
  return  logUnif_y_tilde;
  
}

void CData::set_logUnif_case2(int i_original, int i_tau, double v) {
  logUnif_case2(Original2Faulty[i_original-1],i_tau) = v;
}

// bool CData::ReadDimInfo() { 	// Commented on 02/08/2020
//   FILE *file_import;
//   file_import = fopen(fn_makefilename(data_folder,"dimension.txt").c_str(),"r");
//   fscanf(file_import, "%d", &n_sample);
//   fscanf(file_import, "%d", &n_var);
//   fscanf(file_import, "%d", &n_EditVec);
//   fscanf(file_import, "%d", &n_balance_edit);
//   fscanf(file_import, "%f", &seed_no);
//   fscanf(file_import, "%f", &epsilon);
//   fclose(file_import);
//   n_tau = ipow(2,n_var)-2;
//   return true;
// }

// Matrix CData::ReadMatrix(string file,int row,int col) { 	// Commented on 02/08/2020
//   FILE *file_import;
//   float temp_import;
//   Matrix mat(row,col);
//   file_import = fopen(fn_makefilename(data_folder,file).c_str(),"r");
//   for (int i=1; i<=row; i++){
//     for (int j=1; j<=col; j++){
//       fscanf(file_import, "%f", &temp_import);
//       mat(i,j)=temp_import;
//     }
//   }
//   fclose(file_import);
//   return mat;
// }

// ColumnVector CData::ReadVec(string file,int row) {		// Commented on 02/08/2020
//   FILE *file_import;
//   float temp_import;
//   ColumnVector vec(row);
//   file_import = fopen(fn_makefilename(data_folder,file).c_str(),"r");
//   for (int i=1; i<=row; i++){
//     fscanf(file_import, "%f", &temp_import); vec(i)=temp_import;
//   }
//   fclose(file_import);
//   return vec;
// }

bool CData::is_case(int i_original,int editcase) {//editcase 0,1,2
  return obs_edit_fail(i_original,n_balance_edit+2) == editcase;
}

bool CData::PassEdits(ColumnVector &x) {
  ColumnVector mat_x_q = EditMat * x - EditVec ;
  /*
    if (Debug) {
    cout << "mat_x_q" <<setw(15) << setprecision(10) << scientific<< mat_x_q << endl;
    Debug = false;
    }
  */
  return mat_x_q.maximum()<=0;  
  //return mat_x_q.maximum()<0;  
}

// // CHANGED by Hang, 2014/12/29
//void CData::set_balance_edit_values_for_x_q(ColumnVector &x_q, ColumnVector &item_by_bal) {
//   for (int be = 1; be <= n_balance_edit; be++) {
//    double sum = 0.0;
//    for (int j = 2; j <= BalanceEdits[be-1].nrows(); j++) {
//      sum += x_q(BalanceEdits[be-1](j));
//    }
//    if ( item_by_bal(be) == BalanceEdits[be-1](1) ){  
//        x_q(BalanceEdits[be-1](1)) = sum;
//    } else if ( item_by_bal(be) > 0 ){
//  		  x_q(item_by_bal(be)) = x_q(BalanceEdits[be-1](1))-sum+x_q(item_by_bal(be)) ;
//  	}
//  }
//}

// // CHANGED by Hang, 2014/12/29
void CData::set_balance_edit_values_for_x_q(ColumnVector &s_q, ColumnVector &x_q, ColumnVector &item_by_bal) {

  for (int be = 1; be <= n_balance_edit; be++){

    // MODIFIED 2/2/2015
    
    ColumnVector s_q_bal = subvector_by_index(s_q, BalanceEdits[be-1]);
    
    if ( (s_q_bal.sum()>=1) && (item_by_bal(be)!=BalanceEdits[be-1](1)) ){

      int which_item_by_bal ; 
      double sum = 0.0 ;
      for (int j=2; j<=BalanceEdits[be-1].nrows(); j++){
	if ( item_by_bal(be)==BalanceEdits[be-1](j) ){
	  which_item_by_bal = j ; 
	} else {
	  sum += x_q(BalanceEdits[be-1](j)) * (-1.0*BalanceEdits_coeff[be-1](j)) ; 
	}
      }
      x_q(item_by_bal(be)) = ( x_q(BalanceEdits[be-1](1)) * BalanceEdits_coeff[be-1](1) - sum ) / (-1.0*BalanceEdits_coeff[be-1](which_item_by_bal)) ;

    } // if ( (s_q_bal.sum()>=1) && (item_by_bal(be)!=BalanceEdits[be-1](1)) )
      
  } // for (int be = 1; be <= n_balance_edit; be++)


  for (int be = 1; be <= n_balance_edit; be++){
    double sum = 0.0;
    for (int j = 2; j <= BalanceEdits[be-1].nrows(); j++) {
      sum += x_q(BalanceEdits[be-1](j)) * (-1.0*BalanceEdits_coeff[be-1](j)) ; // MODIFIED 2/2/2015
    }
    ColumnVector s_q_bal = subvector_by_index(s_q, BalanceEdits[be-1]);
    if ( (s_q_bal.sum()>=1) && (item_by_bal(be)==BalanceEdits[be-1](1)) ){
      x_q(BalanceEdits[be-1](1)) = sum ;
    } // if ( s_q_bal.sum() > 1 )
  } // for (int be = 1; be <= n_balance_edit; be++)
  
} // void CData::set_balance_edit_values_for_x_q(ColumnVector &x_q, ColumnVector &item_by_bal)


ColumnVector CData::get_compact_vector(ColumnVector &full) {
  ColumnVector compact = Index_for_compact;
  for (int i = 1; i <= Index_for_compact.nrows(); i++) {
    compact(i) = full(Index_for_compact(i));
  }
  return compact;
}

void CData::update_full_x_for_balance_edit(ColumnVector &x_full) {
  for (int be = 1; be <= n_balance_edit; be++) {
    double sum = 0.0;
    for (int j = 2; j <= BalanceEdits[be-1].nrows(); j++) {
      sum += x_full(BalanceEdits[be-1](j))  * (-1.0*BalanceEdits_coeff[be-1](j)) ; // MODIFIED 2/2/2015
    }
    x_full( BalanceEdits[be-1](1)) = sum;
  }
}

void CData::build_index_for_balance_edit() {
  ColumnVector Index(n_var); Index = 0;
  for (int be = 1; be <= n_balance_edit; be++) {
    Index(BalanceEdits[be-1](1)) = 1; 
  }
  Index_for_compact = ColumnVector(n_var - n_balance_edit);
  for (int i = 1, count =0; i <= n_var; i++) {
    if (Index(i) == 0) {
      Index_for_compact(++count) = i; 
    }
  }
  
  Index = 0;
  Index_for_BE_sum_and_not_in_BE = ColumnVector(n_balance_edit);
  for (int be = 1; be <= n_balance_edit; be++) {
    Index_for_BE_sum_and_not_in_BE(be) = BalanceEdits[be-1](1);
    for (int j = 1; j <= BalanceEdits[be-1].nrows(); j++) {
      Index(BalanceEdits[be-1](j)) = 1;
    }
  }
  int n_temp = Index.sum();
  if (n_temp < n_var) {
    Index_for_var_not_in_balance_edit = ColumnVector(n_var - n_temp);
    for (int i = 1, count =0; count < n_var - n_temp; i++) {
      if (Index(i) < 1) {
        Index_for_var_not_in_balance_edit(++count) = i;
      }
    }
    Index_for_BE_sum_and_not_in_BE &= Index_for_var_not_in_balance_edit;//concatenate
    //cout << Index_for_BE_sum_and_not_in_BE << endl;
  }
}

// whether drawn by item_by_rnorm and balance edits
ColumnVector CData::get_item_by_norm_indicator(ColumnVector &s_q, ColumnVector &item_by_bal) {
  
  ColumnVector item_by_rnorm(n_var) ; item_by_rnorm = 0; 
  item_by_bal = ColumnVector(n_balance_edit); item_by_bal = 0;
    
  for (int be = 1; be <= n_balance_edit; be++) {
    ColumnVector s_q_bal = subvector_by_index(s_q, BalanceEdits[be - 1]);
    if ( s_q_bal.sum() > 1 ){
      int is_balance_decided = 0 ; 
      for (int j_temp=1; j_temp<=BalanceEdits[be - 1].nrows(); j_temp++){
	if ( s_q_bal(j_temp)==1 ){
	  if ( is_balance_decided==0 ){
	    item_by_bal(be) = BalanceEdits[be - 1](j_temp) ; 
	    is_balance_decided = 1 ; 
	  } else {
	    item_by_rnorm(BalanceEdits[be - 1](j_temp))=1 ; 
	  }
	} 
      }
    } else if ( s_q_bal.sum() == 1 ){
      for (int j_temp=1; j_temp<=BalanceEdits[be - 1].nrows(); j_temp++){
	if ( s_q_bal(j_temp)==1 ){
	  item_by_bal(be) = BalanceEdits[be - 1](j_temp) ;  
	}
      }
    } 
  }
    
  // for varaibels not in balance edits
  for (int i = 1; i <= Index_for_var_not_in_balance_edit.nrows();i++) {
    if (s_q(Index_for_var_not_in_balance_edit(i)) == 1) {
      item_by_rnorm(Index_for_var_not_in_balance_edit(i)) = 1;
    } 
  }
  return item_by_rnorm;
    
  // ADDED by Hang, 2014/12/29
  for (int i_var=1; i_var<=n_var; i_var++){
    for (int be = 1; be <= n_balance_edit; be++) {
      if (item_by_rnorm(i_var)==1){
	if (item_by_bal(be)==i_var){
	  item_by_rnorm(i_var) = 0; 
	} // if
      } // if
    } // for (int be = 1; be <= n_balance_edit; be++)
  } // for (int i_var=1; i_var<=n_var; i_var++)
    
} // ColumnVector CData::get_item_by_norm_indicator(ColumnVector &s_q, ColumnVector &item_by_bal)

bool CData::PassStep0(ColumnVector &s, int i_original) {
  for (int be = 1; be <= n_balance_edit; be++) {
    ColumnVector s_bal = subvector_by_index(s, BalanceEdits[be - 1]);
    if ( (obs_edit_fail(i_original,be)==1) && (sum(s_bal)==0) ) return false; 
    if ( (obs_edit_fail(i_original,be)==0) && (s_bal(1)!=(s_bal.rows(2,s_bal.nrows())).maximum()) ) return false;
  }
  ColumnVector s_temp = subvector_by_index(s, Index_for_BE_sum_and_not_in_BE);
  if ( (obs_edit_fail(i_original,n_balance_edit+1)==1) && sum(s_temp)==0 ) return false;
  return true;
}

bool CData::PassStep0_for_init(ColumnVector X_tilde_i, ColumnVector &s, float epsilon) {
  for (int be = 1; be <= n_balance_edit; be++) {
    ColumnVector s_bal = subvector_by_index(s, BalanceEdits[be - 1]);
    if (sum(s_bal) == 1) {
      // double temp = -X_tilde_i(BalanceEdits[be - 1](1));
      double temp = X_tilde_i(BalanceEdits[be-1](1)) * BalanceEdits_coeff[be-1](1) ; // CHANGED 2/2/2015
      for (int j = 2; j <= BalanceEdits[be-1].nrows(); j++){
	temp += X_tilde_i(BalanceEdits[be-1](j)) * BalanceEdits_coeff[be-1](j) ; // CHANGED 2/2/2015
      }
      if (temp < 0) temp = -temp;
      if (temp <= epsilon) return false;
    }
  }
  return true;
}


bool CData::PassStep1(ColumnVector &s, Matrix &A_0, ColumnVector &x_0) {
  ColumnVector which_rows_for0(n_EditVec) ; which_rows_for0 = 1 ; 
  for (int i_row=1; i_row<=n_EditVec; i_row++){
    for (int j_var=1; j_var<=n_var; j_var++){
      if ( (s(j_var)==1) && (EditMat(i_row,j_var)!=0) ){
	which_rows_for0(i_row) = 0 ; // eliminate a row of EditMat which is associated with s_ij=1
      }
    }	
  }
	
  int n_EditVec_small_for0 = which_rows_for0.sum() ; 
  Matrix EditMat_small_for0(n_EditVec_small_for0,x_0.nrows()) ; EditMat_small_for0 = 0 ; 
  ColumnVector EditVec_small_for0(n_EditVec_small_for0) ; EditVec_small_for0 = 0 ; 
	
  int count_row = 0 ; 
  for (int i_row=1; i_row<=n_EditVec; i_row++){
    if (which_rows_for0(i_row)==1){
      count_row++; 
      EditVec_small_for0(count_row) = EditVec(i_row) ; 
      EditMat_small_for0.row(count_row) = A_0.row(i_row) ; 
    }
  } 
  ColumnVector test_Edit_small_for0 = EditMat_small_for0*x_0-EditVec_small_for0; 
  return  (test_Edit_small_for0.maximum() <= 0); 
}

ColumnVector CData::copy_non_balance_edit(ColumnVector &orig) {
  ColumnVector result(n_var); 
  result = 0;
  for (int i = 1; i<= Index_for_compact.nrows();i++) {
    result(Index_for_compact(i)) = orig(Index_for_compact(i));
  }
  return result;
}

string CData::fn_makefilename(string folder, string file) {
  string output = folder;
  return output.append(file);
}

void CData::UpdateCompactMatrix(Matrix &Compact, Matrix &Full) {
  for (int i = 1; i<= Index_for_compact.nrows();i++) {
    Compact.column(i) = Full.column(Index_for_compact(i));
  }
}

void CData::UpdateFullVector(ColumnVector &Compact, ColumnVector &Full) {
  for (int i = 1; i<= Index_for_compact.nrows();i++) {
    Full.row(Index_for_compact(i)) = Compact.row(i);
  }
}

void CData::UpdateCompactVector(ColumnVector &Compact, ColumnVector &Full) {
  for (int i = 1; i<= Index_for_compact.nrows();i++) {
    Compact.row(i) = Full.row(Index_for_compact(i));
  }
}

//add for satisfied balance 
void CData::AddjustCandidate_s_i_ForSatisfiedBalanceEdit(ColumnVector &cand_s, int i_original) {
  for (int be = 1; be <= n_balance_edit; be++) {
    if (obs_edit_fail(i_original,be)==0) {
      ColumnVector s_bal = subvector_by_index(cand_s, BalanceEdits[be - 1]);
      cand_s(BalanceEdits[be - 1](1)) = (s_bal.rows(2,s_bal.nrows())).maximum(); 
    }
  }
}

ColumnVector CData::get_free_s_i(int i_original, ColumnVector &s) {
  ColumnVector free_s_i = s;
  for (int be = 1; be <= n_balance_edit; be++) {
    if (obs_edit_fail(i_original,be)==0) { free_s_i(BalanceEdits[be - 1](1)) = 9;}
  }
  return free_s_i;
}

ColumnVector CData::get_var_by_runif(int i_original, ColumnVector &s) {
  ColumnVector item_by_runif = s; 
  for (int be = 1; be <= n_balance_edit; be++) {
    if (obs_edit_fail(i_original,be)==0) { item_by_runif(BalanceEdits[be - 1](1)) = 0;}
  }
  return item_by_runif;
}

bool CData::InitialRecordValid( ) {
  int NonBalancedEdit = n_EditVec - 2 * n_balance_edit;
  Matrix RatioEditMat = EditMat.rows(1,NonBalancedEdit);
  ColumnVector RatioEditVec = EditVec.rows(1,NonBalancedEdit);
  bool ret = true;
  for (int i_sample=1; i_sample<=n_sample; i_sample++){
    ColumnVector x_i = (D_initial.row(i_sample)).t() ; 
    if (!PassEdits(x_i)) {
      ret = false;
      ColumnVector mat_Ratio_check = RatioEditMat * x_i - RatioEditVec ;  
	  Rprintf( "Error: edit-fail, initilized record\n");
      
	  Rprintf( " i = %d , ratio = %g \n",i_sample,   mat_Ratio_check.maximum());
      //cout << " i = " << i_sample << ", ratio = " << mat_Ratio_check.maximum() << endl ;
      for (int be = 1; be <= n_balance_edit; be++) {
        double temp_diff = x_i(BalanceEdits[be-1](1)) * BalanceEdits_coeff[be-1](1) ; // MODIFIED 2/2/2015
        for (int j = 2; j <= BalanceEdits[be-1].nrows(); j++) {
          temp_diff += x_i(BalanceEdits[be-1](j)) * BalanceEdits_coeff[be-1](j) ; // MODIFIED 2/2/2015
        }
		Rprintf( "diff = %d = %g\n",be,   temp_diff);
      }
      //cout << x_i.t() << endl ; //need to put this line in later on
    }  
  }  
  return ret;
}

int CData::count_Draw_find_s(Matrix &S) {
  if (!has_trueS) return -1; 
  int count = 0;
  for (int i = 1; i <= S.ncols(); i++) {
    if (S.column(i) == True_S_Mat.column(i)) {
      count++;
    }
  }
  return count;
}
int CData::ipow(int base, int exp) {
  int result = 1;
  while (exp) {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
}
