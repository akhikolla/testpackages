
#if !defined(_CDATA_H)
#define _CDATA_H
#include "CHeader.h"

class CData {
  public:
		
	  // CData(string folder); //Constructor   // Commented on 05/21/2015
    
    CData(); //Constructor
	  ~CData(); //destructor
    bool Debug;
    int n_sample, n_var, n_EditVec, n_balance_edit;
    float seed_no, epsilon;
    
    //derived values
    int n_faulty;
    int n_tau; // exclude (00...0) (11...1)
    
    Matrix D_Observed;
    Matrix log_D_Observed; //for optimization only
    
    Matrix EditMat;
    ColumnVector EditVec;
    
    ColumnVector logB_U_L; //B-L, for optimization only
    ColumnVector log_logB_U_L; //for optimization only
    ColumnVector logB_L;
    ColumnVector  logB_U;
    
    bool has_trueS;
    bool has_initialValue;
    
    
    vector<int> Faulty2Original;
    vector<int> Original2Faulty;
    
    Matrix obs_edit_fail;
    
    Matrix logUnif_case2;
    
    Matrix True_S_Mat;  // not used in model fitting, optional 
    
    void init();
    // void init(string folder);	// Commented on 05/21/2015
    int ipow(int base, int exp);
    
    string fn_makefilename(string folder, string file);
    bool InitialRecordValid();
    void UpdateCompactVector(ColumnVector &Compact, ColumnVector &Full);
    void UpdateCompactMatrix(Matrix &Compact, Matrix &Full);
    void UpdateFullVector(ColumnVector &Compact, ColumnVector &Full);
    bool PassEdits(ColumnVector &x);
    bool PassStep0(ColumnVector &s, int i_original);
    bool PassStep0_for_init(ColumnVector X_tilde_i, ColumnVector &s, float epsilon);
    bool PassStep1(ColumnVector &s, Matrix &A_0, ColumnVector &x_0);
    ColumnVector get_var_by_runif(int i_original, ColumnVector &s);
    ColumnVector get_free_s_i(int i_original, ColumnVector &s);
    ColumnVector get_item_by_norm_indicator(ColumnVector &s_q, ColumnVector &item_by_bal);
    // void set_balance_edit_values_for_x_q(ColumnVector &x_q, ColumnVector &item_by_bal); // CHANGED by Hang, 2014/12/29
    void set_balance_edit_values_for_x_q(ColumnVector &s_q, ColumnVector &x_q, ColumnVector &item_by_bal); // CHANGED by Hang, 2014/12/29
    void AddjustCandidate_s_i_ForSatisfiedBalanceEdit(ColumnVector &cand_s, int i_original);
    void update_full_x_for_balance_edit(ColumnVector &x_full);
    ColumnVector copy_non_balance_edit(ColumnVector &orig);
    ColumnVector get_compact_vector(ColumnVector &full);
    
    // Matrix ReadMatrix(string file,int row,int col);  // Commented on 02/08/2020
    // ColumnVector ReadVec(string file,int row);  // Commented on 02/08/2020
    bool is_case(int i_original,int editcase);
    void set_logUnif_case2(int i_original, int i_tau, double v);
    double get_logUnif_y_tilde(ColumnVector &s, int i_tau, int i_original);
    int count_Draw_find_s(Matrix &S);
    int get_feasible_tau(ColumnVector &s);
    Matrix D_initial; // only used to initilize y_in
    Matrix initial_S_Mat; //only used for initilize S_mat
    
    void SetData(Matrix &X_, Matrix &Edit, Matrix &LogB_, int n_blanceedit);
  private:
    string data_folder;
    ColumnVector Index_for_compact;
    ColumnVector Index_for_var_not_in_balance_edit;
    ColumnVector Index_for_BE_sum_and_not_in_BE;
    vector<ColumnVector> BalanceEdits;
    vector<ColumnVector> BalanceEdits_coeff; // ADDED 2/2/2015
    
    // bool ReadData();		// Commented on 05/21/2015
    // bool ReadOptionalData();		// Commented on 05/21/2015
    // bool ReadDimInfo();			// Commented on 02/08/2020
    void build_index_for_balance_edit();
    void initilize_balance_edits();
    void initilize_obs_fail_matrix();
};

#endif  //_CDATA_H
