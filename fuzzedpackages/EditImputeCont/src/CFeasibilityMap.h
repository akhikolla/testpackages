
#if !defined(_CFEASIBILITY_H)
#define _CFEASIBILITY_H

class MapManager {
public:
  MapManager(CData &Data);
  CData& data;
  
  int iter;				//the current iter number, passed to this class so I don't have to change code everywhile  
	vector<int> last_trim; //keep track of iter numbers that have been trimed most recently 
	vector< map<int,int> > map_feasible_tau_t;	// mapping (int,int) to the vector element
	vector<int> last_trim_fn2;
	vector< map<int,int> > map_count_x_out; // for another mapping
		
};
class CFeasibilityMap {
  public:
    bool useMap;
	  Matrix feasibleMap;
    bool Debug;
    
    CFeasibilityMap();
    ~CFeasibilityMap();
    
    void Build(CData &Data);
    void initilize_D_and_S(CData &Data);
    int isCandidateFeasible(ColumnVector &cand_s, int i_faulty);
    static ReturnMatrix tau_to_s_fn( double tau_i, int n_var ) ;
    static ReturnMatrix tau_to_s_fn2( double tau_i, int n_var ) ;
    
    static int s_to_tau_fn( ColumnVector &s_i );
    double Simulate_logUnif_case2(int i_tau, int i_original,int n_simul,Uniform &randUnif);
    int count_x_out_fn(CData &Data,int i_tau, int i_original,int n_simul, Uniform &randUnif);
    int Map_count_x_out_fn(int i_tau, int i_original,int n_simul, Uniform &randUnif);
    void Simulate_logUnif_case2(int n_simul, Uniform &randUnif,CData &Data);
    int EvaluateMove(int i_original, CData &Data, ColumnVector &s_i, int iter, int &what_type_move, 
    double &g_option_q, double &g_mode_q, bool SampleMove);
    int test_feasible_test_fn(CData &Data, ColumnVector &x_tilde_i, ColumnVector &s_i);
    // void test(CData &Data);	// Commented on 05/21/2015
    MapManager* pmm;
    
  private:
    int feasible_test_fn(CData &Data, ColumnVector &x_tilde_i, ColumnVector &s_i, int i_original, bool initD_S,float epsilon,ColumnVector &x);
    int feasible_test_fn(CData &Data, ColumnVector &x_tilde_i, ColumnVector &s_i, int i_original,float epsilon);
  
    bool SolveLP(Matrix & A, ColumnVector &b);
    bool SolveLP(Matrix & A, ColumnVector &b, bool domax);
    bool SolveLP(Matrix &A, ColumnVector &b, ColumnVector &x); //for initilization 
    ColumnVector get_feasible_tau(CData &Data);
    ColumnVector get_order_to_test(int n_tau, int n_var);
    
    
};

#endif  //_CFEASIBILITY_H
