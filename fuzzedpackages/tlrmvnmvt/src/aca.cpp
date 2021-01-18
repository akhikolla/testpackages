#include <vector>
#include "aca.h"

using namespace std;
using namespace Eigen;

/*
	Suppose U, V in tlrNode have extra cols
	Chg (increase) the ncol of U, V only if necessary
*/
int aca(std::function<double(int,int)> f, int sz, double epsl, TLRNode &tlrNode)
{
        vector<VectorXd> u_vec,v_vec;
        VectorXi seq_vec(sz); for(int l = 0; l < sz; l++) seq_vec(l) = l;
        VectorXd R_diag = seq_vec.unaryExpr([&f](int i){return f(i,i);});
        VectorXd R_sec_diag = seq_vec.unaryExpr([&f,&sz](int i){return f(i,sz-1-i);});
        int row_idx,col_idx,k;
        k = 0;
        row_idx = 0;
        // adjust epsilon
        epsl = epsl*R_diag.cwiseAbs().maxCoeff() > 1e-8? epsl*R_diag.cwiseAbs().
                maxCoeff():1e-8;
        function<double(int)> f_row = [&f,&row_idx](int j){return f(row_idx,j);};
        function<double(int)> f_col = [&f,&col_idx](int i){return f(i,col_idx);};
        while(k < sz)
        {
                // crt row
                VectorXd row_i = seq_vec.unaryExpr(f_row);
                for(int l = 0; l < k; l++)
                        row_i.noalias() = row_i - u_vec[l](row_idx) * v_vec[l];
                row_i.cwiseAbs().maxCoeff(&col_idx);
                // max coef idx
                double max_row_coef = row_i[col_idx];
                // stop condition 
                // if crt row is close to 0, check the primary and secondary 
                //      diagonal
                if(abs(max_row_coef) < epsl)
                {
                        double max_diag_coef = R_diag.cwiseAbs().maxCoeff(&row_idx);
                        if(max_diag_coef < epsl)
                        {
                                double max_sec_diag_coef = R_sec_diag.cwiseAbs().
                                        maxCoeff(&row_idx);
                                if(max_sec_diag_coef < epsl)
                                        break;
                        }
                }
                else
                {
                        row_i = row_i / max_row_coef;
                        VectorXd col_j = seq_vec.unaryExpr(f_col);
                        for(int l = 0; l < k; l++)
                                col_j = col_j - v_vec[l](col_idx) * u_vec[l];
                        u_vec.push_back(col_j);
                        v_vec.push_back(row_i);
                        R_diag.noalias() = R_diag - row_i.cwiseProduct(col_j);
                        R_sec_diag.noalias() = R_sec_diag - col_j.cwiseProduct(
                                row_i.reverse());
                        k++;
                }
        }
        if(k==0)
        {
		if(tlrNode.maxColNum < 1)
		{
//	                Rprintf("Warning: the TLRNode needs resizing in TLR "
//				"addition. "
//        	                "Consider enlarge their initial allocations\n");
			tlrNode.U.resize(sz,1);
			tlrNode.V.resize(sz,1);
			tlrNode.maxColNum = 1;
		}
		tlrNode.crtColNum = 1;
		tlrNode.U.col(0).setZero();
		tlrNode.V.col(0).setZero();
        }
        else
        {
 		if(tlrNode.maxColNum < k)
		{
//	                Rprintf("Warning: the TLRNode needs resizing in TLR "
//				"addition. "
//        	                "Consider enlarge their initial allocations\n");
			tlrNode.U.resize(sz,k);
			tlrNode.V.resize(sz,k);
			tlrNode.maxColNum = k;
		}
		tlrNode.crtColNum = k;
                for(int l = 0; l < k; l++)
                {
                        tlrNode.U.col(l) = u_vec[l];
                        tlrNode.V.col(l) = v_vec[l];
                }
        }
        return k;
}

