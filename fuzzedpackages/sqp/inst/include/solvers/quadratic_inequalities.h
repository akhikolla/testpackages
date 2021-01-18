#ifndef __sqp_solvers_quadratic_inequalities_included__        // if include guard for 'solvers/quadratic_inequalities.h' is undefined
#define __sqp_solvers_quadratic_inequalities_included__        // define include guard for 'solvers/quadratic_inequalities.h'  

#include "sqp.h"

namespace sqp {
namespace solvers{



inline sqp::solvers::result inequalities(
    arma::vec                 x,       // Optimization variable (initial value)
    const arma::mat           &Q,      // Quadratic Multiplier
    const arma::mat           &C_eq,   // Equality  Constraint Multiplier
    const arma::mat           &C_ineq, // Inquality Constraint Multiplier
    const arma::vec           &l,      // Linear Multiplier
    const arma::vec           &t_eq,   // Equality   Constraint RHS
    const arma::vec           &t_ineq, // Inequality Constraint RHS (upper bound)
    const double              tol       = 1e-7,
    const unsigned            max_iter  = 500,
    int                       dim_eq    = -1,
    int                       dim_ineq  = -1,
    int                       dim_Q     = -1,
    const unsigned            solver    =  0,
    const bool                fast      = false,
    const bool                debug     = false)
{
  // -------------------------------------------------------
  sqp::misc::section("Inequality constrained sqp (dense)", debug);
  // -------------------------------------------------------
  
  if(dim_eq < 0)
    dim_eq = C_eq.n_rows;
  
  if(dim_ineq < 0)
    dim_ineq = C_ineq.n_rows;
  
  if(dim_Q < 0)
    dim_Q = Q.n_rows;
  
  arma::uvec ineq_active, ineq_inactive;
  
  // Rcpp::Rcout << "t_ineq pre: " << t_ineq(0) << "\n";
  // 
  // t_ineq.replace(-HUGE_VAL,std::numeric_limits<long int>::min());
  // t_ineq.replace(HUGE_VAL,std::numeric_limits<long int>::max());
  // 
  // Rcpp::Rcout << "t_ineq post: " << t_ineq(0) << "\n";
  
  {
    // -------------------------------------------------------
    sqp::misc::section("Check constraints", debug);
    // -------------------------------------------------------
    
    const arma::vec equalities_valid = C_eq * x;
    const arma::vec inequalities_valid = C_ineq * x;
    
    ineq_active = arma::find((inequalities_valid - t_ineq)>=-tol);
    ineq_inactive = arma::find((inequalities_valid - t_ineq)<-tol);
    
    const arma::uvec equalities_invalid = arma::find(arma::abs(equalities_valid - t_eq)>tol);
    const arma::uvec inequalities_invalid = arma::find((inequalities_valid - t_ineq)>tol);
    
    if(equalities_invalid.n_elem > 0)
    {
      std::string invalid_values = "The following equalities are violated by the initial x-Values:\n     ";
      
      for(unsigned k=0;k<equalities_invalid.n_elem;++k)
      {
        invalid_values += std::to_string(arma::as_scalar(equalities_invalid.row(k)));
        
        if(k < equalities_invalid.n_elem-1)
          invalid_values += ", "; 
      }
      
      invalid_values += "\n";
      
      arma::mat err_log(C_eq.n_rows,C_eq.n_cols + 3,arma::fill::zeros);
      err_log(0,0,arma::size(C_eq)) = C_eq;
      err_log(0,C_eq.n_cols,arma::size(equalities_valid)) = equalities_valid;
      err_log(0,C_eq.n_cols+1,arma::size(t_eq)) = t_eq;
      err_log(0,C_eq.n_cols+2,arma::size(t_eq)) = equalities_valid - t_eq;
      
      Rcpp::Rcout << "This are the erroneous equality constraints (last 3 columns: value / target / delta):\n" 
                  << err_log.rows(equalities_invalid) << "\n";
      
      
      Rcpp::stop(invalid_values);
    }
    
    
    if(inequalities_invalid.n_elem > 0)
    {
      arma::mat ineq_viol_output(inequalities_invalid.n_elem,3,arma::fill::zeros);
      
      std::string invalid_values = "The following inequalities are violated by the initial x-Values:\n     ";
      
      for(unsigned k=0;k<inequalities_invalid.n_elem;++k)
      {
        invalid_values += std::to_string(arma::as_scalar(inequalities_invalid.row(k)));
        
        if(k < inequalities_invalid.n_elem-1)
          invalid_values += ", "; 
      }
      
      invalid_values += "\n";
      
      arma::mat err_log(C_ineq.n_rows,C_ineq.n_cols + 3,arma::fill::zeros);
      err_log(0,0,arma::size(C_ineq)) = C_ineq;
      err_log(0,C_ineq.n_cols,arma::size(inequalities_valid)) = inequalities_valid;
      err_log(0,C_ineq.n_cols+1,arma::size(t_ineq)) = t_ineq;
      err_log(0,C_ineq.n_cols+2,arma::size(t_ineq)) = inequalities_valid - t_ineq;
      
      Rcpp::Rcout << "This are the erroneous inequality constraints (last 3 columns: value / target / delta):\n" 
                  << err_log.rows(inequalities_invalid) << "\n";
      
      // Rcpp::Rcout << "-------------------------------------\n" << ineq_viol_output << "\n";
      
      Rcpp::stop(invalid_values);
    }
    
  }
  
  
  // Index vectors
  const arma::uvec Q_active = arma::linspace<arma::uvec>(0,dim_Q-1,dim_Q);
  arma::uvec eq_active(0,0,arma::fill::zeros);
  
  
  if(dim_eq>0) 
    eq_active = arma::linspace<arma::uvec>(0,dim_eq-1,dim_eq) + dim_Q;
  
  const arma::uvec inequalities= arma::linspace<arma::uvec>(0,dim_ineq-1,dim_ineq) + dim_Q + dim_eq;
  
  // Fixed left-hand side (needs subsetting for active constraints)
  arma::mat LHS(dim_Q + dim_eq + dim_ineq,
                dim_Q + dim_eq + dim_ineq,
                arma::fill::zeros);
  
  
  
  LHS(0,0,arma::size(dim_Q,dim_Q)) = Q;
  
  
  if(dim_eq>0)
  {
    LHS(0,dim_Q,arma::size(dim_Q,dim_eq)) = C_eq.t();
    LHS(dim_Q,0,arma::size(dim_eq,dim_Q)) = C_eq;
  }
  
  if(dim_ineq>0)
  {
    LHS(0,dim_Q+dim_eq,arma::size(dim_Q,dim_ineq)) = C_ineq.t();
    LHS(dim_Q+dim_eq,0,arma::size(dim_ineq,dim_Q)) = C_ineq;
  }
  
  
  
  arma::vec lagrange_eq(dim_eq,arma::fill::zeros);
  arma::vec lagrange_ineq(dim_ineq,arma::fill::zeros);
  
  unsigned iterations = 0;
  
  // auto     pre = std::chrono::high_resolution_clock::now();
  // auto     post = std::chrono::high_resolution_clock::now();
  // double solving = 0;
  
  do
  { 
    iterations++;
    // Rcpp::Rcout << "--------------------------------------\nIteration " << iterations << "\n";
    
    // -------------------------------------------------------
    sqp::misc::section("Check KKT conditions", debug);
    // -------------------------------------------------------
    
    const arma::vec eq_values = C_eq * x;
    
    const arma::vec ineq_values = C_ineq * x;
    
    // Derivative of Lagrangean = 0
    const bool d_L = arma::all(arma::abs(Q.t() * x + l + (C_eq.t() * lagrange_eq) + (C_ineq.t() * lagrange_ineq)) <= tol);
    // Primal feasibility
    const bool primal_feasibility = arma::all(eq_values == t_eq) && arma::all(ineq_values <= t_ineq);
    // Dual feasibility
    const bool dual_feasibility = arma::all(lagrange_ineq >= 0);
    
    // Rcpp::Rcout << "lagrange_eq: " << lagrange_eq.n_elem << "\n";
    
    
    // Slackness
    const bool slackness = arma::all(arma::abs(lagrange_ineq%(ineq_values - t_ineq))<= tol);
    
    
    // arma::mat INEQ_TEST(lagrange_ineq.n_elem,5,arma::fill::zeros);
    // INEQ_TEST.col(0) = t_ineq;
    // INEQ_TEST.col(1) = ineq_values;
    // INEQ_TEST.col(2) = ineq_values - t_ineq;
    // INEQ_TEST.col(3) = lagrange_ineq;
    // INEQ_TEST.col(4) = lagrange_ineq%INEQ_TEST.col(2);
    // // 
    // Rcpp::Rcout << "d_L: " << d_L<< "\n";
    // Rcpp::Rcout << "primal_feasibility: " << primal_feasibility<< "\n";
    // Rcpp::Rcout << "dual_feasibility: " << dual_feasibility<< "\n";
    // Rcpp::Rcout << "slackness: " << slackness<< "\n";
    // // 
    // Rcpp::Rcout << "INEQ:\n" << INEQ_TEST << "\n";
    
    if(d_L && primal_feasibility && dual_feasibility && slackness)
    {
      sqp::misc::section("sqp KKT conditions reached", debug);
      
      break;
    }
    
    // const arma::uvec u_ineq = arma::unique(ineq_active);
    // 
    // Rcpp::Rcout << "ineq_active: " << u_ineq.n_elem << " / " << ineq_active.n_elem << "\n" << u_ineq << "\n";
    
    
    const arma::uvec active_entries = arma::join_cols(Q_active,
                                                      arma::join_cols(eq_active,
                                                                      inequalities.elem(ineq_active)));
    
    // -------------------------------------------------------
    sqp::misc::section("Solve sqp with active constraints", debug);
    // -------------------------------------------------------
    
    // Fixed right-hand side
    
    arma::vec rhs = arma::join_cols(-Q.t()*x-l,arma::zeros(dim_eq + ineq_active.n_elem));
    arma::mat lhs_set = LHS.submat(active_entries,active_entries);
    
    arma::vec optim;
    
    arma::uvec activating_constraint, deactivating_constraint;
    
    const arma::vec lhs_min(arma::min(lhs_set,1));
    const arma::vec lhs_max(arma::max(lhs_set,1));
    
    arma::uvec lhs_check = arma::sort(arma::find(lhs_min == 0  && lhs_max == 0));
    // 
    if(lhs_check.n_elem > 0)
    {
      
      rhs.rows(lhs_check).zeros();
      
      for(int i=lhs_check.n_rows-1;i>=0;--i)
      {
        const unsigned pos = arma::as_scalar(lhs_check.row(i));
        
        lhs_set(pos,pos,arma::size(1,1)).ones();
      }
      
      std::string warn = "Got matrix with rows of zeros for quadratic problem with inequalities, which is not solvable!\n";
      Rcpp::Rcout << warn << "\n";
      Rcpp::warning(warn);
    }
    
    
    try{
      if(solver == 0)
      {
        if(fast)
        {
          optim = arma::solve(lhs_set,
                              rhs, 
                              arma::solve_opts::no_approx  +
                                arma::solve_opts::fast 
                                + arma::solve_opts::likely_sympd
          );
          
        } else
        {
          
          optim = arma::solve(lhs_set,
                              rhs, 
                              arma::solve_opts::equilibrate + 
                                arma::solve_opts::likely_sympd
                                + arma::solve_opts::refine
          );
          
        }
      } else if(solver==1)
      {
        arma::mat LU_L,LU_U;
        arma::lu(LU_L,LU_U,lhs_set);
        
        if(fast)
        {
          arma::vec lower_result = arma::solve(LU_L,
                                               rhs,
                                               arma::solve_opts::no_approx +
                                                 arma::solve_opts::fast);
          
          optim = arma::solve(arma::trimatu(LU_U),
                              lower_result,
                              arma::solve_opts::no_approx +
                                arma::solve_opts::fast);
        } else
        {
          arma::vec lower_result = arma::solve(LU_L,
                                               rhs,
                                               arma::solve_opts::equilibrate);
          
          optim = arma::solve(arma::trimatu(LU_U),
                              lower_result);
        }
        
      } else if((solver>1) && (solver<6))
      {
        
#if SQP_USE_VIENNACL == 1
        
        arma::Col<double> rhsx = rhs;
        viennacl::vector<double> rhs_vcl(rhsx.n_rows);
        viennacl::copy(rhsx,rhs_vcl);
        
        arma::SpMat<double> lhs_1x;
        arma::Mat<double> lhs_2x;
        
        if(solver%2==0)
        {
          lhs_1x = arma::SpMat<double>(lhs_set);
        } else
        {
          arma::Mat<double> lhs_1;
          arma::lu(lhs_1,lhs_2x,lhs_set);
          lhs_1x = arma::SpMat<double>(lhs_1);
        }
        
        viennacl::compressed_matrix<double> lhs1_vcl(lhs_1x.n_rows,lhs_1x.n_cols);
        viennacl::matrix<double> lhs2_vcl(lhs_2x.n_rows,lhs_2x.n_cols);
        viennacl::vector<double> optim_vcl(lhs_2x.n_rows);
        
        viennacl::copy(lhs_1x,lhs1_vcl);
        if(lhs_2x.n_rows>0)
          viennacl::copy(lhs_2x,lhs2_vcl);
        
        
        if(solver==2)
        {
          // viennacl::linalg::ilu0_tag ilu0_config;
          // viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<double>> vcl_precond(lhs1_vcl, ilu0_config);
          //
          viennacl::linalg::row_scaling< viennacl::compressed_matrix<double> > vcl_precond(lhs1_vcl, viennacl::linalg::row_scaling_tag());
          
          
          optim_vcl = viennacl::linalg::solve(lhs1_vcl,
                                              rhs_vcl,
                                              viennacl::linalg::cg_tag(),
                                              vcl_precond);
          
        } else if(solver==3)
        {
          // viennacl::linalg::ilu0_tag ilu0_config;
          // viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<double>> vcl_precond(lhs1_vcl, ilu0_config);
          //
          viennacl::linalg::row_scaling< viennacl::compressed_matrix<double> > vcl_precond(lhs1_vcl, viennacl::linalg::row_scaling_tag());
          
          viennacl::vector<double>  lower_result = viennacl::linalg::solve(lhs1_vcl,
                                                                           rhs_vcl,
                                                                           viennacl::linalg::cg_tag(),
                                                                           vcl_precond);
          
          optim_vcl = viennacl::linalg::solve(lhs2_vcl,
                                              lower_result,
                                              viennacl::linalg::cg_tag());
        } else if(solver==4)
        {
          // viennacl::linalg::ilu0_tag ilu0_config;
          // viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<double>> vcl_precond(lhs1_vcl, ilu0_config);
          //
          viennacl::linalg::row_scaling< viennacl::compressed_matrix<double> > vcl_precond(lhs1_vcl, viennacl::linalg::row_scaling_tag());
          
          
          optim_vcl = viennacl::linalg::solve(lhs1_vcl,
                                              rhs_vcl,
                                              viennacl::linalg::gmres_tag(),
                                              vcl_precond);
          
        } else if(solver==5)
        {
          // viennacl::linalg::ilu0_tag ilu0_config;
          // viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<double>> vcl_precond(lhs1_vcl, ilu0_config);
          //
          viennacl::linalg::row_scaling< viennacl::compressed_matrix<double> > vcl_precond(lhs1_vcl, viennacl::linalg::row_scaling_tag());
          
          
          viennacl::vector<double>  lower_result = viennacl::linalg::solve(lhs1_vcl,
                                                                           rhs_vcl,
                                                                           viennacl::linalg::gmres_tag(),
                                                                           vcl_precond);
          
          optim_vcl = viennacl::linalg::solve(lhs2_vcl,
                                              lower_result,
                                              viennacl::linalg::gmres_tag(),
                                              vcl_precond);
        } 
        
        optim.set_size(lhs_1x.n_rows);
        viennacl::copy(optim_vcl,optim);
        
#endif
        
      } else if ((solver > 5) && (solver < 20))
      {
        Eigen::MatrixXd rhs_eig = sqp::misc::convert_eigen(rhs);
        
        Eigen::MatrixXd lhs1_eig;
        Eigen::MatrixXd lhs2_eig;
        Eigen::MatrixXd optim_eig;
        
        if(solver%2==0)
        {
          lhs1_eig = sqp::misc::convert_eigen(lhs_set);
        } else
        {
          arma::mat lhs_1,lhs_2;
          arma::lu(lhs_1,lhs_2,lhs_set);
          lhs1_eig = sqp::misc::convert_eigen(lhs_1);
          lhs2_eig = sqp::misc::convert_eigen(lhs_2);
        }
        
        
        if(solver==6)
        {
          optim_eig = lhs1_eig.householderQr().solve(rhs_eig);
          
        } else if(solver==7)
        {
          Eigen::MatrixXd lower_result = lhs1_eig.householderQr().solve(rhs_eig);
          
          optim_eig = lhs2_eig.householderQr().solve(lower_result);
        } else if(solver==8)
        {
          optim_eig = lhs1_eig.lu().solve(rhs_eig);
          
        } else if(solver==9)
        {
          Eigen::MatrixXd lower_result = lhs1_eig.lu().solve(rhs_eig);
          
          optim_eig = lhs2_eig.lu().solve(lower_result);
        } else if(solver==10)
        {
          
          Eigen::SparseMatrix<double> lhs1_eig_sp = lhs1_eig.sparseView();
          Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
          
          solver.compute(lhs1_eig_sp);
          
          optim_eig = solver.solve(rhs_eig);
          
        } else if(solver==11)
        {
          
          Eigen::SparseMatrix<double> lhs1_eig_sp = lhs1_eig.sparseView();
          Eigen::SparseMatrix<double> lhs2_eig_sp = lhs2_eig.sparseView();
          
          Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
          
          solver.compute(lhs1_eig_sp);
          
          Eigen::MatrixXd lower_result = solver.solve(rhs_eig);
          
          solver.compute(lhs2_eig_sp);
          
          optim_eig = solver.solve(lower_result);
        } else if(solver==12)
        {
          
          Eigen::SparseMatrix<double> lhs1_eig_sp = lhs1_eig.sparseView();
          Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::COLAMDOrdering<int>> solver;
          
          solver.compute(lhs1_eig_sp);
          
          optim_eig = solver.solve(rhs_eig);
          
        } else if(solver==13)
        {
          
          Eigen::SparseMatrix<double> lhs1_eig_sp = lhs1_eig.sparseView();
          Eigen::SparseMatrix<double> lhs2_eig_sp = lhs2_eig.sparseView();
          
          Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::COLAMDOrdering<int> > solver;
          
          solver.compute(lhs1_eig_sp);
          
          Eigen::MatrixXd lower_result = solver.solve(rhs_eig);
          
          solver.compute(lhs2_eig_sp);
          
          optim_eig = solver.solve(lower_result);
        } else if(solver==14)
        {
          
          Eigen::SparseMatrix<double> lhs1_eig_sp = lhs1_eig.sparseView();
          Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >  solver;
          
          solver.compute(lhs1_eig_sp);
          
          optim_eig = solver.solve(rhs_eig);
          
        } else if(solver==15)
        {
          
          Eigen::SparseMatrix<double> lhs1_eig_sp = lhs1_eig.sparseView();
          Eigen::SparseMatrix<double> lhs2_eig_sp = lhs2_eig.sparseView();
          
          Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
          
          solver.compute(lhs1_eig_sp);
          
          Eigen::MatrixXd lower_result = solver.solve(rhs_eig);
          
          solver.compute(lhs2_eig_sp);
          
          optim_eig = solver.solve(lower_result);
        }
        
        optim = sqp::misc::convert_eigen(optim_eig);
      }
    } catch (std::exception& ex)
    {
      Rcpp::Rcout << "Exception in dense sqp with inequalities:\n" << ex.what() << ",\nattempting solution using inverse.\n";
      try{
        optim = arma::inv(LHS.submat(active_entries,active_entries))*rhs;
      } catch (std::exception& ex2)
      {
        Rcpp::stop(ex2.what());
      }
    }
    
    
    
    //post = std::chrono::high_resolution_clock::now();
    
    // solving += std::chrono::duration_cast<std::chrono::milliseconds>(post-pre).count();
    
    lagrange_eq = optim.elem(eq_active);
    
    lagrange_ineq.zeros();
    if(ineq_active.n_elem>0)
    {
      lagrange_ineq.elem(ineq_active) = optim.subvec(dim_Q+dim_eq,arma::size(ineq_active));
    }
    
    arma::vec delta_x = optim.rows(Q_active);
    
    const bool delta_x_zero = arma::all(arma::abs(delta_x) <= tol);
    
    if(!delta_x.is_finite())
    {
      delta_x.replace(arma::datum::nan,0);
    }
    
    
    if(delta_x_zero)
    {
      
      bool mult_ge_zero = true;
      double min_L_ineq = 0;
      
      if(ineq_active.n_elem>0)
      {
        min_L_ineq = arma::min(lagrange_ineq);
        mult_ge_zero = min_L_ineq >= 0;
      }
      
      if(mult_ge_zero)
      {
        // -------------------------------------------------------
        sqp::misc::section("Case a)", debug);
        // -------------------------------------------------------
        // delta_x = 0 and all lagrange multipliers of inequalities >= 0:
        // current value is optimal
        
        break;
      } 
      
      if(ineq_active.n_elem>0)
      {
        if(!mult_ge_zero)
        {
          // -------------------------------------------------------
          sqp::misc::section("Case b)", debug);
          // -------------------------------------------------------
          // delta_x = 0 and some lagrange multipliers of inequalities < 0 
          // reduce active set
          
          GetRNGstate();
          const int rm_active_ineq = arma::as_scalar((arma::shuffle(arma::find(lagrange_ineq.elem(ineq_active) == min_L_ineq))).eval().row(0));
          PutRNGstate();
          
          activating_constraint = ineq_active.row(rm_active_ineq);
          
          ineq_inactive.insert_rows(ineq_inactive.n_rows,activating_constraint);
          ineq_active.shed_row(rm_active_ineq);
        }
      }
    } else
    {
      const arma::vec x_new = x + delta_x;
      
      const bool equality_check = arma::approx_equal(C_eq * x_new,t_eq,"both",tol,tol);
      const bool inequality_check = arma::all((C_ineq * x_new - t_ineq)<tol);
      
      const arma::vec checkx = (C_ineq * x_new - t_ineq);
      
      // Rcpp::Rcout << "check vec: " << checkx.t() << "\n";
      
      if(equality_check & inequality_check)
      {
        // -------------------------------------------------------
        sqp::misc::section("Case c)", debug);
        // -------------------------------------------------------
        // delta_x != 0 and x + delta_x feasible
        // change x
        
        x = x_new;
      } else
      {
        // -------------------------------------------------------
        sqp::misc::section("Case d)", debug);
        // -------------------------------------------------------
        // delta_x != 0 and x + delta_x not feasible
        // change x as far as possible without becoming unfeasible
        // and activate corresponding inequality constraint
        const arma::vec numerator = t_ineq.rows(ineq_inactive) - C_ineq.rows(ineq_inactive) * x;
        const arma::vec denominator = C_ineq.rows(ineq_inactive) * delta_x;
        
        const arma::uvec valid_pos = arma::find(denominator>0);
        
        // Rcpp::Rcout << "Violated lhs:\n" << C_ineq.rows(ineq_inactive.rows(valid_pos));
        // Rcpp::Rcout << "Violated Rhs:\n" << t_ineq.rows(ineq_inactive.rows(valid_pos));
        
        const arma::vec ratio = numerator.rows(valid_pos)/denominator.rows(valid_pos);
        double stepwidth;
        try{
          stepwidth = arma::as_scalar(arma::min(ratio));
        } catch (std::exception& e)
        {
          Rcpp::stop(e.what());
        }
        GetRNGstate();
        const int rm_inactive_ineq = arma::as_scalar((arma::shuffle(valid_pos.rows(arma::find(ratio == stepwidth)))).eval().row(0));
        PutRNGstate();
        
        deactivating_constraint = ineq_inactive.row(rm_inactive_ineq);
        
        ineq_active.insert_rows(ineq_active.n_rows,deactivating_constraint);
        
        ineq_inactive.shed_row(rm_inactive_ineq);
        
        const arma::vec step = stepwidth * delta_x;
        
        x = x + step;
        
        if(arma::all(arma::abs(step) < tol))
          break;
      }
    }
  } while (iterations < max_iter);
  
  sqp::solvers::result output;
  
  output.x = x;
  output.lagrange_eq = lagrange_eq;
  output.lagrange_ineq = lagrange_ineq;
  
  return(output);
}

inline sqp::solvers::result inequalities(
    arma::vec                 x,        // Optimization variable (initial value)
    const arma::sp_mat        &Q,       // Quadratic Multiplier
    const arma::sp_mat        &C_eq,    // Equality  Constraint Multiplier
    arma::sp_mat              C_ineq,   // Inquality Constraint Multiplier
    arma::vec                 l,        // Linear Multiplier
    const arma::vec           &t_eq,    // Equality   Constraint RHS
    const arma::vec           &t_ineq,  // Inequality Constraint RHS (upper bound)
    const double              tol       = 1e-7,
    const unsigned            max_iter  = 500,
    int                       dim_eq    = -1,
    int                       dim_ineq  = -1,
    int                       dim_Q     = -1,
    const bool                debug     = false)
{
  // -------------------------------------------------------
  sqp::misc::section("Inequality constrained sqp (sparse)", debug);
  // -------------------------------------------------------  
  if(dim_eq < 0)
    dim_eq = C_eq.n_rows;
  
  if(dim_ineq < 0)
    dim_ineq = C_ineq.n_rows;
  
  if(dim_Q < 0)
    dim_Q = Q.n_rows;
  
  
  // arma::mat C_ineq_dense(C_ineq);
  
  // Index vectors
  arma::uvec Q_active = arma::linspace<arma::uvec>(0,dim_Q-1,dim_Q);
  arma::uvec eq_active(0,arma::fill::zeros);
  
  if(dim_eq>0)
    eq_active = arma::linspace<arma::uvec>(0,dim_eq-1,dim_eq) + dim_Q;
  
  arma::uvec inequalities= arma::linspace<arma::uvec>(0,dim_ineq-1,dim_ineq) + dim_Q + dim_eq;
  
  // Fixed left-hand side (needs subsetting for active constraints)
  arma::sp_mat LHS(dim_Q + dim_eq + dim_ineq,
                   dim_Q + dim_eq + dim_ineq);
  
  LHS(0,0,arma::size(dim_Q,dim_Q)) = Q;
  if(dim_eq>0)
  {
    LHS(0,dim_Q,arma::size(dim_Q,dim_eq)) = C_eq.t();
    LHS(dim_Q,0,arma::size(dim_eq,dim_Q)) = C_eq;
  }
  
  
  
  unsigned n_active_ineq;
  {
    // -------------------------------------------------------
    sqp::misc::section("Check constraints", debug);
    // -------------------------------------------------------
    const arma::vec equalities_valid = C_eq * x;
    const arma::vec inequalities_valid = C_ineq * x;
    
    // arma::uvec ineq_active = arma::find(inequalities_valid == t_ineq);
    // arma::uvec ineq_inactive = arma::find(inequalities_valid != t_ineq);
    
    arma::uvec ineq_active = arma::find((inequalities_valid - t_ineq)>=-tol);
    arma::uvec ineq_inactive = arma::find((inequalities_valid - t_ineq)<-tol);
    
    // Rcpp::Rcout << "active inequalities:\n " << (inequalities_valid-t_ineq).eval().rows(ineq_active).t() << "\n";
    
    const arma::uvec equalities_invalid = arma::find(arma::abs(equalities_valid - t_eq)>tol);
    const arma::uvec inequalities_invalid = arma::find((inequalities_valid - t_ineq)>tol);
    
    if(equalities_invalid.n_elem)
    {
      std::string invalid_values = "The following equalities are violated by the initial x-Values:\n     ";
      
      for(unsigned k=0;k<equalities_invalid.n_elem;++k)
      {
        invalid_values += std::to_string(arma::as_scalar(equalities_invalid.row(k)));
        
        if(k < equalities_invalid.n_elem-1)
          invalid_values += ", ";
      }
      
      invalid_values += "\n";
      
      arma::mat err_log(C_eq.n_rows,C_eq.n_cols + 3,arma::fill::zeros);
      err_log(0,0,arma::size(C_eq)) = C_eq;
      err_log(0,C_eq.n_cols,arma::size(equalities_valid)) = equalities_valid;
      err_log(0,C_eq.n_cols+1,arma::size(t_eq)) = t_eq;
      err_log(0,C_eq.n_cols+2,arma::size(t_eq)) = equalities_valid - t_eq;
      
      Rcpp::Rcout << "This are the erroneous equality constraints (last 3 columns: value / target / delta):\n" 
                  << err_log.rows(equalities_invalid) << "\n";
      
      Rcpp::stop(invalid_values);
    }
    
    
    if(inequalities_invalid.n_elem)
    {
      
      std::string invalid_values = "The following inequalities are violated by the initial x-Values:\n     ";
      
      for(unsigned k=0;k<inequalities_invalid.n_elem;++k)
      {
        invalid_values += std::to_string(arma::as_scalar(inequalities_invalid.row(k)));
        
        if(k < inequalities_invalid.n_elem-1)
          invalid_values += ", ";
      }
      
      invalid_values += "\n";
      arma::mat err_log(C_ineq.n_rows,C_ineq.n_cols + 2,arma::fill::zeros);
      err_log(0,0,arma::size(C_ineq)) = C_ineq;
      err_log(0,C_ineq.n_cols,arma::size(inequalities_valid)) = inequalities_valid;
      err_log(0,C_ineq.n_cols+1,arma::size(t_ineq)) = t_ineq;
      err_log(0,C_ineq.n_cols+2,arma::size(t_ineq)) = inequalities_valid - t_ineq;
      
      Rcpp::Rcout << "This are the erroneous inequality constraints (last 3 columns: value / target / delta):\n" 
                  << err_log.rows(inequalities_invalid) << "\n";
      
      Rcpp::stop(invalid_values);
    }
    
    // Reduce number of swaps by swapping the lesser of active / inactive constraints
    // in transposed mode (because of internal colwise storage)
    
    arma::sp_mat C_ineq_sorter = C_ineq.t();
    
    if(ineq_inactive.n_elem < ineq_active.n_elem)
    {
      ineq_inactive = arma::sort(ineq_inactive);
      
      for(unsigned i=0;i<ineq_inactive.n_elem;++i)
      {
        const unsigned pos = arma::as_scalar(ineq_inactive(i));
        
        C_ineq_sorter.swap_cols(pos,ineq_active.n_elem+i);
        inequalities.swap_rows(pos,ineq_active.n_elem+i);
      }
    } else
    {
      ineq_active = arma::sort(ineq_active);
      
      for(unsigned i=0;i<ineq_active.n_elem;++i)
      {
        const unsigned pos = arma::as_scalar(ineq_active(i));
        
        C_ineq_sorter.swap_cols(pos,i);
        inequalities.swap_rows(pos,i);
      }
    }
    
    if(dim_ineq>0)
    {
      LHS(0,dim_Q+dim_eq,arma::size(dim_Q,dim_ineq)) = C_ineq_sorter;
      LHS(dim_Q+dim_eq,0,arma::size(dim_ineq,dim_Q)) = C_ineq_sorter.t();
    }
    n_active_ineq = ineq_active.n_elem;
  }
  
  arma::vec lagrange_eq(dim_eq,arma::fill::zeros);
  arma::vec lagrange_ineq(dim_ineq,arma::fill::zeros);
  
  arma::superlu_opts opts;
  
  opts.allow_ugly  = false;
  opts.equilibrate = true;
  opts.symmetric = true;
  
  opts.permutation = arma::superlu_opts::MMD_AT_PLUS_A;
  
  // try({
  // opts.permutation = arma::superlu_opts::COLAMD;
  // }) // not available for older versions of Armadillo
  
  opts.refine = arma::superlu_opts::REF_SINGLE;
  opts.pivot_thresh = 0.001;
  
  unsigned iterations = 0;
  do
  {
    iterations++;
    
    // -------------------------------------------------------
    sqp::misc::section("Check KKT conditions", debug);
    // -------------------------------------------------------
    const arma::vec eq_values = C_eq * x;
    const arma::vec ineq_values = C_ineq * x;
    
    // Derivative of Lagrangean = 0
    const bool d_L = arma::all(arma::abs(Q.t() * x + l + (C_eq.t() * lagrange_eq) + (C_ineq.t() * lagrange_ineq)) <= tol);
    // Primal feasibility
    const bool primal_feasibility = arma::all(eq_values == t_eq) && arma::all(ineq_values <= t_ineq);
    // Dual feasibility
    const bool dual_feasibility = arma::all(lagrange_ineq >= 0);
    // Slackness
    const bool slackness = arma::all(arma::abs(lagrange_ineq%(ineq_values - t_ineq))<= tol);
    
    
    if(d_L && primal_feasibility && dual_feasibility && slackness)
    {
      sqp::misc::section("sqp KKT conditions reached", debug);
      
      break;
    }
    
    // -------------------------------------------------------
    sqp::misc::section("Solve sqp with active constraints", debug);
    // -------------------------------------------------------
    
    // Fixed right-hand side
    arma::vec rhs = arma::join_cols(-Q.t()*x-l,arma::zeros(dim_eq + n_active_ineq));
    
    arma::sp_mat lhs_set = LHS(0,0,arma::size(dim_Q + dim_eq + n_active_ineq,
                                              dim_Q + dim_eq + n_active_ineq));
    
    
    arma::vec optim;
    try{
#if SQP_USE_SUPERLU == 0
      
      optim = arma::spsolve(lhs_set,rhs,"lapack");
      
#elif SQP_USE_SUPERLU == 1
      
      // Avoid singularity due to negligable constrains
      const arma::vec lhs_min(arma::min(lhs_set,1));
      const arma::vec lhs_max(arma::max(lhs_set,1));
      
      arma::uvec lhs_check = arma::sort(arma::find(lhs_min == 0  && lhs_max == 0));
      // 
      if(lhs_check.n_elem > 0)
      {
        rhs.rows(lhs_check).zeros();
        
        for(int i=lhs_check.n_rows-1;i>=0;--i)
        {
          const unsigned pos = arma::as_scalar(lhs_check.row(i));
          
          lhs_set(pos,pos,arma::size(1,1)).ones();
        }
        std::string warn = "Got matrix with rows of zeros for quadratic problem with inequalities, which is not solvable!\n";
        Rcpp::Rcout << warn << "\n";
        Rcpp::warning(warn);
        
      }
      optim = arma::spsolve(lhs_set,rhs,"superlu",opts);
      
#endif // SQP_USE_SUPERLU == 1
    } catch (std::exception& e)
    {
      Rcpp::Rcout << "Exception in Sparse sqp with inequalities:\n" << e.what() << "\n";
      Rcpp::stop(e.what());
    }
    
    lagrange_eq = optim.elem(eq_active);
    
    lagrange_ineq.zeros();
    if(n_active_ineq>0)
    {
      lagrange_ineq.subvec(0,n_active_ineq-1) = 
        optim.subvec(dim_Q+dim_eq,dim_Q+dim_eq+n_active_ineq-1);
    }
    
    arma::vec delta_x = optim.rows(Q_active);
    
    if(!delta_x.is_finite())
    {
      delta_x.replace(arma::datum::nan,0);
      // return(delta_x);
    }
    
    const bool delta_x_zero = arma::all(arma::abs(delta_x) <= tol);
    
    
    if(delta_x_zero)
    {
      bool mult_ge_zero = true;
      double min_L_ineq = 0;
      
      if(n_active_ineq>0)
      {
        min_L_ineq = arma::min(lagrange_ineq);
        mult_ge_zero = min_L_ineq >= 0;
      }
      
      if(mult_ge_zero)
      {
        // -------------------------------------------------------
        sqp::misc::section("Case a)", debug);
        // -------------------------------------------------------
        // delta_x = 0 and all lagrange multipliers of inequalities >= 0:
        // current value is optimal
        
        break;
      }
      
      if(n_active_ineq>0)
      {
        if(!mult_ge_zero)
        {
          // -------------------------------------------------------
          sqp::misc::section("Case b)", debug);
          // -------------------------------------------------------
          // delta_x = 0 and some lagrange multipliers of inequalities < 0
          // reduce active set
          
          GetRNGstate();
          const int rm_active_ineq = arma::as_scalar(arma::shuffle(
            arma::find(lagrange_ineq.subvec(0,n_active_ineq-1) == min_L_ineq)).eval().row(0));
          PutRNGstate();
          
          const unsigned switch_pos_1 = dim_Q+dim_eq+rm_active_ineq;
          const unsigned switch_pos_2 = dim_Q+dim_eq+n_active_ineq-1;
          
          for(int i=0;i<dim_Q;++i)
          {
            const double val_1 = LHS(switch_pos_1,i);
            const double val_2 = LHS(switch_pos_2,i);
            
            LHS(i,switch_pos_1) = val_2;
            LHS(switch_pos_1,i) = val_2;
            
            LHS(i,switch_pos_2) = val_1;
            LHS(switch_pos_2,i) = val_1;
          }
          
          inequalities.swap_rows(rm_active_ineq,n_active_ineq-1);
          lagrange_ineq.swap_rows(rm_active_ineq,n_active_ineq-1);
          n_active_ineq--;
          
          
        }
      }
    } else
    {
      const arma::vec x_new = x + delta_x;
      
      const bool equality_check = arma::approx_equal(C_eq * x_new,t_eq,"both",tol,tol);
      const bool inequality_check = arma::all((C_ineq * x_new - t_ineq)<tol);
      
      if(equality_check & inequality_check)
      {
        // time_c_pre = std::chrono::high_resolution_clock::now();
        
        // -------------------------------------------------------
        sqp::misc::section("Case c)", debug);
        // -------------------------------------------------------
        // delta_x != 0 and x + delta_x feasible
        // change x
        
        x = x_new;
        
        // time_c = std::chrono::high_resolution_clock::now();
        // case_c += std::chrono::duration_cast<std::chrono::milliseconds>(time_c-time_c_pre).count();
        
      } else
      {
        //time_d_pre = std::chrono::high_resolution_clock::now();
        
        // -------------------------------------------------------
        sqp::misc::section("Case d)", debug);
        // -------------------------------------------------------
        // delta_x != 0 and x + delta_x not feasible
        // change x as far as possible without becoming unfeasible
        // and activate corresponding inequality constraint
        
        const arma::uvec ineq_inactive = inequalities.subvec(n_active_ineq,
                                                             inequalities.n_elem-1) - dim_Q - dim_eq;
        
        const arma::vec numerator = t_ineq.rows(ineq_inactive) - 
          LHS(dim_Q+dim_eq+n_active_ineq,0,
              arma::size(dim_ineq-n_active_ineq,dim_Q)) * x;
        
        const arma::vec denominator = LHS(dim_Q+dim_eq+n_active_ineq,0,
                                          arma::size(dim_ineq-n_active_ineq,dim_Q)) * 
                                            delta_x;
        
        const arma::uvec valid_pos = arma::find(denominator>0);
        ;
        const arma::vec ratio = numerator.rows(valid_pos)/denominator.rows(valid_pos);
        
        const double stepwidth = arma::as_scalar(arma::min(ratio));
        
        GetRNGstate();
        const int add_active_ineq =  n_active_ineq + 
          arma::as_scalar(arma::shuffle(valid_pos.rows(arma::find(ratio == stepwidth))).eval().row(0));
        PutRNGstate();
        
        const unsigned switch_pos_1 = dim_Q+dim_eq+add_active_ineq;
        const unsigned switch_pos_2 = dim_Q+dim_eq+n_active_ineq;
        
        for(int i=0;i<dim_Q;++i)
        {
          const double val_1 = LHS(switch_pos_1,i);
          const double val_2 = LHS(switch_pos_2,i);
          
          LHS(i,switch_pos_1) = val_2;
          LHS(switch_pos_1,i) = val_2;
          
          LHS(i,switch_pos_2) = val_1;
          LHS(switch_pos_2,i) = val_1;
        }
        
        lagrange_ineq.swap_rows(add_active_ineq,n_active_ineq);
        
        inequalities.swap_rows(add_active_ineq,n_active_ineq);
        
        n_active_ineq++;
        
        x = x + stepwidth * delta_x;
      }
    }
  } while (iterations < max_iter);
  
  lagrange_ineq.rows(inequalities-dim_Q-dim_eq) = lagrange_ineq;
  
  sqp::solvers::result output;
  
  output.x = x;
  output.lagrange_eq = lagrange_eq;
  output.lagrange_ineq = lagrange_ineq;
  
  return(output);
}



}
}

#endif                                                        // end of include guard for 'solvers/quadratic_inequalities.h'  
