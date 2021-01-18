#include <sundials/sundials_math.h>
#include <sunlinsol_rmumps.h>
// exported functions

// Function to create a new RMUMPS linear solver
SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_RMUMPS(N_Vector y, SUNMatrix A, int  permutation=RMUMPS_PERM_AUTO) {
//Rcout << "call SUNLinSol_RMUMPS\n";
//Rcout << "permutation=" << permutation << "\n";
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_RMUMPS content;
  
  // Check compatibility with supplied SUNMatrix and N_Vector
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE)
    return(NULL);
  if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Columns(A))
    return(NULL);
  if (N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL)
    return(NULL);

  // optimally this function would be replaced with a generic N_Vector routine
  int n = NV_LENGTH_S(y), nz=SM_NNZ_S(A);
  
  // Create linear solver
  S = NULL;
  S = (SUNLinearSolver) malloc(sizeof *S);
  if (S == NULL) return(NULL);
  
  // Create linear solver operation structure
  ops = NULL;
  ops = (SUNLinearSolver_Ops) malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  if (ops == NULL) { free(S); return(NULL); }

  // Attach operations
  ops->gettype           = SUNLinSolGetType_RMUMPS;
  ops->initialize        = SUNLinSolInitialize_RMUMPS;
  ops->setup             = SUNLinSolSetup_RMUMPS;
  ops->solve             = SUNLinSolSolve_RMUMPS;
  ops->lastflag          = NULL;
  ops->space             = NULL;
  ops->free              = SUNLinSolFree_RMUMPS;
  ops->setatimes         = NULL;
  ops->setpreconditioner = NULL;
  ops->setscalingvectors = NULL;
  ops->numiters          = NULL;
  ops->resnorm           = NULL;
  ops->resid             = NULL;

  // Create content
  content = NULL;
  content = (SUNLinearSolverContent_RMUMPS) malloc(sizeof(struct _SUNLinearSolverContent_RMUMPS));
  if (content == NULL) { free(ops); free(S); return(NULL); }

  // Fill content
  content->last_flag = 0;
  if (SUNSparseMatrix_SparseType(A) != CSC_MAT) {
    stop("SUNLinSol_RMUMPS: wrong sparse matrix type, expected CSC_MAT");
  }
  if (n != SM_COLUMNS_S(A))
    stop("SUNLinSol_RMUMPS: ncol(A) (%d) and length(y) (%d) don't concord", SM_COLUMNS_S(A), n);
  if (SM_COLUMNS_S(A) != SM_ROWS_S(A))
    stop("SUNLinSol_RMUMPS: matrix is supposed to be square, instead got %dx%d", SM_ROWS_S(A), SM_COLUMNS_S(A));
  // build jcp array from irp and pc
#if defined(SUNDIALS_INT32_T)
  ivec ir(SM_INDEXVALS_S(A), nz, false);
  ivec pc(SM_INDEXPTRS_S(A), n+1, false);
#else
  ivec ir=conv_to<ivec>::from(Col<sunindextype>(SM_INDEXVALS_S(A), nz, false));
  ivec pc=conv_to<ivec>::from(Col<sunindextype>(SM_INDEXPTRS_S(A), n+1, false));
#endif
//pc.print("pc");
//ir.print("ir");
//vec av(SM_DATA_S(S), nz, false);
//av.print("av");
//printf("s_r: irp=%x; jcp=%x; vp=%x\n", SM_INDEXPTRS_S(A), SM_INDEXVALS_S(A), SM_DATA_S(A));

  content->irp=new Col<MUMPS_INT>(ir+1);
  content->jcp=new Col<MUMPS_INT>(nz);
  
  //XPtr<Rmumps> p(rmumps::Rmumps__ptr_ijv(XPtr<int>(content->irp->begin(), false), XPtr<int>(content->jcp->begin(), false), XPtr<double>(SM_DATA_S(A), false), n, pc[n], 0), false); // 0=non symmetric matrix;
  content->rmu=new XPtr<Rmumps>(rmumps::Rmumps__ptr_ijv(XPtr<int>(content->irp->begin(), false), XPtr<int>(content->jcp->begin(), false), XPtr<double>(SM_DATA_S(A), false), n, pc[n], 0));
  rmumps::Rmumps__set_permutation(*(content->rmu), permutation);
/*
List asl=content->rmu->triplet();
print(wrap("init"));
print(asl["i"]);
print(asl["j"]);
print(asl["v"]);
print(asl["nrow"]);
*/
  // Attach content and ops
  S->content = content;
  S->ops     = ops;

  return(S);
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_RMUMPS(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}


SUNDIALS_EXPORT int SUNLinSolInitialize_RMUMPS(SUNLinearSolver S)
{
//Rcout << "call SUNLinSolInitialize_RMUMPS\n";
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


SUNDIALS_EXPORT int SUNLinSolSetup_RMUMPS(SUNLinearSolver S, SUNMatrix A) {
//Rcout << "call SUNLinSolSetup_RMUMPS\n";
  int n=SM_COLUMNS_S(A); //nzres=SM_NNZ_S(A); // nz reserved, it may be larger than real nz
  // Ensure that A is a sparse matrix
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) {
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return(LASTFLAG(S));
  }
  // update matrix data
#if defined(SUNDIALS_INT32_T)
  ivec pc(SM_INDEXPTRS_S(A), n+1, false);
  int nz=pc[n];
  ivec ir(SM_INDEXVALS_S(A), nz, false);
#else
  ivec pc=conv_to<ivec>::from(Col<sunindextype>(SM_INDEXPTRS_S(A), n+1, false));
  int nz=pc[n];
  ivec ir=conv_to<ivec>::from(Col<sunindextype>(SM_INDEXVALS_S(A), nz, false));
#endif
//ir.print("new ir");
//ir.print("new pc");
  delete RMUMPS_CONTENT(S)->irp;
  RMUMPS_CONTENT(S)->irp=new Col<MUMPS_INT>(ir+1);
  RMUMPS_CONTENT(S)->jcp->resize(nz);
  for (int j=1; j <= n; j++)
    RMUMPS_CONTENT(S)->jcp->subvec(pc[j-1], pc[j]-1).fill(j);
  // if matrix structure changed, reinit Rmumps object
  if (RMUMPS_CONTENT(S)->irp->n_elem != RMU(S).get()->irn.size() || !std::equal(RMUMPS_CONTENT(S)->irp->begin(), RMUMPS_CONTENT(S)->irp->end(), RMU(S).get()->irn.begin())) {
    // new Rmumps
//Rcout << "new Rmumps\n";
    int permutation=rmumps::Rmumps__get_permutation(RMU(S));
    rmumps::Rmumps__del_ptr(RMU(S));
    RMUMPS_CONTENT(S)->irp->resize(nz);
    RMUMPS_CONTENT(S)->irp->subvec(0, nz-1)=ir+1;
//RMUMPS_CONTENT(S)->irp->print("new irp");
//RMUMPS_CONTENT(S)->jcp->print("new jcp");
    RMU(S)=rmumps::Rmumps__ptr_ijv(XPtr<int>((MUMPS_INT *)RMUMPS_CONTENT(S)->irp->begin(), false), XPtr<int>((MUMPS_INT *)RMUMPS_CONTENT(S)->jcp->begin(), false), XPtr<double>(SM_DATA_S(A), false), n, nz, 0); // 0=non symmetric matrix;
    rmumps::Rmumps__set_permutation(RMU(S), permutation);
  } else {
//Rcout << "reset old Rmumps\n";
    rmumps::Rmumps__set_mat_ptr(RMU(S), XPtr<double>(SM_DATA_S(A), false));
  }
/*
List asl=rmumps::Rmumps__triplet(RMU(S));
print(wrap("reset"));
print(asl["i"]);
print(asl["j"]);
print(asl["v"]);
*/
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}

SUNDIALS_EXPORT int SUNLinSolSolve_RMUMPS(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                       N_Vector b, realtype tol) {
//Rcout << "call SUNLinSolSolve_RMUMPS\n";
//static long clk_tck = CLOCKS_PER_SEC;
//clock_t t1, t2;
  //static int ncall=0;
  //ncall++;
  int n=NV_LENGTH_S(x);
  sunindextype *ap=SM_INDEXPTRS_S(A);
  realtype *xdata=N_VGetArrayPointer(x), *bdata=N_VGetArrayPointer(b), *adata=SM_DATA_S(A);
  
  if (xdata == NULL) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(LASTFLAG(S));
  }
  // copy b into x
  if (bdata != xdata)
    std::copy(bdata, bdata+n, xdata);
  // check for identity matrix in A
  if (ap[n] == n && std::count(adata, adata+n, 1.) == n) {
    // identity matrix => we are done
    LASTFLAG(S) = SUNLS_SUCCESS;
    return(LASTFLAG(S));
  }
  

  // Call RMUMPS to solve the linear system
/*// print non identity matrix
List asl=rmumps::Rmumps__triplet(RMU(S));
print(wrap("solve"));
print(asl["i"]);
print(asl["j"]);
print(asl["v"]);
vec(xdata, NV_LENGTH_S(x), false).print("b");
*/
  rmumps::Rmumps__solveptr(RMU(S), XPtr<double>(xdata, false), n, 1);
  /*static mat ad=zeros(3,3);
  static vec xd(3);
  ad(0,0)=adata[0];
  ad(1,0)=adata[1];
  std::copy(adata+2, adata+8, ad.begin()+3);
//t1=clock();
  xd=solve(ad, vec(bdata, 3, false), solve_opts::fast);
//ad.print("ad");
//vec(adata, 8, false).print("adata");
//stop("ad");
//t2=clock();
  std::copy(xd.begin(), xd.end(), xdata);
//printf("%d t solve (s) : %lf \n", ncall, (double)(t2-t1)/(double)clk_tck);
//vec(xdata, NV_LENGTH_S(x), false).print("x");
  */
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}

SUNDIALS_EXPORT long int SUNLinSolLastFlag_RMUMPS(SUNLinearSolver S) {
  // return the stored 'last_flag' value
  return(LASTFLAG(S));
}

SUNDIALS_EXPORT int SUNLinSolFree_RMUMPS(SUNLinearSolver S) {
//Rcout << "call SUNLinSolFree_RMUMPS\n";
  // return with success if already freed
  if (S == NULL)
    return(SUNLS_SUCCESS);
  
  // delete items from the contents structure (if it exists)
  if (S->content) {
    delete RMUMPS_CONTENT(S)->irp;
    delete RMUMPS_CONTENT(S)->jcp;
    rmumps::Rmumps__del_ptr(RMU(S)); // Rmumps destructor
    delete RMUMPS_CONTENT(S)->rmu;
    free(S->content);
    S->content = NULL;
  }
  
  // delete generic structures
  if (S->ops) {
    free(S->ops);  
    S->ops = NULL;
  }
  
  free(S);
  S = NULL;
  return(SUNLS_SUCCESS);
}
