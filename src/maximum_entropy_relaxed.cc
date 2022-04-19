
#include "maximum_entropy_relaxed.hh"


/* TODO: might need to change def of value_t, so to ensure
safe numerical computation  */

Max_ent_solver_return_type /* TODO: specify the type */
max_ent_solver
    ( matrix_t n_ini
    , vector_t y_ini
    , vector_t beta_ini
    , vector_t target_log_vcounts
    , matrix_t f_log_counts
    , matrix_t S
    , matrix_t K
    , vector_t obj_rxn_idx
    )
{
    /* Flip Stoichiometric Matrix */
    matrix_t S_T = S; /* this should do deep copy */
    S.transposeInPlace();
    size_t n_react = S.cols(); /* should be same as numpy.shape(S)[1] */

    /* Set the System Parameters */
    matrix_t FxdM = f_log_counts.reshaped<Eigen::RowMajor>();
    matrix_t K_eq = K.reshaped<Eigen::RowMajor>();

    /* Metabloite params */
    size_t n_M_f = f_log_counts.size();
    size_t n_M_v = n_ini.size();
    size_t n_M = n_M_f + n_M_v;

    /* construct param indices */
    index_list_t react_idx(n_react);
    std::iota(react_idx.begin(), react_idx.end(), 0);

    index_list_t TotM_idx(n_M);
    std::iota(TotM_idx.begin(), TotM_idx.end(), 0);

    index_list_t VarM_idx(n_M_v); /* variable metabolite indices */
    std::iota(VarM_idx.begin(), VarM_idx.end(), 0);

    index_list_t FxdM_idx(n_M - n_M_v); /* fixed metabolite indices */
    std::iota(FxdM_idx.begin(), FxdM_idx.end(), n_M_v);

    /* 
    Split S into the component corresponding to the variable metabolites S_v
    and the component corresponding to the fixed metabolites.

    Reilly: this should do the same thing as numpy.delete
        https://eigen.tuxfamily.org/dox-devel/group__TutorialSlicingIndexing.html */
    matrix_t S_v_T = S_T(Eigen::all, FxdM_idx);
    matrix_t S_f_T = S_T(Eigen::all, VarM_idx);

    matrix_t S_v = S_v_T.transpose();
    matrix_t S_f = S_f_T.transpose();

    /* find a basis for the nullspace of S_v,
        and get the dimension of the null space */
    Eigen::FullPivLU<Eigen::MatrixXf> lu_decomp(S_v);
    matrix_t Sv_N = lu_decomp.kernel();
    size_t dSv_N = lu_decomp.dimensionOfKernel();

    /* precompute SVS */
    matrix_t SvS = (S_v * S_v_T).inverse() * S_v;

    /* Get the system parameters to construct the model */

    /* metabolite parameters */
    
    /* n_M = n_M_f + n_M_v was perfomred above (line 49) */

    /* construct param indices */

    /* react_idx = np.arange(0,n_react) was performed above (lines 52-53) */

    index_list_t beta_idx(dSv_N);
    std::iota(beta_idx.begin(), beta_idx.end(), 0);

    /* Set the initial condition */

    matrix_t b_ini  = (S_v_T * n_ini.reshaped<Eigen::RowMajor>()) 
                    + (S_f_T * FxdM.reshaped<Eigen::RowMajor>() );
    /* Reilly: need to call .eval to avoid memory issues due to aliasing */
    b_ini = b_ini.reshaped<Eigen::RowMajor>().eval();

    auto y_init_array = y_ini.array();

    vector_t h_ini = std::log(2) 
                   - Eigen::log( y_ini_array.abs()
                               + (y_ini_array.pow(2) + 4).sqrt()   
                               );
    if (std::signbit(y_ini)) { h_ini = - h_ini; }

    value_t VarM_lbnd = -300;
    

    /* TODO: IPOPT starting on line 86 */


}
  

matrix_t
rxn_flux
    ( matrix_t v_log_counts
    , matrix_t f_log_counts
    , matrix_t S
    , matrix_t K
    , matrix_t E_regulation
    )
{
    /* TODO */
}