
#include <iostream>
#include <vector>
#include <numeric>
#include <math.h>

#include "Eigen/Core"
#include "Eigen/Dense"

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE = size_t;

typedef long double value_t;
typedef std::vector<value_t> value_list_t;
typedef std::vector<size_t> index_list_t;
typedef Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic> matrix_t;

matrix_t 
Eq_obj
    ( matrix_t beta
    , matrix_t Sv_N
    , matrix_t SvS
    , matrix_t S_v_T
    , matrix_t S_f_T
    , matrix_t alpha
    , matrix_t K
    , matrix_t FxdM 
    )
{
    alpha = alpha.reshaped<Eigen::RowMajor>().eval();
    matrix_t y_e = Sv_N * beta.reshaped<Eigen::RowMajor>();

    auto alpha_a = alpha.array();
    auto K_a = K.array();
    auto y_e_a = y_e.array();
    auto first_part =  ( K_a.log() 
                       - std::log(2) 
                       - alpha_a.log() 
                       + (2 * alpha_a).log() 
                       ).matrix();
    auto second_part = 
        (   ( std::pow(10,10) * y_e_a)
            * Eigen::inverse(y_e_a.abs() * std::pow(10,10) + std::pow(10, -12))
            * ((2*alpha_a).log() - Eigen::log(y_e_a.abs() + Eigen::sqrt( y_e_a.pow(2) + 4 * alpha_a.pow(2))))
        ).matrix();

    matrix_t y_hat_e = first_part + second_part;

    matrix_t b_e = y_hat_e - (S_f_T * FxdM.reshaped<Eigen::RowMajor>());

    matrix_t ret = ((S_v_T * (SvS * b_e) ) - b_e).reshaped<Eigen::RowMajor>();

    return ret;
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
    /* flip Stoichiometric Matrix */
    matrix_t S_T = S; /* S_T is the Stoich matrix with rows as reactions, columns as metabolites */
    S.transposeInPlace(); /* this now is the Stoich matrix with rows metabolites, and columns reactions */

    size_t n_react = S.cols();

    v_log_counts = v_log_counts.reshaped<Eigen::RowMajor>().eval();
    f_log_counts = f_log_counts.reshaped<Eigen::RowMajor>().eval();

    matrix_t tot_log_counts(v_log_counts.size() + f_log_counts.size(), 1);
    tot_log_counts << v_log_counts, f_log_counts;

    K = K.reshaped<Eigen::RowMajor>().eval();

    E_regulation = E_regulation.reshaped<Eigen::RowMajor>().eval();
    
    auto multiplier =   ((-0.25) * S_T * tot_log_counts)
                        .array().exp().pow(4);
    auto forward_odds = K.array() * multiplier;

    auto ret    = E_regulation.array() 
                * (forward_odds - (forward_odds.pow(-1)));

    return ret.matrix();
}

feasible_point_solver_return_type
Feasible_point_solver
    ( matrix_t v_log_counts
    , value_list_t target_log_vcounts
    , matrix_t f_log_counts
    , matrix_t S
    , matrix_t K
    )
{
    /* flip Stoichiometric Matrix */
    matrix_t S_T = S; /* S_T is the Stoich matrix with rows as reactions, columns as metabolites */
    S.transposeInPlace(); /* this now is the Stoich matrix with rows metabolites, and columns reactions */
    size_t n_react = S.cols();

    /* Set the System Parameters */
    matrix_t VarM = v_log_counts;
    matrix_t FxdM = f_log_counts.reshaped<Eigen::RowMajor>();
    matrix_t K_eq = K.reshaped<Eigen::RowMajor>();

    /* Metabloite params */
    size_t n_M_f = f_log_counts.size();
    size_t n_M_v = v_log_counts.size();
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

    /* Find the gradient direction */
    matrix_t b(dSv_N, 1);
    b << Sv_N(Eigen::seq(181, 183), Eigen::all).colwise().sum();

    /*
    Steady State satisfies the following least squares problem
    || S_T*( SVS* y_hat) - y_hat ||^2 = 0
    where SVS is the matrix (S_v * S_v_T)^-1 * S_v
    */

    /* precompute SvS */
    matrix_t SvS = (S_v * S_v_T).inverse() * S_v;

    auto alpha = matrix_t::Constant(n_react, 1, 1.);

    /* Get the initial least squares solution without optimizing regulation */
    matrix_t beta = b;

    /* TODO : least squares (line 107) */
    matrix_t betlsq; /* place holder */

    /* Extract the solution */
    /* TODO: might need to modify depending on output of least squares */
    matrix_t y_sol = Sv_N * betlsq.reshaped<Eigen::RowMajor>();

    auto y_sol_a = y_sol.array();
    matrix_t h_sol  = 
        ( ( std::pow(10, 25) * y_sol_a ) 
        * ( std::pow(10, 25) * y_sol_a.abs() + std::pow(10, -25) ).inverse()
        * ( std::log(2) - Eigen::log(y_sol_a.abs() + Eigen::sqrt(y_sol_a.pow(2) + 4) ) )
        ).matrix();

    matrix_t b_RHS  = K_eq.array().log().matrix()
                    + h_sol
                    - (S_f_T * FxdM.reshaped<Eigen::RowMajor>());

    matrix_t n_sol = SvS * b_RHS;

    matrix_t b_sol = 
        ( ( S_v_T * n_sol.reshaped<Eigen::RowMajor>() )
        + ( S_f_T * FxdM.reshaped<Eigen::RowMajor>() )
        ).reshaped<Eigen::RowMajor>();

    matrix_t beta_sol = beta;


    /* get the system parameters to contruct a model */

    /* metabolite params */
    n_M = n_M_f + n_M_v; /* total number of metabolites */

    /* contruct param indices */
    /* react_idx contructed above (line 118 - 119) */

    matrix_t n_ini    = n_sol.reshaped<Eigen::RowMajor>();
    matrix_t b_ini    = b_sol.reshaped<Eigen::RowMajor>();
    matrix_t y_ini    = y_sol.reshaped<Eigen::RowMajor>();
    matrix_t beta_ini = beta_sol.reshaped<Eigen::RowMajor>();
    matrix_t h_ini    = h_sol.reshaped<Eigen::RowMajor>();

    /* Set optimization parameters */
    value_t VarM_bind = -300;

    /* TODO: optimizaton (lines 148 - 386) */

    return;
}





/*
Feasible Point Solver to generate sampling of initial conditions
*/

matrix_t
Eq_lobj
    ( matrix_t beta
    , matrix_t Sv_N
    , matrix_t SvS
    , matrix_t S_v_T
    , matrix_t S_f_T
    , matrix_t lalpha
    , matrix_t K
    , matrix_t FxdM
    )
{
    lalpha = lalpha.reshaped<Eigen::RowMajor>().eval();

    matrix_t y_e = Sv_N * beta.reshaped<Eigen::RowMajor>();

    auto K_a = K.array();
    auto y_e_a = y_e.array();
    auto lalpha_a = lalpha.array();
    matrix_t y_hat_e = 
        (   K_a.log()
        +   (   ( std::pow(10, 10) * y_e_a ) 
            *   ( std::pow(10,10) * y_e_a.abs() + std::pow(10, -12) )
                .inverse()
            *   ( std::log(2) 
                + lalpha_a 
                - Eigen::log( y_e_a.abs() 
                            + Eigen::sqrt   (  y_e_a.pow(2)
                                            + 4 * (2*lalpha_a).exp()
                                            )
                            )
                )
            )
        ).matrix();

    matrix_t b_e = y_hat_e - (S_f_T * FxdM.reshaped<Eigen::RowMajor>() );

    matrix_t ret = ((S_v_T * (SvS * b_e)) - b_e).reshaped<Eigen::RowMajor>();

    return ret;
}

Random_Feasible_Point_solver_return_type
Random_Feasible_Point_solver
    ( matrix_t v_log_counts
    , value_list_t target_log_vcounts
    , matrix_t f_log_counts
    , matrix_t S
    , matrix_t K
    )
{
    /* flip Stoichiometric Matrix */
    matrix_t S_T = S; /* S_T is the Stoich matrix with rows as reactions, columns as metabolites */
    S.transposeInPlace(); /* this now is the Stoich matrix with rows metabolites, and columns reactions */
    size_t n_react = S.cols();

    /* Set the System Parameters */
    matrix_t VarM = v_log_counts;
    matrix_t FxdM = f_log_counts.reshaped<Eigen::RowMajor>();
    matrix_t K_eq = K.reshaped<Eigen::RowMajor>();

    /* Metabloite params */
    size_t n_M_f = f_log_counts.size();
    size_t n_M_v = v_log_counts.size();
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
        https://eigen.tuxfamily.org/dox-devel/group__TutorialSlicingIndexing.html 
    */
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

    /* Get the initial least squares solution with optimizing requlation */
    auto beta = matrix_t::Random(dSv_N, 1); /* check: this should be uniform double [-1,1] */
    auto alpha = matrix_t::Constant(n_react, 1, 1.);

    /* TODO : least squares (line 564) */
    matrix_t betlsq; /* place holder */

    /* extract the solution */
    matrix_t y_sol = Sv_N * betlsq.reshaped<Eigen::RowMajor>();

    /* Set the initional condition */

    auto y_sol_a = y_sol.array();
    auto h_sol =
        (   y_sol_a * std::pow(10, 25) 
        *   Eigen::inverse(y_sol_a.abs() * std::pow(10, 25) + std::pow(10, -25))
        *   (std::log(2) - Eigen::log( y_sol_a.abs() + Eigen::sqrt(y_sol_a.pow(2) + 4))) 
        );

    auto K_eq_a = K_eq.array();
    matrix_t b_RHS =
        (   K_eq_a + h_sol 
        -   (S_f_T * FxdM.reshaped<Eigen::RowMajor>()).array()
        )   .matrix();

    matrix_t n_sol = SvS * b_RHS;

    matrix_t b_sol = 
        (   S_v_T * n_sol.reshaped<Eigen::RowMajor>()
        +   S_f_T * FxdM.reshaped<Eigen::RowMajor>()
        )   .reshaped<Eigen::RowMajor>();

    matrix_t beta_sol = beta;


    /* Get the system parameters to contruct the model */

    /* Metabolite params */
    n_M = n_M_f + n_M_v;

    /* construct param indices */
    index_list_t react_idx(n_react);
    std::iota(react_idx.begin(), react_idx.end(), 0);

    index_list_t beta_idx(dSv_N);
    std::iota(beta_idx.begin(), beta_idx.end(), 0);

    /* set the optimization parameters */
    value_t VarM_lbnd = -300l;

    /* TODO: model defintion (lines 607-845) */

    return;
}
