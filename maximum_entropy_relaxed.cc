
#include <iostream>
#include <vector>
#include <numeric>

#include "Eigen/Dense"

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE = size_t;

typedef double value_t;
typedef long integer_t;
typedef std::vector<value_t> value_list_t;
typedef std::vector<size_t> index_list_t;
typedef Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic> matrix_t;


matrix_t remove_indices(
    const matrix_t& M, 
    const index_list_t indxs,
    const size_t axis)
{
}

Max_ent_solver_return_type /* TODO: specify the type */
max_ent_solver( 
    value_list_t n_ini,  /* TODO: figure out types */
    value_t y_ini,
    value_t beta_ini,
    size_t target_log_vcounts,
    value_list_t f_log_counts,
    matrix_t S,
    matrix_t K,
    size_t obj_rxn_idx
)
{
    /* Flip Stoichiometric Matrix */
    matrix_t S_T = S; /* this should do deep copy */
    S.transposeInPlace();
    size_t n_react = S.cols(); /* should be same as numpy.shape(S)[1] */

    /* Set the System Parameters */
    /* TODO? FxdM = np.reshape(f_log_counts,(len(f_log_counts),1) ) */
    matrix_t K_eq = K;
    K_eq.resize(n_react, 1); /* need to look into behavior of resize vs np.reshape */

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

    /* Split S into the component corresponding to the variable metabolites S_v
    and the component corresponding to the fixed metabolites */
    matrix_t S_v_T = remove_indices(S_T, FxdM_idx, 1);
    /* this should do the same thing as numpy.delete
            might not need to do anything,
            see page on indexing 
                https://eigen.tuxfamily.org/dox-devel/group__TutorialSlicingIndexing.html */






}
  

