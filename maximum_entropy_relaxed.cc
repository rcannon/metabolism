
#include <iostream>
#include <vector>

#include "Eigen/Dense"

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE = unsigned

typedef double value_t;
typedef long integer_t;
typedef Eigen::Index index_t;
typedef std::vector<value_t> value_list_t;
typedef Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic> matrix_t;

Max_ent_solver_return_type /* TODO: specify the type */
max_ent_solver( 
    value_list_t n_ini,  /* TODO: figure out types */
    value_t y_ini,
    value_t beta_ini,
    index_t target_log_vcounts,
    value_list_t f_log_counts,
    matrix_t S,
    matrix_t K,
    index_t obj_rxn_idx
)
{
    /* Flip Stoichiometric Matrix */
    matrix_t S_T = S; /* this should do deep copy */
    S.transposeInPlace()
    index_t n_react = S.columns() /* should be same as numpy.shape(S)[1] */

    /* Set the System Parameters */
    /* TODO? FxdM = np.reshape(f_log_counts,(len(f_log_counts),1) ) */
    matrix_t K_eq = K;
    K_eq.resize(n_react, 1); /* need to look into behavior of resize vs np.reshape */

    /* Metabloite params */
    index_t n_M_f = f_log_counts.size();
    index_t n_M_v = n_ini.size();
    index_t n_M = n_M_f + n_M_v;



}
  





)