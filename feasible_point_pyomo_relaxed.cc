
#include <iostream>
#include <vector>
#include <numeric>
#include <math.h>

#include "Eigen/Core"
#include "Eigen/Dense"

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE = size_t;

typedef double value_t;
typedef std::vector<value_t> value_list_t;
typedef std::vector<size_t> index_list_t;
typedef Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic> matrix_t;

matrix_t 
Eq_obj( matrix_t beta,
        matrix_t Sv_N,
        matrix_t SvS,
        matrix_t S_v_T,
        matrix_t S_f_T,
        matrix_t alpha,
        matrix_t K,
        matrix_t FxdM )
{
    alpha = alpha.reshaped<Eigen::RowMajor>().eval();
    matrix_t y_e = Sv_N * beta.reshaped<Eigen::RowMajor>();

    auto alpha_a = alpha.array();
    auto K_a = K.array();
    auto y_e_a = y_e.array();
    matrix_t y_hat_e = ( K_a.log() 
                       - std::log(2) 
                       - alpha_a.log() 
                       + (2 * alpha_a).log() 
                       ).matrix();

    auto the_rest = ( std::pow(10,10) * y_e_a).matrix() 
                    * Eigen::inverse(y_e_a.abs() * std::pow(10,10) + std::pow(10, -12)).matrix()
                    * ((2*alpha_a).log() - Eigen::log(y_e_a.abs() + Eigen::sqrt( y_e_a.pow(2) + 4 * alpha_a.pow(2)))).matrix();

    y_hat_e = (y_hat_e + the_rest).eval();

    matrix_t b_e = y_hat_e - (S_f_T * FxdM.reshaped<Eigen::RowMajor>());

    matrix_t ret = ((S_v_T * (SvS * b_e) ) - b_e).reshaped<Eigen::RowMajor>();

    return ret;
}

matrix_t
rxn_flux( matrix_t v_log_counts,
          matrix_t f_log_counts,
          matrix_t S,
          matrix_t K,
          matrix_t E_regulation)
{
    matrix_t S_T = S;
    S.transposeInPlace();

    size_t n_react = S.cols();

    v_log_counts = v_log_counts.reshaped<Eigen::RowMajor>().eval();
    f_log_counts = f_log_counts.reshaped<Eigen::RowMajor>().eval();

    



}