
#include "includes_and_types.hh"
#include "classes_maximum_entropy.hh"
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include <ifopt/test_vars_constr_cost.h>

void
max_ent_solver
    ( matrix_t n_ini
    , vector_t y_ini
    , vector_t beta_ini
    , vector_t target_log_vcounts
    , matrix_t f_log_counts
    , matrix_t S
    , matrix_t K
    , vector_t obj_rxn_idx
    );

vector_t
rxn_flux
    ( vector_t v_log_counts
    , vector_t f_log_counts
    , matrix_t S
    , matrix_t K
    , matrix_t E_regulation
    );

