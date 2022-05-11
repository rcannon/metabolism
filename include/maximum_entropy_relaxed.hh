
#include "includes_and_types.hh"
#include "ifopt_variable_class.hh"
#include "ifopt_constraint_class.hh"
#include "ifopt_cost_class.hh"
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include <ifopt/test_vars_constr_cost.h>

#pragma once

std::tuple<vector_t, vector_t, vector_t, vector_t, vector_t>
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
reaction_flux
    ( vector_t v_log_counts
    , vector_t f_log_counts
    , matrix_t S
    , matrix_t K
    , matrix_t E_regulation
    );

