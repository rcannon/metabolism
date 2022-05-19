
#include "includes_and_types.hh"
#include "ifopt_variable_class.hh"
#include "ifopt_constraint_classes.hh"
#include "ifopt_cost_class.hh"
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include <ifopt/test_vars_constr_cost.h>

#pragma once

std::tuple<vector_t, vector_t, vector_t, vector_t, vector_t, vector_t>
maximum_entropy_solver
    ( const vector_t& variable_metabolites_init // aka n
    , const vector_t& fixed_metabolites
    , const vector_t& flux_variables_init // aka y
    , const vector_t& beta_init // \beta
    , const vector_t& target_log_variable_metabolites_counts
    , const matrix_t& stoich_matrix
    , const vector_t& equilibrium_constants // aka K
    , const index_list_t& objective_reaction_indices
    );

vector_t
reaction_flux
    ( const vector_t& variable_metabolites_log_counts
    , const vector_t& fixed_metabolites_log_counts
    , const matrix_t& stoich_matrix
    , const matrix_t& equilibrium_constants
    , const matrix_t& E_regulation
    );

