
#include "includes_and_types.hh"
#include "ifopt_variable_class.hh"
#include "ifopt_constraint_class.hh"
#include "ifopt_cost_class.hh"
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include <ifopt/test_vars_constr_cost.h>

#pragma once

std::tuple<vector_t, vector_t, vector_t, vector_t, vector_t>
maximum_entropy_solver
    ( vector_t variable_metabolites_init // aka n
    , vector_t fixed_metabolites
    , vector_t flux_variables_init // aka y
    , vector_t beta_init // \beta
    , vector_t target_log_variable_metabolites_counts
    , matrix_t stoichiometric_matrix
    , vector_t equilibrium_constants // aka K
    , index_list_t objective_reaction_indices
    );

vector_t
reaction_flux
    ( vector_t variable_metabolites_log_counts
    , vector_t fixed_metabolites_log_counts
    , matrix_t stochiometric_matrix
    , matrix_t equilibrium_constants
    , matrix_t E_regulation
    );

