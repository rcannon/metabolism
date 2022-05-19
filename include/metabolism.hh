
#include "includes_and_types.hh"
#include "maximum_entropy_relaxed.hh"

#pragma once

std::tuple<vector_t, vector_t, vector_t>
initialize 
    ( const vector_t& variable_metabolites
    , const vector_t& fixed_metabolites
    , const vector_t& target_log_counts
    , const matrix_t& stoichiometric_matrix
    , const vector_t& equilibrium_constants
    );

std::tuple<vector_t, vector_t, vector_t, vector_t>
optimize
    ( const vector_t& variable_metabolites_init
    , const vector_t& flux_variables_init
    , const vector_t& beta_init
    , const vector_t& fixed_metabolite_log_counts
    , const vector_t& target_log_variable_metabolites_counts
    , const matrix_t& stoichiometric_matrix
    , const vector_t& equilibrium_constants
    , const index_list_t& objective_reaction_indices
    );

vector_t
odds_diff
    ( const vector_t& variable_metabolites_counts
    , const vector_t& fixed_metabolites_counts
    , const double mu0
    , const matrix_t& stoichiometric_matrix
    , const matrix_t& negative_stoich_matrix
    , const matrix_t& positive_stoich_matrix
    , const vector_t& delta
    , const vector_t& equilibrium_constants
    , const vector_t& e_regulation
    );

vector_t
calc_odds
    ( const vector_t& log_counts
    , const matrix_t& negative_stoich_matrix
    , const matrix_t& positive_stoich_matrix
    , const vector_t& delta
    , const vector_t& equilibrium_constants
    , const double direction = 1.0
    );

std::tuple<vector_t, vector_t>
run ( const vector_t& variable_metabolites
    , const vector_t& fixed_metabolites
    , const vector_t& target_log_variable_metabolites_count
    , const matrix_t& stoichiometric_matrix
    , const vector_t& equilibrium_constants
    , const int i_uptake
    , const int Vmax
    , const double Km
    , const double s
    , const index_list_t& objective_reaction_indices
    );