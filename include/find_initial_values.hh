
#include "includes_and_types.hh"
#include "read_data_files.hh"

#pragma once

std::tuple<vector_t, vector_t, vector_t>
find_initial_values
    ( const vector_t& variable_metabolites
    , const vector_t& fixed_metabolites
    , const vector_t& target_log_counts
    , const matrix_t& stoichiometric_matrix
    , const vector_t& equilibrium_constants
    );