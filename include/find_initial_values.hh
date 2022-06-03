
#include "includes_and_types.hh"
#include "read_data_files.hh"

#pragma once

std::tuple<vector_t, vector_t, vector_t>
find_initial_values
    ( const matrix_t& stoichiometric_matrix
    );