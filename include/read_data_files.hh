
#include "includes_and_types.hh"

#pragma once

using namespace Eigen;

std::tuple<size_t,vector_t>
load_concentrations( const std::string& path);

std::tuple<size_t,size_t, index_list_t, vector_t>
load_equilibrium_constants( const std::string& path);

matrix_t
load_stochiometric_matrix( const std::string& path);