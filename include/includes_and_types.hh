
#pragma once

/* general */
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <numeric>
#include <math.h>

/* IFOPT */
#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>

/* Eigen */
#include "Eigen/Core"
#include "Eigen/Dense"

/* Type definitions and aliases */
// #define EIGEN_DEFAULT_DENSE_INDEX_TYPE = int;
typedef double value_t;
typedef std::vector<value_t> value_list_t;
typedef std::vector<int> index_list_t;
typedef Eigen::MatrixXd matrix_t;
typedef Eigen::VectorXd vector_t;