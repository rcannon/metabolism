
#pragma once

/* general */
#include <iostream>
#include <vector>
#include <string>
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
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE = int;
typedef double value_t;
typedef std::vector<value_t> value_list_t;
typedef std::vector<int> index_list_t;
typedef Eigen::MatrixXd matrix_t;
typedef Eigen::VectorXd vector_t;

namespace ifopt {

/* generic variable class */
class Variables : public VariableSet {
public:

    Variables( const std::string& name
             , const int num_variables
             );
    Variables( const std::string& name
             , const int num_variables
             , const vector_t& init_values
             )
    void SetVariables( const vector_t& new_vars ) override;
    vector_t GetValues() const override;
    int GetNumVariables() const;

private:
    const int num_variables_;
    vector_t variables_;

} /* end variables */

/* arguement struct for initializing constraint class */
struct ConstraintArguements {
    std::string& name;
    int n_metabolites;
    int n_variable_metabolites;
    vector_t& fixed_metabolites;
    std::string& variable_metabolites_name;
    std::string& flux_variables_name;
    std::string& steady_state_variables_name; // S^T * eta 
    std::string& h_variables_name;
    std::string& u_variable_name;
    std::string& beta_variables_name;
    matrix_t& null_space_matrix;
    int null_space_dimension;
    matrix_t& stochiometric_matrix;
    double big_M_value;
    vector_t& K_vector;
    double variable_metabolites_upper_bound;
    double variable_metabolites_lower_bound;
}; /* end arguement struct */

/* specific constraint_class that holds all constraints for the problem */
class Constraints : public ConstraintSet {
public:

    Constraints( const ConstraintArguements& args );
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock( std::string var_set
                          , Jacobian& jac_block
                          ) const override;

private:

    vector_t GetVariablesVectorByName( const std::string& variables_name) const;
    vector_t GetAllMetabolites() const;
    vector_t CalculateNullSpaceRepresentationConstraint() const;
    vector_t CalculateSteadyStateConstraint() const;
    vector_t CalculateSmoothConstraint() const;
    vector_t CalculateRelaxedFluxUpperConstraint() const;
    vector_t CalculateRelaxedFluxLowerConstraint() const;
    vector_t CalculateSignConstraint() const;
    vector_t CalculateRelaxedFluxSignConstraint() const;
    vector_t CalculateMetabolitesUpperBoundConstraint() const;
    vector_t CalculateMetaboliteLowerBoundConstraint() const;
    vector_t CalculateSmoothConstraintGradientFluxVariables() const;
    vector_t CalculateSignConstraintGadientFluxVariables() const;
    vector_t CalculateSignConstraintGradientSteadyStateVariables() const;

    const int n_metabolites_;
    const int n_variable_metabolites_;
    const vector_t fixed_metabolites_; // no reference
    const std::string variable_metabolites_name_;
    const std::string flux_variables_name_;
    const std::string steady_state_variables_name_; // S^T * eta 
    const std::string h_variables_name_;
    const std::string u_variables_name_;
    const std::string beta_variables_name_;
    const matrix_t null_space_matrix_;
    const double null_space_dimension_;
    const matrix_t stochiometric_matrix_;
    const double big_M_value_;
    const vector_t K_vector_;
    const double variable_metabolites_upper_bound_;
    const double variable_metabolites_lower_bound_;
}


/* specific cost function class */
class Cost : public CostTerm {
public:

    MetabCost( const std::string& name 
             , const std::string& reaction_vars_name 
             );
    double GetCost() const override;
    void FillJacobianBlock( std::string var_set
                          , Jacobian& jac
                          ) const override;
private:
    const std::string reaction_vars_name_;
};

} /* ifopt namespace */



