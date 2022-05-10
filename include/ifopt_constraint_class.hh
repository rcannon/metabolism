
#include "includes_and_types.hh"

#pragma once

namespace ifopt {

/* specific constraint_class that holds all constraints for the problem */
class Constraints : public ConstraintSet {
public:

    Constraints ( const std::string& name
                , const int n_metabolites
                , const int n_variable_metabolites
                , const vector_t& fixed_metabolites
                , const std::string& variable_metabolites_name
                , const std::string& flux_variables_name
                , const std::string& steady_state_variables_name // S^T * eta 
                , const std::string& h_variables_name
                , const std::string& u_variable_name
                , const std::string& beta_variables_name
                , const matrix_t& null_space_matrix
                , const int null_space_dimension
                , const matrix_t& stochiometric_matrix
                , const double big_M_value
                , const vector_t& equilibrium_constants
                , const double variable_metabolites_upper_bound
                , const double variable_metabolites_lower_bound
                ) 
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
    const int n_reactions_;
    const vector_t fixed_metabolites_; // no reference
    const std::string variable_metabolites_name_;
    const std::string flux_variables_name_;
    const std::string steady_state_variables_name_; // S^T * eta 
    const std::string h_variables_name_;
    const std::string u_variables_name_;
    const std::string beta_variables_name_;
    const matrix_t null_space_matrix_;
    const double null_space_dimension_;
    const matrix_t stochiometric_matrix_T_;
    const double big_M_value_;
    const vector_t equilibrium_constants_;
    const vector_t variable_metabolites_upper_bound_;
    const double variable_metabolites_lower_bound_;
};

} // namespace