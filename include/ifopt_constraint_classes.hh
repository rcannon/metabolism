
#include "includes_and_types.hh"

#pragma once

namespace ifopt {


class NullSpaceRepresentationConstraint : public ConstraintSet {
    // MEPPF 94
    NullSpaceRepresentationConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    const int n_reactions_;
    const std::string beta_variables_name_;
    const std::string flux_variables_name_;
    const matrix_t null_space_matrix;
    const int dim_null_space;
}


class SteadyStateConstraint : public ConstraintSet {
    // MEPPF 95
public:
    SteadyStateConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    const int n_reactions_;
    const int n_metabolites_;
    const int n_variable_metabolites_;
    const vector_t fixed_metabolites_;
    const std::string variable_metabolites_name_;
    const std::string stead_state_variables_name_;
    const matrix_t variable_metabolite_stoich_matrix_T
    const matrix_t fixed_metabolite_stoich_matrix_T
}


class SmoothConstraint : public ConstraintSet {
    // MEPPF 96
public:
    SmoothConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    vector_t CalculateSmoothConstraintGradientFluxVariables() const;
    const int n_reactions_;
    const std::string flux_variables_name_;
    const std::string h_variables_name_;
}


class RelaxedFluxUpperConstraint : public ConstraintSet {
    // MEPPF 97
public:
    RelaxedFluxUpperConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    const int n_reactions_;
    const std::string steady_state_variables_name_;
    const std::string h_variables_name_;
    const std::string u_variables_name_;
    const double big_M_value_;
    const vector_t equilibrium_constants_;
}


class RelaxedFluxLowerConstraint : public ConstraintSet {
    // MEPPF 98
public:
    RelaxedFluxLowerConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    const int n_reactions_;
    const std::string steady_state_variables_name_;
    const std::string h_variables_name_;
    const std::string u_variables_name_;
    const double big_M_value_;
    const vector_t equilibrium_constants_;
}


class SignConstraint : public ConstraintSet {
    // MEPPF 99
public:
    SignConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    vector_t CalculateSignConstraintGadientFluxVariables() const;
    vector_t CalculateSignConstraintGradientSteadyStateVariables() const;
    const int n_reactions_;
    const std::string steady_state_variable_names_;
    const std::string flux_variable_names_;
    const vector_t equilibrium_constants_;
}


class RelaxedFluxSignConstraint : public ConstraintSet {
    // MEPPF 100
public:
    RelaxedFluxSignConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    const int n_reactions_;
    const std::string flux_variables_name_;
    const std::string u_variables_name_;
}


class MetabolitesUpperBoundConstraint : public ConstraintSet {
    // MEPPF 101 upper
public:
    MetabolitesUpperBoundConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    const int n_variable_metabolites_;
    const std::string variables_metabolites_name_;
    const vector_t variable_metabolites_upper_bound_;
}


class MetabolitesLowerBoundConstraint : public ConstraintSet {
    // MEPPF 101 lower
public:
    MetabolitesLowerBoundConstraint();
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    const int n_variable_metabolites_;
    const std::string variables_metabolites_name_;
    const vector_t variable_metabolites_lower_bound_;
}

} // namespace ifopt