
#include "includes_and_types.hh"

#pragma once

namespace ifopt {


class NullSpaceRepresentationConstraint : public ConstraintSet {
    // MEPPF 94
public:
    NullSpaceRepresentationConstraint
        ( const std::string& name
        , const int n_reactions
        , const std::string& flux_variables_name
        , const std::string& beta_variables_name
        , const matrix_t& null_space_matrix
        , const int dim_null_space
        );
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    const int n_reactions_;
    const std::string beta_variables_name_;
    const std::string flux_variables_name_;
    const matrix_t null_space_matrix_;
    const int dim_null_space_;
};


class SteadyStateConstraint : public ConstraintSet {
    // MEPPF 95
public:
    SteadyStateConstraint
        ( const std::string& name
        , const int n_reactions
        , const int n_metabolites
        , const int n_variable_metabolites
        , const vector_t& fixed_metabolites
        , const std::string& variable_metabolites_name
        , const std::string& stead_state_variables_name
        , const matrix_t& variable_metabolite_stoich_matrix_T
        , const matrix_t& fixed_metabolite_stoich_matrix_T
        );
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
    const std::string steady_state_variables_name_;
    const matrix_t variable_metabolite_stoich_matrix_T_;
    const matrix_t fixed_metabolite_stoich_matrix_T_;
};


class SmoothConstraint : public ConstraintSet {
    // MEPPF 96
public:
    SmoothConstraint
        ( const std::string& name
        , const int n_reactions
        , const std::string& flux_variables_name
        , const std::string& h_variables_name
        );
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
};


class RelaxedRegulationUpperConstraint : public ConstraintSet {
    // MEPPF 97
public:
    RelaxedRegulationUpperConstraint
        ( const std::string& name
        , const int n_reactions
        , const std::string& steady_state_variables_name
        , const std::string& h_variables_name
        , const std::string& u_variables_name
        , const double big_M_value
        , const vector_t& equilibrium_constants
        );
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
};


class RelaxedRegulationLowerConstraint : public ConstraintSet {
    // MEPPF 98
public:
    RelaxedRegulationLowerConstraint
        ( const std::string& name
        , const int n_reactions
        , const std::string& steady_state_variables_name
        , const std::string& h_variables_name
        , const std::string& u_variables_name
        , const double big_M_value
        , const vector_t& equilibrium_constants
        );
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
};


class SignConstraint : public ConstraintSet {
    // MEPPF 99
public:
    SignConstraint
        ( const std::string& name
        , const int n_reactions
        , const std::string& steady_state_variables_name
        , const std::string& flux_variables_name
        , const vector_t equilibrium_constants
        );
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    vector_t CalculateSignConstraintGadientFluxVariables() const;
    vector_t CalculateSignConstraintGradientSteadyStateVariables() const;
    const int n_reactions_;
    const std::string steady_state_variables_name_;
    const std::string flux_variables_name_;
    const vector_t equilibrium_constants_;
};


class RelaxedRegulationSignConstraint : public ConstraintSet {
    // MEPPF 100
public:
    RelaxedRegulationSignConstraint
        ( const std::string& name
        , const int n_reactions
        , const std::string& flux_variables_name
        , const std::string& u_variables_name
        );
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    const int n_reactions_;
    const std::string flux_variables_name_;
    const std::string u_variables_name_;
};


class MetabolitesUpperBoundConstraint : public ConstraintSet {
    // MEPPF 101 upper
public:
    MetabolitesUpperBoundConstraint
        ( const std::string& name
        , const int n_variable_metabolites
        , const std::string& variable_metabolites_name
        , const vector_t& variable_metabolites_upper_bound
        );
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    const int n_variable_metabolites_;
    const std::string variable_metabolites_name_;
    const vector_t variable_metabolites_upper_bound_;
};


class MetabolitesLowerBoundConstraint : public ConstraintSet {
    // MEPPF 101 lower
public:
    MetabolitesLowerBoundConstraint
        ( const std::string& name
        , const int n_variable_metabolites
        , const std::string& variables_metabolites_name
        , const double variable_metabolites_lower_bound
        );
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    void FillJacobianBlock  ( std::string var_set
                            , Jacobian& jac_block
                            ) const override;
private:
    const int n_variable_metabolites_;
    const std::string variable_metabolites_name_;
    const double variable_metabolites_lower_bound_;
};

} // namespace ifopt