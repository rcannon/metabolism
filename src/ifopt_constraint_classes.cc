
#include "ifopt_constraint_classes.hh"

namespace ifopt {

//
// NullSpaceRepresentationConstraint
//

NullSpaceRepresentationConstraint::NullSpaceRepresentationConstraint
    ( const std::string& name
    , const int n_reactions
    , const std::string& flux_variables_name
    , const std::string& beta_variables_name
    , const matrix_t& null_space_matrix
    , const int dim_null_space
    )
    : ConstraintSet(n_reactions, name)
    , n_reactions_(n_reactions)
    , flux_variables_name_(flux_variables_name)
    , beta_variables_name_(beta_variables_name)
    , null_space_matrix_(null_space_matrix)
    , dim_null_space_(dim_null_space)
{
    assert(null_space_matrix_.rows() == n_reactions_);
    assert(null_space_matrix_.cols() == dim_null_space_);
}

vector_t 
NullSpaceRepresentationConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 94

    const vector_t beta_vec = GetVariables()->GetComponent(beta_variables_name_)->GetValues();
    const vector_t fluxes = GetVariables()->GetComponent(flux_variables_name_)->GetValues();

    const vector_t result = ((null_space_matrix_ * beta_vec) - fluxes);

    return result;
}

ConstraintSet::VecBound
NullSpaceRepresentationConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(0.0, 0.0));

    return b;
}

void
NullSpaceRepresentationConstraint::FillJacobianBlock
    ( std::string var_set
    , Jacobian& jac_block
    )
const
{
    int row_end = n_reactions_;

    if (var_set == flux_variables_name_) {
        // aka y in MEPPF
        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef(diag, diag) = -1.0;
        }
    }
    else if (var_set == beta_variables_name_) {
        // aka \beta in MEPPF
        for (int row = 0; row < row_end; row++) {
            for (int col = 0; col < dim_null_space_; col++) {
                jac_block.coeffRef(row, col) = null_space_matrix_(row,col);
            }
        }
    }
}

//
// End NullSpaceRepresentationConstraint
//

//
// SteadyStateConstraint
//

SteadyStateConstraint::SteadyStateConstraint
    ( const std::string& name
    , const int n_reactions
    , const int n_metabolites
    , const int n_variable_metabolites
    , const vector_t& fixed_metabolites
    , const std::string& variable_metabolites_name
    , const std::string& stead_state_variables_name
    , const matrix_t& variable_metabolite_stoich_matrix_T
    , const matrix_t& fixed_metabolite_stoich_matrix_T
    )
    : ConstraintSet(n_reactions, name )
    , n_reactions_(n_reactions)
    , n_metabolites_(n_metabolites)
    , n_variable_metabolites_(n_variable_metabolites)
    , fixed_metabolites_(fixed_metabolites)
    , variable_metabolites_name_(variable_metabolites_name)
    , steady_state_variables_name_(stead_state_variables_name)
    , variable_metabolite_stoich_matrix_T_(variable_metabolite_stoich_matrix_T)
    , fixed_metabolite_stoich_matrix_T_(fixed_metabolite_stoich_matrix_T)
{
    assert(fixed_metabolites_.size() == n_metabolites_ - n_variable_metabolites_);
    assert(variable_metabolite_stoich_matrix_T_.rows() == n_reactions_);
    assert(variable_metabolite_stoich_matrix_T_.cols() == n_variable_metabolites_);
    assert(fixed_metabolite_stoich_matrix_T_.rows() == n_reactions_);
    assert(fixed_metabolite_stoich_matrix_T_.cols() == fixed_metabolites_.size());
}

vector_t 
SteadyStateConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 95

    const vector_t steady_state_variables = GetVariables()->GetComponent(steady_state_variables_name_)->GetValues();
    const vector_t variable_metabolites = GetVariables()->GetComponent(variable_metabolites_name_)->GetValues();

    const vector_t result = (variable_metabolite_stoich_matrix_T_ * variable_metabolites) 
                    + (fixed_metabolite_stoich_matrix_T_ * fixed_metabolites_)
                    - steady_state_variables;

    return result;
}

ConstraintSet::VecBound
SteadyStateConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(0.0, 0.0));

    return b;
}

void
SteadyStateConstraint::FillJacobianBlock
    ( std::string var_set
    , Jacobian& jac_block
    )
const
{
    int row_end = n_reactions_;

    if (var_set == variable_metabolites_name_) {
        // aka \eta in MEPPF
        for (int row = 0; row < row_end; row++) {
            for (int col = 0; col < n_variable_metabolites_; col++) {
                jac_block.coeffRef(row, col) = variable_metabolite_stoich_matrix_T_(row, col);
            }
        }
    }
    else if (var_set == steady_state_variables_name_) {
        // aka g in MEPPF
        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef( diag, diag ) = -1.0;
        }
    }
}

//
// End SteadyStateConstraint
//

//
// SmoothConstraint
//

SmoothConstraint::SmoothConstraint
    ( const std::string& name
    , const int n_reactions
    , const std::string& flux_variables_name
    , const std::string& h_variables_name
    )
    : ConstraintSet(n_reactions, name)
    , n_reactions_(n_reactions)
    , flux_variables_name_(flux_variables_name)
    , h_variables_name_(h_variables_name)
{}

vector_t
SmoothConstraint::GetValues()
const
{
    const vector_t fluxes = GetVariables()->GetComponent(flux_variables_name_)->GetValues();
    const vector_t coeff = 
        ( fluxes.array().sign()
        *   ( std::log(2) 
            - Eigen::log( fluxes.array().abs()
                        + Eigen::sqrt( fluxes.array().pow(2) + 4 )
                        ) 
            )
        ).matrix();
    const vector_t h_variables = GetVariables()->GetComponent(h_variables_name_)->GetValues();
    const vector_t result = coeff - h_variables;

    return result;
}

ConstraintSet::VecBound
SmoothConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(0.0, 0.0));

    return b;
}

void
SmoothConstraint::FillJacobianBlock
    ( std::string var_set
    , Jacobian& jac_block
    )
const
{
    int row_end = n_reactions_;

    if (var_set == flux_variables_name_) {
        // aka y in MEPPF
        vector_t gradient = CalculateSmoothConstraintGradientFluxVariables().eval();
        for (int diag = 0; diag < row_end; diag ++) {
            jac_block.coeffRef( diag , diag ) = gradient(diag);
        }
    }
    else if (var_set == h_variables_name_) {
        // aka h in MEPPF
        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef(diag, diag) = -1.0;
        }
    }
}

vector_t
SmoothConstraint::CalculateSmoothConstraintGradientFluxVariables()
const
{ 
    const vector_t fluxes = GetVariables()->GetComponent(flux_variables_name_)->GetValues().array().matrix();
    
    const vector_t result = 
        (
            -1.0 / (fluxes.array().pow(2) + 4).sqrt()
        ).matrix();

    return result;
}

//
// End SmoothConstraint
//

//
// RelaxedRegulationUpperConstraint
//

RelaxedRegulationUpperConstraint::RelaxedRegulationUpperConstraint
    ( const std::string& name
    , const int n_reactions
    , const std::string& steady_state_variables_name
    , const std::string& h_variables_name
    , const std::string& u_variables_name
    , const double big_M_value
    , const vector_t& equilibrium_constants
    )
    : ConstraintSet( n_reactions, name)
    , n_reactions_(n_reactions)
    , steady_state_variables_name_(steady_state_variables_name)
    , h_variables_name_(h_variables_name)
    , u_variables_name_(u_variables_name)
    , big_M_value_(big_M_value)
    , equilibrium_constants_(equilibrium_constants)
{
    assert(equilibrium_constants_.size() == n_reactions_);
}

vector_t
RelaxedRegulationUpperConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 97
    
    const vector_t h_variables = GetVariables()->GetComponent(h_variables_name_)->GetValues();
    const vector_t u_variables = GetVariables()->GetComponent(u_variables_name_)->GetValues();
    const vector_t steady_state_variables = GetVariables()->GetComponent(steady_state_variables_name_)->GetValues();
    const vector_t log_equilibrium_constants = equilibrium_constants_.array().log().matrix();

    const vector_t result = h_variables - (big_M_value_ * u_variables) + log_equilibrium_constants - steady_state_variables;

    return result;
}

ConstraintSet::VecBound
RelaxedRegulationUpperConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(-inf, 0.0));

    return b;
}

void
RelaxedRegulationUpperConstraint::FillJacobianBlock
    ( std::string var_set
    , Jacobian& jac_block
    )
const
{
    int row_end = n_reactions_;

    if (var_set == steady_state_variables_name_) {
        // aka g in MEPPF
        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef( diag, diag) = -1.0;
        }

    }
    else if (var_set == h_variables_name_) {
        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef( diag, diag) = 1.0;
        }
    }
    else if (var_set == u_variables_name_) {
        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef( diag, diag ) = -big_M_value_;
        }
    }
}

//
// End RelaxedRegulationUpperConstraint
//

//
// RelaxedRegulationLowerConstraint
//

RelaxedRegulationLowerConstraint::RelaxedRegulationLowerConstraint
    ( const std::string& name
    , const int n_reactions
    , const std::string& steady_state_variables_name
    , const std::string& h_variables_name
    , const std::string& u_variables_name
    , const double big_M_value
    , const vector_t& equilibrium_constants
    )
    : ConstraintSet( n_reactions, name)
    , n_reactions_(n_reactions)
    , steady_state_variables_name_(steady_state_variables_name)
    , h_variables_name_(h_variables_name)
    , u_variables_name_(u_variables_name)
    , big_M_value_(big_M_value)
    , equilibrium_constants_(equilibrium_constants)
{
    assert(equilibrium_constants_.size() == n_reactions_);
}

vector_t
RelaxedRegulationLowerConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 98
    
    const vector_t h_variables = GetVariables()->GetComponent(h_variables_name_)->GetValues();
    const vector_t u_variables = GetVariables()->GetComponent(u_variables_name_)->GetValues();
    const vector_t one_minus_u = (1 - u_variables.array()).matrix();
    const vector_t steady_state_variables = GetVariables()->GetComponent(steady_state_variables_name_)->GetValues();
    const vector_t log_equilibrium_constants = equilibrium_constants_.array().log().matrix();

    const vector_t result = h_variables + (big_M_value_ * one_minus_u) + log_equilibrium_constants - steady_state_variables;

    return result;
}

ConstraintSet::VecBound
RelaxedRegulationLowerConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(0.0, +inf));

    return b;
}

void
RelaxedRegulationLowerConstraint::FillJacobianBlock
    ( std::string var_set
    , Jacobian& jac_block
    )
const
{
    int row_end = n_reactions_;

    if (var_set == steady_state_variables_name_) {
        // aka g in MEPPF
        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef( diag, diag) = -1.0;
        }

    }
    else if (var_set == h_variables_name_) {
        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef( diag, diag) = 1.0;
        }
    }
    else if (var_set == u_variables_name_) {
        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef( diag, diag ) = -big_M_value_;
        }
    }
}

//
// End RelaxedRegulationLowerConstraint
//

//
// SignConstraint
//

SignConstraint::SignConstraint
    ( const std::string& name
    , const int n_reactions
    , const std::string& steady_state_variables_name
    , const std::string& flux_variables_name
    , const vector_t equilibrium_constants
    )
    : ConstraintSet( n_reactions, name)
    , n_reactions_(n_reactions)
    , steady_state_variables_name_(steady_state_variables_name)
    , flux_variables_name_(flux_variables_name)
    , equilibrium_constants_(equilibrium_constants)
{}

vector_t
SignConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 99

    const vector_t fluxes = GetVariables()->GetComponent(flux_variables_name_)->GetValues();
    const vector_t log_equilibrium_constants = equilibrium_constants_.array().log().matrix();
    const vector_t steady_state_variables = GetVariables()->GetComponent(steady_state_variables_name_)->GetValues();

    const vector_t result = 
        (
            (log_equilibrium_constants.array() - steady_state_variables.array()) * fluxes.array()
        ).matrix();

    return result.matrix();
}

ConstraintSet::VecBound
SignConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(0.0, +inf));

    return b;
}

void
SignConstraint::FillJacobianBlock
    ( std::string var_set
    , Jacobian& jac_block
    )
const
{
    int row_end = n_reactions_;

    if (var_set == flux_variables_name_) {
        // aka y in MEPPF
        const vector_t gradient = CalculateSignConstraintGadientFluxVariables().eval();
        for (int diag = 0; diag < row_end; diag ++) {
            jac_block.coeffRef( diag, diag ) = gradient(diag);
        }

    }
    else if (var_set == steady_state_variables_name_) {
        // aka g in MEPPF
        const vector_t gradient = CalculateSignConstraintGradientSteadyStateVariables().eval();
        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef( diag, diag ) = gradient(diag);
        }
    }
}

vector_t
SignConstraint::CalculateSignConstraintGadientFluxVariables() 
const
{
    const vector_t steady_state_variables = GetVariables()->GetComponent(steady_state_variables_name_)->GetValues();
    const vector_t result = equilibrium_constants_.array().log().matrix() - steady_state_variables;

    return result;
}

vector_t
SignConstraint::CalculateSignConstraintGradientSteadyStateVariables()
const
{
    const vector_t flux_variables = GetVariables()->GetComponent(flux_variables_name_)->GetValues();
    const vector_t result = -1.0 * flux_variables;

    return result;
}

//
// End SignConstraint
//

//
// RelaxedRegulationSignConstraint
//

RelaxedRegulationSignConstraint::RelaxedRegulationSignConstraint
    ( const std::string& name
    , const int n_reactions
    , const std::string& flux_variables_name
    , const std::string& u_variables_name
    )
    : ConstraintSet( n_reactions, name)
    , n_reactions_(n_reactions)
    , flux_variables_name_(flux_variables_name)
    , u_variables_name_(u_variables_name)
{}

vector_t
RelaxedRegulationSignConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 100

    const vector_t fluxes = GetVariables()->GetComponent(flux_variables_name_)->GetValues();
    const vector_t sign_fluxes = 
        (
            fluxes.array() / (fluxes.array().abs() + std::pow(10, -50))
        ).matrix();
    const vector_t u_variables = GetVariables()->GetComponent(u_variables_name_)->GetValues();

    const vector_t result = 
        (
            (2 * u_variables.array()) - 1 - sign_fluxes.array()
        ).matrix();

    return result;
}

ConstraintSet::VecBound
RelaxedRegulationSignConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(0.0, 0.0));

    return b;
}

void
RelaxedRegulationSignConstraint::FillJacobianBlock
    ( std::string var_set
    , Jacobian& jac_block
    )
const
{
    
    int row_end = n_reactions_;

    const vector_t fluxes = GetVariables()->GetComponent(flux_variables_name_)->GetValues();
    const vector_t gradient = 
        (
            (-1.0 * std::pow(10, -50)) / (fluxes.array().abs() + std::pow(10, -50)).pow(2)
        ).matrix();
    if (var_set == flux_variables_name_){
        for (int diag =0; diag < row_end; diag++) {
            jac_block.coeffRef(diag, diag) = gradient(diag);
        }
    }
    else if (var_set == u_variables_name_){
        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef( diag, diag) = 2.0;
        }
    }
    
}

//
// End RelaxedRegulationSignConstraint
//

//
// MetabolitesUpperBoundConstraint
//

MetabolitesUpperBoundConstraint::MetabolitesUpperBoundConstraint
    ( const std::string& name
    , const int n_variable_metabolites
    , const std::string& variable_metabolites_name
    , const vector_t& variable_metabolites_upper_bound
    )
    : ConstraintSet(n_variable_metabolites, name)
    , n_variable_metabolites_(n_variable_metabolites)
    , variable_metabolites_name_(variable_metabolites_name)
    , variable_metabolites_upper_bound_(variable_metabolites_upper_bound)
{
    assert(variable_metabolites_upper_bound_.size() == n_variable_metabolites_);
}

vector_t
MetabolitesUpperBoundConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 101 upper bounds

    const vector_t variable_metabolites = GetVariables()->GetComponent(variable_metabolites_name_)->GetValues();

    const vector_t result = 
        (
            variable_metabolites_upper_bound_.array() - variable_metabolites.array()
        ).matrix();

    return result;
}

ConstraintSet::VecBound
MetabolitesUpperBoundConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(0.0, +inf));

    return b;
}

void
MetabolitesUpperBoundConstraint::FillJacobianBlock
    ( std::string var_set
    , Jacobian& jac_block
    )
const
{
    int row_end = n_variable_metabolites_;

    if (var_set == variable_metabolites_name_){

        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef( diag, diag ) = -1.0;
        }
    }
}

//
// End MetabolitesUpperBoundConstraint
//

//
// MetaboliteLowerBoundConstraint
//

MetabolitesLowerBoundConstraint::MetabolitesLowerBoundConstraint
    ( const std::string& name
    , const int n_variable_metabolites
    , const std::string& variable_metabolites_name
    , const double variable_metabolites_lower_bound
    )
    : ConstraintSet(n_variable_metabolites, name)
    , n_variable_metabolites_(n_variable_metabolites)
    , variable_metabolites_name_(variable_metabolites_name)
    , variable_metabolites_lower_bound_(variable_metabolites_lower_bound)
{}

vector_t
MetabolitesLowerBoundConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 101 lower bounds

    const vector_t variable_metabolites = GetVariables()->GetComponent(variable_metabolites_name_)->GetValues();

    const vector_t result = 
        (
            variable_metabolites_lower_bound_ - variable_metabolites.array()
        ).matrix();

    return result;
}

ConstraintSet::VecBound
MetabolitesLowerBoundConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(-inf, 0.0));

    return b;
}

void
MetabolitesLowerBoundConstraint::FillJacobianBlock
    ( std::string var_set
    , Jacobian& jac_block
    )
const
{
    int row_end = n_variable_metabolites_;

    if (var_set == variable_metabolites_name_){

        for (int diag = 0; diag < row_end; diag++) {
            jac_block.coeffRef( diag, diag ) = -1.0;
        }
    }
}

//
// End MetaboliteLowerBoundConstraint
//

} // namespace ifopt