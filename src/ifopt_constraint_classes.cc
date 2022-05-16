
#include "ifopt_constraint_classes.hh"

namespace ifopt {

//
// NullSpaceRepresentationConstraint
//

NullSpaceRepresentationConstraint::NullSpaceRepresentationConstraint
    ( const string& name
    , const int n_reactions
    , const std::string& beta_variables_name
    , const std::string& flux_variables_name
    , const matrix_t& null_space_matrix
    , const int dim_null_space
    )
    : ConstraintSet(n_reactions, name)
    , n_reactions_(n_reactions)
    , beta_variables_name_(beta_variables_name)
    , flux_variables_name_(flux_variables_name)
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

    vector_t beta_vec = GetVariables()->GetComponent(beta_variables_name_)->GetValues();
    vector_t fluxes = GetVariables()->GetComponent(flux_variables_name_)->GetValues();

    vector_t result = ((null_space_matrix_ * beta_vec) - fluxes);

    return result;
}

ConstraintSet::VecBound
NullSpaceRepresentationConstraint()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(0.0, 0.0));

    return b;
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
    , const matrix_t& stoichiometric_matrix_T
    )
    : ConstraintSet(n_reactions, name )
    , n_reactions_(n_reactions)
    , n_metabolites_(n_metabolites)
    , n_variable_metabolites_(n_variable_metabolites)
    , fixed_metabolites_(fixed_metabolites)
    , variable_metabolites_name_(variable_metabolites_name)
    , steady_state_variables_name_(stead_state_variables_name)
    , stoichiometric_matrix_T_(stoichiometric_matrix_T)
{
    assert(fixed_metabolites_.size() == n_metabolites_ - n_variable_metabolites_);
    assert(stoichiometric_matrix_T_.rows() == n_reactions_);
    assert(stoichiometric_matrix_T_.cols() == n_metabolites_);
}

vector_t 
SteadyStateConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 95

    vector_t steady_state_variables = GetVariables()->GetComponent(steady_state_variables_name_)->GetValues();

    vector_t variable_metabolites = GetVariables()->GetComponent(variable_metabolites_name_)->GetValues();
    vector_t all_metabolites( n_metabolites_ );
    metabolites << variable_metabolites, fixed_metabolites_;

    vector_t result = ((stoichiometric_matrix_T_ * metabolites) - steady_state_variables);

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
    // Maximum Entropy Production Problem Formulation 96

    auto h_variables = GetVariables()->GetComponent(h_variables_name_)->GetValues().array();
    auto fluxes_array = GetVariables()->GetComponent(flux_variables_name_)->GetValues().array();
    auto sign_fluxes = fluxes_array.sign();
    auto abs_fluxes = fluxes_array.abs();
    auto square_fluxes = fluxes_array.pow(2);

    auto result = sign_fluxes
                * ( std::log(2) 
                  - Eigen::log( abs_fluxes
                              + Eigen::sqrt( square_fluxes + 4 )
                              ) 
                  )
                - h_variables;

    return result.matrix();
}

ConstraintSet::VecBound
SmoothConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(0.0, 0.0));

    return b;
}

//
// End SmoothConstraint
//

//
// RelaxedFluxUpperConstraint
//

RelaxedFluxUpperConstraint::RelaxedFluxUpperConstraint
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
    assert(equilibrium_constants_.size() == n_reactions_)
}

vector_t
RelaxedFluxUpperConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 97
    
    vector_t h_variables = GetVariables()->GetComponent(h_variables_name_)->GetValues();
    vector_t u_variables = GetVariables()->GetComponent(u_variables_name_)->GetValues();
    vector_t steady_state_variables = GetVariables()->GetComponent()(steady_state_variables_name_)->GetValues();
    vector_t log_equilibrium_constants = equilibrium_constants_.array().log().matrix();

    vector_t result = h_variables - (big_M_value_ * u_variables) + log_equilibrium_constants - steady_state_variables;

    return result;
}

ConstraintSet::VecBound
RelaxedFluxUpperConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(-inf, 0.0));

    return b;
}

//
// End RelaxedFluxUpperConstraint
//

//
// RelaxedFLuxLowerConstraint
//

RelaxedFluxLowerConstraint:RelaxedFluxLowerConstraint
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
    assert(equilibrium_constants_.size() == n_reactions_)
}

vector_t
RelaxedFluxLowerConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 98
    
    vector_t h_variables = GetVariables()->GetComponent(h_variables_name_)->GetValues();
    vector_t u_variables = GetVariables()->GetComponent(u_variables_name_)->GetValues();
    vector_t one_minus_u = (1 - u_variables.array()).matrix();
    vector_t steady_state_variables = GetVariables()->GetComponent(steady_state_variables_name_)->GetValues();
    vector_t log_equilibrium_constants = equilibrium_constants_.array().log().matrix();

    vector_t result = h_variables + (big_M_value_ * one_minus_u) + log_equilibrium_constants - steady_state_variables;

    return result;
}

ConstraintSet::VecBound
RelaxedFluxLowerConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(0.0, +inf));

    return b;
}

//
// End RelaxedFLuxLowerConstraint
//

//
// SignConstraint
//

SignConstraint::SignConstraint
    ( const std::string& name
    , const int n_reactions
    , const std::string& steady_state_variable_names
    , const std::string& flux_variable_names
    , const vector_t equilibrium_constants
    )
    : ConstraintSet( n_reactions, name)
    , n_reactions_(n_reactions)
    , steady_state_variable_names_(steady_state_variable_names)
    , flux_variable_names_(flux_variable_names)
    , equilibrium_constants_(equilibrium_constants)
{}

vector_t
SignConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 99

    auto fluxes_array = GetVariables()->GetComponent(flux_variables_name_)->GetValues().array();
    auto log_equilibrium_constants_array = equilibrium_constants_.array().log();
    auto steady_state_variables_array = GetVariables()->GetComponent(steady_state_variables_name_)->GetValues().array();

    auto result = (log_equilibrium_constants_array - steady_state_variables_array) *  fluxes_array;

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

//
// End SignConstraint
//

//
// RelaxedFluxSignConstraint
//

RelaxedFluxSignConstraint::RelaxedFluxSignConstraint
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
RelaxedFluxSignConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 100

    auto fluxes_array = GetVariables()->GetComponent(flux_variables_name_)->GetValues().array();
    auto sign_fluxes = fluxes_array.sign();
    auto u_variables_array = GetVariables()->GetComponent(u_variables_name_)->GetValues().array();

    auto result = (2 * u_variables_array) - 1 - sign_fluxes;

    return result.matrix();
}

ConstraintSet::VecBound
RelaxedFluxSignConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(0.0, 0.0));

    return b;
}

//
// End RelaxedFluxSignConstraint
//

//
// MetabolitesUpperBoundConstraint
//

MetabolitesUpperBoundConstraint::MetabolitesUpperBoundConstraint
    ( const std::string& name
    , const int n_variable_metabolites
    , const std::string& variables_metabolites_name
    , const vector_t& variable_metabolites_upper_bound
    )
    : ConstraintSet(n_variable_metabolites, name)
    , n_variable_metabolites_(n_variable_metabolites)
    , variables_metabolites_name_(variables_metabolites_name)
    , variable_metabolites_upper_bound_(variable_metabolites_upper_bound)
{}

vector_t
MetabolitesUpperBoundConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 101 upper bounds

    vector_t variable_metabolites = GetVariables()->GetComponent(variable_metabolites_name_)->GetValues();
    auto metabolites_array = variable_metabolites.array();

    auto result = variable_metabolites_upper_bound_.array() - metabolites_array;

    return result.matrix();
}

ConstraintSet::VecBound
MetabolitesUpperBoundConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(0.0, +inf));

    return b;
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
    , const std::string& variables_metabolites_name
    , const vector_t& variable_metabolites_lower_bound
    )
    : ConstraintSet(n_variable_metabolites, name)
    , n_variable_metabolites_(n_variable_metabolites)
    , variables_metabolites_name_(variables_metabolites_name)
    , variable_metabolites_lower_bound_(variable_metabolites_upper_bound)
{}

vector_t
MetabolitesLowerBoundConstraint::GetValues()
const
{
    // Maximum Entropy Production Problem Formulation 101 lower bounds

    vector_t variable_metabolites = GetVariables()->GetComponent(variable_metabolites_name_)->GetValues();
    auto metabolites_array = variable_metabolites.array();

    auto result = variable_metabolites_lower_bound_ - metabolites_array;

    return result.matrix();
}

ConstraintSet::VecBound
MetabolitesLowerBoundConstraint::GetBounds()
const
{
    VecBound b(GetRows());
    std::fill(b.begin(), b.end(), Bounds(-inf, 0.0));

    return b;
}

//
// End MetaboliteLowerBoundConstraint
//

} // namespace ifopt