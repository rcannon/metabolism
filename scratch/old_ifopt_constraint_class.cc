
#include "ifopt_constraint_class.hh"

namespace ifopt {

Constraints::Constraints( const std::string& name
                        , const int n_metabolites
                        , const int n_variable_metabolites
                        , const int n_reactions
                        , const vector_t& fixed_metabolites
                        , const std::string& variable_metabolites_name
                        , const std::string& flux_variables_name
                        , const std::string& steady_state_variables_name // S^T * eta 
                        , const std::string& h_variables_name
                        , const std::string& u_variable_name
                        , const std::string& beta_variables_name
                        , const matrix_t& null_space_matrix
                        , const int null_space_dimension
                        , const matrix_t& stoichiometric_matrix_T
                        , const double big_M_value
                        , const vector_t& equilibrium_constants
                        , const vector_t& variable_metabolites_upper_bound
                        , const double variable_metabolites_lower_bound
                        )
                        : ConstraintSet( 7 * n_reactions + 2 * n_variable_metabolites
                                       , name)
                        , n_metabolites_(n_metabolites)
                        , n_variable_metabolites_(n_variable_metabolites)
                        , n_reactions_(n_reactions)
                        , fixed_metabolites_(fixed_metabolites)
                        , variable_metabolites_name_(variable_metabolites_name)
                        , flux_variables_name_(flux_variables_name)
                        , steady_state_variables_name_(steady_state_variables_name)
                        , h_variables_name_(h_variables_name)
                        , u_variables_name_(u_variable_name)
                        , beta_variables_name_(beta_variables_name)
                        , null_space_matrix_(null_space_matrix)
                        , null_space_dimension_(null_space_dimension)
                        , stoichiometric_matrix_T_(stoichiometric_matrix_T)
                        , big_M_value_(big_M_value)
                        , equilibrium_constants_(equilibrium_constants)
                        , variable_metabolites_upper_bound_(variable_metabolites_upper_bound)
                        , variable_metabolites_lower_bound_(variable_metabolites_lower_bound)
{
    assert(fixed_metabolites_.size() == n_metabolites_ - n_variable_metabolites_);
    assert(null_space_matrix_.rows() == n_reactions_);
    assert(null_space_matrix_.cols() == null_space_dimension_);
    assert(stoichiometric_matrix_T_.rows() == n_reactions_);
    assert(stoichiometric_matrix_T_.cols() == n_metabolites_);
    assert(equilibrium_constants_.size() == n_reactions_);
    assert(variable_metabolites_upper_bound_.minCoeff() > variable_metabolites_lower_bound_);
}

vector_t 
Constraints::GetValues() const
{
    /* 
    Equation references:
        MEPPF = Maximum Entropy Production Problem Formulation
        RMP   = Regulation Methods Paper
                ("An Approach to Maximize Growth and Entropy Production Rates in Metabolism")

    NOTE: If the numbers do not match, the constraints are 
        stll coded in the order that they appear in MEPPF.
    */

    vector_t result(GetRows()); // result vector to be returned

    std::cout << "\nhere get vals begin\n";

    result << CalculateNullSpaceRepresentationConstraint() // MEPPF 94
            , CalculateSteadyStateConstraint()             // MEPPF 95
            , CalculateSmoothConstraint()                  // MEPPF 96
            , CalculateRelaxedFluxUpperConstraint()        // MEPPF 97
            , CalculateRelaxedFluxLowerConstraint()        // MEPPF 98
            , CalculateSignConstraint()                    // MEPPF 99
            , CalculateRelaxedFluxSignConstraint()         // MEPPF 100
            , CalculateMetabolitesUpperBoundConstraint()   // MEPPF 101 upper
            , CalculateMetaboliteLowerBoundConstraint();   // MEPPF 101 lower

    std::cout << "\nhere get vals end\n";

    return result;
}

Constraints::VecBound 
Constraints::GetBounds() 
const
{
    std::cout << "\nhere get bounds begin\n";
    VecBound b(GetRows());

    assert(b.end() == (b.begin() + 7 * n_reactions_ + 2 * n_variable_metabolites_));

    // MEPPF 94 bounds
    auto next_begin = b.begin();
    auto next_end   = b.begin() + n_reactions_;
    std::fill( next_begin, next_end, Bounds(0.0, 0.0) ); 

    // MEPPF 95 bounds
    next_begin = b.begin() + n_reactions_;
    next_end   = b.begin() + 2 * n_reactions_;
    std::fill( next_begin, next_end, Bounds(0.0, 0.0) );

    // MEPPF 96 bounds
    next_begin = b.begin() + 2 * n_reactions_;
    next_end =   b.begin() + 3 * n_reactions_;
    std::fill( next_begin , next_end , Bounds(0.0, 0.0) ); 

    // MEPPF 97 bounds
    next_begin = b.begin() + 3 * n_reactions_;
    next_end   = b.begin() + 4 * n_reactions_;
    std::fill( next_begin , next_end , Bounds(-inf, 0.0) ); 

    // MEPPF 98 bounds
    next_begin = b.begin() + 4 * n_reactions_;
    next_end   = b.begin() + 5 * n_reactions_;
    std::fill( next_begin , next_end , Bounds(0.0, +inf) ); 

    // MEPPF 99 bounds
    next_begin = b.begin() + 5 * n_reactions_;
    next_end   = b.begin() + 6 * n_reactions_;
    std::fill( next_begin , next_end , Bounds(0.0, +inf) ); 

    // MEPPF 100 bounds
    next_begin = b.begin() + 6 * n_reactions_;
    next_end   = b.begin() + 7 * n_reactions_;
    std::fill( next_begin , next_end , Bounds(0.0, 0.0) ); 

    // MEPPF 101 upper bounds
    next_begin = b.begin() + 7 * n_reactions_;
    next_end   = b.begin() + 7 * n_reactions_ + n_variable_metabolites_;
    std::fill( next_begin , next_end , Bounds(0.0, +inf) ); 

    // MEPPF 101 lower bounds
    next_begin = b.begin() + 7 * n_reactions_ + n_variable_metabolites_;
    next_end   = b.begin() + 7 * n_reactions_ + 2 * n_variable_metabolites_;
    assert(next_end == b.end());
    std::fill( next_begin , next_end , Bounds(-inf, 0.0) ); 

    std::cout << "\nhere get bounds end\n";

    return b;
}

void 
Constraints::FillJacobianBlock( std::string var_set
                              , Jacobian& jac_block
                              ) 
const
{
    std::cout << "\nhere fill jac begin\n";
    int row;
    int col;
    int diag;

    const int start = 0;
    int end = n_reactions_;

    // index offests by constraint
    int meppf94 = 0;
    int meppf95 = 1 * n_reactions_;
    int meppf96 = 2 * n_reactions_;
    int meppf97 = 3 * n_reactions_;
    int meppf98 = 4 * n_reactions_;
    int meppf99 = 5 * n_reactions_;
    int meppf100 = 6 * n_reactions_;
    int meppf101_upper = 7 * n_reactions_;
    int meppf101_lower = 7 * n_reactions_ + n_variable_metabolites_;
    int row_offset; // generally equal to constraint_offset * n_metabolites_
    int col_offset; // increases by n_metabolities_ after each constraint jacobian fill block within each variable if-statement,
                    // since jac_block ionly contains columns relevent to particular variable

    vector_t gradient;

    if (var_set == variable_metabolites_name_) {
        // aka \eta in MEPPF
        end = n_variable_metabolites_;
        std::cout << "\nhere fill jac end\n";

        // MEPPF 95
        col_offset = 0;
        row_offset = meppf95;
        for (row = start; row < end; row ++){
            for (col = start; col < end; col++){
                jac_block.coeffRef( row + row_offset
                                  , col + col_offset
                                  ) 
                                  = stoichiometric_matrix_T_(row, col);
            }
        }

        // MEPPF 101 upper
        col_offset += end;
        row_offset = meppf101_upper;
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( row + row_offset
                              , col + col_offset
                              ) 
                              = -1.0;
        }

        // MEPPF 101 lower
        col_offset += end;
        row_offset = meppf101_lower;
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( row + row_offset
                              , col + col_offset
                              ) 
                              = -1.0;
        }

    } else if (var_set == beta_variables_name_) {
        // \beta in MEPPF

        // MEPPF 94
        col_offset = 0;
        row_offset = meppf94;
        for (row = start; row < end; row++) {
            for (col = 0; col < null_space_dimension_; col++) {
                jac_block.coeffRef( row + row_offset
                                  , col + col_offset
                                  ) 
                                  = null_space_matrix_(row,col);
            }
        }

    } else if (var_set == steady_state_variables_name_) {
        // aka g in MEPPF

        // MEPPF 95
        col_offset = 0;
        row_offset = meppf95;
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = -1.0;
        }

        // MEPPF 97
        col_offset += end;
        row_offset = meppf97;
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = -1.0;
        }

        // MEPPF 98
        col_offset += end;
        row_offset = meppf98;
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = -1.0;
        }

        // MEPPF 99
        col_offset += end;
        row_offset = meppf99;
        gradient = CalculateSignConstraintGradientSteadyStateVariables();
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = gradient(diag);
        }

    } else if (var_set == flux_variables_name_) {
        // aka y in MEPPF

        // MEPPF 94
        col_offset = 0;
        row_offset = meppf94;
        for (diag = start; diag < end ; row++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = -1.0;
        }

        // MEPPF 96
        col_offset += end;
        row_offset = meppf96;
        gradient = CalculateSmoothConstraintGradientFluxVariables().eval();
        for (diag = start; diag < end; diag ++) {
            jac_block.coeffRef( diag + row_offset
                              , diag + col_offset
                              )
                              = gradient(diag);
        }

        // MEPPF 99
        col_offset += end;
        row_offset = meppf99;
        gradient = CalculateSignConstraintGadientFluxVariables().eval();
        for (diag = start; diag < end; diag ++) {
            jac_block.coeffRef( diag + row_offset
                              , diag + col_offset
                              )
                              = gradient(diag);
        }

        // MEPPF 100
        // see MEPPF, gradient is 0

    } else if (var_set == h_variables_name_) { 
        // aka h in MEPPF

        // MEPPF 96
        col_offset = 0;
        row_offset = meppf96;
        for (diag = start; diag < end; diag ++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = -1.0;
        }

        // MEPPF 97
        col_offset += end;
        row_offset = meppf97;
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = 1.0;
        }

        // MEPPF 98
        col_offset += end;
        row_offset = meppf98;
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = 1.0;
        }
        
    } else if (var_set == u_variables_name_) {
        // aka u in MEPPF

        // MEPPF 97
        col_offset = 0;
        row_offset = meppf97;
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = -big_M_value_;
        }

        // MEPPF 98
        col_offset += end;
        row_offset = meppf98;
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = -big_M_value_;
        }

        // MEPPF 100
        col_offset += end;
        row_offset = meppf100;
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = 2.0;
        }
    }

}

vector_t 
Constraints::GetVariablesVectorByName( const std::string& variables_name )
const
{
    return GetVariables()->GetComponent(variables_name)->GetValues();
}

vector_t
Constraints::GetAllMetabolites()
const
{
    vector_t variable_metabolites = GetVariablesVectorByName(variable_metabolites_name_);
    vector_t all_metabolites( n_metabolites_ );
    all_metabolites << variable_metabolites, fixed_metabolites_;

    return all_metabolites;
}

vector_t
Constraints::CalculateNullSpaceRepresentationConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 94

    vector_t beta_vec = GetVariablesVectorByName(beta_variables_name_);
    vector_t fluxes = GetVariablesVectorByName(flux_variables_name_);

    vector_t result = ((null_space_matrix_ * beta_vec) - fluxes);

    return result;
}

vector_t
Constraints::CalculateSteadyStateConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 95

    vector_t steady_state_variables = GetVariablesVectorByName(steady_state_variables_name_);
    vector_t metabolites = GetAllMetabolites();

    vector_t result = ((stoichiometric_matrix_T_ * metabolites) - steady_state_variables);

    return result;
}

vector_t
Constraints::CalculateSmoothConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 96

    auto h_variables = GetVariablesVectorByName(h_variables_name_).array();
    auto fluxes_array = GetVariablesVectorByName(flux_variables_name_).array();
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

vector_t
Constraints::CalculateRelaxedFluxUpperConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 97
    
    vector_t h_variables = GetVariablesVectorByName(h_variables_name_);
    vector_t u_variables = GetVariablesVectorByName(u_variables_name_);
    vector_t steady_state_variables = GetVariablesVectorByName(steady_state_variables_name_);
    vector_t log_equilibrium_constants = equilibrium_constants_.array().log().matrix();

    vector_t result = h_variables - (big_M_value_ * u_variables) + log_equilibrium_constants - steady_state_variables;

    return result;
}

vector_t
Constraints::CalculateRelaxedFluxLowerConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 98
    
    vector_t h_variables = GetVariablesVectorByName(h_variables_name_);
    vector_t u_variables = GetVariablesVectorByName(u_variables_name_);
    vector_t one_minus_u = (1 - u_variables.array()).matrix();
    vector_t steady_state_variables = GetVariablesVectorByName(steady_state_variables_name_);
    vector_t log_equilibrium_constants = equilibrium_constants_.array().log().matrix();

    vector_t result = h_variables + (big_M_value_ * one_minus_u) + log_equilibrium_constants - steady_state_variables;

    return result;
}

vector_t
Constraints::CalculateSignConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 99

    auto fluxes_array = GetVariablesVectorByName(flux_variables_name_).array();
    auto log_equilibrium_constants_array = equilibrium_constants_.array().log();
    auto steady_state_variables_array = GetVariablesVectorByName(steady_state_variables_name_).array();

    auto result = (log_equilibrium_constants_array - steady_state_variables_array) *  fluxes_array;

    return result.matrix();
}

vector_t
Constraints::CalculateRelaxedFluxSignConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 100

    auto fluxes_array = GetVariablesVectorByName(flux_variables_name_).array();
    auto sign_fluxes = fluxes_array.sign();
    auto u_variables_array = GetVariablesVectorByName(u_variables_name_).array();

    auto result = (2 * u_variables_array) - 1 - sign_fluxes;

    return result.matrix();
}

vector_t
Constraints::CalculateMetabolitesUpperBoundConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 101 upper bounds

    auto metabolites_array = GetAllMetabolites().array();
    auto result = variable_metabolites_upper_bound_.array() - metabolites_array;
    return result.matrix();
}

vector_t
Constraints::CalculateMetaboliteLowerBoundConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 101 lower bounds

    auto metabolites_array = GetAllMetabolites().array();
    auto result = variable_metabolites_lower_bound_ - metabolites_array;

    return result.matrix();
}

vector_t
Constraints::CalculateSmoothConstraintGradientFluxVariables()
const 
{
    auto flux_array_squared = GetVariablesVectorByName(flux_variables_name_).array().pow(2);
    auto result = (flux_array_squared + 4).sqrt().inverse();

    return result.matrix();
}

vector_t
Constraints::CalculateSignConstraintGadientFluxVariables() 
const
{
    auto steady_state_variables = GetVariablesVectorByName(steady_state_variables_name_);
    vector_t result = equilibrium_constants_.array().log().matrix() - steady_state_variables;

    return result;
}

vector_t
Constraints::CalculateSignConstraintGradientSteadyStateVariables()
const
{
    vector_t flux_variables = GetVariablesVectorByName(flux_variables_name_);
    vector_t result = -1.0 * flux_variables;

    return result;
}

} // namespace