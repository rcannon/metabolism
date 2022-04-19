
#include "classes_maximum_entropy.hh"

namespace ifopt {

/*
VARIABLES CLASS
*/

Variables::Variables( const std::string& name
                    , const int num_variables
                    )
                    : VariableSet(num_variables, name)
                    , num_variables_(num_variables)
                    , variables_(vector_t::Zero(nun_variables))
{}

Variables::Variables( const std::string& name
                    , const int num_variables
                    , const vector_t& init_values
                    )
                    : VariableSet(num_variables, name)
                    , num_variables_(num_variables)
                    , variables_(init_values)
{}

void 
Variables::SetVariables( const vector_t& new_vars ) 
override
{
    assert(n_variables_ == new_vars.size());
    variables_ = new_vars;
} 

vector_t 
Variables::GetValues() 
const override
{
    return variables_;
}

int 
Variables::GetNumVariables() 
const
{
    return n_variables_;
}

/*
END VARIABLES CLASS
*/


/*
CONSTRAINTS CLASS
*/

Constraints::Constraints(const ConstraintArguements& args) 
                        : ConstraintSet( 7 * args.n_metabolites + 2 * args.n_variable_metabolites 
                                       , args.name)
                        , n_metabolites_(args.n_metabolites)
                        , n_variable_metabolites_(args.n_variable_metabolites)
                        , fixed_metabolites_(args.fixed_metabolites_)
                        , variable_metabolites_name_(args.variable_metabolites_name)
                        , flux_variables_name_(args.flux_variables_name)
                        , steady_state_variables_name_(args.steady_state_variables_name)
                        , h_variables_name_(args.h_variables_name)
                        , u_variables_name_(args.u_variable_name)
                        , beta_variables_name_(args.beta_variables_name)
                        , null_space_matrix_(args.null_space_matrix)
                        , null_space_dimension_(args.null_space_dimension)
                        , stochiometric_matrix_(args.stochiometric_matrix)
                        , big_M_value_(args.big_M)
                        , K_vector_(args.K_vector)
                        , variable_metabolites_upper_bound_(args.variable_metabolites_upper_bound)
                        , variable_metabolites_lower_bound_(args.variable_metabolites_lower_bound)
{
    assert(n_variable_metabolites_ == GetVariablesVectorByName(variable_metabolites_name_).size());
    assert(fixed_metabolites_.size() == n_metabolites_ - n_variable_metabolites_);
    assert(GetVariablesVectorByName(flux_variables_name_).size() == n_metabolites_);
    assert(GetVariablesVectorByName(steady_state_variables_name_).size() == n_metabolites_);
    assert(GetVariablesVectorByName(h_variables_name_).size() == n_metabolites_);
    assert(GetVariablesVectorByName(u_variables_name_).size() == n_metabolites_);
    assert(null_space_matrix_.rows() == n_metabolites_);
    assert(null_space_dimension_ == null_space_matrix_.cols());
    assert(null_space_dimension_ == GetVariablesVectorByName(beta_variables_name_).size());
    assert(stochiometric_matrix_.rows() == n_metabolites_);
    assert(stochiometric_matrix_.cols() == n_metabolites_);
    assert(K_vector_.size() == n_metabolites_);
    assert(variable_metabolites_upper_bound_ > variable_metabolites_lower_bound_);
}

vector_t 
Constraints::GetValues() const override
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

    result << CalculateNullSpaceRepresentationConstraint() // MEPPF 94
            , CalculateSteadyStateConstraint()             // MEPPF 95
            , CalculateSmoothConstraint()                  // MEPPF 96
            , CalculateRelaxedFluxUpperConstraint()        // MEPPF 97
            , CalculateRelexedFluxLowerConstraint()        // MEPPF 98
            , CalculateSignConstraint()                    // MEPPF 99
            , CalculateRelaxedFluxSignConstraint()         // MEPPF 100
            , CalculateMetabolitesUpperBoundConstraint()   // MEPPF 101
            , CalculateMetaboliteLowerBoundConstraint()    // MEPPF 101 lower

    return result;
}

VecBound 
GetBounds() 
const override 
{
    VecBound b(GetRows());

    assert(b.end() == (b.begin() + 7 * n_metabolites_ + 2 * n_variable_metabolites_));

    std::fill( b.begin()                      , b.begin() +     n_metabolites_ , Bounds(0.0, 0.0)  ); // MEPPF 94 bounds
    std::fill( b.begin() +     n_metabolites_ , b.begin() + 2 * n_metabolites_ , Bounds(0.0, 0.0)  ); // MEPPF 95 bounds
    std::fill( b.begin() + 2 * n_metabolites_ , b.begin() + 3 * n_metabolites_ , Bounds(0.0, 0.0)  ); // MEPPF 96 bounds
    std::fill( b.begin() + 3 * n_metabolites_ , b.begin() + 4 * n_metabolites_ , Bounds(-inf, 0.0) ); // MEPPF 97 bounds
    std::fill( b.begin() + 4 * n_metabolites_ , b.begin() + 5 * n_metabolites_ , Bounds(0.0, +inf) ); // MEPPF 98 bounds
    std::fill( b.begin() + 5 * n_metabolites_ , b.begin() + 6 * n_metabolites_ , Bounds(0.0, +inf) ); // MEPPF 99 bounds
    std::fill( b.begin() + 6 * n_metabolites_ , b.begin() + 7 * n_metabolites_ , Bounds(0.0, 0.0)  ); // MEPPF 100 bounds
    std::fill( b.begin() + 7 * n_metabolites_ , b.begin() + 7 * n_metabolites_ 
                                                          + n_variable_metabolites_ , Bounds(0.0, +inf) ); // MEPPF 101 upper bounds
    std::fill( b.begin() + 7 * n_metabolites_ + n_variable_metabolites_ , b.end()   , Bounds(-inf, 0.0) ); // MEPPF 101 lower bounds

    return b;
}

void 
Constraints::FillJacobianBlock( std::string var_set
                              , Jacobian& jac_block
                              ) 
const override 
{

    int row;
    int col;
    int diag;

    const int start = 0;
    int end = n_metabolites_;

    // index offests by constraint
    int meppf94 = 0;
    int meppf95 = 1 * n_metabolites_;
    int meppf96 = 2 * n_metabolites_;
    int meppf97 = 3 * n_metabolites_;
    int meppf98 = 4 * n_metabolites_;
    int meppf99 = 5 * n_metabolites_;
    int meppf100 = 6 * n_metabolites_;
    int meppf101 = 7 * n_metabolites_;
    int meppf102 = 8 * n_metabolites_;
    int row_offset; // generally equal to constraint_offset * n_metabolites_
    int col_offset; // increases by n_metabolities_ after each constraint jacobian fill block within each variable if-statement,
                    // since jac_block ionly contains columns relevent to particular variable

    vector_t gradient;

    if (var_set == variable_metabolites_name_) {
        // aka \eta in MEPPF
        end = n_variable_metabolites_;

        // MEPPF 95
        col_offest = 0;
        row_offset = meppf95;
        for (row = start; row < end; row ++){
            for (col = start; col < end; col++){
                jac_block.coeffRef( row + row_offset
                                  , col + col_offset
                                  ) 
                                  = stochiometric_matrix(row, col);
            }
        }

        // MEPPF 101
        col_offest = += end;
        row_offset = meppf101;
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( row + row_offset
                              , col + col_offset
                              ) 
                              = -1.0;
            }
        }

        // MEPPF 102
        col_offset += end;
        row_offset = meppf102;
        for (diag = start; diag < end; diag++) {
            jac_block.coeffRef( row + row_offset
                              , col + col_offset
                              ) 
                              = -1.0;
            }
        }

    } else if (var_set == beta_variables_name_) {
        // \beta in MEPPF

        // MEPPF 94
        col_offest = 0;
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
        col_offest = 0;
        row_offset = meppf95;
        for (diag = start; diag < diag_end; diag++) {
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
        for (diag = start; diag < diag_end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = gradient(diag);
        }

    } else if (var_set == flux_variables_name_) {
        // aka y in MEPPF

        // MEPPF 94
        col_offest = 0;
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
        for (diag = start; diag < diag_end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = 1.0;
        }

        // MEPPF 98
        col_offset += end;
        row_offset = meppf98;
        for (diag = start; diag < diag_end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = 1.0;
        }
        
    } else if (var_set == u_variables_name_) {
        // aka u in MEPPF

        // MEPPF 97
        col_offest = 0;
        row_offset = meppf97;
        for (diag = start; diag < diag_end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = -big_M_value_;
        }

        // MEPPF 98
        col_offset += end;
        row_offset = meppf98;
        for (diag = start; diag < diag_end; diag++) {
            jac_block.coeffRef( diag + row_offset 
                              , diag + col_offset
                              ) 
                              = -big_M_value_;
        }

        // MEPPF 100
        col_offset += end;
        row_offset = meppf100;
        for (diag = start; diag < diag_end; diag++) {
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

    return ((null_space_matrix_ * beta_vec) - fluxes);
}

vector_t
Constraints::CalculateSteadyStateConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 95

    vector_t steady_state_variables = GetVariablesVectorByName(steady_state_variables_name_);
    vector_t metabolites = GetAllMetabolites();

    return ((stochiomtric_matrix_ * metabolites) - steady_state_variables);
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
Constriaint::CalculateRelaxedFluxUpperConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 97
    
    vector_t h_variables = GetVariablesVectorByName(h_variables_name_);
    vector_t u_variables = GetVariablesVectorByName(u_variables_name_);
    vector_t steady_state_variables = GetVariablesVectorByName(steady_state_variables_name_);
    vector_t log_K_vector = K_vector_.array().log().matrix();

    return h_variables - (big_M_value_ * u_variables) + log_K_vector - steady_state_variables;
}

vector_t
Constriaint::CalculateRelaxedFluxLowerConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 98
    
    vector_t h_variables = GetVariablesVectorByName(h_variables_name_);
    vector_t u_variables = GetVariablesVectorByName(u_variables_name_);
    vector_t one_minus_u = (1 - u_variables.array()).matrix()
    vector_t steady_state_variables = GetVariablesVectorByName(steady_state_variables_name_);
    vector_t log_K_vector = K_vector_.array().log().matrix();

    return h_variables + (big_M_value_ * one_minus_u) + log_K_vector - steady_state_variables;
}

vector_t
Constraints::CalculateSignConstraint()
const
{
    // Maximum Entropy Production Problem Formulation 99

    auto fluxes_array = GetVariablesVectorByName(flux_variables_name_).array();
    auto log_K_array = K_vector_.array().log();
    auto steady_state_variables_array = GetVariablesVectorByName(steady_state_variables_name_).array();

    auto result = (log_K_array - steady_state_variables_array) *  fluxes_array;

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
    auto result = variable_metabolites_upper_bound_ - metabolites_array;
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
    vector_t result = K_vector_.array().log().matrix() - steady_state_variables;

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

/*
END CONSTRAINTS CLASS
*/

/*
COST FUNCTION CLASS
*/

MetabConst::MetabCost( const std::string& name 
                     , const std::string& variable_metabolites_name 
                     ) 
                     : CostTerm(name)
                     , variable_metabolites_name_(variable_metabolites_name) 
{}


double 
MetabCost::GetCost() 
const override
{
    vector_t x = GetVariables()->GetComponent(variable_metabolites_name_)->GetValues();
    return -x.sum();
}


void 
MetabCost::FillJacobianBlock( std::string var_set
                            , Jacobian& jac
                            ) 
const override
{
    if (var_set == variable_metabolites_name_) {
        vector_t x = GetVariables()->GetComponent(variable_metabolites_name_)->GetValues();
        int num_vars = GetVariables()->GetComponent(variable_metabolites_name_)->GetNumVariables();

        for (int i = 0; i < num_vars; i++) {
            /*
            Derivative of cost wrt to react variable i.
            Note: IFOPT will put these in correct location in overall jacobian matrix.
            */
            jac.coeffRef(0, i) = -1.0;
        }
}

/*
END COST FUNCTION CLASS
*/

} // end ifopt namespace