
#include "ifopt_cost_class.hh"

namespace ifopt {


Cost::Cost  ( const std::string& name 
            , const std::string& flux_variables_name
            , const index_list_t& objective_reaction_indices
            ) 
            : CostTerm(name)
            , flux_variables_name_(flux_variables_name) 
            , objective_reaction_indices_(objective_reaction_indices)
{}


double 
Cost::GetCost() 
const
{
    vector_t x = GetVariables()->GetComponent(flux_variables_name_)->GetValues()(objective_reaction_indices_);
    return -x.sum();
}


void 
Cost::FillJacobianBlock ( std::string var_set
                        , Jacobian& jac
                        ) 
const
{
    /*std::cout << var_set << "\n";
    for (int row =0; row < 1000; row++){
        for (int col=0; col < 1000; col++){
            std::cout << "jac " << jac.coeffRef(row, col);
            std::cout << " row " << row << " col " << col << "\n";
        }
    }*/
    
    if (var_set == flux_variables_name_) {
        // vector_t x = GetVariables()->GetComponent(variable_metabolites_name_)->GetValues()(objective_reaction_indices);
        //int num_vars = GetVariables()->GetComponent(variable_metabolites_name_)->GetNumVariables();

        for (const auto idx : objective_reaction_indices_) {
            /*
            Derivative of cost wrt to each reaction variable i.
            Note: IFOPT will put these in correct location in overall jacobian matrix.
            */
            jac.coeffRef(0, idx) = -1.0;
        }
    }
}


} // namespace