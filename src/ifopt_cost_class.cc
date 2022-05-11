
#include "ifopt_cost_class.hh"

namespace ifopt {


Cost::Cost  ( const std::string& name 
            , const std::string& variable_metabolites_name 
            , const index_list_t& objective_reaction_indices
            ) 
            : CostTerm(name)
            , variable_metabolites_name_(variable_metabolites_name) 
            , objective_reaction_indices_(objective_reaction_indices)
{}


double 
Cost::GetCost() 
const override
{
    vector_t x = GetVariables()->GetComponent(variable_metabolites_name_)->GetValues()(objective_reaction_indices);
    return -x.sum();
}


void 
Cost::FillJacobianBlock ( std::string var_set
                        , Jacobian& jac
                        ) 
const override
{
    if (var_set == variable_metabolites_name_) {
        // vector_t x = GetVariables()->GetComponent(variable_metabolites_name_)->GetValues()(objective_reaction_indices);
        //int num_vars = GetVariables()->GetComponent(variable_metabolites_name_)->GetNumVariables();

        for (int& idx : objective_reaction_indices_) {
            /*
            Derivative of cost wrt to each reaction variable i.
            Note: IFOPT will put these in correct location in overall jacobian matrix.
            */
            jac.coeffRef(0, idx) = -1.0;
        }
}


} // namespace