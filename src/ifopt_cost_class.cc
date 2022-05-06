
#include "ifopt_cost_class.hh"

namespace ifopt {


Cost::Cost  ( const std::string& name 
            , const std::string& variable_metabolites_name 
            ) 
            : CostTerm(name)
            , variable_metabolites_name_(variable_metabolites_name) 
{}


double 
Cost::GetCost() 
const override
{
    vector_t x = GetVariables()->GetComponent(variable_metabolites_name_)->GetValues();
    return -x.sum();
}


void 
Cost::FillJacobianBlock ( std::string var_set
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


} // namespace