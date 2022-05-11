
#include "includes_and_types.hh"

#pragma once 

namespace ifopt {

/* specific cost function class */
class Cost : public CostTerm {
public:

    Cost( const std::string& name 
        , const std::string& reaction_vars_name 
        , const index_list_t& objective_reaction_indices
        );
    double GetCost() const override;
    void FillJacobianBlock( std::string var_set
                          , Jacobian& jac
                          ) const override;
private:
    const std::string reaction_vars_name_;
    index_list_t objective_reaction_indices_;
};

} // namespace