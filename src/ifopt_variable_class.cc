
#include "ifopt_variable_class.hh"

namespace ifopt {

Variables::Variables( const std::string& name
                    , const int num_variables
                    )
                    : VariableSet(num_variables, name)
                    , num_variables_(num_variables)
                    , variables_(vector_t::Zero(num_variables))
{}

Variables::Variables( const std::string& name
                    , const int num_variables
                    , const vector_t& init_values
                    )
                    : VariableSet(num_variables, name)
                    , num_variables_(num_variables)
                    , variables_(init_values)
{
    assert(num_variables_ == variables_.size());
}

void 
Variables::SetVariables( const vector_t& new_vars ) 
{
    assert(num_variables_ == new_vars.size());
    variables_ = new_vars;
} 

vector_t 
Variables::GetValues() 
const
{
    return variables_;
}

int 
Variables::GetNumVariables() 
const
{
    return num_variables_;
}

Variables::VecBound
Variables::GetBounds() // metabolite bounds moved to constraints class
const
{
    VecBound bounds(GetRows());
    for (int i = 0; i < num_variables_; i++){
        bounds.at(i) = NoBound; 
    }
    return bounds;
}

} // ifopt

