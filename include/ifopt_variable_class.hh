
#include "includes_and_types.hh"

#pragma once

namespace ifopt {

/* generic variable class */
class Variables : public VariableSet {
public:

    Variables( const std::string& name
             , const int num_variables
             , const bool is_u_variables
             );
    Variables( const std::string& name
             , const int num_variables
             , const bool is_u_variables
             , const vector_t& init_values
             );
    void SetVariables( const vector_t& new_vars ) override;
    vector_t GetValues() const override;
    VecBound GetBounds() const override;
    int GetNumVariables() const;

private:
    const int num_variables_;
    const bool is_u_variables_;
    vector_t variables_;

};

} // namespace
