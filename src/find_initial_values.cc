
#include "find_initial_values.hh"

std::tuple<vector_t, vector_t, vector_t>
find_initial_values
    ( const vector_t& variable_metabolites
    , const vector_t& fixed_metabolites
    , const vector_t& target_log_counts
    , const matrix_t& stoichiometric_matrix
    , const vector_t& equilibrium_constants
    )
{
    //
    // Replace met with something better
    //

    int n_reactions = stoichiometric_matrix.rows();
    vector_t variable_metabolites_init = variable_metabolites;
    auto beta_init = vector_t::Random(n_reactions);
    auto flux_init = vector_t::Random(n_reactions);
    
    return { beta_init, flux_init, variable_metabolites_init };
}