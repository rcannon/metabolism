
#include "metabolism.hh"

std::tuple<vector_t, vector_t>
run ( vector_t variable_metabolites
    , vector_t fixed_metabolites
    , vector_t target_log_variable_metabolites_count
    , matrix_t stoichiometric_matrix
    , vector_t equilibrium_constants
    , int i_uptake
    , int Vmax
    , double Km
    , double s
    , index_list_t objective_reaction_indices
    )
{
    // get initial values
    auto initial_results = initialize
        ( variable_metabolites
        , fixed_metabolites
        , target_log_variable_metabolites_count
        , stoichiometric_matrix
        , equilibrium_constants
        );
    vector_t beta_init = std::get<0>(initial_results);
    vector_t flux_values_init = std::get<1>(initial_results);
    vector_t variable_metabolites_init = std::get<1>(initial_results);

    // get solution
    auto optimize_results = optimize
        ( variable_metabolites_init
        , flux_values_init
        , beta_init
        , fixed_metabolites
        , target_log_variable_metabolites_count
        , stoichiometric_matrix
        , equilibrium_constants
        , objective_reaction_indices
        );
    vector_t steady_state_sol = std::get<0>(optimize_results);
    vector_t flux_sol = std::get<1>(optimize_results);
    vector_t variable_metabolites_sol = std::get<2>(optimize_results);
    vector_t e_regulation = std::get<3>(optimize_results);


    // run odds difference
    matrix_t positive_stoich_matrix = stoichiometric_matrix;
    for ( row : positive_stoich_matrix.rowwise()) {
        for (auto elem : row) {
            if (elem < 0.0) {
                elem = 0.0;
            }
        }
    }
    matrix_t negative_stoich_matrix = stoichiometric_matrix;
    for ( row : negative_stoich_matrix.rowwise()) {
        for (auto elem : row) {
            if (elem > 0.0) {
                elem = 0.0;
            }
        }
    }
    int num_total_metabolites = variable_metabolites.size() + fixed_metabolites.size();
    auto delta = vector_t::Constant(num_total_metabolites, 1.0);
    vector_t reaction_flux = odds_difference
        ( variable_metabolites
        , fixed_metabolites
        , 1.0
        , stoichiometric_matrix
        , negative_stoich_matrix
        , positive_stoich_matrix
        , delta
        , equilibrium_constants
        , e_regulation
        ;)

    return { reaction_flux, variable_metabolites_sol };
}