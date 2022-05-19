
#include "metabolism.hh"

std::tuple<vector_t, vector_t, vector_t>
initialize 
    ( const vector_t& variable_metabolites
    , const vector_t& fixed_metabolites
    , const vector_t& target_log_counts
    , const matrix_t& stoichiometric_matrix
    , const vector_t& equilibrium_constants
    )
{
    //
    // modify me for better initialization
    //
    auto initial_values = find_initial_values
        ( variable_metabolites
        , target_log_counts
        , fixed_metabolites 
        , stoichiometric_matrix
        , equilibrium_constants
        );

    vector_t beta_init = std::get<0>(initial_values);
    vector_t flux_init = std::get<1>(initial_values);
    vector_t variable_metabolites_init = std::get<2>(initial_values);
    return { beta_init, flux_init, variable_metabolites_init };
}

std::tuple<vector_t, vector_t, vector_t, vector_t>
optimize
    ( const vector_t& variable_metabolites_init
    , const vector_t& flux_variables_init
    , const vector_t& beta_init
    , const vector_t& fixed_metabolite_log_counts
    , const vector_t& target_log_variable_metabolites_counts
    , const matrix_t& stoichiometric_matrix
    , const vector_t& equilibrium_constants
    , const index_list_t& objective_reaction_indices
    )
{
    auto results_tuple = maximum_entropy_solver 
        ( variable_metabolites_init
        , fixed_metabolite_log_counts
        , flux_variables_init
        , beta_init
        , target_log_variable_metabolites_counts
        , stoichiometric_matrix
        , equilibrium_constants
        , objective_reaction_indices
        );
    vector_t steady_state_sol = std::get<0>(results_tuple);
    vector_t flux_sol = std::get<1>(results_tuple);
    vector_t alpha_sol = std::get<2>(results_tuple);
    vector_t h_sol = std::get<3>(results_tuple);
    vector_t beta_sol = std::get<4>(results_tuple);
    vector_t variable_metabolite_sol = std::get<5>(results_tuple);

    vector_t e_regulation = alpha_sol;
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
    int num_total_metabolites = variable_metabolites_counts.size() + fixed_metabolites_counts.size();
    vector_t all_metabolites(num_total_metabolites) << variable_metabolites_counts, fixed_metabolites_counts;

    auto delta = vector_t::Zero(num_total_metabolites);

    double mu0 = 1.0;
    reaction_flux = odds_diff
        ( variable_metabolite_sol
        , fixed_metabolite_log_counts
        , mu0
        , stoichiometric_matrix
        , negative_stoich_matrix
        , positive_stoich_matrix
        , delta
        , equilibrium_constants
        , e_regulation
        );

    return { steady_state_sol, flux_sol, variable_metabolite_sol, alpha_sol};
}

vector_t
odds_diff
    ( const vector_t& variable_metabolites_counts
    , const vector_t& fixed_metabolites_counts
    , const double mu0
    , const matrix_t& stoichiometric_matrix
    , const matrix_t& negative_stoich_matrix
    , const matrix_t& positive_stoich_matrix
    , const vector_t& delta
    , const vector_t& equilibrium_constants
    , const vector_t& e_regulation
    )
{
    int num_total_metabolites = variable_metabolites_counts.size() + fixed_metabolites_counts.size();
    vector_t all_metabolites(num_total_metabolites) << variable_metabolites_counts, fixed_metabolites_counts;
    vector_t KQ_f = calc_odds
        ( all_metabolites
        , negative_stoich_matrix
        , positive_stoich_matrix
        , delta
        , equilibrium_constants
        );
    vector_t equilibrium_constants_inverse = equilibrium_constants_inverse.array().inverse().matrix();
    double direction = -1.0;
    vector_t KQ_r = calc_odds
        ( all_metabolites
        , negative_stoich_matrix
        , positive_stoich_matrix
        , delta
        , equilibrium_constants_inverse
        , direction
        );
    vector_t KQ_diff = ( e_regulation.array() * (KQ_f - KQ_r).array() ).matrix();
    
    return KQ_diff;
}

vector_t
calc_odds
    ( const vector_t& log_counts
    , const matrix_t& negative_stoich_matrix
    , const matrix_t& positive_stoich_matrix
    , const vector_t& delta
    , const vector_t& equilibrium_constants
    , const double direction = 1.0
    )
{
    vector_t counts = log_counts.array().exp().matrix();
    vector_t delta_counts = counts + delta;
    vector_t log_delta = delta_counts.array().log().matrix();
    vector_t Q_inv = 
        ( -direction * (R * log_counts)
        + (P * log_delta)
        ).array().exp().matrix();
    vector_t KQ = (equilibrium_constants.array() * Q_inv.array()).matrix();
    return KQ;
}

std::tuple<vector_t, vector_t>
run ( const vector_t& variable_metabolites
    , const vector_t& fixed_metabolites
    , const vector_t& target_log_variable_metabolites_count
    , const matrix_t& stoichiometric_matrix
    , const vector_t& equilibrium_constants
    , const int i_uptake
    , const int Vmax
    , const double Km
    , const double s
    , const index_list_t& objective_reaction_indices
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
    auto delta = vector_t::Zero(num_total_metabolites);
    double mu0 = 1.0
    vector_t reaction_flux = odds_diff
        ( variable_metabolites
        , fixed_metabolites
        , mu0
        , stoichiometric_matrix
        , negative_stoich_matrix
        , positive_stoich_matrix
        , delta
        , equilibrium_constants
        , e_regulation
        ;)

    return { reaction_flux, variable_metabolites_sol };
}