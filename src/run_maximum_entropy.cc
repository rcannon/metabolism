
#include "run_maximum_entropy.hh"

void
run_maximum_entropy ()
{
    // read the concentrations data from the file
    auto concentrations_tuple = load_concentrations('./data/concentrations.csv');
    num_variable_metabolites = std::get<0>(concentrations_tuple);
    vector_t all_metabolites = std::get<1>(concentrations_tuple);
    vector_t variable_metabolites = all_metabolites(Eigen::seq(0, num_variable_metabolites));
    vector_t fixed_metabolites = all_metabolites(Eigen::seg(num_variable_metabolites, all_metabolites.size()));

    // read the equilibrium constant data from the file
    auto equlibrium_constants_tuple = load_equilibrium_constants("./data/EquilibriumConstants.csv");
    int num_uptake = std::get<0>(equlibrium_constants_tuple);
    int num_output = std::get<1>(equlibrium_constants_tuple);
    std::vector<int> objective_reaction_indices = std::get<2>(equlibrium_constants_tuple);
    vector_t equilibrium_constants = std::get<3>(equlibrium_constants_tuple);

    // read the stochiometric matrix from the file
    matrix_t stochiometric_matrix = load_stochiometric_matrix("./data/StochiometricMatrix.csv");

    int num_reations = stochiometric_matrix.rows();

    // get the target log counts for the variable metabolites
    value_t n_Avagadro = std::pow(6.022140857, 23);
    value_t cell_volume = std::pow(10, -15);
    value_t concentration_to_count = n_Avagadro * coll_volume;
    value_t concentration_increment = std::pow(concentration_to_count, -1);
    vector_t target_log_variable_metabolites_count = 
        ( vector_t::Constant(num_reactions, 1.0) 
        * std::pow(10, -3)
        * concentration_to_count
        ).array().log().matrix();

    // set the initial flux_variables
    flux_variables = vector_t::Constant(num_reactions, 1.0);

    // set the initial beta_variables
    flux_variables = vector_t::Constant(num_reactions, 1.0);

    // run maximum_entropy solver
    auto results_tuple = maximum_entropy_solver ( variable_metabolites
                                                , fixed_metabolites 
                                                , flux_variables
                                                , beta_variables
                                                , target_log_variable_metabolites_count
                                                , stoichiometric_matrix
                                                , equilibrium_constants
                                                , objective_reaction_indices
                                                );

    vector_t steady_state_sol = std::get<0>(results_tuple);
    vector_t flux_sol = std::get<1>(results_tuple);
    vector_t alpha_sol = std::get<2>(results_tuple);
    vector_t beta_sol = std::get<3>(results_tuple);
    vector_t metabolite_sol = std::get<4>(results_tuple);

    std::cout << "optimsation done\n";
}