
#include "read_data_files.hh"
#include "ifopt_constraint_classes.cc"
#include "find_initial_values.hh"

int main()
{
    // read the concentrations data from the file
    auto concentrations_tuple = load_concentrations("../data/concentrations.csv");
    int num_variable_metabolites = std::get<0>(concentrations_tuple);
    vector_t all_metabolites = std::get<1>(concentrations_tuple);
    vector_t temp_variable_metabolites = all_metabolites(Eigen::seq(0, num_variable_metabolites-1));
    //
    vector_t fixed_metabolites = all_metabolites(Eigen::seq(num_variable_metabolites, all_metabolites.size()-1));

    // read the equilibrium constant data from the file
    auto equlibrium_constants_tuple = load_equilibrium_constants("../data/EquilibriumConstants.csv");
    int num_uptake = std::get<0>(equlibrium_constants_tuple);
    int num_output = std::get<1>(equlibrium_constants_tuple);
    std::vector<int> objective_reaction_indices = std::get<2>(equlibrium_constants_tuple);
    //
    vector_t equilibrium_constants = std::get<3>(equlibrium_constants_tuple);

    // read the stochiometric matrix from the file
    matrix_t stoichiometric_matrix = read_to_matrix("../data/StoichiometricMatrix.csv");

    // get the target log counts for the variable metabolites
    value_t n_Avagadro = 6.022140857 * std::pow(10,23);
    value_t cell_volume = std::pow(10, -15);
    value_t concentration_to_count = n_Avagadro * cell_volume;
    value_t concentration_increment = std::pow(concentration_to_count, -1);
    double target_val = std::log( concentration_to_count * std::pow(10, -3) );
    //
    vector_t target_log_variable_metabolites_count = vector_t::Constant(num_variable_metabolites, target_val);

    // read in initial values
    auto initial_values = find_initial_values
        ( temp_variable_metabolites
        , target_log_counts
        , fixed_metabolites_log_counts
        , stoichiometric_matrix
        , equilibrium_constants
        );

    vector_t beta_init = std::get<0>(initial_values);
    vector_t flux_init = std::get<1>(initial_values);
    vector_t variable_metabolites_init = std::get<2>(initial_values);

    


}