
#include "includes_and_types.hh"
#include "maximum_entropy_relaxed.hh"
#include "read_data_files.hh"
#include "metabolism.hh"

int
main()
{
    // read the concentrations data from the file
    auto concentrations_tuple = load_concentrations("../data/concentrations.csv");
    int num_variable_metabolites = std::get<0>(concentrations_tuple);
    vector_t all_metabolites = std::get<1>(concentrations_tuple);
    // Eigen::seq is upper inclusive
    vector_t variable_metabolites = all_metabolites(Eigen::seq(0, num_variable_metabolites-1));
    vector_t fixed_metabolites = all_metabolites(Eigen::seq(num_variable_metabolites, all_metabolites.size()-1));

    // read the equilibrium constant data from the file
    auto equlibrium_constants_tuple = load_equilibrium_constants("../data/EquilibriumConstants.csv");
    int num_uptake = std::get<0>(equlibrium_constants_tuple);
    int num_output = std::get<1>(equlibrium_constants_tuple);
    std::vector<int> objective_reaction_indices = std::get<2>(equlibrium_constants_tuple);
    vector_t equilibrium_constants = std::get<3>(equlibrium_constants_tuple);


    // read the stochiometric matrix from the file
    matrix_t stoichiometric_matrix = read_to_matrix("../data/StoichiometricMatrix.csv");

    int num_reactions = stoichiometric_matrix.rows();

    // get the target log counts for the variable metabolites
    value_t n_Avagadro = 6.022140857 * std::pow(10,23);
    value_t cell_volume = std::pow(10, -15);
    value_t concentration_to_count = n_Avagadro * cell_volume;
    value_t concentration_increment = std::pow(concentration_to_count, -1);
    double target_val = std::log( concentration_to_count * std::pow(10, -3) );
    vector_t target_log_variable_metabolites_count = vector_t::Constant(num_variable_metabolites, target_val);
    
    int Vmax = 1000;
    double s = 0.085;
    double Km = 0.5*s;
    int iuptake = 32;
    int ioutput = 32;
    auto run_results = run_metabolism
        ( variable_metabolites
        , fixed_metabolites
        , target_log_variable_metabolites_count
        , stoichiometric_matrix
        , equilibrium_constants
        , objective_reaction_indices
        );
    vector_t reaction_flux = std::get<0>(run_results);
    variable_metabolites = std::get<1>(run_results);

    double nutrient_ratio = 1;
    double s_new;
    vector_t scaled_flux2;
    for (int step = 0; step < 5; step++){
        std::cout << "\nStep " << step << "\n" << std::endl;
        s_new = s + 0.04 * s;
        nutrient_ratio = nutrient_ratio * (s_new/s);

        if ((nutrient_ratio > 1.05) || (nutrient_ratio < 0.95)){
            /*auto run_results = run_metabolism
                ( variable_metabolites
                , fixed_metabolites
                , target_log_variable_metabolites_count
                , stoichiometric_matrix
                , equilibrium_constants
                , objective_reaction_indices
                );
            vector_t reaction_flux = std::get<0>(run_results);
            variable_metabolites = std::get<1>(run_results);*/

            nutrient_ratio = 1;
        }
        scaled_flux2 = scale_flux(reaction_flux, iuptake, Vmax, Km, s);

        std::cout << "scale_flux2: " << scaled_flux2(iuptake) << "\n" <<std::endl;
    }
    return 0;
}