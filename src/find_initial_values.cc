
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
    // Gets initial values from the csv files
    // located in data/python_feasible_point.
    // These files are generated from running
    // the python code:
    //
    // $ python3 metabolism_driver.py
    //

    vector_t variable_metabolites_init = read_to_vector("../data/python_feasible_point/variable_metabolites_log_counts.csv");
    vector_t beta_init = read_to_vector("../data/python_feasible_point/null_space_variables.csv");
    vector_t flux_init = read_to_vector("../data/python_feasible_point/flux_variables.csv");

    assert(variable_metabolites_init.size() == 182);
    assert(flux_init.size() == 198);
    assert(beta_init.size() == 16);

    assert(variable_metabolites_init.size() <= stoichiometric_matrix.cols());
    assert(flux_init.size() == stoichiometric_matrix.rows()); // n_reaction
    
    return { beta_init, flux_init, variable_metabolites_init };
}

vector_t
read_to_vector(const std::string path)
{
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<value_t> values;

    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
    }
    vector_t res = Eigen::Map<vector_t>(values.data(), values.size(), 1);

    return res;
}
