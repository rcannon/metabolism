
#include "includes_and_types.hh"

using namespace Eigen;

/*
These functions are adapted from the stack overflow answer:
https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix
*/

std::tuple<size_t,vector_t>
load_concentrations( const std::string& path)
{
    /* 
    This function will read in the concentration.csv file
    and produce a pair containing
        (1) the number of variable concentrations and
        (2) a vector of the concentrations.
    
    The file should be formatted as follows:
    write the number of variable concentrations on the first line,
    write the total number of concentrations on the second line,
    and on the third line write the values that will comprise the vector.

    Here is an example:

    ```concentrations.csv
    2
    3
    -1.0,2.0,-3.0
    ```
    */

    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<value_t> values;

    //vector_t res_vector;
    size_t num_variable;
    size_t num_total;
    size_t rows = 0;

    // read numvar
    std::getline(indata, line);
    num_variable = std::stoi(line);

    // read numtot
    std::getline(indata, line);
    num_total = std::stoi(line);

    // read concentrations vector
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }

    auto res_vector = Map<vector_t>(values.data(), values.size(), 1);

    return { num_variable, res_vector };
}

std::tuple<size_t,size_t, index_list_t, vector_t>
load_equilibrium_constants( const std::string& path)
{
    /* 
    This function will read in the concentration.csv file
    and produce a pair containing
        (1) the number of uptake reactions
        (2) the number of output reactions
        (3) an std::vector of the obj_idxs
        (4) a vector of the equilibrium constants
    
    The file should be formatted as follows:
    write the number of uptake reactions on the first line,
    write the number of output reactions on the second line,
    list the obj_idxs on the third line (separated by commas)
    and on the fourth line write the values that will comprise the vector.

    Here is an example:

    ```EquilibriumConstants.csv
    2
    3
    0,1
    -1.0,2.0,-3.0
    ```
    */

    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<value_t> values;

    vector_t res_vector;
    index_list_t obj_idxs;
    size_t num_uptake;
    size_t num_output;
    size_t rows = 0;

    // read num_uptake
    std::getline(indata, line);
    num_uptake = std::stoi(line);

    // read num_output
    std::getline(indata, line);
    num_output = std::stoi(line);

    // read obj_idxs
    std::getline(indata, line);
    std::stringstream lineStream(line);
    std::string cell;
    while (std::getline(lineStream, cell, ',')) {
        obj_idxs.push_back(std::stoi(cell));
    }

    // read equilibrium constants vector
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    res_vector = Map<vector_t>(values.data(), values.size(), 1);

    return { num_uptake, num_output, obj_idxs, res_vector };
}

matrix_t
load_stochiometric_matrix( const std::string& path)
{
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;

    int rows = 0;
    // read matrix
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    // Matrix<double, Dynamic, Dynamic, RowMajor>
    return Map<matrix_t>(values.data(), rows, values.size() / rows);
}