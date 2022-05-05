
#include "maximum_entropy_relaxed.hh"


std::tuple<vector_t, vector_t, vector_t, vector_t, vector_t>
max_ent_solver
    ( vector_t variable_metabolites_ini // aka n
    , vector_t fixed_metabolites
    , vector_t flux_vars_ini // aka y
    , vector_t beta_ini // \beta
    , vector_t target_log_variable_metabolites_counts
    , matrix_t stoichiometric_matrix
    , vector_t equilibrium_constants // aka K
    , index_list_t obj_rxn_idx
    )
{

    TODO: 
    If n_react <= n_total_metabolites:
        trim stoich matrix, flux, steady_state, h_vars to only include ones less than n_react
    where n_react in line 37
    
    //
    // setup 
    //
    
    size_t num_fixed_metabolites = fixed_metabolites.size();
    size_t num_variable_metabolites = variable_metabolites_ini.size();
    size_t num_total_metabolies = num_fixed_metabolites + num_variable_metabolites;

    value_t variable_metabolite_lower_bound = -300;
    value_t big_M_value = 1000;
    vector_t flux_variables_ini_signs = flux_variables_ini.array().sign().matrix();

    // Flip Stoichiometric Matrix
    matrix_t stoichiometric_matrix_T = stoichiometric_matrix; // this should do deep copy 
    stoichiometric_matrix.transposeInPlace();
    size_t n_react = stoichiometric_matrix.cols(); // should be same as numpy.shape(S)[1]

    // Split S into the component corresponding to the variable metabolites S_v
    // and the component corresponding to the fixed metabolites. 
        // Reilly: below should do the same thing as numpy.delete
        // https://eigen.tuxfamily.org/dox-devel/group__TutorialSlicingIndexing.html
    index_list_t variable_metabolite_indices(num_variable_metabolites); /* variable metabolite indices */
    std::iota(variable_metabolite_indexes.begin(), variable_metabolite_indexes.end(), 0);
    index_list_t fixed_metabolite_indices(num_fixed_metabolites); /* fixed metabolite indices */
    std::iota(fixed_metabolite_indices.begin(), fixed_metabolite_indices.end(), num_fixed_metabolites);

    matrix_t stoich_matrix_variable_metab_section_T = stoichiometric_matrix_T(Eigen::all, variable_metabolite_indices);
    matrix_t stoich_matrix_fixed_metab_section_T = stoichiometric_matrix_T(Eigen::all, fixed_metabolite_indices);

    matrix_t stoich_matrix_variable_metab_section = stoich_matrix_variable_metab_section_T.transpose();
    matrix_t stoich_matrix_fixed_metab_section = stoich_matrix_fixed_metab_section_T.transpose();

    // find a basis for the nullspace of S_v,
    // and get the dimension of the null space 
    Eigen::FullPivLU<Eigen::MatrixXf> lu_decomp(stoich_matrix_variable_metab_section);
    matrix_t null_space_stoich_matrix_variable_metab_section = lu_decomp.kernel();
    size_t dim_null_space = lu_decomp.dimensionOfKernel();

    // precompute SVS
    matrix_t svs_matrix =   ( stoich_matrix_variable_metab_section 
                            * stoich_matrix_variable_metab_section_T
                            ).inverse() 
                            * stoich_matrix_variable_metab_section;

    //
    // initialize the ifopt problem
    //

    ifopt::Problem nlp;

    // add variables classes

    // add metabolite variables
    std::string variable_metabolites_variables_name = "log_variable_metabolites";
    nlp.AddVariableSet(std::make_shared<Variables>  ( variable_metabolites_variable_name
                                                    , num_variable_metabolites
                                                    , variable_metabolites_ini
                                                    ));

    // add flux variables (aka y)
    std::string flux_variables_name = "flux_variables";
    size_t num_flux_variables = flux_vars_ini.size();
    nlp.AddVariableSet(std::make_shared<Variables>  ( flux_variables_name
                                                    , num_flux_variables
                                                    , flux_vars_ini
                                                    ));

    // add steady state variables (aka g in paper, b in python version)
    std::string steady_state_variables_name = "steady_state_variables";
    vector_t steady_state_vars_ini  
        = (S_v_T * variable_metabolites_ini) 
        + (S_f_T * fixed_metabolites );
    size_t num_steady_state_variables = steady_state_variables_ini.size();
    nlp.AddVariableSet(std::make_shared<Variables>  ( steady_state_variables_name
                                                    , num_steady_state_variables
                                                    , steady_state_vars_ini
                                                    ));

    // add beta variables
    std::string beta_variables_name = "beta_variables";
    nlp.AddVariableSet(std::make_shared<Variables>  ( beta_variables_name
                                                    , dSv_N
                                                    , beta_ini
                                                    ));

    // add h variables
    std::string h_variables_name = "h_variables";
    auto fvia = flux_vars_ini.array();
    vector_t h_ini = 
        (   flux_variables_ini_signs.array()
        *   (   std::log(2) 
            -   Eigen::log  (   fvia.abs()
                            +   (fvia.pow(2) + 4).sqrt()   
                            )
            )
        ).matrix();
    size_t h_ini_size = h_ini.size()
    nlp.AddVariableSet(std::make_shared<Variables>  ( h_variables_name
                                                    , h_ini_size
                                                    , h_ini
                                                    ));

    // add u variables
    std::string u_variable_name = "u_variables";
    vector_t u_ini = 0.5 + 0.5 * flux_variables_ini_signs;
    size_t u_ini_size = u_ini.size()
    nlp.AddVariableSet(std::make_shared<Variables>  ( u_variable_name
                                                    , u_ini_size
                                                    , u_ini
                                                    ));

    //
    // add cost function_class
    //
    std::string cost_class_name = "metab_cost"; // ifopt stuff really likes it when you name it
    nlp.AddCostSet(std::make_shared<Cost>   ( cost_class_name
                                            , metabolites_variable_name
                                            ));
    
    //
    // add constraint class
    //
    std::string cost_class_name = "metab_constraints"; // ifopt stuff really likes it when you name it
    nlp.AddConstraintSet(std::make_shared<Constraint>   ( cost_class_name
                                                        , num_total_metabolies
                                                        , num_variable_metabolites
                                                        , fixed_metabolites
                                                        , variable_metabolites_variables_name
                                                        , flux_variables_name
                                                        , steady_state_variables_name
                                                        , h_variables_name
                                                        , u_variable_name
                                                        , beta_variables_name
                                                        , null_space_stoich_matrix_variable_metab_section
                                                        , dim_null_space
                                                        , stoichiometric_matrix
                                                        , big_M_value
                                                        , equilibrium_constants
                                                        , target_log_variable_metabolites_counts
                                                        , variable_metabolite_lower_bound
                                                        ));

    //
    // run the solver
    //
    
    // print out initial state
    nlp.PrintCurrent();

    // initialize solver
    ifopt::IpoptSolver ipopt;
    ipopt.SetOption("linear_solver", "ma57"); // change if desired
    ipopt.SetOption("jacobian_approximation", "exact"); // keep this because we specified jacobians

    // solve the problem
    ipopt.Solve(nlp);

    // get the found variable solutions
    vector_t metabolites_sol = nlp.GetOptVariables()->GetValues(variable_metabolites_variables_name);
    vector_t steady_state_sol = nlp.GetOptVariables()->GetValues(steady_state_variables_name);
    vector_t flux_sol = nlp.GetOptVariables()->GetValues(flux_variables_name);
    vector_t beta_sol = nlp.GetOptVariables()->GetValues(beta_variables_name);
    vector_t h_sol = nlp.GetOptVariables()->GetValues(h_variables_name);

    vector_t E_regulation = vector_t::Constant(flux_sol.size(), 1.0);
    vector_t unreq_rxn_flux = rxn_flux  ( metabolites_sol
                                        , fixed_metabolites
                                        , stoichiometric_matrix_T
                                        , equilibrium_constants
                                        , E_regulation
                                        );

    vector_t alpha_sol = 
        ( flux_sol.array() / unreq_rxn_flux.array() ).matrix();

    return { steady_state_sol, flux_sol, alpha_sol, beta_sol, metabolite_sol}; 

}
  

vector_t
rxn_flux
    ( vector_t v_log_counts
    , vector_t f_log_counts
    , matrix_t S
    , matrix_t K
    , matrix_t E_regulation
    )
{
    /* Flip Stoichiometric Matrix */
    matrix_t S_T = S; /* this should do deep copy */
    S.transposeInPlace();

    size_t n_react = S.cols();

    v_log_counts = v_log_counts.reshaped<Eigen::RowMajor>().eval();
    f_log_counts = f_log_counts.reshaped<Eigen::RowMajor>().eval();

    vector_t tot_log_count(v_log_counts.size() + f_log_counts.ize());
    tot_log_counts << v_log_counts, f_log_counts;

    K = K.reshaped<Eigen::RowMajor>().eval();
    E_regulation = E_regulation.reshaped<Eigen::RowMajor>().eval();

    auto coeff = (-0.25 * S_T * tot_log_counts).array().exp().pow(4)
    auto forward_odds = K.array() * coeff;
    auto reverse_odds = K.array().inverse() * coeff;

    vector_t result = 
        ( E_regulation.array()
        * ( forward_odds - reverse_odds )
        ).matrix()

    return result;
}