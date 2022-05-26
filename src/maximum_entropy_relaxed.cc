
#include "maximum_entropy_relaxed.hh"

using namespace ifopt;


std::tuple<vector_t, vector_t, vector_t, vector_t, vector_t, vector_t>
maximum_entropy_solver
    ( const vector_t& variable_metabolites_init // aka n
    , const vector_t& fixed_metabolites
    , const vector_t& flux_variables_init // aka y
    , const vector_t& beta_init // \beta
    , const vector_t& target_log_variable_metabolites_counts
    , const matrix_t& stoich_matrix
    , const vector_t& equilibrium_constants // aka K
    , const index_list_t& objective_reaction_indices
    )
{
    //
    // setup 
    //
    
    size_t num_fixed_metabolites = fixed_metabolites.size();
    size_t num_variable_metabolites = variable_metabolites_init.size();
    size_t num_total_metabolies = num_fixed_metabolites + num_variable_metabolites;

    value_t variable_metabolite_lower_bound = -300;
    value_t big_M_value = 1000;

    matrix_t stoichiometric_matrix = stoich_matrix; // copying because bad practive but easy, also see line 314
    
    // Stoichiometric Matrix - rows are reactions, cols are metabolites
    matrix_t stoichiometric_matrix_T = stoichiometric_matrix; // this should do deep copy 
    stoichiometric_matrix.transposeInPlace(); // now rows are metabolites, cols are reactions
    size_t n_reactions = stoichiometric_matrix.cols(); // should be same as numpy.shape(S)[1]

    // Split S into the component corresponding to the variable metabolites S_v
    // and the component corresponding to the fixed metabolites. 
        // Reilly: below should do the same thing as numpy.delete
        // https://eigen.tuxfamily.org/dox-devel/group__TutorialSlicingIndexing.html

    matrix_t stoich_matrix_variable_metab_section_T 
        = stoichiometric_matrix_T(Eigen::all, Eigen::seq(0, num_variable_metabolites-1));
    matrix_t stoich_matrix_fixed_metab_section_T 
        = stoichiometric_matrix_T(Eigen::all, Eigen::seq(num_variable_metabolites, stoichiometric_matrix_T.cols()-1));

    matrix_t stoich_matrix_variable_metab_section = stoich_matrix_variable_metab_section_T.transpose();
    matrix_t stoich_matrix_fixed_metab_section = stoich_matrix_fixed_metab_section_T.transpose();

    // find a basis for the nullspace of S_v,
    // and get the dimension of the null space 
    Eigen::JacobiSVD<matrix_t> svd(stoich_matrix_variable_metab_section, Eigen::ComputeFullV);
    int rnk = svd.rank();
    matrix_t null_space_stoich_matrix_variable_metab_section 
        = svd.matrixV()(Eigen::seq(rnk, Eigen::last), Eigen::all).transpose().conjugate();
    
    //matrix_t null_space_stoich_matrix_variable_metab_section = read_to_matrix("../data/python_feasible_point/null_space_matrix.csv");
    size_t dim_null_space = null_space_stoich_matrix_variable_metab_section.cols();
    

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
    bool is_u_variables;

    std::cout << "\ntarget metabolites " << std::endl;
    for (int i = 0; i < 10; i++   ) {
        std::cout << target_log_variable_metabolites_counts(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "\nfixed_metabolites " << std::endl;
    for (int i = 0; i < 10; i++   ) {
        std::cout << fixed_metabolites(i) << " ";
    }
    std::cout << std::endl;

    // add metabolite variables
    is_u_variables = false;
    std::string variable_metabolites_variables_name = "log_variable_metabolites";
    nlp.AddVariableSet(std::make_shared<Variables>  ( variable_metabolites_variables_name
                                                    , num_variable_metabolites
                                                    , is_u_variables
                                                    , variable_metabolites_init
                                                    ));
    std::cout << "\nvariable_metabolites " << std::endl;
    for (int i = 0; i < 10; i++   ) {
        std::cout << variable_metabolites_init(i) << " ";
    }
    std::cout << std::endl;


    // add beta variables
    std::string beta_variables_name = "beta_variables";
    size_t num_beta_variables = dim_null_space;
    vector_t beta_variables_for_use = beta_init(Eigen::seq(0, num_beta_variables-1));
    is_u_variables = false;
    nlp.AddVariableSet(std::make_shared<Variables>  ( beta_variables_name
                                                    , num_beta_variables
                                                    , is_u_variables
                                                    , beta_variables_for_use
                                                    ));
    std::cout << "\nbeta " << std::endl;
    for (int i = 0; i < 10; i++   ) {
        std::cout << beta_variables_for_use(i) << " ";
    }
    std::cout << std::endl;

    // add flux variables (aka y)
    std::string flux_variables_name = "flux_variables";
    size_t num_flux_variables = n_reactions;
    is_u_variables = false;
    //vector_t flux_variables_init = null_space_stoich_matrix_variable_metab_section * beta_variables_for_use;
    nlp.AddVariableSet(std::make_shared<Variables>  ( flux_variables_name
                                                    , num_flux_variables
                                                    , is_u_variables
                                                    , flux_variables_init
                                                    ));
    std::cout << "\nflux " << std::endl;
    for (int i = 0; i < 10; i++   ) {
        std::cout << flux_variables_init(i) << " ";
    }
    std::cout << std::endl;

    // add steady state variables (aka g in paper, b in python version)
    std::string steady_state_variables_name = "steady_state_variables";
    vector_t steady_state_vars_init
        = (stoich_matrix_variable_metab_section_T * variable_metabolites_init) 
        + (stoich_matrix_fixed_metab_section_T * fixed_metabolites );
    size_t num_steady_state_variables = n_reactions;
    is_u_variables = false;
    nlp.AddVariableSet(std::make_shared<Variables>  ( steady_state_variables_name
                                                    , num_steady_state_variables
                                                    , is_u_variables
                                                    , steady_state_vars_init
                                                    ));
    std::cout << "\nsteady state " << std::endl;
    for (int i = 0; i < 10; i++   ) {
        std::cout << steady_state_vars_init(i) << " ";
    }
    std::cout << std::endl;

    /*std::cout << "\nvariables" << std::endl;
    for ( auto row : stoich_matrix_fixed_metab_section_T.rowwise()) {
        for (auto elem : row) {
            std::cout << elem << "\n";
        }
        std::cout << "\n" << std::endl;
    }*/
    /*    std::cout << "\nvariables" << std::endl;
    for (auto elem : flux_variables_init_signs) {
        std::cout << elem << std::endl;
    }*/

    // add h variables
    std::string h_variables_name = "h_variables";
    vector_t flux_variables_init_signs = flux_variables_init.array().sign().matrix();
    vector_t fvia = flux_variables_init.array();
    vector_t h_init = 
        (   flux_variables_init_signs.array()
        *   (   std::log(2) 
            -   Eigen::log  (   fvia.array().abs()
                            + Eigen::sqrt(fvia.array().pow(2) + 4)  
                            )
            )
        ).matrix();
    size_t num_h_variables = n_reactions;
    is_u_variables = false;
    //std::cout << "\nmax ent h " << h_init(0) << std::endl;
    //std::cout << "\nmax ent flux " << flux_variables_init(0) << std::endl;
    nlp.AddVariableSet(std::make_shared<Variables>  ( h_variables_name
                                                    , num_h_variables
                                                    , is_u_variables
                                                    , h_init
                                                    ));
    std::cout << "\nh vars " << std::endl;
    for (int i = 0; i < 10; i++   ) {
        std::cout << h_init(i) << " ";
    }
    std::cout << std::endl;

    // add u variables
    std::string u_variable_name = "u_variables";
    vector_t u_init = 
        ( 0.5 
        + (0.5 * flux_variables_init_signs.array())
        ).matrix();
    size_t num_u_variables = n_reactions;
    is_u_variables = true;
    nlp.AddVariableSet(std::make_shared<Variables>  ( u_variable_name
                                                    , num_u_variables
                                                    , is_u_variables
                                                    , u_init
                                                    ));
    std::cout << "\nu vars " << std::endl;
    for (int i = 0; i < 10; i++   ) {
        std::cout << u_init(i) << " ";
    }
    std::cout << std::endl;

    //
    // add cost function_class
    //
    std::string cost_class_name = "metab_cost";
    nlp.AddCostSet(std::make_shared<Cost>   ( cost_class_name
                                            , flux_variables_name
                                            , objective_reaction_indices
                                            ));
    std::cout << "\nobjective reaction indices " << std::endl;
    for (int i = 0; i < 3; i++   ) {
        std::cout << objective_reaction_indices.at(i) << " ";
    }
    std::cout << std::endl;
    
    //
    // add constraint classes
    //

    // maximum entropy problem formulation 94
    std::string null_space_constraint_class_name = "null_space_constraints";
    nlp.AddConstraintSet(std::make_shared<NullSpaceRepresentationConstraint>  
        ( null_space_constraint_class_name
        , n_reactions
        , flux_variables_name
        , beta_variables_name
        , null_space_stoich_matrix_variable_metab_section
        , dim_null_space
        ));

    // maximum entropy problem formulation 95
    std::string steady_state_constraint_class_name = "steady_state_constraints";
    nlp.AddConstraintSet(std::make_shared<SteadyStateConstraint>  
        ( steady_state_constraint_class_name
        , n_reactions
        , num_total_metabolies
        , num_variable_metabolites
        , fixed_metabolites
        , variable_metabolites_variables_name
        , steady_state_variables_name
        , stoich_matrix_variable_metab_section_T
        , stoich_matrix_fixed_metab_section_T
        ));

    // maximum entropy problem formulation 96
    std::string smooth_constraint_class_name = "smooth_constraints";
    nlp.AddConstraintSet(std::make_shared<SmoothConstraint>  
        ( smooth_constraint_class_name
        , n_reactions
        , flux_variables_name
        , h_variables_name
        ));

    // maximum entropy problem formulation 97
    std::string relaxed_flux_upper_constraint_class_name = "relaxed_regulation_upper_constraints";
    nlp.AddConstraintSet(std::make_shared<RelaxedRegulationUpperConstraint>  
        ( relaxed_flux_upper_constraint_class_name
        , n_reactions
        , steady_state_variables_name
        , h_variables_name
        , u_variable_name
        , big_M_value
        , equilibrium_constants
        ));

    // maximum entropy problem formulation 98
    std::string relaxed_flux_lower_constraint_class_name = "relaxed_regulation_lower_constraints";
    nlp.AddConstraintSet(std::make_shared<RelaxedRegulationLowerConstraint>  
        ( relaxed_flux_lower_constraint_class_name
        , n_reactions
        , steady_state_variables_name
        , h_variables_name
        , u_variable_name
        , big_M_value
        , equilibrium_constants
        ));

    // maximum entropy problem formulation 99
    std::string sign_constraint_class_name = "sign_constraints";
    nlp.AddConstraintSet(std::make_shared<SignConstraint>  
        ( sign_constraint_class_name
        , n_reactions
        , steady_state_variables_name
        , flux_variables_name
        , equilibrium_constants
        ));

    // maximum entropy problem formulation 100
    std::string relaxed_regulation_sign_constraint_class_name = "relaxed_regulation_sign_constraints";
    nlp.AddConstraintSet(std::make_shared<RelaxedRegulationSignConstraint>  
        ( relaxed_regulation_sign_constraint_class_name
        , n_reactions
        , flux_variables_name
        , u_variable_name
        ));

    // maximum entropy problem formulation 101 upper
    std::string metabolites_upper_bound_constraint_class_name = "metabolites_upper_bound_constraints";
    nlp.AddConstraintSet(std::make_shared<MetabolitesUpperBoundConstraint>  
        ( metabolites_upper_bound_constraint_class_name
        , num_variable_metabolites
        , variable_metabolites_variables_name
        , target_log_variable_metabolites_counts
        ));

    // maximum entropy problem formulation 101 lower
    std::string metabolites_lower_bound_constraint_class_name = "metabolites_lower_bound_constraints";
    nlp.AddConstraintSet(std::make_shared<MetabolitesLowerBoundConstraint>  
        ( metabolites_lower_bound_constraint_class_name
        , num_variable_metabolites
        , variable_metabolites_variables_name
        , variable_metabolite_lower_bound
        ));
    
    //
    // run the solver
    //

    // print out initial state
    nlp.PrintCurrent();

    // initialize solver
    ifopt::IpoptSolver ipopt;
    ipopt.SetOption("linear_solver", "mumps");
    //ipopt.SetOption("jacobian_approximation", "exact"); // keep this because we specified jacobians
    ipopt.SetOption("jacobian_approximation", "finite-difference-values"); 
    ipopt.SetOption("max_cpu_time", 800000.0);
    ipopt.SetOption("max_iter", 10000);
    //ipopt.SetOption("print_level", 7);
    //ipopt.SetOption("tol", 0.00001)

    // solve the problem
    ipopt.Solve(nlp);
    nlp.PrintCurrent();

    // get the found variable solutions
    // reilly: make sure these are in the same order as added to problem (nlp.AddVariableSet).
    // hopefully github issue allows for better access
    //
    // note: Eigen::seq is inclusive on both ends
    vector_t all_sols = nlp.GetOptVariables()->GetValues();
    int begin = 0;
    int end = num_variable_metabolites -1;
    vector_t metabolites_sol = all_sols(Eigen::seq(begin, end));
    begin = end+1;
    end = begin + num_flux_variables -1;
    vector_t flux_sol = all_sols(Eigen::seq(begin, end));
    begin = end+1;
    end = begin + num_steady_state_variables -1;
    vector_t steady_state_sol = all_sols(Eigen::seq(begin, end));
    begin = end + 1;
    end = begin + num_beta_variables -1;
    vector_t beta_sol = all_sols(Eigen::seq(begin, end));
    begin = end +1;
    end = begin + num_h_variables-1;
    vector_t h_sol = all_sols(Eigen::seq(begin, end));
    begin = end + 1;
    end = begin + num_u_variables-1;
    assert( end + 1 == all_sols.size() );
    // not needed : vector_t u_sol = all_sols(Eigen::seq(begin, end));

    vector_t E_regulation = vector_t::Constant(flux_sol.size(), 1.0);
    vector_t unreq_rxn_flux = reaction_flux ( metabolites_sol
                                            , fixed_metabolites
                                            , stoichiometric_matrix_T
                                            , equilibrium_constants
                                            , E_regulation
                                            );

    vector_t alpha_sol = 
        ( flux_sol.array() / unreq_rxn_flux.array() ).matrix();

    return { steady_state_sol, flux_sol, alpha_sol, h_sol, beta_sol, metabolites_sol}; 

}
  

vector_t
reaction_flux
    ( const vector_t& variable_metabolites_log_counts
    , const vector_t& fixed_metabolites_log_counts
    , const matrix_t& stoich_matrix
    , const matrix_t& equilibrium_constants
    , const matrix_t& E_regulation
    )
{

    matrix_t stoichiometric_matrix = stoich_matrix;
    matrix_t stoichiometric_matrix_T = stoichiometric_matrix;
    stoichiometric_matrix.transposeInPlace();

    //size_t n_reactions = stochiometric_matrix.cols();

    vector_t total_log_counts(variable_metabolites_log_counts.size() + fixed_metabolites_log_counts.size());
    total_log_counts << variable_metabolites_log_counts, fixed_metabolites_log_counts;

    vector_t coeff = (-0.25 * stoichiometric_matrix_T * total_log_counts).array().exp().pow(4).matrix();
    vector_t forward_odds = 
        (
            equilibrium_constants.array() * coeff.array()
        ).matrix();
    vector_t reverse_odds = 
        ( 
            equilibrium_constants.array().inverse() * coeff.array()
        ).matrix();

    vector_t result = 
        ( E_regulation.array()
        * ( forward_odds.array() - reverse_odds.array() )
        ).matrix();

    return result;
}