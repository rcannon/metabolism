
import numpy as np
import maximum_entropy_pyomo_relaxed as mep
from metabolism_driver import read_input

def test():

    N_avogadro = 6.022140857e+23
    VolCell = 1.0e-15
    Concentration2Count = N_avogadro * VolCell
    concentration_increment = 1/(N_avogadro*VolCell)

    conc_file = 'concentrations.csv'
    K_file = 'EquilibriumConstants.csv'
    S_file = 'StoichiometricMatrix.csv'

    v_log_counts, f_log_counts, obj_rxn_idx, K, S = read_input(conc_file, K_file, S_file)
    target_log_vcounts = np.log(np.ones(np.shape(v_log_counts)) *1.0e-03*Concentration2Count)
    
    num_reactions = np.shape(S)[0]

    n_init = v_log_counts
    beta_init = np.ones(num_reactions)
    y_init = np.ones(num_reactions)

    [b_sol, y_sol, alpha_sol, h_sol, beta_sol, n_sol] = mep.Max_ent_solver(n_init,y_init,beta_init,target_log_vcounts, f_log_counts, S, K, obj_rxn_idx)

    return 0
