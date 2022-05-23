import numpy as np
import pyomo
import pyomo.environ as pe
import itertools
import numpy.random as nprd
from scipy.linalg import norm
from scipy.optimize import least_squares
import scipy.linalg as spL
import pyutilib.services
from pyomo.opt import TerminationCondition
from importlib import reload
import timeit

import maximum_entropy_pyomo_relaxed as mep
import feasible_point_pyomo_relaxed as fep
mep = reload(mep)


def initialize(v_log_counts,f_log_counts,target_log_vcounts, S, K):

  ### SOLVE FOR AN INITIAL FEASIBLE POINT
  [solved,feas_obj,bi_sol, yi_sol, alphai_sol, hi_sol, betai_sol, ni_sol] = fep.Feasible_point_solver(v_log_counts, target_log_vcounts, f_log_counts, S, K)                            

  v_log_counts = ni_sol

  # find out if it terminated on optimal solution
  print(feas_obj)
  if solved == 0:
    print("solver failed to converge")
  else: 
    print("Solver Converged")

  #set an objective tolerance to call solution feasible
  feasible_tol = 1e-12

  if feas_obj > feasible_tol:
    print("Solver did not find a feasible point")
    print("Feasibility tolerance = ", feasible_tol, "\tObjective tolerance = ", feas_obj)    

  # set the initial condition as the solution from the feasible point solve
  n_ini = np.reshape(ni_sol,(len(ni_sol),) )
  b_ini = np.reshape(bi_sol,(len(bi_sol),) )
  y_ini = np.reshape(yi_sol,(len(yi_sol),) )
  beta_ini = np.reshape(betai_sol,(len(betai_sol),) )
  h_ini = np.reshape(hi_sol,(len(hi_sol),) )

  # save for c++
  path = "../data/python_feasible_point/"
  np.savetxt(path + "variable_metabolites_log_counts.csv",n_ini, delimiter=',')
  np.savetxt(path + "flux_variables.csv", y_ini, delimiter=',')
  np.savetxt(path + "null_space_variables.csv",beta_ini, delimiter=',')

  return(beta_ini, y_ini,n_ini)

def optimize(n_init, y_init, beta_init, f_log_counts,target_log_vcounts, S, K,obj_rxn_idx, max_iter = 10000):

  start_t = timeit.timeit()
  [b_sol, y_sol, alpha_sol, h_sol, beta_sol, n_sol] = mep.Max_ent_solver(n_init,y_init,beta_init,target_log_vcounts, f_log_counts, S, K, obj_rxn_idx)
  end_t = timeit.timeit()
  print('time Max_ent_solver =', end_t - start_t)

  E_regulation = alpha_sol
  v_log_counts = n_sol
  P = np.where(S>0,S,0)
  R = np.where(S<0, S, 0)
  log_counts = x = np.concatenate((v_log_counts, f_log_counts))
  delta = delta = np.zeros(log_counts.size)
  rxn_flux = oddsDiff(v_log_counts, f_log_counts, 1, S, R, P, delta, K, E_regulation)
  dndt = S.T.dot(rxn_flux)
  #print('Rxn Flux:', rxn_flux)
  #print('Metabolite Derivatives:',dndt)

  return (b_sol, y_sol, n_sol, alpha_sol)

def oddsDiff(vcounts,fcounts,mu0,S, R, P, delta,Keq,E_Regulation):
  metabolites = np.append(vcounts,fcounts)
  KQ_f = odds(metabolites,mu0,S, R, P, delta,Keq);
  Keq_inverse = np.power(Keq,-1);
  KQ_r = odds(metabolites,mu0,-S, P, R, delta,Keq_inverse,-1);
  
  #WARNING: Multiply regulation here, not on individual Keq values.
  KQdiff =  E_Regulation * (KQ_f - KQ_r);
  return(KQdiff)

def odds(log_counts,mu0,S, R, P, delta, K, direction = 1):
  counts = np.exp(log_counts)
  delta_counts = counts+delta;
  log_delta = np.log(delta_counts);
  Q_inv = np.exp(-direction*(R.dot(log_counts) + P.dot(log_delta)))
  KQ = np.multiply(K,Q_inv);
  return(KQ)

def scale_flux(rxn_flux,iglucose, Vmax, Km, s):
#  inputs:
#  Kmichealis constant, Km
#  Vmax, Vmax
# External glucose concentration, s
# index of glucose uptake reaction, iglucose
# reaction flux from optimization, rxn_flux

  V = Vmax * s/(Km + s)
  scaled_flux = V/rxn_flux[iglucose] * rxn_flux
  return(scaled_flux)

def run(v_log_counts,f_log_counts,target_log_vcounts, S, K, iuptake, Vmax, Km, s, obj_rxn_idx):

  beta_init, y_init,n_init = initialize(v_log_counts,f_log_counts,target_log_vcounts, S, K)
  b_sol, y_sol, n_sol, alpha_sol = optimize(n_init, y_init, beta_init, f_log_counts,target_log_vcounts, S, K,obj_rxn_idx)

  E_regulation = alpha_sol
  v_log_counts = n_sol
  P = np.where(S>0,S,0)
  R = np.where(S<0, S, 0)
  log_counts = x = np.concatenate((v_log_counts, f_log_counts))
  delta = delta = np.zeros(log_counts.size)
  rxn_flux = oddsDiff(v_log_counts, f_log_counts, 1, S, R, P, delta, K, E_regulation)
  #dndt = S.T.dot(rxn_flux)

  #scaled_flux = scale_flux(rxn_flux,iuptake, Vmax, Km, s)
  return(rxn_flux,n_sol)

  

#def report():
#  inputs:
#  Kmichealis
#  Vmax
#  index of cell wall
#  index of glucose
#  alpha_sol - regulation
#  n_sol - variable log counts of metabolites 
#  f_log_counts - fixed log counts of metabolites
#  K, S - equilibrium constants and stoichiometric matrix

#  compute P, R, delta
#  E_regulation_pyomo = alpha_sol
#  rxn_flux_ip = oddsDiff(n_sol, f_log_counts, mu0, S, R, P, delta, K, E_regulation_pyomo)
#  print('Metabolite Derivatives:')
#   V = V_max*s/(km + s)
#  dndt = S_active.T.dot(rxn_flux_ip)
#  print(dndt)

  #output:
#  return(V, glucose_uptake, cell_wall_output)
