#!/Users/d3k137/.pyenv/shims/python
# -*- coding: utf-8 -*-

from os.path import exists
import numpy as np
import metabolism
from importlib import reload
N_avogadro = 6.022140857e+23
VolCell = 1.0e-15

Concentration2Count = N_avogadro * VolCell
concentration_increment = 1/(N_avogadro*VolCell)
#print("exists")
#print(exists('../data/concentrations.csv'))

def read_input(conc_file, K_file, S_file):
  with open('../data/concentrations.csv', 'r', encoding='UTF8') as f:
      lines = f.read().splitlines()
      x = lines[0].split()
      nvar = int(x[0])
      x = np.asarray(lines[2].split(","),dtype=np.float64)
      v_log_counts = x[:nvar]
      f_log_counts = x[nvar:]
    
  with open('../data/EquilibriumConstants.csv', 'r', encoding='UTF8') as f:
      lines = f.read().splitlines()
      x = lines[0].split()
      iuptake = int(x[0])
      x = lines[1].split()
      ioutput = int(x[0])
      #x = lines[2].split(" ",1)
      #y = (x[1].replace('[','').replace(']','').split(", "))
      obj_rxn_idx = [int(i) for i in lines[2].split(',')]
      #print(obj_rxn_idx)

      K = np.asarray(lines[3].split(","),dtype=np.float64)
    
  with open('../data/StoichiometricMatrix.csv', 'r', encoding='UTF8') as f:
      S = np.loadtxt(f, delimiter=",")

  return(v_log_counts, f_log_counts, obj_rxn_idx, K, S)

def metabolism_driver():
  conc_file = '../data/concentrations.csv'
  K_file = '../data/EquilibriumConstants.csv'
  S_file = '../data/StoichiometricMatrix.csv'

  Vmax = 1000 #int
  s = 0.085 #float
  Km = 0.5*s  #float
  iuptake = 32
  ioutput = 32

  print('Running initial metabolism model')
  v_log_counts, f_log_counts, obj_rxn_idx, K, S = read_input(conc_file, K_file, S_file)
  target_log_vcounts = np.log(np.ones(np.shape(v_log_counts)) *1.0e-03*Concentration2Count) # 1-D array of floats  
  # Run initial metabolism:
  rxn_flux,v_log_counts = metabolism.run(v_log_counts,f_log_counts,target_log_vcounts, S, K, iuptake, Vmax, Km, s, obj_rxn_idx)
  #rxn_flux = rxn_flux_ip

  nutrient_ratio = 1
  # Mimic a run of the fungal model with updates to external nutrient s:
  for i in range(5):
    print('Step ', i)
    s_new = s + 0.04*s

    nutrient_ratio = nutrient_ratio * s_new/s #Whether or not metabolism is recalculated depends on whether s is significantly different,
                         # not on whether V = Vmax*s/(Km + s) is significantly different.
    if ((nutrient_ratio > 1.05) or (nutrient_ratio < 0.95)):
      #call detailed metabolism
      #rxn_flux, v_log_counts = metabolism.run(v_log_counts,f_log_counts,target_log_vcounts, S, K,iuptake, Vmax, Km, s, obj_rxn_idx)
      nutrient_ratio = 1

    scaled_flux2 = metabolism.scale_flux(rxn_flux,iuptake, Vmax, Km, s)

    print(scaled_flux2[iuptake])

  return()

metabolism_driver()
