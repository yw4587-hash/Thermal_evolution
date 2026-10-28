"""Shock transition (from shock_transition.ipynb). Used for Fig. 2."""
import numpy as np

from cooling import gamma, kb, tcool_P


def shock_tran(M_1, rho_1, T_1, m):
  # presrure and bulk velocity
  # m here is the average mass of the gases rho_tot = n_tot * m
  P_1 = rho_1* kb * T_1 /m
  u_1 = M_1 * (gamma * P_1 /rho_1)**.5 # M = u / v_sound
  # after shock density, tem, and bulk velocity
  rho_2 = rho_1 * ((gamma+1) * M_1**2)/((gamma - 1) * M_1**2 + 2)
  T_2 = ((gamma-1)* M_1** 2 + 2)*(2* gamma* M_1**2 - (gamma - 1))/((gamma + 1)**2 * M_1**2)*T_1
  u_2 = u_1 * rho_1/rho_2
  # return the ratio
  return u_1/u_2, rho_1/rho_2, T_2/T_1


def time_diff(M_1, T_1, n_1, m, metallicity):
  """
  cooling time diff: log10(t_2/t_1)
  M_1 has no unit
  T_1 in unit of k
  n_1 in unit of cm **-3
  """
  rho_1 = n_1 * m
  u_r, rho_r, T_r = shock_tran(M_1, rho_1, T_1, m)
  T_2 = T_r * T_1
  n_2 = n_1 * ((gamma+1) * M_1**2)/((gamma - 1) * M_1**2 + 2)
  tcool1 = tcool_P(T_1, n_1, metallicity)
  tcool2 = tcool_P(T_2, n_2, metallicity)
  diff = tcool2 - tcool1
  return diff
