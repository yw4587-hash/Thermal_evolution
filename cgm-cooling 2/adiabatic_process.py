"""Adiabatic transition (from adiabatic_process.ipynb). Used for Figs. 1 and 3."""
import numpy as np

from cooling import gamma, kb, mH, mu, tcool_P


# shock transition
def shock_tran(M_1, n_1, T_1):
  # presrure and bulk velocity
  P_1 = n_1* kb * T_1
  rho_1 = n_1*mH*mu
  u_1 = M_1 * (gamma * P_1 /(n_1*mH*mu))**.5 # M = u / v_sound
  # after shock density, tem, and bulk velocity
  rho_2 = rho_1 * ((gamma+1) * M_1**2)/((gamma - 1) * M_1**2 + 2)
  T_2 = ((gamma-1)* M_1** 2 + 2)*(2* gamma* M_1**2 - (gamma - 1))/((gamma + 1)**2 * M_1**2)*T_1
  u_2 = u_1 * rho_1/rho_2
  # return the ratio
  return u_2/u_1, rho_2/rho_1, T_2/T_1


# adiabatic transition
def adi_tran(n1, T1, ka):
  """
  ka is p2/p1
  """
  p1 = kb * n1 * T1
  T2 = T1 * ka ** ((gamma-1)/gamma)
  n2 = n1 * ka ** (1/gamma)
  return n2, T2


def time_diff(n1, T1, kai, metallicity):
  """
  cooling time diff: log10(t_2/t_1)
  n in unit of cm **-3
  T_1 in unit of k
  """
  n2, T2 = adi_tran(n1, T1, kai)
  tcool1 = tcool_P(T1, n1, metallicity)
  tcool2 = tcool_P(T2, n2, metallicity)
  diff = tcool2 - tcool1
  return diff


def kai_func(T, n):
  return np.log(T/n**(2/3))
