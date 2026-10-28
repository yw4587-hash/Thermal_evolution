"""Continuous compression with radiative cooling (from infall_galaxy.ipynb). Used for Figs. 4 and 5."""
import numpy as np
from scipy.integrate import solve_ivp

from cooling import Lambda_z0, gamma, kb, mu, muH


def cont_radi(t, n, P0, kai, period, metallicity=0.3):
  """ODE of changing density with compression and radiative cooling: dn/dt = func"""
  P = P0 + P0* (kai - 1)/period * t
  T = P/(kb * n)
  nH = min((mu/ muH) * n, 1)
  result = n * (kai - 1)* P0 /(gamma * period * P) + (gamma - 1)* n* Lambda_z0((np.log10(nH), np.log10(T), metallicity)) * nH**2 / (gamma* (P))
  return result


def cool_radiation(n0, T0, kai, metallicity, period):
  """
  The func that gives log10 of the final density for a given period with radiation
  """
  P0 = kb * n0 * T0
  t_span = np.linspace(0, period, 1000)
  result = solve_ivp(lambda t, n: cont_radi(t, n, P0, kai, period, metallicity), [t_span[0], t_span[-1]], [n0], t_eval= t_span)
  n_end = result.y[0, :][-1]
  return np.log10(n_end)


def entropy_diff(n0, T0, kai, metallicity, period):
  """
  The function that gives the entropy ratio log10(K_start/K_end) for a given period with radiation
  """
  P0 = kb * n0 * T0
  P_end = kai * P0
  n_end = 10** cool_radiation(n0, T0, kai, metallicity, period)
  T_end = P_end/(kb * n_end)
  K_end = T_end / n_end**(2/3)
  K_start = T0 / n0**(2/3)
  return np.log10(K_start/K_end)


def cool_radiation_single(n0, T0, kai, metallicity, period):
  """
  Single-point evolution: returns log10(n(t)) and log10(T(t))
  """
  P0 = kb * n0 * T0
  t_span = np.linspace(0, period, 10000)
  result = solve_ivp(lambda t, n: cont_radi(t, n, P0, kai, period, metallicity), [t_span[0], t_span[-1]], [n0], t_eval= t_span)
  P_end = kai * P0
  p_span = P0 + (- P0 + P_end)/ period * result.t
  n_span = result.y
  T_span = p_span / (kb* n_span)
  return np.log10(n_span), np.log10(T_span)
