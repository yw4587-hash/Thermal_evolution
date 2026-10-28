"""
Constants, cooling function, and cooling time (shared by all three studies).

Define of all the unit:
CGS system
P in units of erg * cm **-3
T in unit of K
m in unit of g
rho in unit of g cm **-3
metallicity in units of solar metallicity
n is the TOTAL number density; nH = n * mu / muH
"""
import os

import h5py
import numpy as np
from scipy import interpolate

"""
Define of all constant:
"""
mu = 0.62
muH = 1/0.75
kb = 1.380649e-16
gamma = 5/3
Z_solar = 0.02
mH = 1.673536 * 10**(-24)

"""
set up: cooling function lambda
"""
_here = os.path.dirname(os.path.abspath(__file__))
Cooling_File = os.environ.get("CGM_COOLING_TABLE", os.path.join(_here, "data", "z_0.000.hdf5"))  ### From Wiersma et al. (2009) appropriate for z=0 UVB
f = h5py.File(Cooling_File, 'r')
i_X_He = -3
Metal_free   = f.get('Metal_free')
Total_Metals = f.get('Total_Metals')
log_Tbins    = np.array(np.log10(Metal_free['Temperature_bins']))
log_nHbins   = np.array(np.log10(Metal_free['Hydrogen_density_bins']))
Cooling_Metal_free   = np.array(Metal_free['Net_Cooling'])[i_X_He]
Cooling_Total_Metals = np.array(Total_Metals['Net_cooling'])
HHeCooling           = interpolate.RectBivariateSpline(log_Tbins,log_nHbins, Cooling_Metal_free)
ZCooling             = interpolate.RectBivariateSpline(log_Tbins,log_nHbins, Cooling_Total_Metals)
Zs          = np.logspace(-2,1,31)
f.close()
Lambda_tab  = np.array([[[HHeCooling.ev(lT,ln)+Z*ZCooling.ev(lT,ln) for Z in Zs] for lT in log_Tbins] for ln in log_nHbins])
# cooling time is dependent on initial temp and the number density.
Lambda_z0 = interpolate.RegularGridInterpolator((log_nHbins,log_Tbins,Zs), Lambda_tab, bounds_error=False, fill_value=-1e-30)


def tcool_P(T,n, metallicity):
  """
  cooling time function
  T in units of K
  n total number density in cm**-3
  metallicity in units of solar metallicity
  returning time is log10(s)
  """
  T = np.where(T>10**8.98, 10**8.98, T)
  T = np.where(T<10**2, 10**2, T)
  nH_actual = np.clip(n*(mu/muH), None, 1)
  # n here is the total number density as well, but we convert total number density to hydrogen number density
  nH = np.where(nH_actual>100, 100, nH_actual)
  nH = np.where(nH<10**-8, 10**-8, nH)
  return np.log10((1./(gamma-1.)) * (muH/mu) * kb * T / ( nH_actual * Lambda_z0((np.log10(nH),np.log10(T), metallicity)).squeeze()))
