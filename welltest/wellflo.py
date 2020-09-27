"""
Modeling and Analysis of Well Tests

@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

def radius_dimensionless(re, rw):
  """Calculate dimensionless radius (rD)"""
    return re / rw

def time_dimensionless(perm, t, poro, mu, ct, rw):
  """Calculate dimensionless time (tD)"""
  return (.0002637 * perm * t) / (poro * mu * ct * (rw**2))

def pressure_multirate(pD, delta_q, pi, B, mu, perm, h):
  """Calculate Flowing Pressure as Sum of Constant Rates"""
  import numpy as np
  return pi - ((B * mu / (.007082 * perm * h)) * (np.sum(pD * delta_q)))

def rate_multipressure(qD, delta_p, B, mu, perm, h):
  """Calculate Rate as Sum of Constant Flowing Pressures"""
  import numpy as np
  return ((.007082 * perm * h) / (B * mu)) * (np.sum(qD * delta_p))

def time_finite_acting(perm, poro, mu, ct, rw, re):
  """Calculate time at flow starts behaving infinite-acting"""
  r_D = re / rw
  t_Dw = 0.25 * r_D**2
  return (poro * mu * ct * (rw**2) * t_Dw) / (.0002637 * perm)
