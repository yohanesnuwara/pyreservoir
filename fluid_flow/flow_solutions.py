def time_to_dimensionless(time, rw, poro, mu_oil, ct, k):
  return (0.000263 * k * time) / (poro * mu_oil * ct * (rw**2))

def dimensionless_radius(re, rw): return re / rw

def time_finite_acting(re, rw, poro, mu_oil, ct, k):
  r_eD = re / rw # dimensionless radius
  t_Dw = 0.25 * r_eD**2 # dimensionless time
  return (poro * mu_oil * ct * (rw**2) * t_Dw) / (0.0002637 * k)

def constant_terminal_rate(time, distance, re, rw, pi, q, poro, ct, k, h, mu_oil, Bo):
  """
  Constant Terminal Rate Solution (Approximation Method)

  INPUT:

  time: Time at which flow is evaluated, hour
  distance: Distance from the wellbore, ft (NOT distance from centre of wellbore)
  re: Reservoir extent, ft
  rw: Wellbore radius, ft
  pi: Initial reservoir pressure, psia
  q: Wellbore flowing rate, STB/D
  poro: Porosity
  ct: Total compressibility, sip
  k: Permeability, md
  h: Reservoir net thickness, ft
  mu_oil: Oil viscosity, cp
  Bo: Oil FVF, RB/STB

  OUTPUT:

  td: Dimensionless time
  pd: Dimensionless pressure
  pwf: Wellbore flowing pressure (psia)
  """
  import numpy as np
  from scipy.special import expi
  
  # Access Ei-function table
  Ei_table = np.loadtxt("https://raw.githubusercontent.com/yohanesnuwara/pyreservoir/master/fluid_flow/Ei_table.txt")
  
  r = rw + distance
  t_finite_acting = time_finite_acting(re, rw, poro, mu_oil, ct, k)

  if time > 0 and time < t_finite_acting:
    """Time behaving infinite acting"""
    td = time_to_dimensionless(time, rw, poro, mu_oil, ct, k)
    if r==rw:
      # Your distance is at the wellbore
      if td > 100:
        # Eq 6.20
        pd = 0.5 * (((np.log(td)) + 0.80907)) 
        pwf = pi - ((pd * q * Bo * mu_oil) / (0.007082 * k * h))
      if td < 100:
        # No solution
        pd = np.nan
        pwf = np.nan

    if r>rw:
      # Your distance is away from the wellbore, in the reservoir
      td = time_to_dimensionless(time, r, poro, mu_oil, ct, k)
      if td > 12.5:
        # pd can be approximated using Eq 6.28

        pd = 0.5 * (np.log(td) + 0.80907)
        # pd_arr.append(float(pd))
        
        "Calculate pwf after n hours"
        pwf = pi - ((pd * q * Bo * mu_oil) / (0.007082 * k * h))  

      if td < 12.5:
        # pd calculated using Eq 6.26. Find the value of integral exponent function -Ei(-x) using tabulation    

        x = 0.25 * (1 / td)

        if x >= 0 and x <= 0.209:
          x_new = round(x, 3)

          # "Tabulation value finder"
          index = np.where(Ei_table[:,0] == x_new)
          index = np.array((index)[0])
          index = int(index)
          minusEi = Ei_table[index, 1]

        if x > 0.209 and x <= 2.09:
          x_new = round(x, 2) 

          # "Tabulation value finder"
          index = np.where(Ei_table[:,0] == x_new)
          index = np.array((index)[0])
          index = int(index)
          minusEi = Ei_table[index, 1]

        if x > 2.09 and x <= 10.9:
          x_new = round(x, 1) 

          # "Tabulation value finder"
          index = np.where(Ei_table[:,0] == x_new)
          index = np.array((index)[0])
          index = int(index)
          minusEi = Ei_table[index, 1]

        if x > 10.9:
          # if x above 10.9, meaning Table A-1 can't be used because it's limited to only x below 10.9. so use scipy
          x_new = x
          minusEi = -expi(-x) # from scipy.expi
        
        "Calculate pd"
        pd = 0.5 * minusEi

        "Calculate pwf after n hours"
        pwf = pi - ((pd * q * Bo * mu_oil) / (0.007082 * k * h))
      # pd, pwf = np.nan, np.nan

  elif time == 0:
    """Time at start of flow"""
    td = time_to_dimensionless(time, rw, poro, mu_oil, ct, k)
    pd = np.nan
    pwf = pi
    
  elif time >= t_finite_acting:
    """Time behaving finite acting"""
    td = time_to_dimensionless(time, rw, poro, mu_oil, ct, k)
    if r==rw:
      # Your distance is at wellbore
      if td > 25:
        # Calculate dimensionless radius at reservoir outer boundary (r=re)
        r_eD = dimensionless_radius(re, rw)
        # Eq 6.19
        pd = (2 * td / (r_eD**2)) + np.log(r_eD) - 0.75 # Eq 6.19
        pwf = pi - ((pd * q * Bo * mu_oil) / (0.007082 * k * h))

      if td < 25:
        # No solution
        pd = np.nan
        pwf = np.nan
    
    if r>rw:
      # Your distance is outside the wellbore, inside the reservoir
      pd, pwf = np.nan, np.nan # No solution

  return td, pd, pwf
