"""
Modeling Single-Phase Flow in a Well Test

@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com

Note: 
The finite-acting time part of constant pressure solution is still not solved
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

def pressure_dimensionless(rD, tD):
    """
    Calculate Dimensionless Pressure from Constant Rate Flow
    """
    import numpy as np
    if tD < (0.25 * rD**2):
        # Infinite-acting solution for constant-rate (Towler, Eq. 6.20; from Lee, 1982)
        pD = 0.5 * (np.log(tD) + .80907)
    if tD > (0.25 * rD**2):
        # Finite-acting solution for constant-rate (Towler, Eq. 6.19; from Lee, 1982)
        pD = (2 * tD / rD**2) + np.log(rD) - .75
    return pD

def rate_dimensionless(rD, tD):
    """
    Calculate Dimensionless Rate from Constant Pressure Flow
    """
    import numpy as np
    import pandas as pd
    from scipy.interpolate import griddata
    
    if tD < (0.25 * rD**2):
        # Infinite-acting solution for constant-rate (Towler, Eq. 6.42, 6.43; from Edwardson et al, 1962)
        if tD > 0.01 and tD < 200:
          # Eq. 6.42
          qD = (26.7544 + (45.5537 * np.sqrt(tD)) + (13.3813 * tD) + (0.492949 * tD * np.sqrt(tD))) / ((47.4210 * np.sqrt(tD)) + (35.5372 * tD) + (2.60967 * tD * np.sqrt(tD)))
        if tD >= 200:
          # Eq. 6.43
          qD = ((2.02623 * tD * (np.log(tD) - 1)) + 3.90086) / (tD * ((np.log(tD))**2))

    if tD > (0.25 * rD**2):
        # Finite-acting solution for constant-rate (Towler, Eq. 7.32; from Jacob and Lohman, 1952)
        qD = np.nan
        qD = 2 / (np.log(tD) + .80907)
#         columns = ['rd', 'td', 'qd']
#         veh = pd.read_csv('/content/pyreservoir/welltest/Appendix A-4.csv', names=columns)
#         rd = veh['rd'].values
#         td = veh['td'].values
#         qd = veh['qd'].values

#         ## gridding and interpolation
#         data = np.stack((rd, td), axis=1)
#         qD = griddata(data, qd, [rD, tD], method='linear')

    return qD
  
def check_validity(solver='constant_rate', time='infinite', tmin=0.1, rw=0.5, re=1000, perm=100, poro=0.2, mu=2, ct=3E-6):
  """Check validity of using the Approaches to Flow Solutions"""
  import numpy as np
  if solver == 'constant_rate':
      if time == 'infinite':
          # Infinite-acting solution for constant-rate (Towler, Eq. 6.20; from Lee, 1982)
          rw_lim = np.sqrt((.0002637 * perm * tmin) / (100 * poro * mu * ct))
          if rw < rw_lim:
              print('valid')
          else:
              print('invalid')
      if time == 'finite':
          # Finite-acting solution for constant-rate (Towler, Eq. 6.19; from Lee, 1982)
          rw_lim = np.sqrt((.0002637 * perm * tmin) / (25 * poro * mu * ct))
          rD2 = (re / rw)**2
          if rw < rw_lim and rD2 > 1:
              print('valid')
          else:
              print('invalid')     
  if solver == 'constant_pressure':
      if time == 'infinite':
          # Infinite-acting solution for constant-pressure (Towler, Eq. 6.42, 6.43; from Edwardson et al, 1962)
          rw_lim = np.sqrt((.0002637 * perm * tmin) / (.01 * poro * mu * ct))
          if rw < rw_lim:
              print('valid')
          else:
              print('invalid')
      if time == 'finite':
          # Finite-acting solution for constant-pressure (Towler, Appendix A-4)
          print('valid')      

def simulate_multirate_test(p_initial, t_step, t_change, q_change,
                            re, rw, perm, poro, mu, ct, Bo, h):
  """
  Simulate the Multiple Constant Rate Test Started from 0th Hour 
  Based on Superposition Principle
  """
  import numpy as np
  import matplotlib.pyplot as plt
  import matplotlib.patches as mpl_patches
  
  # calculate finite-acting time
  t_finite_acting = time_finite_acting(perm, poro, mu, ct, rw, re)

  # produce time array
  t_end = t_change[-1]
  time = np.arange(0, t_end+1, t_step)

  # calculate dimensionless radius
  rD = re / rw

  # calculate delta rate (Δq)
  t_change = np.append(0, t_change)
  delta_q = [j-i for i, j in zip(q_change[:-1], q_change[1:])]
  delta_q = np.concatenate((np.array([0, q_change[0]]), delta_q))

  # create rate step profile
  tmax = t_change[-1] + 1
  t = []
  q = []
  pwf = []

  for i in range(len(time)):  
      for j in range(0, len(t_change)-1):
          if time[i] > t_change[j] and time[i] <= t_change[j+1]:
              # produce t and q profile
              t.append(time[i])
              q.append(q_change[j])
              
              # calculate dimensionless time tD (tD1, tD2, ..., tDn) at each time
              tn = time[i] - t_change[:j+1] # is an array   
              tD = time_dimensionless(perm, tn, poro, mu, ct, rw)
              
              # calculate dimensionless pressure pD at each time
              pD = []
              for k in range(len(tD)):
                  _ = pressure_dimensionless(rD, tD[k])
                  # _ = pd(rD, tD[k])
                  pD.append(_)
              
              # calculate final pressure after superposition
              delta_qn = delta_q[1:j+2] # is an array 
              
              pwf_ = pressure_multirate(pD, delta_qn, p_initial, Bo, mu, perm, h)
              pwf.append(pwf_)       

  t, q, pwf = np.append(0, t), np.append(q_change[0], q), np.append(p_initial, pwf)

  # plot well rate and flowing pressure profile
  plt.figure(figsize=(17,5))

  ## output the finite-acting time into the plot
  labels = []
  labels.append("Time @ Finite-acting = {} hours".format(np.round(t_finite_acting, 2)))

  handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                  lw=0, alpha=0)] * 1

  ## plot rate
  plt.subplot(1,2,1)
  plt.step(t, q, color='blue')
  plt.title('Well Rate Profile', size=20, pad=15)
  plt.xlim(0, t_end)
  plt.ylim(ymax=max(q)+200)
  plt.xlabel('Time (hours)'); plt.ylabel('Rate (STB/D)')

  plt.legend(handles, labels, loc='upper right', fontsize=12, 
              fancybox=True, framealpha=0.7, 
              handlelength=0, handletextpad=0) 

  ## plot BHFP
  plt.subplot(1,2,2)
  # t = np.arange(len(pwf))
  plt.plot(t, pwf, color='red')
  plt.title('Well Flowing Pressure Profile', size=20, pad=15)
  plt.xlim(0, t_end)
  plt.xlabel('Time (hours)'); plt.ylabel('BHFP (psia)')

  plt.show()
  
def simulate_multipressure_test(p_initial, t_step, t_change, p_change,
                                re, rw, perm, poro, mu, ct, Bo, h):
  """
  Simulate the Multiple Constant Borehole Flowing Pressure (BHFP) Test 
  Based on Superposition Principle
  """
  import numpy as np
  import matplotlib.pyplot as plt
  import matplotlib.patches as mpl_patches
  
  # calculate finite-acting time
  t_finite_acting = time_finite_acting(perm, poro, mu, ct, rw, re)  

  # produce time array
  t_end = t_change[-1]
  time = np.arange(0, t_end+1, t_step)

  # calculate dimensionless radius
  rD = re / rw

  # calculate delta rate (Δq)
  t_change = np.append(0, t_change)
  pi_min_p0 = p_initial - p_change[0]
  delta_p = [i-j for i, j in zip(p_change[:-1], p_change[1:])]
  delta_p = np.concatenate((np.array([0, pi_min_p0]), delta_p))

  # create rate step profile
  tmax = t_change[-1] + 1
  t = []
  pwf = []
  q = []

  for i in range(len(time)):  
      for j in range(0, len(t_change)-1):
          if time[i] > t_change[j] and time[i] <= t_change[j+1]:
              # produce t and p profile
              t.append(time[i])
              pwf.append(p_change[j])
              
              # calculate dimensionless time tD (tD1, tD2, ..., tDn) at each time
              tn = time[i] - t_change[:j+1] # is an array   
              tD = time_dimensionless(perm, tn, poro, mu, ct, rw)
              
              # calculate dimensionless rate qD at each time
              qD = []
              for k in range(len(tD)):
                  _ = rate_dimensionless(rD, tD[k])
                  # _ = qd(rD, tD[k])
                  qD.append(_)
              
              # calculate final rate after superposition
              delta_pn = delta_p[1:j+2] # is an array 
              
              q_ = rate_multipressure(qD, delta_pn, Bo, mu, perm, h)
              q.append(q_)     

  # plot flowing pressure and well rate profile
  plt.figure(figsize=(17,5))

  ## output the finite-acting time into the plot
  labels = []
  labels.append("Time @ Finite-acting = {} hours".format(np.round(t_finite_acting, 2)))

  handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                  lw=0, alpha=0)] * 1

  ## plot BHFP
  plt.subplot(1,2,1)
  plt.step(t, pwf, color='red')
  plt.title('Well Flowing Pressure Profile', size=20, pad=15)
  plt.xlim(0, t_end)
  plt.ylim(ymax=max(pwf)+200)
  plt.xlabel('Time (hours)'); plt.ylabel('Pressure (psia)')

  plt.legend(handles, labels, loc='upper right', fontsize=12, 
              fancybox=True, framealpha=0.7, 
              handlelength=0, handletextpad=0) 

  ## plot rate
  plt.subplot(1,2,2)
  # t = np.arange(len(pwf))
  plt.plot(t, q, color='blue')
  plt.title('Well Rate Profile', size=20, pad=15)
  plt.xlim(0, t_end)
  plt.xlabel('Time (hours)'); plt.ylabel('Rate (STB/D)')

  plt.show()
  
            
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

  # Access to Ei-function table
  Ei_table = Ei_table = np.loadtxt("https://raw.githubusercontent.com/yohanesnuwara/reservoir-engineering/master/Appendix%20A.%20Values%20and%20Tabulations/Table%20A-1-Values%20for%20Exponential%20Integral.txt") 

  r = rw + distance
  t_finite_acting = time_finite_acting(re, rw, poro, mu_oil, ct, k)

  if time > 0 and time < t_finite_acting:
    """Time behaving infinite acting"""
    td = time_dimensionless(time, rw, poro, mu_oil, ct, k)
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
      td = time_dimensionless(time, r, poro, mu_oil, ct, k)
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
    td = time_dimensionless(time, rw, poro, mu_oil, ct, k)
    pd = np.nan
    pwf = pi
    
  elif time >= t_finite_acting:
    """Time behaving finite acting"""
    td = time_dimensionless(time, rw, poro, mu_oil, ct, k)
    if r==rw:
      # Your distance is at wellbore
      if td > 25:
        # Calculate dimensionless radius at reservoir outer boundary (r=re)
        r_eD = radius_dimensionless(re, rw)
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
