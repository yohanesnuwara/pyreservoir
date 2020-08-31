"""
Code to Calculate Water Influx (Schilthuis, VEH, Fetkovich, and Material Balance Method)

@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

def fetkovich_initial_encroachable_water(pi, ct, r_R, r_aq, h_aq, poro, theta):
  "calculate initial encroachable water"
  import numpy as np
  # r_R: reservoir size (radius of cylindrical-assumed reservoir), in ft
  # r_aq: aquifer size, in ft
  # theta: for full circle cylindrical, theta=360. if half-circle, theta=180
  Wei = (pi * ct * np.pi * ((r_aq**2) - (r_R**2)) * h_aq * poro * theta) / (5.61458 * 360)
  return(Wei)

def fetkovich_productivity_index(perm, h_aq, mu_w, r_aq, r_R, theta, flow='constant'):
  "calculate productivity index"
  import numpy as np

  if flow == 'constant':
      # mu_w: water viscosity
      J = (0.007082 * perm * h_aq * theta) / ((mu_w * (np.log(r_aq / r_R)) * 360))
      return(J)

  if flow == 'no flow':
      # mu_w: water viscosity
      J = (0.007082 * perm * h_aq * theta) / ((mu_w * (np.log(r_aq / r_R) - 0.75) * 360))
      return(J)

def fetkovich(datetime, pressure, Wei, J):
    """
    Calculate aquifer influx (We) using Fetkovich Pseudo-steady Method
    """
    import numpy as np

    "Subtracting datetimes to get time differences (how many days) IN INTEGER"
    diff = [j-i for i, j in zip(datetime[:-1], datetime[1:])]
    diff = np.array(diff)

    # convert datetime format to integer
    diffr_arr = []
    for k in range(len(diff)):
      diffr = diff[k] / np.timedelta64(1,'D')
      diffr_arr.append(float(diffr))

    # append 0 to the first index of numpy
    diffr_arr = np.append(0, diffr_arr) # now diff has same dimension with time data (a)
    delta_time = diffr_arr

    "Initial conditions"
    We = 0 # We at initial production date (NOTE: different from Wei, initial encroachable water)
    pi = pressure[0]
    pRn_min_one = pn_min_one = pi

    "Calculate aquifer influx"
    We_fetkovich = []

    for i in range(len(datetime)):
        # calculate p_Rn average, Eq 8.29
        p_Rn = 0.5 * (pRn_min_one + pressure[i])

        # update value of pRn-1 equals to current pressure
        pRn_min_one = pressure[i]

        # calculate (p_n-1 - p_Rn average), Eq 8.30
        pn_min_prn = pn_min_one - p_Rn

        # calculate delta Wen, Eq 8.30
        delta_We = (Wei / pi) * pn_min_prn * (1 - np.exp(-(J * pi * delta_time[i]) / (Wei)))

        # calculate We, Eq 8.31
        We = We + delta_We

        # update p_n-1 for the next timestep, Eq 8.32
        pn_min_one = pi * (1 - (We / Wei))

        We_fetkovich.append(We)

    return We_fetkovich
