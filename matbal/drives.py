"""
Program for Drive Indices Calculation and Energy Plot

@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

class saturated_nonvolatile_totaloil():
  """
  Saturated Non-Volatile Oil Reservoir Drive Indices Calculation
  (for PVT that only have Bg, Bto, and Rsi; No Rs)
  """
  def calculate_params(self, p, Bg, Bto, Rsi, Np, Gp, cf, cw, swi):
    """Material Balance parameters"""
    pi, Boi, Bgi = p[0], Bto[0], Bg[0]

    # Formation expansion factor
    Efw = ((cf + cw * swi) / (1 - swi)) * (pi - p)

    # Oil expansion factor
    Eo = Bto - Boi

    # Gas expansion factor
    Btg = Bg  
    Eg = Btg - Bgi

    # Reservoir voidage
    F = (Np * (Bto - (Rsi * Bg))) + (Gp * Bg)
    return F, Efw, Eo, Eg

  def indices(self, F, Efw, Eo, Eg, Nfoi, Gfgi, Boi, Bgi, We, Bw, Wp, Wi):
    """Calculate Drive Indices"""
    import numpy as np
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)    

    # Depletion Drive Index
    Idd = (Nfoi * Eo) / F 
    # Segregation Drive Index
    Isd = (Gfgi * Eg) / F 
    # Formation Drive Index
    Ifd = (((Nfoi * Boi) + (Gfgi * Bgi)) * Efw) / F 
    # Water Drive Index
    deltaW = We - (Bw * Wp)
    Iwd = deltaW / F
    # Water Injection Index
    Iwi = Wi * Bw / F

    return Idd, Isd, Ifd, Iwd, Iwi

def energy_plot(t, Idd, Isd, Ifd, Iwi):
  """Energy Plot of the Drive Indices"""
  import matplotlib.pyplot as plt
  t = t[1:]
  Idd_curve = Idd[1:]
  Isd_curve = Idd[1:] + Isd[1:]
  Ifd_curve = Idd[1:] + Isd[1:] + Ifd[1:]
  Iwi_curve = Idd[1:] + Isd[1:] + Ifd[1:] + Iwi[1:]
  
  plt.figure(figsize=(12,7))

  plt.title('Energy Plot', size=20, pad=15)
  plt.plot(t, Idd_curve, color='green', linewidth=1)
  plt.plot(t, Isd_curve, color='orange', linewidth=1)
  plt.plot(t, Iwi_curve, color='grey', linewidth=1)
  plt.plot(t, Ifd_curve, color='brown', linewidth=1)

  y0 = np.full(len(t), 0)
  y1 = np.full(len(t), 1)

  wd = plt.fill_between(t, y1, Iwi_curve, color='blue', 
                        alpha=.8) 
  wi = plt.fill_between(t, Iwi_curve, Ifd_curve, color='grey', 
                        alpha=.8) 
  fd = plt.fill_between(t, Ifd_curve, Isd_curve, color='brown', 
                        alpha=.8) 
  sd = plt.fill_between(t, Isd_curve, Idd_curve, color='orange', 
                        alpha=.8) 
  dd = plt.fill_between(t, Idd_curve, y0, color='green', 
                        alpha=.8) 
  plt.legend(handles=[wd, wi, fd, sd, dd], 
             labels=['Water Drive', 'Water Injection', 'Formation Drive',
                     'Segregation Drive', 'Depletion Drive'])

  plt.xlim(min(t), max(t))
  plt.ylim(0,1) 
  plt.xlabel('Year')
  plt.ylabel('Index') 

  plt.show()
