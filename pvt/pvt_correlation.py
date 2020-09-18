"""
Codes for gas, oil, and water PVT correlations
@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

"""
GAS
"""

def gas_pseudoprops(temp, pressure, sg, x_h2s, x_co2):
  """
  Calculate Gas Pseudo-critical and Pseudo-reduced Pressure and Temperature
  * Pseudo-critical properties
    For range: 0.57 < sg < 1.68
    (Sutton, 1985)
  * Pseudo-reduced properties
    For range: x_h2s (mol%) < 0.738; x_co2 (mol%) < 0.544; 154 < p (psia) < 7026; 40 < temp (°F) < 300 (error 0.97%)
    (Wichert and Aziz, 1972)
  """
  import numpy as np
  temp = temp + 459.67 # convert to Rankine

  # calculate pseudocritical properties (Sutton, valid for 0.57<sg<1.68)
  P_pc = 756.8 - (131.07 * sg) - (3.6 * sg**2)
  T_pc = 169.2 + (349.50 * sg) - (74 * sg**2) # in Rankine

  # calculate adjustment to pseudocritical properties for sour gas (Wiechert-Aziz, valid for x_co2<0.544 and x_h2s<0.738)
  e = (120 * (((x_h2s + x_co2)**0.9) - ((x_h2s + x_co2)**1.6))) + (15 * (x_h2s**0.5 - x_h2s**4))
  T_pc = T_pc - e # corrected T_pc
  P_pc = (P_pc * T_pc) / (T_pc - x_h2s * e * (1-x_h2s))

  # calculate pseudoreduced properties
  P_pr = pressure / P_pc
  T_pr = temp / T_pc

  return(P_pc, T_pc, P_pr, T_pr)

def gas_zfactor(T_pr, P_pr):
  """
  Calculate Gas Compressibility Factor
  For range: 0.2 < P_pr < 30; 1 < T_pr < 3 (error 0.486%)
  (Dranchuk and Aboukassem, 1975)
  """
  # T_pr : calculated pseudoreduced temperature
  # P_pr : calculated pseudoreduced pressure   
  from scipy.optimize import fsolve # non-linear solver
  import numpy as np

  a1 = 0.3265; a2 = -1.0700; a3 = -0.5339; a4 = 0.01569; a5 = -0.05165; a6 = 0.5475
  a7 = -0.7361; a8 = 0.1844; a9 = 0.1056; a10 = 0.6134; a11 = 0.7210

  def f(y):
    rho_pr, z = y
    c1 = a1 + (a2/T_pr) + (a3/(T_pr**3))+ (a4/(T_pr**4))+ (a5/(T_pr**5))
    c2 = a6 + (a7/T_pr) + (a8/(T_pr**2))
    c3 = a9*((a7/T_pr) + (a8/(T_pr**2)))
    c4 = (a10)*(1+(a11*(rho_pr**2)))*((rho_pr**2)/(T_pr**3))*(np.exp(-a11*(rho_pr**2)))

    f1 = z + (c3*(rho_pr**5)) - (c2*(rho_pr**2)) - (c1*(rho_pr**1)) - c4 - 1
    f2 = rho_pr - ((0.27 * P_pr) / (z * T_pr))
    return[f1, f2]

  solve = fsolve(f, [1, 1]) # initial guess
  return(solve[0], solve[1]) # result is density, z-factor

def gas_density(temp, pressure, sg, z):
  temp = temp + 459.67
  R = 10.732 # gas constant in (ft3*psi)/(lb-mol*R) 
  rhogas = (28.97 * sg * pressure) / (z * R * temp)
  return rhogas  

def gas_fvf(z, temp, pressure):
  """
  Calculate Gas FVF
  For range: this is not a correlation, so valid for infinite intervals
  """
  temp = temp + 459.67
  Bg = 0.0282793 * z * temp / pressure 
  return(Bg)

def gas_fvf2(unit='unit1', z=0.8, temp=186, pressure=2000):
  """
  Gas FVF calculated in other units
  unit: choice of units (unit1: RB/scf, unit2: res m3/std m3)
  for unit1, inputs temp in Rankine (Fahrenheit + 460), pressure in psia or psig
  for unit2, inputs temp in Kelvin, pressure in psia or psig
  """
  if unit == 'unit1':
    return(0.00503676 * z * temp / pressure) 
  if unit == 'unit2':
    return(0.350958 * z * temp / pressure)

def gas_mu(temp, rhogas, sg):
  """
  Calculate Gas Viscosity 
  For gas with CO2 and N2 composition
  For range: 100 < temp (°F) < 340; 0.9 < x_CO2 (mol%) < 3.2; x_N2 (mol%) < 4.8 (std 2.7-9.0%)
  (Lee et al, 1996)
  """
  import numpy as np
  temp = temp + 459.67
  Mg = 28.97 * sg
  rhogas_lee = rhogas * 0.0160185 # lbm/ft3 converted to gas density unit of Lee et al (g/cm3)
  K = ((0.00094 + 2E-06)*(temp**1.5)) / (209 + 19*Mg + temp)
  x = 3.5 + (986 / temp) + (0.01 * Mg)
  y = 2.4 - 0.2*x  
  viscogas = K * np.exp(x * (rhogas_lee**y))
  return viscogas

def gas_compressibility(T_pr, P_pr, rho_pr, z, P_pc):
  """
  Calculate Gas Isothermal Compressibility
  For range: unspecified
  (Trube, 1957; Mattar, 1975)
  """
  import numpy as np

  a1 = 0.3265; a2 = -1.0700; a3 = -0.5339; a4 = 0.01569; a5 = -0.05165; a6 = 0.5475
  a7 = -0.7361; a8 = 0.1844; a9 = 0.1056; a10 = 0.6134; a11 = 0.7210

  do = ((a1 + (a2/T_pr) + (a3/T_pr**3) +(a4/T_pr**4) + (a5/T_pr**5)) * rho_pr) + \
      (2 * ((a6 + (a7/T_pr) + (a8/T_pr**2))) * rho_pr**2) - \
      (5 * a9 * (((a7/T_pr) + (a8/T_pr**2))) * rho_pr**4) + (1 + (a11 * rho_pr**2) - (a11 * rho_pr**2)**2) \
      * ((2 * a10 * rho_pr / T_pr**3)*np.exp(-a11 * rho_pr**2))

  c_pr_analytical = (1 / P_pr) - ((0.27 / (z**2 * T_pr)) * (do / (1 + ((rho_pr / z) * do))))
  cgas_analytical = c_pr_analytical / P_pc
  return(cgas_analytical)           

"""
OIL
"""

def oil_pbubble(Rsb, sg2, api, temp2):
  """
  Calculate Oil Bubble-Point Pressure
  For range: 20 < Rsb (scf/STB) < 2,070; 0.56 < sg < 1.18; 16 < api < 58; 70 < temp (°F) < 295 (err=0.7%)
  (Vazquez and Beggs, 1980)
  """
  import numpy as np
  # c1, c2, c3 coefficient from Vazquez-Beggs
  if api <=30:
    c1 = 0.0362
    c2 = 1.0937
    c3 = 25.7240
  if api > 30:
    c1 = 0.0178
    c2 = 1.187
    c3 = 23.9310

  P_bubble_vaz = (Rsb / (c1 * sg2 * np.exp((c3 * api)/(temp2 + 459.67))))**(1 / c2) # convert temp to Rankine
  return(P_bubble_vaz)

def oil_fvf(P_bubble, api, Rsb, sg2, temp2, pressure2):
  """
  Calculate Oil FVF
  * Above bubble-point pressure
    For range: unspecified
    (Vazquez and Beggs, 1980)
  * At and bubble-point pressure
    For range: unspecified
    (Levitan and Murtha, 1999)
  """

  import numpy as np
  # FVF of oil at bubblepoint pressure using Levitan-Murtha
  so = 141.5 / (api + 131.5)
  Bo_bubble = 1 + ((0.0005 * Rsb) * ((sg2 / so)**0.25)) + ((0.0004*(temp2- 60)) / (so * sg2)) # temp in def F

  Bo_array = []

  if pressure2 < P_bubble: # use Vazquez-Beggs
    if api <= 30:
      # use Vazquez-Beggs 
      c1 = 0.0362
      c2 = 1.0937
      c3 = 25.7240
      c4 = 4.677E-4
      c5 = 1.751E-5
      c6 = -1.811E-8
    if api > 30:
      c1 = 0.0178
      c2 = 1.187
      c3 = 23.9310
      c4 = 4.670E-4
      c5 = 1.100E-5
      c6 = 1.337E-9
    Rsc = (pressure2**c2) * c1 * sg2 * np.exp((c3 * api) / (temp2 + 459.67))
    Bo = 1 + (c4 * Rsc) + (c5 * (temp2 - 60) * (api / sg2)) + (c6 * Rsc *(temp2 - 60) * (api / sg2)) # temp in deg F
  if pressure2 == P_bubble:
    # use Levitan-Murtha
    Bo = Bo_bubble
  if pressure2 > P_bubble:
    # Calculate oil compressibility first using Levitan-Murtha
    coil = ((5 * Rsb) + (17.2 * temp2) - (1180 * sg2) + (12.61 * api) - 1433) / (1E+05 * pressure2)
    # Calculate Bo using Levitan-Murtha
    Bo = Bo_bubble * np.exp(coil * (P_bubble - pressure2))

  return Bo
  
def oil_mu(pressure2, P_bubble, sg2, api, temp2, Rsb):
  """
  Calculate Oil Viscosity
  * Below and at bubble-point pressure
    For range: 0 < p (psia) < 5,250; range sg unspecified; 16 < api < 58; 70 < temp (°F) < 295; 20 < Rs (scf/STB) < 2,070 (err=1.83%)
    (Beggs and Robinson, 1975; Chew and Connally, 1959)
  * Above bubble-point pressure
    For range: 126 < p (psia) < 9,500; 0.511 < sg < 1.351; 15.3 < api < 59.5; range temp unspecified; 9.3 < Rs (scf/STB) < 2199 (err=7.54%)
    (Vazquez and Beggs, 1980)
  """
  # Calculate viscosity of oil
  import numpy as np

  mu_oil_array = []

  if pressure2 <= P_bubble:
    if api <=30:
      c1 = 0.0362
      c2 = 1.0937
      c3 = 25.7240
    if api > 30:
      c1 = 0.0178
      c2 = 1.187
      c3 = 23.9310

    # use Beggs and Robinson
    # valid for: 0 < pressure < 5250 psig, 70 < temp < 295 F, 20 < Rs < 2070 scf/STB, 16 < api < 58 API 
    x = (temp2**(-1.163)) * np.exp(6.9824 - (0.04658 * api))
    mu_dead_oil = 10**x - 1
    Rs = (pressure2**c2) * c1 * sg2 * np.exp((c3 * api) / (temp2 + 459.67)) # gas-oil ratio at any pressure BELOW BUBBLEPOINT using Vazquez-Beggs
    a = 10.715 * ((Rs + 100)**(-0.515))
    b = 5.44 * ((Rs + 150)**(-0.338))
    mu_live_oil = a * (mu_dead_oil**b)

  if pressure2 > P_bubble:
    if api <=30:
      c1 = 0.0362
      c2 = 1.0937
      c3 = 25.7240
    if api > 30:
      c1 = 0.0178
      c2 = 1.187
      c3 = 23.9310

    # use Vazquez and Beggs
    # valid for: 126 < pressure < 9500 psig, 9.3 < Rs < 2199 scf/STB, 15.3 < api < 59.5 API, 0.511 < sg < 1.351 

    # compute oil viscosity at bubblepoint first
    x_bubble = (temp2**(-1.163)) * np.exp(6.9824 - (0.04658 * api))
    mu_dead_oil_bubble = 10**x_bubble - 1
    
    Rsb = (P_bubble**c2) * c1 * sg2 * np.exp((c3 * api) / (temp2 + 459.67)) # gas-oil ratio at any pressure BELOW BUBBLEPOINT using Vazquez-Beggs
    a_bubble = 10.715 * ((Rsb + 100)**(-0.515))
    b_bubble = 5.44 * ((Rsb + 150)**(-0.338))
    
    mu_live_oil_bubble = a_bubble * (mu_dead_oil_bubble**b_bubble)

    m = 2.6 * (pressure2**1.187) * np.exp(-11.513 - (8.98E-05 * pressure2))
    mu_live_oil = mu_live_oil_bubble * ((pressure2 / P_bubble)**m)

  return mu_live_oil

def oil_compressibility(pressure2, P_bubble, temp2, api, Rsb, sg2):
  """
  Calculate Oil Isothermal Compressibility
  * Below bubble-point pressure
    For range: unspecified
    (McCain, 1988)
  * Above and at bubble-point pressure
    For range: unspecified
    (Vazquez and Beggs, 1980)
  """
  import numpy as np
  from math import e

  # oil isothermal compressibility

  coil_array = []

  if pressure2 < P_bubble:
    # use McCain
    ln_coil = -7.573 - (1.45 * np.log(pressure2)) - (0.383 * np.log(P_bubble)) + (1.402 * np.log(temp2)) + (0.256 * np.log(api)) + (0.449 * np.log(Rsb))  
    coil = np.exp(ln_coil)
  if pressure2 >= P_bubble:
    # use Vazquez-Beggs
    coil = ((5 * Rsb) + (17.2 * temp2) - (1180 * sg2) + (12.61 * api) - 1433) / (1E+05 * pressure2)

  return coil


def gasoilratio(pressure2, P_bubble, sg2, api, temp2):
  """
  Calculate Solution Gas-Oil Ratio in Oil Phase
  * Below Bubble-Point
    For range: unspecified
    (Vazquez and Beggs, 1980)
  * At and Above Bubble-Point 
    is constant 
    Approach used here: Using Vazquez and Beggs (1980) to calculate Rs at bubble-point pressure
  """
  import numpy as np
  Rs_array = []

  if pressure2 < P_bubble:
    # Using Vazquez and Beggs
    if api <=30:
      c1 = 0.0362
      c2 = 1.0937
      c3 = 25.7240
    if api > 30:
      c1 = 0.0178
      c2 = 1.187
      c3 = 23.9310
    Rs = (pressure2**c2) * c1 * sg2 * np.exp((c3 * api) / (temp2 + 459.67)) 
    
  if pressure2 >= P_bubble:
    # Because Rs will be constant above BB
    pressure2 = P_bubble
    # Calculate with Vazquez and Beggs
    if api <=30:
      c1 = 0.0362
      c2 = 1.0937
      c3 = 25.7240
    if api > 30:
      c1 = 0.0178
      c2 = 1.187
      c3 = 23.9310
    Rs = (pressure2**c2) * c1 * sg2 * np.exp((c3 * api) / (temp2 + 459.67))     

  return Rs

"""
WATER
"""

def waterfvf(temp, p):
  "Water FVF (Bw)"
  # temp in Fahrenheit
  # p pressure in psia
  Vwp = (-1.95301E-9 * p * temp) - (1.72834E-13 * (p**2) * temp) - (3.588922E-7 * p) - (2.25341E-10 * p**2)
  Vwt = (-1.001E-2) + (1.33391E-4 * temp) + (5.50654E-7 * temp**2)
  Bw = (1 + Vwt) * (1 + Vwp)
  return(Bw)
