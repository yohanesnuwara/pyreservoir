"""
Correlation Models for Gas Properties
@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

from scipy.optimize import fsolve  # non-linear solver
import numpy as np

def density():
    return()

def pseudoprops(pressure, temp, sg, x_h2s, x_co2):
    """
    Pseudoproperties of gas
    """

    temp = temp + 459.67  # convert to Rankine

    # calculate pseudocritical properties (Sutton, valid for 0.57<sg<1.68)
    P_pc = 756.8 - (131.07 * sg) - (3.6 * sg ** 2)
    T_pc = 169.2 + (349.50 * sg) - (74 * sg ** 2)  # in Rankine

    # calculate adjustment to pseudocritical properties for sour gas (Wiechert-Aziz, valid for x_co2<0.544 and x_h2s<0.738)
    e = (120 * (((x_h2s + x_co2) ** 0.9) - ((x_h2s + x_co2) ** 1.6))) + (15 * (x_h2s ** 0.5 - x_h2s ** 4))
    T_pc_corr = T_pc - e  # corrected T_pc
    P_pc_corr = (P_pc * T_pc_corr) / (T_pc - x_h2s * e * (1 - x_h2s))

    # calculate pseudoreduced properties
    P_pr = pressure / P_pc_corr
    T_pr = temp / T_pc_corr

    return(P_pr, T_pr)

class compressibility_factor():
    """
    Compressibility Factor for Gas
    """
    def dranchuk_aboukassem(self, P_pr, T_pr):
        """
        Dranchuk and Aboukassem (1975)

        Input:
        T_pr : calculated pseudoreduced temperature
        P_pr : calculated pseudoreduced pressure
        """

        a1 = 0.3265;
        a2 = -1.0700;
        a3 = -0.5339;
        a4 = 0.01569;
        a5 = -0.05165;
        a6 = 0.5475
        a7 = -0.7361;
        a8 = 0.1844;
        a9 = 0.1056;
        a10 = 0.6134;
        a11 = 0.7210

        def f(y):
            rho_pr, z = y
            c1 = a1 + (a2 / T_pr) + (a3 / (T_pr ** 3)) + (a4 / (T_pr ** 4)) + (a5 / (T_pr ** 5))
            c2 = a6 + (a7 / T_pr) + (a8 / (T_pr ** 2))
            c3 = a9 * ((a7 / T_pr) + (a8 / (T_pr ** 2)))
            c4 = (a10) * (1 + (a11 * (rho_pr ** 2))) * ((rho_pr ** 2) / (T_pr ** 3)) * (np.exp(-a11 * (rho_pr ** 2)))

            f1 = z + (c3 * (rho_pr ** 5)) - (c2 * (rho_pr ** 2)) - (c1 * (rho_pr ** 1)) - c4 - 1
            f2 = rho_pr - ((0.27 * P_pr) / (z * T_pr))
            return [f1, f2]

        solve = fsolve(f, [1, 1])  # initial guess

        rho_pr = solve[0]
        z = solve[1]

        return(z)

# # testing
# pressure = 2010 # in psia
# temp = 75 # in deg F
# sg = 0.7 # specific gravity
# x_h2s = 0.07 # mole fraction
# x_co2 = 0.1

# P_pr, T_pr = pseudoprops(pressure, temp, sg, x_h2s, x_co2)
# z = compressibility_factor.dranchuk_aboukassem(compressibility_factor(), P_pr, T_pr)
# print(z)

class fvf():
    def 

class viscosity():
    def

class isocompressibility():
    def
   
