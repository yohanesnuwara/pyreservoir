"""
Gas Condensate Material Balance

@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

import numpy as np
import matplotlib.pyplot as plt

class mbplot():
    """
    Gas Condensate Material Balance Plot
    Ideal Data
    """
    def plot1(self, Bg, Bo, Np, Gp, Gi, Rs, Rv, output=None):
        """
        Plot 1: F vs Eg

        OUTPUT:
        
        array Eg (x-axis values), array F (y-axis values), and plot
        """
        # initial conditions
        Rvi = Rv[0]
        Rsi = Rs[0]
        Bgi = Bg[0]

        Btg = ((Bg * (1 - (Rs * Rvi))) + (Bo * (Rvi - Rv))) / (1 - (Rv * Rs))  # in RB/STB
        Bto = ((Bo * (1 - (Rv * Rsi))) + (Bg * (Rsi - Rs))) / (1 - (Rv * Rs))  # in RB/scf

        F = (Np * ((Bo - (Rs * Bg)) / (1 - (Rv * Rs)))) + ((Gp - Gi) * ((Bg - (Rv * Bo)) / (1 - (Rv * Rs))))
        Eg = Btg - Bgi

        if output == 'allparams':
            return(Btg, Bto, Eg, F)

        if output == None:
            # plot
            plt.plot(Eg, F, '.')
            plt.title('Plot 1: F vs Eg')
            plt.xlabel('Eg (RB/scf)');
            plt.ylabel('F (res bbl)')
            plt.xlim(xmin=0);
            plt.ylim(ymin=0)
            plt.show()
            return(Eg, F)

    def plot2(self, Gp, p, z, output=None):
        """
        Plot 2: p/z vs Gp

        OUTPUT:
        
        array Gp (x-axis values), array p_z (y-axis values), plot
        """
        # calculate plotting parameters
        p_z = p / z

        if output == 'allparams':
            return(p_z)

        if output == None:
            # plotting
            plt.plot(Gp, p_z, '.')
            plt.title('Plot 2: p/z vs Gp')
            plt.xlabel('Gp (scf)');
            plt.ylabel('p/z (psia)')
            plt.xlim(xmin=0);
            plt.ylim(ymin=0)
            plt.show()
            return(Gp, p_z)

    def plot3(self, Bg, Bo, Np, Gp, Gi, Rs, Rv, output=None):
        """
        Plot 3: F/Eg vs Gp

        OUTPUT:
        
        array Gp (x-axis values), array F_Eg (y-axis values), and plot
        """
        # initial conditions
        Rvi = Rv[0]
        Rsi = Rs[0]
        Bgi = Bg[0]

        Btg = ((Bg * (1 - (Rs * Rvi))) + (Bo * (Rvi - Rv))) / (1 - (Rv * Rs))  # in RB/STB
        Bto = ((Bo * (1 - (Rv * Rsi))) + (Bg * (Rsi - Rs))) / (1 - (Rv * Rs))  # in RB/scf

        F = (Np * ((Bo - (Rs * Bg)) / (1 - (Rv * Rs)))) + ((Gp - Gi) * ((Bg - (Rv * Bo)) / (1 - (Rv * Rs))))
        Eg = Btg - Bgi

        N = F / Eg

        if output == 'allparams':
            return(Btg, Bto, Eg, F)
        if output == None:
            # plot
            plt.plot(Gp, N, '.')
            plt.title('Plot 3: F/Eg vs Gp')
            plt.xlabel('Gp (scf)');
            plt.ylabel('F/Eg (scf)')
            plt.xlim(xmin=0);
            plt.ylim(ymin=0)
            plt.show()
            return(Gp, N)

    def plot6(self, p, Bg, Bo, Np, Gp, Gi, cf, cw, swi, Rs, Rv, output=None):
        """
        Plot 6: F vs (Eg+Bgi*Efw)

        OUTPUT:
        
        array Eg_Bgi_Efw (x-axis values), array F (y-axis values), and plot
        """
        
        # initial conditions
        Rvi = Rv[0]
        Rsi = Rs[0]
        Bgi = Bg[0]
        pi = p[0]

        Btg = ((Bg * (1 - (Rs * Rvi))) + (Bo * (Rvi - Rv))) / (1 - (Rv * Rs))  # in RB/STB
        Bto = ((Bo * (1 - (Rv * Rsi))) + (Bg * (Rsi - Rs))) / (1 - (Rv * Rs))  # in RB/scf

        F = (Np * ((Bo - (Rs * Bg)) / (1 - (Rv * Rs)))) + ((Gp - Gi) * ((Bg - (Rv * Bo)) / (1 - (Rv * Rs))))
        Eg = Btg - Bgi

        Efw = ((cf + cw * swi) / (1 - swi)) * (pi - p)
        Eg_Bgi_Efw = Eg + Bgi * Efw


        if output == 'allparams':
            return(Efw, Btg, Bto, Eg, F)
        if output == None:
            # plot
            plt.plot(Eg_Bgi_Efw, F, '.')
            plt.title('Plot 6: Eg+Bgi*Efw vs F')
            plt.xlabel('Eg+Bgi*Efw (res ft3/scf)');
            plt.ylabel('F (res ft3)')
            plt.xlim(xmin=0);
            plt.ylim(ymin=0)
            plt.show()
            return(Eg_Bgi_Efw, F)

    def plot7(self, p, z, Gp, cf, cw, swi, output=None):
        """
        Plot 7 (p/z*(1-Efw)) vs Gp

        OUTPUT:
        
        array Gp (x-axis values), array p_z_Efw (y-axis values), plot
        """
        # initial conditions
        pi = p[0]

        p_z = p / z
        Efw = ((cf + cw * swi) / (1 - swi)) * (pi - p)

        p_z_Efw = (p / z) * (1 - Efw)

        if output == 'allparams':
            return(Efw, p_z)
        if output == None:
            # plot
            plt.plot(Gp, p_z_Efw, '.')
            plt.title('Plot 7: (p/z)*(1-Efw) vs Gp')
            plt.xlabel('Gp (scf)');
            plt.ylabel('(p/z)*(1-Efw) (psia)')
            plt.xlim(xmin=0);
            plt.ylim(ymin=0)
            plt.show()
            return(Gp, p_z_Efw)
