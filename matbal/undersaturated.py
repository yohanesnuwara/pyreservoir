"""
Undersaturated (Volatile and Non-volatile) Oil Material Balance

@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

import numpy as np
import matplotlib.pyplot as plt

class mbplot():
    """
    Undersaturated Oil Material Balance Plot
    Ideal Data
    """
    def plot1(self, p, Bg, Bo, Np, Gp, Gi, cf, cw, swi, Rs, Rv, output=None):
        """
        Plot 1: F vs Eo+(Bti* Efw)
        """
        # initial conditions
        pi = p[0]
        Boi = Bo[0]
        Rsi = Rs[0]

        Bto = []
        F = []
        for i in range(len(p)):

            if Rv[i] == 0:
                # reservoir is non-volatile undersaturated
                Bto_ = Bo[i] + Bg[i] * (Rsi - Rs[i])
                F_ = Np[i](Bo[i] - Rs[i] * Bg[i]) + ((Gp[i] - Gi[i]) * Bg[i])
            Bto.append(Bto_)
            F.append(F_)

            if Rv[i] != 0:
                # reservoir is volatile undersaturated
                Bto_ = ((Bo[i] * (1 - (Rv[i] * Rsi))) + (Bg[i] * (Rsi - Rs[i]))) / (1 - (Rv[i] * Rs[i]))
                F_ = (Np * ((Bo - (Rs * Bg)) / (1 - (Rv * Rs)))) + ((Gp[i] - Gi[i]) * ((Bg - (Rv * Bo)) / (1 - (Rv * Rs))))
            Bto.append(Bto_)
            F.append(F_)

        Bto = np.array(Bto)
        F = np.array(F)

        # calculate Eo+(Boi*Efw)
        Efw = ((cf + (cw * swi)) / (1 - swi)) * (pi - p)
        Eo = Bto - Boi

        Eo_Boi_Efw = Eo + Boi * Efw

        if output == 'allparams':
            return(Bto, Efw, Eo, F)
        if output == None:
            # plot
            plt.plot(Eo_Boi_Efw, F, '.')
            plt.title('Plot 1: F vs Eo+(Bti* Efw)')
            plt.xlim(xmin=0);
            plt.ylim(ymin=0)
            plt.xlabel('Eo+(Bti*Efw) (RB/STB)')
            plt.ylabel('F (res bbl)')
            plt.show()
            return(Eo_Boi_Efw, F)

    def plot2(self, p, Bg, Bo, Np, Gp, cf, cw, swi, Rs, Rv, output=None):
        """
        Plot 2: F/Eo+(Bti* Efw) vs Np
        """
        # initial conditions
        pi = p[0]
        Boi = Bo[0]
        Rsi = Rs[0]

        Bto = []
        F = []
        for i in range(len(p)):

            if Rv[i] == 0:
                # reservoir is non-volatile undersaturated
                Bto_ = Bo[i] + Bg[i] * (Rsi - Rs[i])
                F_ = Np[i](Bo[i] - Rs[i] * Bg[i]) + ((Gp[i] - Gi[i]) * Bg[i])
            Bto.append(Bto_)
            F.append(F_)

            if Rv[i] != 0:
                # reservoir is volatile undersaturated
                Bto_ = ((Bo[i] * (1 - (Rv[i] * Rsi))) + (Bg[i] * (Rsi - Rs[i]))) / (1 - (Rv[i] * Rs[i]))
                F_ = (Np * ((Bo - (Rs * Bg)) / (1 - (Rv * Rs)))) + ((Gp[i] - Gi[i]) * ((Bg - (Rv * Bo)) / (1 - (Rv * Rs))))
            Bto.append(Bto_)
            F.append(F_)

        Bto = np.array(Bto)
        F = np.array(F)

        # calculate Eo+(Boi*Efw)
        Efw = ((cf + (cw * swi)) / (1 - swi)) * (pi - p)
        Eo = Bto - Boi

        Eo_Boi_Efw = Eo + Boi * Efw

        N = F / Eo_Boi_Efw

        if output == 'allparams':
            return(Bto, Efw, Eo, F)
        if output == None:
            # plot
            plt.plot(Np, N, '.')
            plt.title('Plot 2: F/Eo+(Bti* Efw) vs Np')
            plt.xlim(xmin=0);
            plt.ylim(ymin=0)
            plt.xlabel('Np (STB)')
            plt.ylabel('F/Eo+(Bti* Efw) (STB)')
            plt.show()
            return(Np, N)
