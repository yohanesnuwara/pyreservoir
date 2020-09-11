"""
Material Balance Plots
@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

class drygas():
    """
    Dry-Gas Material Balance Plot
    """
    def calculate_params(self, p, Bg, Gp, cf, cw, swi):
        """Calculate Material Balance Paramaters for Dry-Gas Reservoir"""
        import numpy as np

        pi = p[0]
        Bgi = Bg[0]

        # total gas FVF equals the gas FVF itself (for dry-gas)
        Btg = Bg

        # calculate Efw
        Efw = ((cf + cw * swi) / (1 - swi)) * (pi - p)

        F = []; Eg = []
        for i in range(len(p)):
            F_ = Bg[i] * Gp[i]
            Eg_ = Btg[i] - Bgi
            F.append(F_); Eg.append(Eg_)

        F = np.array(F); Eg = np.array(Eg)
        return F, Btg, Efw, Eg

    def plot(self, p, z, Gp, F, Btg, Efw, Eg):
        """Create Material Balance Plots for Dry-Gas Reservoir"""
        import matplotlib.pyplot as plt
        from scipy.optimize import curve_fit

        # plot attributes
        title_size = 12
        title_pad = 10

        # linear function for curve-fit
        def linear_zero_intercept(x, m):
            y = m * x
            return y

        def linear_with_intercept(x, m, c):
            y = m * x + c
            return y

        # Plot 1: F vs Eg
        plt.subplot(3,2,1)
        x1, y1 = Eg, F
        plt.plot(x1, y1, '.-')
        plt.title('Plot 1: F vs Eg', size=title_size, pad=title_pad)
        plt.xlabel('Eg (RB/scf)')
        plt.ylabel('F (res ft3)')

        ## curve-fitting to calculate the slope as OGIP
        x1_norm = x1 / max(x1) # normalize x
        y1_norm = y1 / max(y1) # normalize y
        popt, pcov = curve_fit(linear_zero_intercept, x1_norm, y1_norm)

        m = popt[0]
        Gfgi = m * max(y1) / max(x1) # denormalize the slope, hence the OGIP

        ## plot the regression line
        x1_fit = np.linspace(min(x1), max(x1), 5)
        y1_fit = linear_zero_intercept(x1_fit, Gfgi)
        plt.plot(x1_fit, y1_fit, label='{} MMSCF'.format(np.round(Gfgi * 1E-6, 3)))
        plt.legend()

        # Plot 2: p/z vs Gp
        plt.subplot(3,2,2)
        x2, y2 = Gp, (p / z)
        plt.plot(x2, y2, '.-')
        plt.title('Plot 2: p/z vs Gp', size=title_size, pad=title_pad)
        plt.xlabel('Gp (scf)')
        plt.ylabel('p/z (psia)')

        ## curve-fitting to calculate the slope as OGIP
        x2_norm = x2 / max(x2) # normalize x
        y2_norm = y2 / max(y2) # normalize y
        popt, pcov = curve_fit(linear_with_intercept, x2_norm, y2_norm)

        m, c = popt[0], popt[1]
        Gfgi = (-c / m) * max(x2) # OGIP is the intercept at x-axis, and denormalized
        m = m * max(y2) / max(x2) # denormalize the slope
        c = c * max(y2) # denormalize the intercept

        ## plot the regression line
        x2_fit = np.linspace(min(x2), max(x2), 5)
        y2_fit = linear_with_intercept(x2_fit, m, c)
        plt.plot(x2_fit, y2_fit, label='{} MMSCF'.format(np.round(Gfgi * 1E-6, 3)))
        plt.legend()        

        # Plot 3: F/Eg vs Gp
        plt.subplot(3,2,3)
        x3, y3 = Gp, (F / Eg)
        plt.plot(x3, y3, '.-')
        plt.title('Plot 3: Waterdrive Diagnostic Plot', size=title_size, pad=title_pad)
        plt.xlabel('Gp (scf)')
        plt.ylabel('F/Eg (scf)')

        ## curve-fitting to calculate the slope as OGIP, here [1:] because NaN is removed
        x3_norm = x3[1:] / max(x3[1:]) # normalize x
        y3_norm = y3[1:] / max(y3[1:]) # normalize y
        popt, pcov = curve_fit(linear_with_intercept, x3_norm, y3_norm)

        m, c = popt[0], popt[1]
        m = m * max(y3[1:]) / max(x3[1:]) # denormalize the slope
        Gfgi = c * max(y3[1:]) # denormalize the intercept, hence the OGIP

        ## plot the regression line
        x3_fit = np.linspace(min(x3[1:]), max(x3[1:]), 5)
        y3_fit = linear_with_intercept(x3_fit, m, Gfgi)
        plt.plot(x3_fit, y3_fit, label='{} MMSCF'.format(np.round(Gfgi * 1E-6, 3)))
        plt.legend()          

        # Plot 6: F vs (Eg+Bgi*Efw)
        plt.subplot(3,2,4)
        Bgi = Btg[0]
        x6, y6 = (Eg + Bgi * Efw), F
        plt.plot(x6, y6, '.-')
        plt.title('Plot 6: F vs (Eg+Bgi*Efw)', size=title_size, pad=title_pad)
        plt.xlabel('Eg+Bgi*Efw (res ft3/scf)')
        plt.ylabel('F (res ft3)')

        ## curve-fitting to calculate the slope as OGIP
        x6_norm = x6 / max(x6) # normalize x
        y6_norm = y6 / max(y6) # normalize y
        popt, pcov = curve_fit(linear_zero_intercept, x6_norm, y6_norm)

        m = popt[0]
        Gfgi = m * max(y6) / max(x6) # denormalize the slope, hence the OGIP

        ## plot the regression line
        x6_fit = np.linspace(min(x6), max(x6), 5)
        y6_fit = linear_zero_intercept(x6_fit, Gfgi)
        plt.plot(x6_fit, y6_fit, label='{} MMSCF'.format(np.round(Gfgi * 1E-6, 3)))
        plt.legend()        

        # Plot 7: ((p/z)*(1-Efw)) vs Gp
        plt.subplot(3,2,5)
        x7, y7 = Gp, ((p / z) * (1 - Efw))
        plt.plot(x7, y7, '.-')
        plt.title('Plot 7: ((p/z)*(1-Efw)) vs Gp', size=title_size, pad=title_pad)
        plt.xlabel('Gp (scf)')
        plt.ylabel('(p/z)*(1-Efw) (psia)')

        ## curve-fitting to calculate the slope as OGIP
        x7_norm = x7 / max(x7) # normalize x
        y7_norm = y7 / max(y7) # normalize y
        popt, pcov = curve_fit(linear_with_intercept, x7_norm, y7_norm)

        m, c = popt[0], popt[1]
        Gfgi = (-c / m) * max(x7) # OGIP is the intercept at x-axis, and denormalized
        m = m * max(y7) / max(x7) # denormalize the slope
        c = c * max(y7) # denormalize the intercept

        ## plot the regression line
        x7_fit = np.linspace(min(x7), max(x7), 5)
        y7_fit = linear_with_intercept(x7_fit, m, c)
        plt.plot(x7_fit, y7_fit, label='{} MMSCF'.format(np.round(Gfgi * 1E-6, 3)))
        plt.legend()  

        plt.tight_layout(pad=1.5)
        plt.show()

        return F, Eg, Efw
