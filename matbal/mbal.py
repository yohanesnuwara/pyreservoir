"""
Material Balance Plots
@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

def initial_hydrocarbon_in_place(Nfoi, Gfgi, Rv, Rs):
  """
  Calculate OOIP and OGIP from Nfoi and Gfgi
  And output the result to labels in the plot
  """
  Rvi, Rsi = Rv[0], Rs[0]
  OOIP = Nfoi + Gfgi * Rvi
  OGIP = Gfgi + Nfoi * Rsi

  labels = []
  labels.append("Nfoi = {0:.4g} STB".format(Nfoi))
  labels.append("Gfgi = {0:.4g} SCF".format(Gfgi))
  labels.append("OOIP = {0:.4g} STB".format(OOIP))
  labels.append("OGIP = {0:.4g} SCF".format(OGIP))

  handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                  lw=0, alpha=0)] * 4
  return labels, handles, OOIP, OGIP   

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
        import numpy as np
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

class gascondensate():
    """
    Gas-Condensate Material Balance Plot
    """
    def calculate_params(self, p, pdew, Bg, Bo, Np, Gp, Gi, cf, cw, swi, Rs, Rv):
        """Calculate Material Balance Paramaters for Gas-Condensate Reservoir"""
        import numpy as np
        pi = p[0]
        Rvi = Rv[0]
        Bgi = Bg[0]

        # calculate Efw
        Efw = ((cf + cw * swi) / (1 - swi)) * (pi - p)

        # calculate F and Btg
        F = []; Btg = []; Eg = []
        for i in range(len(p)):
            if p[i] >= pdew:
                # gas-condensate above dewpoint pressure
                F_ = Bg[i] * Gp[i]
                Btg_ = Bg[i]
                Eg_ = Btg_ - Bgi

            if p[i] < pdew:
                # gas-condensate below dewpoint pressure
                F_ = (Np[i] * ((Bo[i] - (Rs[i] * Bg[i])) / (1 - (Rv[i] * Rs[i])))) + ((Gp[i] - Gi[i]) * ((Bg[i] - (Rv[i] * Bo[i])) / (1 - (Rv[i] * Rs[i]))))
                Btg_ = ((Bg[i] * (1 - (Rs[i] * Rvi))) + (Bo[i] * (Rvi - Rv[i]))) / (1 - (Rv[i] * Rs[i]))  # in RB/STB
                Eg_ = Btg_ - Bgi

            F.append(F_); Btg.append(Btg_); Eg.append(Eg_)

        F, Btg, Eg = np.array(F), np.array(Btg), np.array(Eg)

        return F, Btg, Efw, Eg

    def plot(self, p, z, Gp, F, Btg, Efw, Eg, Rv):
        """Create Material Balance Plots for Dry-Gas Reservoir"""
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy.optimize import curve_fit

        def calculate_condensate_inplace(Gfgi, Rv):
            """Calculate initial condensate-in-place from the calculated OGIP"""
            Rvi = Rv[0]
            condensate_inplace = Rvi * Gfgi # in STB
            return condensate_inplace        

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

        ## calculate condensate-in-place
        condensate_inplace = calculate_condensate_inplace(Gfgi, Rv)

        ## plot the regression line
        x1_fit = np.linspace(min(x1), max(x1), 5)
        y1_fit = linear_zero_intercept(x1_fit, Gfgi)
        plt.plot(x1_fit, y1_fit, label='(G) {} MMSCF (C) {} MSTB'.format(np.round(Gfgi * 1E-6, 3), np.round(condensate_inplace * 1E-3, 3)))
        plt.legend()

        # Plot 2: p/z vs Gp
        plt.subplot(3,2,2)
        plt.title('Plot 2: p/z vs Gp', size=title_size, pad=title_pad)
        plt.xlabel('Gp (scf)')
        plt.ylabel('p/z (psia)')

        if np.all(z==0) == False:        
          x2, y2 = Gp, (p / z)
          plt.plot(x2, y2, '.-')

          ## curve-fitting to calculate the slope as OGIP
          x2_norm = x2 / max(x2) # normalize x
          y2_norm = y2 / max(y2) # normalize y
          popt, pcov = curve_fit(linear_with_intercept, x2_norm, y2_norm)

          m, c = popt[0], popt[1]
          Gfgi = (-c / m) * max(x2) # OGIP is the intercept at x-axis, and denormalized
          m = m * max(y2) / max(x2) # denormalize the slope
          c = c * max(y2) # denormalize the intercept

          ## calculate condensate-in-place
          condensate_inplace = calculate_condensate_inplace(Gfgi, Rv)          

          ## plot the regression line
          x2_fit = np.linspace(min(x2), max(x2), 5)
          y2_fit = linear_with_intercept(x2_fit, m, c)
          plt.plot(x2_fit, y2_fit, label='(G) {} MMSCF (C) {} MSTB'.format(np.round(Gfgi * 1E-6, 3), np.round(condensate_inplace * 1E-3, 3)))
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

        ## calculate condensate-in-place
        condensate_inplace = calculate_condensate_inplace(Gfgi, Rv)        

        ## plot the regression line
        x3_fit = np.linspace(min(x3[1:]), max(x3[1:]), 5)
        y3_fit = linear_with_intercept(x3_fit, m, Gfgi)
        plt.plot(x3_fit, y3_fit, label='(G) {} MMSCF (C) {} MSTB'.format(np.round(Gfgi * 1E-6, 3), np.round(condensate_inplace * 1E-3, 3)))
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

        ## calculate condensate-in-place
        condensate_inplace = calculate_condensate_inplace(Gfgi, Rv)        

        ## plot the regression line
        x6_fit = np.linspace(min(x6), max(x6), 5)
        y6_fit = linear_zero_intercept(x6_fit, Gfgi)
        plt.plot(x6_fit, y6_fit, label='(G) {} MMSCF (C) {} MSTB'.format(np.round(Gfgi * 1E-6, 3), np.round(condensate_inplace * 1E-3, 3)))
        plt.legend()        

        # Plot 7: ((p/z)*(1-Efw)) vs Gp
        plt.subplot(3,2,5)
        plt.title('Plot 7: ((p/z)*(1-Efw)) vs Gp', size=title_size, pad=title_pad)
        plt.xlabel('Gp (scf)')
        plt.ylabel('(p/z)*(1-Efw) (psia)')        

        if np.all(z==0) == False:
          x7, y7 = Gp, ((p / z) * (1 - Efw))
          plt.plot(x7, y7, '.-')

          ## curve-fitting to calculate the slope as OGIP
          x7_norm = x7 / max(x7) # normalize x
          y7_norm = y7 / max(y7) # normalize y
          popt, pcov = curve_fit(linear_with_intercept, x7_norm, y7_norm)

          m, c = popt[0], popt[1]
          Gfgi = (-c / m) * max(x7) # OGIP is the intercept at x-axis, and denormalized
          m = m * max(y7) / max(x7) # denormalize the slope
          c = c * max(y7) # denormalize the intercept

          ## calculate condensate-in-place
          condensate_inplace = calculate_condensate_inplace(Gfgi, Rv)          

          ## plot the regression line
          x7_fit = np.linspace(min(x7), max(x7), 5)
          y7_fit = linear_with_intercept(x7_fit, m, c)
          plt.plot(x7_fit, y7_fit, label='(G) {} MMSCF (C) {} MSTB'.format(np.round(Gfgi * 1E-6, 3), np.round(condensate_inplace * 1E-3, 3)))
          plt.legend()  

        plt.tight_layout(pad=1.5)
        plt.show()

        return F, Eg, Efw  

class oil():
    """
    Oil (Undersaturated and saturated; Volatile and Non-volatile) Material Balance Plot
    """
    def calculate_params(self, p, Bo, Bg, Rv, Rs, Np, Gp, Gi, cf, cw, swi):
        """
        Calculate Material Balance Paramaters for Oil Reservoir
        """
        pi = p[0]
        Rsi = Rs[0]
        Rvi = Rv[0]
        Boi = Bo[0]
        Bgi = Bg[0]

        # calculate Efw
        Efw = ((cf + cw * swi) / (1 - swi)) * (pi - p)

        # calculate F, Bto, and Btg
        F = (Np * ((Bo - (Rs * Bg)) / (1 - (Rv * Rs)))) + ((Gp - Gi) * ((Bg - (Rv * Bo)) / (1 - (Rv * Rs))))
        Btg = ((Bg * (1 - (Rs * Rvi))) + (Bo * (Rvi - Rv))) / (1 - (Rv * Rs))  # in RB/STB
        Bto = ((Bo * (1 - (Rv * Rsi))) + (Bg * (Rsi - Rs))) / (1 - (Rv * Rs))  # in RB/scf

        # calculate Eo and Eg
        Eo = Bto - Boi
        Eg = Btg - Bgi

        return F, Bto, Btg, Efw, Eo, Eg

    def gascap(self, Gfgi, Nfoi, Bg, Bo):
      """
      Calculate Total Oil+Gas Expansion Factor from known Gas Cap ratio
      Gfgi and Nfoi known from volumetrics
      """
      Bgi, Boi = Bg[0], Bo[0]

      m = (Gfgi * Bgi) / (Nfoi * Boi)
      return m

    def plot(self, oil_type, F, Bto, Btg, Efw, Eo, Eg, Bo, Rs, Rv, figsize=(10,5)):
      """
      Create Material Balance Plots for Oil Reservoir
      
      Input:
      oil_type: 'undersaturated' or 'saturated'
      """
      import numpy as np
      import matplotlib.pyplot as plt
      from scipy.optimize import curve_fit
      import matplotlib.patches as mpl_patches

      # plot attributes
      title_size = 15
      title_pad = 14

      # linear function for curve-fit
      def linear_zero_intercept(x, m):
          y = m * x
          return y

      def linear_with_intercept(x, m, c):
          y = m * x + c
          return y

      if oil_type == 'undersaturated':

        plt.figure(figsize=figsize)

        " Plot 1: F vs (Eg+Boi*Efw) "

        plt.subplot(1,2,1)
        Boi = Bo[0]
        x1, y1 = (Eg + Boi * Efw), F
        plt.plot(x1, y1, '.-')
        plt.title(r'Plot 1: $F$ vs $(E_o+B_{oi}*E_{fw})$', size=title_size, pad=title_pad)
        plt.xlabel(r'$E_o+B_{oi}E_{fw}$ (RB/STB)', size=15)
        plt.ylabel(r'$F$ (res bbl)', size=15)

        ## curve-fitting to calculate the slope as OOIP
        x1_norm = x1 / max(x1) # normalize x
        y1_norm = y1 / max(y1) # normalize y
        popt, pcov = curve_fit(linear_zero_intercept, x1_norm, y1_norm)

        m = popt[0]
        Nfoi = m * max(y1) / max(x1) # denormalize the slope, hence the OGIP

        ## Calculate OOIP and OGIP from Nfoi
        Rsi = Rs[0]
        Gfgi = 0 # no free gas phase in undersaturated oil
        OOIP = Nfoi
        OGIP = Nfoi * Rsi

        ## Output results into text in plot
        labels, handles, OOIP, OGIP = initial_hydrocarbon_in_place(Nfoi, Gfgi, Rv, Rs)    

        ## plot the regression line
        x1_fit = np.linspace(min(x1), max(x1), 5)
        y1_fit = linear_zero_intercept(x1_fit, Nfoi)
        plt.plot(x1_fit, y1_fit, label='{} MMSTB'.format(np.round(Nfoi * 1E-6, 3)))

        plt.legend(handles, labels, loc='best', fontsize='small', 
                   fancybox=True, framealpha=0.7, 
                   handlelength=0, handletextpad=0) 

        " Plot 2: F/(Eg+Boi*Efw) vs Np (Waterdrive Diagnostic Plot) "

        plt.subplot(1,2,2)
        x2, y2 = Np, F / (Eg + Boi * Efw)
        plt.plot(x2, y2, '.-')
        plt.title('Plot 2: Waterdrive Diagnostic Plot', size=title_size, pad=title_pad)
        plt.xlabel(r'$N_p$ (STB)', size=15)
        plt.ylabel(r'$\frac{F}{(E_o+B_{oi}E_{fw})}$ (STB)', size=15)

        ## curve-fitting to calculate the slope as OOIP, here [1:] because NaN is removed
        x2_norm = x2[1:] / max(x2[1:]) # normalize x
        y2_norm = y2[1:] / max(y2[1:]) # normalize y
        popt, pcov = curve_fit(linear_with_intercept, x2_norm, y2_norm)

        m, c = popt[0], popt[1]
        m = m * max(y2[1:]) / max(x2[1:]) # denormalize the slope
        Nfoi = c * max(y2[1:]) # denormalize the intercept, hence the OGIP

        ## Calculate OOIP and OGIP from Nfoi
        Rsi = Rs[0]
        Gfgi = 0 # no free gas phase in undersaturated oil
        OOIP = Nfoi
        OGIP = Nfoi * Rsi

        ## Output results into text in plot
        labels, handles, OOIP, OGIP = initial_hydrocarbon_in_place(Nfoi, Gfgi, Rv, Rs)           

        ## plot the regression line
        x2_fit = np.linspace(min(x2[1:]), max(x2[1:]), 5)
        y2_fit = linear_with_intercept(x2_fit, m, Nfoi)
        plt.plot(x2_fit, y2_fit, label='{} MMSTB'.format(np.round(Nfoi * 1E-6, 3)))

        plt.legend(handles, labels, loc='best', fontsize='small', 
                   fancybox=True, framealpha=0.7, 
                   handlelength=0, handletextpad=0)  
        
        plt.tight_layout(1)
        plt.show()

      if oil_type == 'saturated':

        plt.figure(figsize=figsize)

        " Plot 1: F/Eo vs Eg/Eo "

        plt.subplot(1,3,1)
        x1, y1 = (Eg / Eo), (F / Eo)
        plt.plot(x1, y1, '.-')
        plt.title('Plot 1: F/Eo vs Eg/Eo', size=title_size, pad=title_pad)
        plt.xlabel(r'$\frac{Eg}{Eo}$ (STB/scf)', size=15)
        plt.ylabel(r'$\frac{F}{Eo}$ (STB)', size=15)

        ## curve-fitting to calculate the slope as Gfgi, intercept as Nfoi
        x1_norm = x1[1:] / max(x1[1:]) # normalize x
        y1_norm = y1[1:] / max(y1[1:]) # normalize y
        popt, pcov = curve_fit(linear_with_intercept, x1_norm, y1_norm)

        m, c = popt[0], popt[1]
        Gfgi = m = m * max(y1[1:]) / max(x1[1:]) # denormalize the slope
        Nfoi = c = c * max(y1[1:]) # denormalize the intercept

        ## calculate OOIP and OGIP from Nfoi and Gfgi
        Rsi, Rvi = Rs[0], Rv[0]
        OOIP = Nfoi + Gfgi * Rvi
        OGIP = Gfgi + Nfoi * Rsi

        ## Output results into text in plot
        labels, handles, OOIP, OGIP = initial_hydrocarbon_in_place(Nfoi, Gfgi, Rv, Rs)    

        ## plot the regression line
        x1_fit = np.linspace(min(x1[1:]), max(x1[1:]), 5)
        y1_fit = linear_with_intercept(x1_fit, m, c)
        plt.plot(x1_fit, y1_fit)

        plt.legend(handles, labels, loc='best', fontsize='small', 
                   fancybox=True, framealpha=0.7, 
                   handlelength=0, handletextpad=0)

        " Plot 2: p/z vs Gp "

        plt.subplot(1,3,2)
        x2, y2 =  (Eo / Eg), (F / Eg)
        plt.plot(x2, y2, '.-')
        plt.title('Plot 2: F/Eg vs Eo/Eg', size=title_size, pad=title_pad)
        plt.xlabel(r'$\frac{Eo}{Eg}$ (scf/STB)', size=15)
        plt.ylabel(r'$\frac{F}{Eg}$ (scf)', size=15)

        ## curve-fitting to calculate the slope as Nfoi, intercept as Gfgi
        x2_norm = x2[1:] / max(x2[1:]) # normalize x
        y2_norm = y2[1:] / max(y2[1:]) # normalize y
        popt, pcov = curve_fit(linear_with_intercept, x2_norm, y2_norm)

        m, c = popt[0], popt[1]
        Nfoi = m = m * max(y2[1:]) / max(x2[1:]) # denormalize the slope
        Gfgi = c = c * max(y2[1:]) # denormalize the intercept

        ## calculate OOIP and OGIP from Nfoi and Gfgi
        Rsi, Rvi = Rs[0], Rv[0]
        OOIP = Nfoi + Gfgi * Rvi
        OGIP = Gfgi + Nfoi * Rsi

        ## Output results into text in plot
        labels, handles, OOIP, OGIP = initial_hydrocarbon_in_place(Nfoi, Gfgi, Rv, Rs)    

        ## plot the regression line
        x2_fit = np.linspace(min(x2[1:]), max(x2[1:]), 5)
        y2_fit = linear_with_intercept(x2_fit, m, c)
        plt.plot(x2_fit, y2_fit)

        plt.legend(handles, labels, loc='best', fontsize='small', 
                   fancybox=True, framealpha=0.7, 
                   handlelength=0, handletextpad=0) 

        plt.tight_layout(1)                 

        plt.show()    
