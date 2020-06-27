"""
Material Balance Plots
@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

import numpy as np
import matplotlib.pyplot as plt

def Efw(cf, cw, swi, p, pi):
    """
    Calculate formation expansion factor
    """
    Efw = ((cf + cw * swi) / (1 - swi)) * (pi - p)
    return (Efw)


def condensate_belowdew(Rs, Rv, Rsi, Rvi, Bo, Bg, Bgi, Np, Gp):
    """
    Calculate the parameters for material balance plot of gas-condensate reservoirs
    below dewpoint pressure
    Input:
    Rs: array
    Rv: array
    Rsi: initial Rs, float (NOTE: if data doesn't provide, calculate it with calculate_condensate_params function)
    Rvi: initial Rv, float (from data Rv)
    Bo: array
    Bg: array
    Np: array
    Gp: array
    Material balance plots:
    * Plot 10.1: F vs Eg
    Output:
    F: array
    Eg: array
    """
    Btg = ((Bg * (1 - (Rs * Rvi))) + (Bo * (Rvi - Rv))) / (1 - (Rv * Rs))  # in RB/STB
    Bto = ((Bo * (1 - (Rv * Rsi))) + (Bg * (Rsi - Rs))) / (1 - (Rv * Rs))  # in RB/scf

    Gi = 0
    F = (Np * ((Bo - (Rs * Bg)) / (1 - (Rv * Rs)))) + ((Gp - Gi) * ((Bg - (Rv * Bo)) / (1 - (Rv * Rs))))
    Eg = Btg - Bgi
    return (F, Eg)


def condensate_abovedew(Bg, Bgi, Gp, Gpi):
    """
    Calculate the parameters for material balance plot of gas-condensate reservoirs
    above dewpoint pressure
    Input:
    Bg: array
    Bgi: initial Bg, float
    Gp: array
    Gpi: initial Gp, float
    Material balance plots:
    * Plot 10.1: F vs Eg
    Output:
    F: array
    Eg: array
    """
    Eg = Bg - Bgi
    F = Bg * (Gp - Gpi)
    return (F, Eg)


class condensate():
    """
    Gas Condensate Material Balance Plot
    """

    def plot1(self, Pdp, p, Bg, Bgi, Bo, Np, Gp, Gpi, Rv, Rvi, Rs=None, Rsi=None):
        """
        Plot 1: F vs Eg

        Input:
        array: p, Bg, Bo, Np, Gp, Rv, Rs
        float: Pdp (dewpoint-pressure), Bgi (initial Bg), Gpi (initial Gp), Rvi (initial Rv), Rsi (initial Rs)

        For Rs and Rsi, specify as None, if there's no data. Through this function, both will be calculated theoretically
        If there are data of Rs and Rsi, specify the variable.

        Output:
        array: Eg (x-axis values), F (y-axis values), and plot
        """

        # for above dewpoint pressure
        id_above = np.where(p > Pdp)[0]

        Bg_above = [];
        Gp_above = []

        for i in id_above:
            Bg_above_ = Bg[i]
            Gp_above_ = Gp[i]
            Bg_above.append(Bg_above_);
            Gp_above.append(Gp_above_)

        F_above, Eg_above = condensate_abovedew(np.array(Bg_above), Bgi, np.array(Gp_above), Gpi)

        # for below dewpoint pressure
        id_below = np.where(p <= Pdp)[0]
        if Rs == None and Rsi == None:
            Rs = 1 / Rv
            Rsi = Rs[0]
        else:
            Rs = Rs
            Rsi = Rs[0]

        Rs_below = [];
        Rv_below = [];
        Bo_below = [];
        Bg_below = []
        Np_below = [];
        Gp_below = []

        for i in id_below:
            Rs_below_ = Rs[i]
            Rv_below_ = Rv[i]
            Bo_below_ = Bo[i]
            Bg_below_ = Bg[i]
            Np_below_ = Np[i]
            Gp_below_ = Gp[i]

            Rs_below.append(Rs_below_);
            Rv_below.append(Rv_below_);
            Bo_below.append(Bo_below_)
            Bg_below.append(Bg_below_);
            Np_below.append(Np_below_);
            Gp_below.append(Gp_below_)

        F_below, Eg_below = condensate_belowdew(np.array(Rs_below), np.array(Rv_below),
                                                Rsi, Rvi, np.array(Bo_below), np.array(Bg_below),
                                                Bgi, np.array(Np_below), np.array(Gp_below))

        # append F and Eg of the below- and above-dewpoint pressure
        F = np.append(F_above, F_below)
        Eg = np.append(Eg_above, Eg_below)

        # plot
        plt.plot(Eg, F, '.')
        plt.title('Plot 1: F vs Eg')
        plt.xlabel('Eg (RB/scf)');
        plt.ylabel('F (res bbl)')
        plt.xlim(xmin=0);
        plt.ylim(ymin=0)
        plt.show()

        return (F, Eg)

    def plot2(self, Gp, p, z):
        """
        Plot 2: p/z vs Gp

        Input:
        array: Gp, p, z

        Output:
        array: Gp (x-axis values), p_z (y-axis values), plot
        """
        # calculate plotting parameters
        p_z = p / z
        Gp = Gp

        # plotting
        plt.plot(Gp, p_z, '.')
        plt.title('Plot 2: p/z vs Gp')
        plt.xlabel('Gp (scf)');
        plt.ylabel('p/z (psia)')
        plt.xlim(xmin=0);
        plt.ylim(ymin=0)
        plt.show()
        return (Gp, p_z)

    def plot3(self, Pdp, p, Bg, Bgi, Bo, Np, Gp, Gpi, Rv, Rvi, Rs=None, Rsi=None):
        """
        Plot 1: F/Eg vs Gp

        Input:
        array: p, Bg, Bo, Np, Gp, Rv, Rs
        float: Pdp (dewpoint-pressure), Bgi (initial Bg), Gpi (initial Gp), Rvi (initial Rv), Rsi (initial Rs)

        For Rs and Rsi, specify as None, if there's no data. Through this function, both will be calculated theoretically
        If there are data of Rs and Rsi, specify the variable.

        Output:
        array: Gp (x-axis values), F_Eg (y-axis values), and plot
        """
        # for above dewpoint pressure
        id_above = np.where(p > Pdp)[0]

        Bg_above = [];
        Gp_above = []

        for i in id_above:
            Bg_above_ = Bg[i]
            Gp_above_ = Gp[i]
            Bg_above.append(Bg_above_);
            Gp_above.append(Gp_above_)

        F_above, Eg_above = condensate_abovedew(np.array(Bg_above), Bgi, np.array(Gp_above), Gpi)

        # for below dewpoint pressure
        id_below = np.where(p <= Pdp)[0]
        if Rs == None and Rsi == None:
            Rs = 1 / Rv
            Rsi = Rs[0]
        else:
            Rs = Rs
            Rsi = Rs[0]

        Rs_below = [];
        Rv_below = [];
        Bo_below = [];
        Bg_below = []
        Np_below = [];
        Gp_below = []

        for i in id_below:
            Rs_below_ = Rs[i]
            Rv_below_ = Rv[i]
            Bo_below_ = Bo[i]
            Bg_below_ = Bg[i]
            Np_below_ = Np[i]
            Gp_below_ = Gp[i]

            Rs_below.append(Rs_below_);
            Rv_below.append(Rv_below_);
            Bo_below.append(Bo_below_)
            Bg_below.append(Bg_below_);
            Np_below.append(Np_below_);
            Gp_below.append(Gp_below_)

        F_below, Eg_below = condensate_belowdew(np.array(Rs_below), np.array(Rv_below),
                                                Rsi, Rvi, np.array(Bo_below), np.array(Bg_below),
                                                Bgi, np.array(Np_below), np.array(Gp_below))

        # append F and Eg of the below- and above-dewpoint pressure
        F = np.append(F_above, F_below)
        Eg = np.append(Eg_above, Eg_below)

        # calculate parameters for plotting
        F_Eg = F / Eg

        plt.plot(Gp, F_Eg)
        plt.show()
        return (Gp, F_Eg)

    def plot6(self, Pdp, p, pi, Bg, Bgi, Bo, Np, Gp, Gpi, cf, cw, swi, Rv, Rvi, Rs=None, Rsi=None):
        """
        Plot 6: F vs (Eg+Bgi*Efw)

        Input:
        array: p, Bg, Bo, Np, Gp, Rv, Rs
        float: Pdp (dewpoint-pressure), Bgi (initial Bg), Gpi (initial Gp), Rvi (initial Rv), Rsi (initial Rs), pi (initial p), swi (initial sw), cf, cw

        For Rs and Rsi, specify as None, if there's no data. Through this function, both will be calculated theoretically
        If there are data of Rs and Rsi, specify the variable.

        Output:
        array: Eg_Bgi_Efw (x-axis values), F (y-axis values), and plot
        """
        # for above dewpoint pressure
        id_above = np.where(p > Pdp)[0]

        Bg_above = [];
        Gp_above = []

        for i in id_above:
            Bg_above_ = Bg[i]
            Gp_above_ = Gp[i]
            Bg_above.append(Bg_above_);
            Gp_above.append(Gp_above_)

        F_above, Eg_above = condensate_abovedew(np.array(Bg_above), Bgi, np.array(Gp_above), Gpi)

        # for below dewpoint pressure
        id_below = np.where(p <= Pdp)[0]
        if Rs == None and Rsi == None:
            Rs = 1 / Rv
            Rsi = Rs[0]
        else:
            Rs = Rs
            Rsi = Rs[0]

        Rs_below = [];
        Rv_below = [];
        Bo_below = [];
        Bg_below = []
        Np_below = [];
        Gp_below = []

        for i in id_below:
            Rs_below_ = Rs[i]
            Rv_below_ = Rv[i]
            Bo_below_ = Bo[i]
            Bg_below_ = Bg[i]
            Np_below_ = Np[i]
            Gp_below_ = Gp[i]

            Rs_below.append(Rs_below_);
            Rv_below.append(Rv_below_);
            Bo_below.append(Bo_below_)
            Bg_below.append(Bg_below_);
            Np_below.append(Np_below_);
            Gp_below.append(Gp_below_)

        F_below, Eg_below = condensate_belowdew(np.array(Rs_below), np.array(Rv_below),
                                                Rsi, Rvi, np.array(Bo_below), np.array(Bg_below),
                                                Bgi, np.array(Np_below), np.array(Gp_below))

        # append F and Eg of the below- and above-dewpoint pressure
        F = np.append(F_above, F_below)
        Eg = np.append(Eg_above, Eg_below)

        # calculate parameters for plotting
        Efw = Efw(cf, cw, swi, p, pi)
        Eg_Bgi_Efw = Eg + Bgi * Efw

        # plotting
        plt.plot(Eg_Bgi_Efw, F)
        plt.show()
        return (Eg_Bgi_Efw, F)

    def plot7(self, Gp, p, pi, z, cf, cw, swi):
        """
        Plot 7 (p/z*(1-Efw)) vs Gp

        Input:
        array: Gp, p, z
        float: pi (initial p), swi (initial sw), cf, cw

        Output:
        array: Gp (x-axis values), p_z_Efw (y-axis values), plot
        """
        # calculate plotting parameters
        Efw = Efw(cf, cw, swi, p, pi)
        p_z_Efw = (p / z) * Efw

        # plotting
        plt.plot(Gp, p_z_Efw)
        plt.show()
        return (Gp, p_z_Efw)


class drygas():
    """
    Dry Gas Material Balance Plot
    """

    def plot1(self):
        plt.plot(Eg, F)
        plt.show()
        return (Eg, F)

    def plot2(self):
        plt.plot(Gp, p_z)
        plt.show()
        return ((Gp, p_z))

    def plot3(self):
        plt.plot(Gp, F_eg)
        plt.show()
        return (Gp, F_eg)

    def plot6(self):
        plt.plot(Eg_Bgi_Efw, F)
        plt.show()
        return (Eg_Bgi_Efw, F)

    def plot7(self):
        plt.plot(Gp, p_z_Efw)
        plt.show()
        return (Gp, p_z_Efw)


def regression(x, y):
    import numpy as np
    x = np.array(x);
    y = np.array(y)
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    err_x = x - mean_x
    err_y = y - mean_y
    err_mult = err_x * err_y
    numerator = np.sum(err_mult)
    err_x_squared = err_x ** 2
    denominator = np.sum(err_x_squared)
    slope = numerator / denominator
    intercept = mean_y - B1 * mean_x
    return (intercept, slope)
