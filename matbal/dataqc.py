import numpy as np

def dataqc(data, Pdb, reservoir='gas', cleansing='Yes', plot='Yes'):

    # check NaN values in Gp, Np, and Rv
    if (data['Gp'].isnull().values.any()) == True or (data['Np'].isnull().values.any()) == True or (data['Rv'].isnull().values.any()) == True:
        if cleansing == 'Yes':

            # drop rows with missing values
            data.dropna(inplace=True)

            # re-index dataframe (in order, started from 0)
            data = data.reset_index(drop=True)
    else:
        data = data

    p = data['p'].values

    # separate the data for P above bubblepoint/dewpoint
    above = data[p >= Pdb]

    p_above = above['p'].values
    Rs_above = above['Rs'].values
    Bo_above = above['Bo'].values
    Bg_above = above['Bg'].values
    Rv_above = above['Rv'].values
    Rvi = Rv_above[0]

    Rs_above_qc = []
    Bo_above_qc = []
    Bg_above_qc = []

    if reservoir == 'gas':
        for i in range(len(p_above)):
            if Rs_above[i] == np.nan:
                Rs_above_ = 1 / Rvi
            elif Rs_above[i] != np.nan:
                Rs_above_ = Rs_above[i]
            if Bo_above[i] == np.nan:
                Bo_above_ = Bg_above[i] * Rv_above[i]
            elif Bo_above[i] != np.nan:
                Bo_above_ = Bo_above[i]
            if Bg_above[i] == np.nan:
                Bg_above_ = gasfvf()
            elif Bg_above[i] != np.nan:
                Bg_above_ = Bg_above[i]

            Rs_above_qc.append(Rs_above_)
            Bo_above_qc.append(Bo_above_)
            Bg_above_qc.append(Bg_above_)

    # separate the data for P below bubblepoint/dewpoint

    below = data[p < Pdb]

    p_below = below['p'].values
    Rs_below = below['Rs'].values
    Bo_below = below['Bo'].values
    Bg_below = below['Bg'].values
    Rv_below = below['Rv'].values
    Rvi = Rv_below[0]

    Rs_below_qc = []
    Bo_below_qc = []
    Bg_below_qc = []

    if reservoir == 'gas':
        for i in range(len(p_below)):
            if Rs_below[i] == np.nan:
                Rs_below_ = gasoilratio()
            elif Rs_below[i] != np.nan:
                Rs_below_ = Rs_below[i]
            if Bo_below[i] == np.nan:
                Bo_below_ = oilfvf()
            elif Bo_below[i] != np.nan:
                Bo_below_ = Bo_below[i]
            if Bg_below[i] == np.nan:
                Bg_below_ = gasfvf()
            elif Bg_below[i] != np.nan:
                Bg_below_ = Bg_below[i]

            Rs_below_qc.append(Rs_below_)
            Bo_below_qc.append(Bo_below_)
            Bg_below_qc.append(Bg_below_)
