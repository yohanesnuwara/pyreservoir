import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from scipy.optimize import curve_fit
from sklearn.covariance import EllipticEnvelope

def elliptic_outlier(df, porosity_column, permeability_column, 
                     contamination=0.27):
  """
  Remove Outliers in Porosity and Log(Permeability) Space
  """
  # DataFrame
  df_complete = df.copy().reset_index(drop=True)
  df = df[[porosity_column, permeability_column]].reset_index(drop=True)

  # Convert k to log k
  df[permeability_column] = np.log(df[permeability_column].values)

  # Outlier detection
  model = EllipticEnvelope(contamination=contamination) 
  out = model.fit_predict(df)

  # Plot outlier vs. non outlier
  id_not_out = list(np.where(out==1)[0])
  dfx = df.loc[id_not_out,:]

  plt.scatter(df[porosity_column], df[permeability_column], 
              c='r', label='Outlier')
  plt.scatter(dfx[porosity_column], dfx[permeability_column], 
              c='k', label='Not Outlier')
  plt.xlabel('Porosity')
  plt.ylabel('Log k [log.mD]')
  plt.legend()

  # Output dataframe (complete)
  df_out = df_complete.loc[id_not_out,:]

  return df_out

def timur(poro, swi, a, b, c):
  return a * poro**b / swi**c

def logk(poro, swirr, a, b, c):
  return np.exp(1)**((a* np.log(poro)) + (b * np.log(swirr)) + c)

def fit_timur(df, porosity_column, permeability_column):
  # Fit with Timur equation
  [swi, a, b, c], _ = curve_fit(timur, df[porosity_column], df[permeability_column])

  # Plot Timur line
  poro = np.linspace(0, 1, 100)
  k = timur(poro, swi, a, b, c)
  plt.plot(poro, k, c='r', lw=3, label=f'Swirr={swi:.3f}\na={a:.3f}\nb={b:.3f}\nc={c:.3f}')

  # Plot multiple lines
  swirr = np.arange(0.1, 1, 0.1)
  for i in swirr:
    k = timur(poro, i, a, b, c)
    plt.plot(poro, k, lw=0.7, c='k')

  # Plot points
  plt.scatter(df[porosity_column], df[permeability_column], c=df['Measured Depth'])
  plt.colorbar()

  # Axis transformation
  plt.yscale('log')
  plt.xlabel('Porosity [fraction]')
  plt.ylabel('Permeability [mD]')  
  plt.legend(loc='lower right')
