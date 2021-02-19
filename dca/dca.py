def arps_fit(t, q):
  """
  Arps Decline Curve Analysis using Non-Linear Curve-Fitting
  
  Input:
  t = time array (in numpy datetime64)
  q = production rate array (unit: STB/day, or SCF/day)
  Output:
  qi = initial production rate (unit: STB/day, or SCF/day)
  di = initial decline rate (unit: STB/day, or SCF/day)
  b = decline exponent 
  """
  import numpy as np
  import datetime
  from scipy.optimize import curve_fit

  def hyperbolic(t, qi, di, b):
    return qi / (np.abs((1 + b * di * t))**(1/b))
  
  def rmse(y, yfit):
    N = len(y)
    return np.sqrt(np.sum(y-yfit)**2 / N)

  # subtract one datetime to another datetime
  timedelta = [j-i for i, j in zip(t[:-1], t[1:])]
  timedelta = np.array(timedelta)
  timedelta = timedelta / datetime.timedelta(days=1)

  # take cumulative sum over timedeltas
  t = np.cumsum(timedelta)
  t = np.append(0, t)
  t = t.astype(float)

  # normalize the time and rate data
  t_normalized = t / max(t)
  q_normalized = q / max(q)  

  # fitting the data with the hyperbolic function
  popt, pcov = curve_fit(hyperbolic, t_normalized, q_normalized)
  qi, di, b = popt

  # RMSE is calculated on the normalized variables
  qfit_normalized = hyperbolic(t_normalized, qi, di, b)
  RMSE = rmse(q_normalized, qfit_normalized)

  # De-normalize qi and di
  qi = qi * max(q)
  di = di / max(t)

  return qi, di, b, RMSE
