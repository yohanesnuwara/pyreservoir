import matplotlib.pyplot as plt

class drygas():
  def plot1(self):
      plt.plot(Eg, F)
      plt.show()
      return(Eg, F)
  def plot2(self):
      plt.plot(Gp, p_z)
      plt.show()
      return((Gp, p_z))
  def plot3(self):
      plt.plot(Gp, F_eg)
      plt.show()
      return(Gp, F_eg)
  def plot6(self):
      plt.plot(Eg_Bgi_Efw, F)
      plt.show()
      return(Eg_Bgi_Efw, F)
  def plot7(self):
      plt.plot(Gp, p_z_Efw)
      plt.show()
      return(Gp, p_z_Efw)

class condensategas():
  def plot1(self):
      plt.plot(Eg, F)
      plt.show()
      return(Eg, F)
  def plot2(self):
      plt.plot(Gp, p_z)
      plt.show()
      return((Gp, p_z))
  def plot3(self):
      plt.plot(Gp, F_eg)
      plt.show()
      return(Gp, F_eg)
  def plot6(self):
      plt.plot(Eg_Bgi_Efw, F)
      plt.show()
      return(Eg_Bgi_Efw, F)
  def plot7(self):
      plt.plot(Gp, p_z_Efw)
      plt.show()
      return(Gp, p_z_Efw)

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
