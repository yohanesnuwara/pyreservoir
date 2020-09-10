"""
Volumetrics calculation

@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

def get_contours(figs, xi, yi, plot='Yes'):
  # get the contour lines using "allsegs" method
  import numpy as np
  import matplotlib.pyplot as plt
  
  lines = figs.allsegs

  contour_all = []
  for i in range(len(lines)):
    contour_each = []
    for j in range(len(lines[i])):
      for k in range(len(lines[i][j])):
        _ = list(lines[i][j][k])
        contour_each.append(_)  
    contour_all.append(contour_each)
  
  if plot=='Yes':
    # plot the individual contour lines
    plt.figure(figsize=(17,17))
    contours = np.linspace(0, 1, len(lines))

    for i in range(len(contour_all)):
      plt.subplot(4,3,i+1)
      seg = np.array(contour_all[i])
      for j in range(len(seg)):
        seg = np.array(contour_all[i])
        x_contour = seg[:,0]
        y_contour = seg[:,1]  
        plt.plot(x_contour, y_contour, '.')
      plt.title('Contour {}'.format(np.round(contours[i], 1)))    
      plt.xlim(np.min(xi), np.max(xi)); plt.ylim(np.min(yi), np.max(yi))  
  else:
    pass

  return contour_all  

def compute_area(contour_all):
  # Use Green's theorem to compute the area
  # enclosed by the given contour.
  def area(vs):
      a = 0
      x0,y0 = vs[0]
      for [x1,y1] in vs[1:]:
          dx = x1-x0
          dy = y1-y0
          a += 0.5*(y0*dx - x0*dy)
          x0 = x1
          y0 = y1
      return a

  contour_area = []
  for i in range(len(contour_all)):
      seg = np.array(contour_all[i])
      # Compute area enclosed by vertices.
      a = area(seg)
      contour_area.append(float(np.abs(a))) 

  return contour_area

def trapezoid(y_vals, h):

  # simplification: V = h/2*(a1+an) + h/2(2*a2+2*a3+2*an-1) as V = first_term + second_term
  first_term = 0.5 * h * (y_vals[0] + y_vals[-1])

  # delete a0 and an
  y_trapez = np.delete(y_vals, 0)
  y_trapez = np.delete(y_trapez, -1)

  second_term = np.sum(2 * y_trapez)
  second_term = 0.5 * h * second_term
  oip_trapezo = first_term + second_term
  
  return oip_trapezo

def pyramidal(y_vals, h):
  import numpy as np

  # sort out a_0 and a_n, calculate the first term, 1/3*delta_z*(a_0+a_n)
  minus_a0 = np.delete(y_vals, 0)
  minus_a0_an = np.delete(minus_a0, -1)
  first_term = y_vals[0] + y_vals[-1]

  # multiply the the first part of second term from sorted list: second_term_1 = (delta_z/3) * (2a1+2a2+...+2an-1)
  second_term_1 = 2 * (np.sum(minus_a0_an))

  # process the second part of second term from sorted list: second_term_2 = (delta_z/3) * (sqrt(a1*a2)+sqrt(a2*a3)+...+sqrt(an-2*an-1))
  second_term_2 = [np.sqrt(j*i) for i, j in zip(y_vals[:-1], y_vals[1:])]
  second_term_2 = np.sum(second_term_2)

  oip_pyramidal = (h / 3) * (first_term + second_term_1 + second_term_2)
  return oip_pyramidal

def simpson(y_vals, h):
    i = 1
    total = y_vals[0] + y_vals[-1]
    for y in y_vals[1:-1]:
        if i % 2 == 0:
            total += 2 * y
        else:
            total += 4 * y
        i += 1
    oip_simpson = total * (h / 3.0)
    return oip_simpson
