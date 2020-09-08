"""
Volumetrics calculation

@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

def get_areas(figs, xi, yi, plot='Yes'):
  # get the contour lines using "allsegs" method
  import numpy as np
  
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
