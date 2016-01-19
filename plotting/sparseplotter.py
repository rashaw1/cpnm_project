# Script to generate a plot of the sparsity of the matrix
# based on density data outputted from the main program. 

import matplotlib.pyplot as plt
from matplotlib import mpl
from matplotlib.colors import LogNorm
import numpy as np

SIZE = 27
FILENAME = 'thc-8_dif.sparse'

densities = np.zeros([SIZE, SIZE])  

# Read data from file
z = np.loadtxt(FILENAME)

for i in range(len(z)):
    densities[z[i][0]][z[i][1]] = 1 + z[i][2]

cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', ['white', 'black'], 256)

img = plt.imshow(densities, interpolation='nearest', cmap = cmap)

img.axes.get_xaxis().set_ticks([])
img.axes.get_yaxis().set_ticks([])

plt.show()
        
