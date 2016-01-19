# Script to generate a plot of the sparsity of the matrix
# based on density data outputted from the main program. 

import matplotlib.pyplot as plt
from matplotlib import mpl
from matplotlib.colors import LogNorm
import numpy as np

S1 = 27
F1 = 'thc-4.sparse'
S2 = 23
F2 = 'taxol-4.sparse'
S3 = 21
F3 = 'pin1ppiase-4.sparse'
S4 = 21
F4 = 'kinase-4.sparse'

d1 = np.zeros([S1, S1])  
d2 = np.zeros([S2, S2])
d3 = np.zeros([S3, S3])
d4 = np.zeros([S4, S4])
# Read data from file
z1 = np.loadtxt(F1)
z2 = np.loadtxt(F2)
z3 = np.loadtxt(F3)
z4 = np.loadtxt(F4)

for i in range(len(z1)):
    d1[z1[i][0]][z1[i][1]] = z1[i][2]

for i in range(len(z2)):
    d2[z2[i][0]][z2[i][1]] = z2[i][2]

for i in range(len(z3)):
    d3[z3[i][0]][z3[i][1]] = z3[i][2]

for i in range(len(z4)):
    d4[z4[i][0]][z4[i][1]] = z4[i][2]

cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', ['white', 'black'], 256)

f, ((img1, img2), (img3, img4)) = plt.subplots(2, 2)

img1.imshow(d1, interpolation='nearest', cmap = cmap)
img1.axes.get_xaxis().set_ticks([])
img1.axes.get_yaxis().set_ticks([])

img2.imshow(d2, interpolation='nearest', cmap = cmap)
img2.axes.get_xaxis().set_ticks([])
img2.axes.get_yaxis().set_ticks([])

img3.imshow(d3, interpolation='nearest', cmap = cmap)
img3.axes.get_xaxis().set_ticks([])
img3.axes.get_yaxis().set_ticks([])

img4.imshow(d4, interpolation='nearest', cmap = cmap)
img4.axes.get_xaxis().set_ticks([])
img4.axes.get_yaxis().set_ticks([])

plt.show()
        
