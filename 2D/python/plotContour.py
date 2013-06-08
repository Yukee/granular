from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt

# loads data file into a 2D array (arr[i][j] is the element at line i, column j in the tsv file)

x = np.arange(0,1,0.1)
y = np.arange(0,1,0.1)

mdata = np.loadtxt("../Results/2DNSContour_0.tsv")

# builds coordinate matrices from x and y coordinate vectors
xx, yy = np.meshgrid(x,y)

# loads info file containing dimensions of the finite-volume cells, timestep and current time
# arrInfo = np.loadtxt("info.tsv",usecols=(1,2,3))

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

ax.contourf(xx,yy,mdata)

plt.show()
