from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import numpy as np
import matplotlib.pyplot as plt

arr = np.loadtxt("test.tsv")
x = arr[:,0]
y = arr[:,1]
u = arr[:,2]
#x, y = np.meshgrid(x, y)

fig = plt.figure()
axes = fig.gca(projection='3d')

z = np.sqrt(x+y)
print(z)

surf = axes.plot_surface(x,y,z)#,rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)

 #c = plt.contourf(x,y,u)

plt.show()
