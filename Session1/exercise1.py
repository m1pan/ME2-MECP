import numpy as np
import math as mt
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#Task A
a = np.arange(-5,-2.5,0.5)
b = np.arange(-2.05,3,0.05)
c = np.arange(3,5.5,0.5)
x = np.concatenate((a,b,c))

f = np.sin(x)
print(f[21])
plt.scatter(x,f)

#Task B
x = np.linspace(-2*np.pi, 2*np.pi, 41)
y = np.linspace(-1*np.pi, 2.1*np.pi, 31)
Nx = len(x)
Ny = len(y)

f = np.ndarray((Ny,Nx))
for i in range(0,Ny):
    for j in range(0,Nx):
        f[i,j] = mt.sin(x[j]) * mt.cos(y[i])

g = np.ndarray((Ny,Nx))
for i in range(0,Ny):
    for j in range(0,Nx):
        g[i,j] = mt.sin(y[i]) * mt.cos(x[j])

Xg,Yg = np.meshgrid(x,y)

p = f * g
s = f + g

# task c
ax = plt.axes(projection='3d')
# ax.plot_surface(Xg,Yg,s)

t = np.arange(0,10.05,0.05)
r = f * np.exp(-0.5*t[100])
plt.contour(Xg,Yg,r)
plt.show()
