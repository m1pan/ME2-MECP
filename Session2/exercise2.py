import numpy as np
import matplotlib.pyplot as plt

def traprule(x,y):
    """trapezium rule for numerical integration with non-uniform intervals"""
    return np.sum((x[1:] - x[:-1]) * (y[1:] + y[:-1]) / 2)

def test(x):
    return 1/sqrt(x**17.1 +2023)

# x = np.zeros()
# for i in range(0,4):
#     x[i] = np.linspace(0.1,1,10**i)

with open('Thames.txt','r') as f:
    data = f.readlines()
    data = [line.split(',')  for line in data]

north = np.zeros((len(data),2))
south = np.zeros((len(data),2))
for i in range(0,len(data)):
    north[i] = float(data[i][0]),float(data[i][1])
    south[i] = float(data[i][2]),float(data[i][3])

plt.plot(north[:,0],north[:,1],color='red')
plt.plot(south[:,0],south[:,1],color='blue')
plt.axis('equal')
plt.show()

area = traprule(north[:,0],north[:,1]) - traprule(south[:,0],south[:,1])
