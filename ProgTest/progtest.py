import math
import numpy as np
import matplotlib.pyplot as plt

# define a mathematical function
def func(x):
    return (1/x**0.5) + x**0.5

# define a function to integrate using trapezium rule with step size
def trapezium(xn,yn):
    # initialise the sum
    sum = 0
    # loop through the list of x values
    for i in range(len(xn)-1):
        # add the area of the trapezium to the sum
        sum += (xn[i+1]-xn[i])*(yn[i+1]+yn[i])/2
    # return the sum
    return sum

# x = np.arange(2,(23**0.5)+0.1,0.1)
# y = func(x)
# print(trapezium(x,y))
def doub(x,y):
    return np.exp(-x)+np.cos(y)
x = np.arange(0,1.01,0.01)
y = np.arange(0,2.01,0.01)
G = np.zeros(len(x))
for i in range(len(x)):
    G[i] = trapezium(y,doub(x[i],y))
print(trapezium(x,G))


# function to get kth derivative from x, y and k
def fwddiff(x,y,k):
    n = len(x)-1
    h = (x[-1]-x[0])/n
    df=np.zeros(n+1)
    for yn in range(n-k+1):
        for i in range(k+1):
            delta = ((-1)**i)*math.factorial(k)*y[yn+k-i]/(math.factorial(i)*math.factorial(k-i))
            df[yn] += delta/(h**k)
    return df
 
def newfunc(t):
    return np.exp(-t)+np.cos(t)**2

# t = np.arange(1,3.1,0.1)
# y = newfunc(t)
# print(fwddiff(t,y,4)[-11])

def lagrangian(j,xp,xn):
    '''takes j,xp and array xn and returns Lj'''
    length = len(xn)
    lj = 1
    for i in range(0,length):
        if i != j and xn[j] != xn[i]:
            lj = lj * (xp-xn[i])/(xn[j]-xn[i])
    return lj

def LagrInterp(xn,yn,x):
    '''receives sets of known yn, xn, and interpolates x'''
    order = len(xn)
    length = len(x)
    y = np.empty(length)
    for i in range(length):
        tmp = 0
        for j in range(order):
            tmp += yn[j] * lagrangian(j,x[i],xn)
        y[i] = tmp
    return y
