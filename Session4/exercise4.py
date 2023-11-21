import numpy as np
import matplotlib.pyplot as plt

def Simpson(x,y):
    a=x[0]
    b=x[-1]
    n=len(x)-1
    h = (b-a) / n
    return (h / 3) * (y[0] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-1:2]) + y[-1])

x = np.arange(0,2.01,0.01)
y = np.power(4,x-1)
# print(Simpson(x,y))

def Simpson2(f,a,b,e):
    n=2
    h = (b-a) / n
    S1 = (h / 3) * (f(a) + 4 * f((a+b)/2) + f(b))
    while True:
        n *= 2
        h = (b-a) / n
        S2 = S1
        x = np.linspace(a,b,n+1)
        S1 = (h/3)*(f(a) + 4*np.sum(f(x[1:-1:2])) + 2*np.sum(f(x[2:-1:2])) + f(b))
        if (1/15)*np.abs(S1 - S2) < e:
            break
    return S1,n+1
# print(Simpson2(lambda x: 1/(1+(x**2)),0,2,0.0001))

with open('Rocket.txt','r') as f:
    data = f.readlines()
    data = [line.strip().split() for line in data]
    data = np.array(data,dtype=float)