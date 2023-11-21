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

interval = 100
t = np.arange(0,1300+interval,interval)
# plt.plot(t,data)
# plt.show()

# function to get kth derivative from x, y and k
def fwddiff(x,y,k):
    n = len(x)-1
    h = (x[-1]-x[0])/n
    df=np.zeros(n+1)
    for yn in range(len(y)):
        for i in range(k):
            delta = ((-1)**i)*np.math.factorial(k)*y[yn+k-i-1]/(np.math.factorial(i)*np.math.factorial(k-i))
            df[yn] = delta/(h**k)
    return df

tspline = np.linspace(0,1300,140)

def splines(xn,yn,yprime,x):
    splined = np.zeros(len(x))
    for i in range(len(x)):
        for j in range(len(xn)-1):
            if x[i] >= xn[j] and x[i] <= xn[j+1]:
                splined[i] = yn[j] + yprime[j]*(x[i]-xn[j]) + (3*(yn[j+1]-yn[j])/(xn[j+1]-xn[j])**2 - (yprime[j+1]+2*yprime[j])/(xn[j+1]-xn[j]))*(x[i]-xn[j])**2 + (-2*(yn[j+1]-yn[j])/(xn[j+1]-xn[j])**3 + (yprime[j+1]+yprime[j])/(xn[j+1]-xn[j])**2)*(x[i]-xn[j])**3
    return splined

plt.plot(tspline,splines(t,data,fwddiff(t,data,1),tspline))
plt.plot(t,data)
plt.show()

# print(splines(t,data,fwddiff(t,data,1),tspline)[121])
