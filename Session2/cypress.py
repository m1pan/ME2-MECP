import math as mt
import numpy as np
import matplotlib.pyplot as plt

#Task A
#Q1
def f(x):
    y = 1/mt.sqrt(pow(x,17.1)+2023)
    return y

def trapzeqd(x,y):
    dy = len(y)-1
    dx=x[1]-x[0]
    sum = 0
    for n in range(dy):
        sum += (y[n]+y[n+1])*(dx)/2
    return sum

b=2
x = np.linspace(0,b,10)

y=[]
for n in x:
    y.append(f(n))

print(x)
print(trapzeqd(x,y))

#Q2
b = [10,100,1000,10000]
I = []

for i in b:
    x = np.linspace(0,i,5)

    y=[]
    for n in x:
        y.append(f(n))
    
    I.append(trapzeqd(x,y))

#plt.scatter(b,I)
#plt.show()

#Q3
b = [10,100,1000,10000]
I = []

for i in b:
    h=0.5
    x = np.arange(0,i+h,h)

    y=[]
    for n in x:
        y.append(f(n))
    
    I.append(trapzeqd(x,y))

#plt.scatter(b,I)
#plt.show()

#Task B
#Q4
def g(x):
    y = 1/mt.sqrt(pow(x,1.1)+2023)
    return y
b = [10,100,1000,10000]
I = []

for i in b:
    x = np.linspace(0,i,5)

    y=[]
    for n in x:
        y.append(g(n))
    
    I.append(trapzeqd(x,y))

#plt.scatter(b,I)
#plt.show()


b = [10,100,1000,10000]
I = []

for i in b:
    h=0.5
    x = np.arange(0,i+h,h)

    y=[]
    for n in x:
        y.append(g(n))
    
    I.append(trapzeqd(x,y))

#plt.scatter(b,I)
#plt.show()

#Task C
def trapz(x,y):
    ln = len(x)-1
    I = 0
    for n in range(ln):
        I += (y[n+1]+y[n])*(x[n+1]-x[n])/2
    return I

#Task D
f = open("Thames.txt","r")
a = f.readlines()

xn = []
yn = []
xs = []
ys = []

for n in a:
    x = n.split(",")
    ct = 0
    for n in x:
        if ct == 0:
            xn.append(float(n))
            ct += 1
        elif ct == 1:
            yn.append(float(n))
            ct += 1
        elif ct == 2:
         xs.append(float(n))
         ct += 1
        else:
            ys.append(float(n))
            ct = 0

#plt.plot(xn,yn)
#plt.plot(xs,ys)
#plt.axis('equal')
#plt.show()

print(trapz(xn,yn)-trapz(xs,ys))

#Task E
# def z(x,y):
#     z = 
#     for n in x:
#         for m in y:
#         z = mt.sqrt(pow(25,2)*(1-pow(x/67,2)-pow(y/56,2)))
#     return z

# h=0.5
# x = np.arange(-67,67+h,h)
# y = np.arange(-56,56+h,h)

# xx,yy = np.meshgrid(x,y)
# z = z(xx,yy)