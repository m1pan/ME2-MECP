import numpy as np
import matplotlib.pyplot as plt

'''Task A'''
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
    y = []
    for i in range(length):
        tmp = 0
        for j in range(order):
            tmp += yn[j] * lagrangian(j,x[i],xn)
        y.append(tmp)
    return y

x = np.arange(0,3.05,0.05)
f = np.sin(x)

xlin = np.linspace(1,2,2)
xquad = np.linspace(1,2,3)
xcub = np.linspace(1,2,4)

linterp = LagrInterp(xlin,np.sin(xlin),x)
quadterp = LagrInterp(xquad,np.sin(xquad),x)
cubterp = LagrInterp(xcub,np.sin(xcub),x)
# plt.plot(x,f)
# plt.plot(x,quadterp,color='r')
# plt.show()
# print(x[16])
# print(cubterp[16])

'''Task B'''
def NewtDivDiff(xn,yn):
    size=len(yn)
    if size==1:
        return yn
    ynn = []
    for i in range(size-1):
        ynn.append((yn[i+1]-yn[i])/(xn[i+1+len(xn)-size]-xn[i]))
    return NewtDivDiff(xn,ynn)