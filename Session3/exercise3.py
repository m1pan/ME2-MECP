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
    y = np.empty(length)
    for i in range(length):
        tmp = 0
        for j in range(order):
            tmp += yn[j] * lagrangian(j,x[i],xn)
        y[i] = tmp
    return y

domain = np.arange(0,3.05,0.05)
f = np.sin(domain)

xlin = np.linspace(1,2,2)
xquad = np.linspace(1,2,3)
xcub = np.linspace(1,2,4)

linterp = LagrInterp(xlin,np.sin(xlin),domain)
quadterp = LagrInterp(xquad,np.sin(xquad),domain)
cubterp = LagrInterp(xcub,np.sin(xcub),domain)
# plt.plot(domain,f)
# plt.plot(domain,quadterp,color='r')
# plt.show()
# print(domain[16])
# print(cubterp[16])

'''Task B'''
def NewtDivDiff(xn,yn):
    '''computes value of Newton's Divided difference from lists of xn and yn'''
    size=len(yn)
    if size==1:
        return yn[0]
    ynn = []
    for i in range(size-1):
        ynn.append((yn[i+1]-yn[i])/(xn[i+1+len(xn)-size]-xn[i]))
    return NewtDivDiff(xn,ynn)

def NewtonInterp(xn,yn,x):
    '''takes known xn,yn and returns interpolated y based on x'''
    order = len(xn)
    length = len(x)
    y = np.empty(length)
    for m in range(length):
        tmp = yn[0]
        for i in range(1,order):
            div = NewtDivDiff(xn[:i+1],yn[:i+1])
            for j in range(i):
                div *= (x[m]-xn[j])
            tmp += div
        y[m]=tmp
    return y

newtlin = NewtonInterp(xlin,np.sin(xlin),domain)
newtquad = NewtonInterp(xquad,np.sin(xquad),domain)
newtcub = NewtonInterp(xcub,np.sin(xcub),domain)

# plt.plot(domain,f)
# plt.plot(domain,newtquad,color='r')
# plt.show()
order11 = np.linspace(-1, 1, 12)
test = [1/(1+25*x**2) for x in order11]
# print(NewtonInterp(order11,test,[0.8]))


'''Task C'''
def Splines(xn,yn,x,ya,yb):
    """takes known xn,yn, points to be interpolated x, and boundary conditions y'(a) and y'(b). Returns interpolated y"""
    