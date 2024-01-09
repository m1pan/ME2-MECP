import numpy as np
import matplotlib.pyplot as plt

def TaskA():
    # myodebc takes boundaries of domain a,b, the value of solutions at these points ya, yb, number N of desired intervals, and function func of the external function
    # it returns the grid points xi and the solution yi at these points
    def myodebc(a,b,ya,yb,N,func):
        xi = np.linspace(a,b,N+1)
        yi = np.zeros(N+1)
        yi[0] = ya
        yi[-1] = yb
        h = (b-a)/N
        # initialise the matrix lineq of linear equations
        lineq = np.zeros((N+1,N+1))
        (lineq[0,0],lineq[-1,-1]) = (1,1)

    def funcA(x,y):
        yn = np.zeros(3)
    return 0

if __name__ == '__main__':
    TaskA()