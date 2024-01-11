import numpy as np
import matplotlib.pyplot as plt

def TaskA():
    # myodebc takes boundaries of domain a,b, the value of solutions at these points ya, yb, number N of desired intervals, and function func of the external function
    # it returns the grid points xi and the solution yi at these points
    def myodebc(a,b,ya,yb,N):
        xi = np.linspace(a,b,N+1)
        p = np.zeros(N+1)
        p[0] = ya
        p[-1] = yb
        h = (b-a)/N
        # initialise the matrix lineq of linear equations
        lineq = np.zeros((N+1,N+1))
        (lineq[0,0],lineq[-1,-1]) = (1,1)
        # fill matrix with unknowns and yi with knowns
        for i in range(1,N):
            lineq[i,i-1] = 1/h**2 - funcA(xi[i])[0]/(2*h)
            lineq[i,i] = funcA(xi[i])[1] - 2/h**2
            lineq[i,i+1] = 1/h**2 + funcA(xi[i])[0]/(2*h)
            p[i] = funcA(xi[i])[2]
        y = np.linalg.solve(lineq,p)
        return xi,y

    def funcA(x):
        f = 2*x
        g = 2
        p = np.cos(3*x)
        return f,g,p

    x,y = myodebc(0,np.pi,1.5,0,10)
    print(x[4],y[4])
    plt.plot(x,y)
    plt.show()
    return 0

if __name__ == '__main__':
    TaskA()
