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


def TaskB():
    # myodebc takes boundaries of domain a,b, boundary conditions bca and bcb, number N of desired intervals, and array R specifying type of BCs
    # it returns the grid points xi and the solution yi at these points
    def myodebc(a,b,bca,bcb,N,R):
        xi = np.linspace(a,b,N+1)
        h = (b-a)/N
        p = np.zeros(N+1)
        p[0] = bca
        p[-1] = bcb
        # initialise the matrix lineq of linear equations
        lineq = np.zeros((N+1,N+1))
        lineq[0,0],lineq[0,1] = R[1] - R[0]/h, R[0]/h
        lineq[-1,-2],lineq[-1,-1] = -R[2]/h, R[2]/h + R[3]

        # fill matrix with unknowns and yi with knowns
        for i in range(1,N):
            lineq[i,i-1] = 1/h**2 - funcB(xi[i])[0]/(2*h)
            lineq[i,i] = funcB(xi[i])[1] - 2/h**2
            lineq[i,i+1] = 1/h**2 + funcB(xi[i])[0]/(2*h)
            p[i] = funcB(xi[i])[2]
        y = np.linalg.solve(lineq,p)
        return xi,y
    
    def funcB(x):
        f = x
        g = 1
        p = 5*x
        return f,g,p
    
    # conditions b2
    R = np.array([1,0,0,1])
    x,y = myodebc(0,2,0,5,50,R)
    print(x[3],y[3])

    return 0

def TaskC():
    # myodebc takes boundaries of domain a,b, boundary conditions bca and bcb, number N of desired intervals, and array R specifying type of BCs
    # it returns the grid points xi and the solution yi at these points
    def myodebc(a,b,bca,bcb,N,R):
        xi = np.linspace(a,b,N+1)
        h = (b-a)/N
        p = np.zeros(N+1)
        p[0] = bca
        p[-1] = bcb
        # initialise the matrix lineq of linear equations
        lineq = np.zeros((N+1,N+1))
        lineq[0,0],lineq[0,1] = R[1] - R[0]/h, R[0]/h
        lineq[-1,-2],lineq[-1,-1] = -R[2]/h, R[2]/h + R[3]

        # fill matrix with unknowns and yi with knowns
        for i in range(1,N):
            f,g,p[i] = funcC(xi[i]) 
            lineq[i,i-1] = 1/h**2 - f/(2*h)
            lineq[i,i] = g - 2/h**2
            lineq[i,i+1] = 1/h**2 + f/(2*h)

        y = np.linalg.solve(lineq,p)
        return xi,y
    
    #constant values
    rad = 15e-3
    w = 3e-3
    k = 16.75
    Tw = 490

    def funcC(x):
        f = 1/x
        g = 0
        p = -(10**8)*np.exp(-x/rad)/(x*k)
        return f,g,p

    x,y = myodebc(rad,rad+w,-6.32*10**5/k,Tw,50,np.array([1,0,k/6e4,1]))
    print(x[40],y[40])
    plt.plot(x,y)
    plt.show()

def TaskD():
    
if __name__ == '__main__':
    TaskC()