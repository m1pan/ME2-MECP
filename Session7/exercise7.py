import numpy as np
import matplotlib.pyplot as plt


'''task a'''
def TaskA():
    ALPHA = 1.172e-5

    # pde takes in initial BCs at t=0, boundary conditions at domain a,b, initial uniform temp T
    # step size of dx and dt, and final time t_end in seconds
    # it returns the grid points xi and the solution T(x,t) at these points
    def pde(Ta,Tb,T0,dx,dt,a,b,t_end):
        # number of grid points
        N = int((b-a)/dx) + 1
        # time interval spacing
        M = int(t_end/dt) + 1
        # initialise the grid
        xi = np.arange(a,b+dx,dx)
        ti = np.arange(0,t_end+dt,dt)
        # initialise the solution matrix
        T = np.zeros((M,N))
        # set initial conditions
        T[0,:]= T0
        T[0,0]=Ta
        T[0,-1]=Tb
        for p in range(1,M):
            T[p,0]=Ta
            T[p,-1]=Tb
            for i in range(1,N-1):
                T[p,i] = (ALPHA * dt/(dx**2)) * (T[p-1,i+1] + T[p-1,i-1]) + (1 - 2*ALPHA*dt/(dx**2))*T[p-1,i]
        return xi,ti,T
 
    x,t,T = pde(50,70,10,0.01,1,0,0.5,3600)
    print(T[-1,10])
    for i in range(0,3600):
        plt.plot(x,T[i,:])
    plt.show()
    return 0

'''task b'''
def TaskB():
    ALPHA = 1.172e-5
    h = 500
    k = 40
    # pde takes in boundary conditions at centre, initial uniform temp T0, water temp Tw
    # step size of dx and dt, and final time t_end in seconds
    # it returns the grid points xi and the solution T(x,t) at these points
    def pde(Tmid,T0,Tw,dx,dt,a,b,t_end):
        # number of grid points
        N = int((b-a)/dx) + 1
        # time interval spacing
        M = int(t_end/dt) + 1
        # initialise the grid
        xi = np.arange(a,b+dx,dx)
        ti = np.arange(0,t_end+dt,dt)
        # initialise the solution matrix
        T = np.zeros((M,N))
        # set initial conditions
        T[0,:]= T0
        for p in range(1,M):
            T[p,0]=(h*Tw + k*T[p-1,1]/dx)/(h+k/dx)
            T[p,-1]=(h*Tw + k*T[p-1,-2]/dx)/(h+k/dx)
            for i in range(1,N-1):
                if i != int(N/2):
                    T[p,i] = (ALPHA * dt/(dx**2)) * (T[p-1,i+1] + T[p-1,i-1]) + (1 - 2*ALPHA*dt/(dx**2))*T[p-1,i]
            T[p,int(N/2)] = Tmid
        return xi,ti,T

    x,t,T = pde(100,10,5,0.01,1,0,0.5,1200)
    for i in range(0,1200):
        if int(T[i,44])==30:
            print(t[i])
            return 0
    
    plt.plot(t,T[:,44])
    plt.show()


if __name__ == "__main__":
    TaskB()
