import numpy as np
import matplotlib.pyplot as plt

def TaskA():
    # Task A and B
    # define a differential equation
    def func(y,t):
        return -2*y*t-3*t**3

    #define a function to solve a general ODE by adopting a forward euler numerical scheme
    # def FwdEuler(t0,y0,tend,h):
    #     t = np.arange(t0,tend+h,h)
    #     y = np.zeros(len(t))
    #     y[0] = y0
    #     for i in range(len(t)-1):
    #         y[i+1] = y[i] + h*func(y[i],t[i])
    #     return t,y

    def FwdEuler(t0,y0,tend,h):
        t = np.arange(t0,tend+h,h)
        y = np.zeros(len(t))
        y[0] = y0
        for i in range(len(t)-1):
            try:
                y[i+1] = y[i] + h*func(y[i],t[i])
            except OverflowError:
                y[i+1] = np.inf
        return t,y

    # define a function to solve a general ODE by adopting a 4th order Runge-Kutta numerical scheme
    def ORDER4(t0,y0,tend,h):
        t = np.arange(t0,tend+h,h)
        y = np.zeros(len(t))
        y[0] = y0
        for i in range(len(t)-1):
            k1=h*func(y[i],t[i])
            k2=h*func(y[i]+k1/2,t[i]+h/2)
            k3=h*func(y[i]+k2/2,t[i]+h/2)
            k4=h*func(y[i]+k3,t[i]+h)
            y[i+1]=y[i]+(k1+2*k2+2*k3+k4)/6
        return t,y

    #define a function to solve a general ODE by adopting a backward euler numerical scheme
    def BwEuler(t0,y0,tend,h):
        t = np.arange(t0,tend+h,h)
        y = np.zeros(len(t))
        y[0] = y0
        for i in range(len(t)-1):
            y[i+1] = (y[i]-2*h*t[i+1]**3)/(1+2*h*t[i+1])
        return t,y 


    # print(FwdEuler(0,800,1,0.1)[1][-1])
    # print(ORDER4(0,800,1,0.1)[1][-1])
    # print(BwEuler(0,800,1,0.1)[1][-1])
    # plt.plot(FwdEuler(0,-10,100,0.1)[0],FwdEuler(0,-10,100,0.1)[1],label='Forward Euler')
    # plt.plot(ORDER4(0,-10,100,0.1)[0],ORDER4(0,-10,100,0.1)[1],label='4th Order Runge-Kutta')
    # plt.plot(BwEuler(0,-10,100,0.1)[0],BwEuler(0,-10,100,0.1)[1],label='Backward Euler')
    # plt.show()
    return 0

def TaskC():
    # Task C2
    # define a forward euler function to solve a system of ODEs
    def FwEulerN(t0,tend,y0,h,func):
        t = np.arange(t0,tend+h,h)
        y = np.zeros((len(t),len(y0)))
        y[0] = y0
        for i in range(len(t)-1):
            y[i+1] = y[i] + h*func(y[i])
        return t,y

    def housemarket(y):
        yn=np.zeros(len(y))
        yn[0]=0.3*y[0]*y[1]-0.8*y[0]
        yn[1]=1.1*y[1]-y[0]*y[1]
        return yn

    # initial conditions
    incon = np.array([0.8,7])
    t,y = FwEulerN(0,40,incon,0.005,housemarket)
    print(y[3999])
    plt.plot(t,y[:,0],label='House Price')
    plt.plot(t,y[:,1],label='Houses Sold')
    plt.legend()
    plt.show()
    return 0


# Task D1
def TaskD():
    def pendulum(y,t):
        a
