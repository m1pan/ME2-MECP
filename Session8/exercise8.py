import numpy as np
import matplotlib.pyplot as plt

def TaskA():
    # plot cosine wave and lagged wave
    amplitude = 10
    frequency = 0.5*2*np.pi #amgular frequency
    phi = np.pi/4
    t = np.arange(0,np.pi+0.1,0.1)
    y1 = amplitude*np.cos(frequency*t)
    y2 =amplitude*np.cos(2*frequency*t-phi)
    y12 = y1+y2
    
    print(y12[1])
    # plt.plot(t,y1,label='y1')
    # plt.plot(t,y2,label='y2')
    # plt.plot(t,y12,label='y1+y2')
    # plt.show()

def TaskB():
    # analogue filters and bode plot
    # task i

    # w = np.linspace(0,10000,1000)
    # H = 1/(1+0.1j*w)
    # dB = 20*np.log10(abs(H))
    # plt.plot(w,dB)
    # plt.show()

    # task ii
    w = 0.1
    zh = 1/((-1j/(w*1e-3))+1000)
    zv = 1/((-1j/(w*2e-3))+2000)
    H = 1/(1+(zh/zv))
    print(20*np.log10(abs(H)))

def TaskC():
    #task ii
    t = 0.2
    N = 8
    T = 5
    y = 0 
    for n in range(N):
        y+= (4/np.pi) * np.sin(2*n*np.pi*t/T)/N
    print(y)

def TaskD():
    def DFT(yn):
        '''Discrete Fourier Transform, receives a set of numerical values yn
        and return the Discrete Fourier Transform of the input, FT'''
        N = len(yn)
        FT = np.zeros(N,dtype=complex)
        for k in range(N):
            for n in range(N):
                FT[k] += yn[n]*np.exp(-2j*np.pi*k*n/N)
        return FT
    
    def DFTInv(FT):
        '''Inverse Discrete Fourier Transform, receives a set of numerical values FT
        and return the Inverse Discrete Fourier Transform of the input, yn'''
        N = len(FT)
        yn = np.zeros(N,dtype=complex)
        for n in range(N):
            for k in range(N):
                yn[n] += (1/N)*FT[k]*np.exp(2j*np.pi*k*n/N)
        return yn

    # task iii (e)
    dt = 0.1
    t = np.arange(0, 6*np.pi+dt, dt)
    N = len(t)
    y = np.exp(-0.25*(t-5)**2)
    df = 1/(N*dt)
    f = np.arange(0,1/dt,df)
    FT = DFT(y)
    # index for f = 0.2105 is 4
    print(abs(FT[4]))

def TaskE():
    with open ('vibration.txt','r') as file:
        vib = [float(i) for i in file.read().split()]
    
    def DFT(yn):
        '''Discrete Fourier Transform, receives a set of numerical values yn
        and return the Discrete Fourier Transform of the input, FT'''
        N = len(yn)
        FT = np.zeros(N,dtype=complex)
        for k in range(N):
            for n in range(N):
                FT[k] += yn[n]*np.exp(-2j*np.pi*k*n/N)
        return FT
    
    def DFTInv(FT):
        '''Inverse Discrete Fourier Transform, receives a set of numerical values FT
        and return the Inverse Discrete Fourier Transform of the input, yn'''
        N = len(FT)
        yn = np.zeros(N,dtype=complex)
        for n in range(N):
            for k in range(N):
                yn[n] += (1/N)*FT[k]*np.exp(2j*np.pi*k*n/N)
        return yn
    # N = len(vib)
    # dt = 0.01
    # df = 1/(N*dt)
    # f = np.arange(0,1/dt,df)
    # t = np.arange(0,N*0.01,0.01)
    # ft = DFT(vib)
    # plt.plot(f[:int(N/2)],np.abs(ft[:int(N/2)]))
    # plt.show()

    with open('noisy.txt','r') as file:
        noisy = [float(i) for i in file]
    
    N = len(noisy)
    dt = 1/20
    df = 1/(N*dt)
    f = np.arange(0,N*df,df)
    t = np.arange(0,N*dt,dt)
    ft = DFT(noisy)
    fc = 0.5 #cut off frequency
    H = 1/(1+1j*f/fc)
    invft = ft*H
    ynew = DFTInv(invft)

    print(abs(ynew[100]))
    # plt.plot(t,noisy)
    # # plt.show()
    # # plt.plot(f[:int(N/2)],np.abs(ft[:int(N/2)]))
    # # plt.show()
    # plt.plot(t,ynew)
    # plt.show()

    

if __name__ == "__main__":
    TaskB()
