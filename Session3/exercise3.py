import numpy as np
import matplotlib.pyplot as plt

def lagrangian(j,xp,xn):
    length = len(xp)
    lj = (xp - xn[0])/(xn[j]-xn[0])
    for i in range(1,length):
        if i != j:
            lj = lj * (xp-xn[i])/(xn[j]-xn[i])
    return lj
