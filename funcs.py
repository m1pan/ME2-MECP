import numpy as np
import matplotlib.pyplot as plt

# function trapzeqd: compute numerical integration with trapezium rule with equidistant nodes
def trapzeqd(x,y):
    # get the interval h: distance between any two consecutives nodes
    h = x[1] - x[0]
    # compute the integral
    I = h * (y[0]/2 + np.sum(y[1:-1]) + y[-1]/2 )
    
    return I

# plot surface
ax = plt.axes(projection='3d')
ax.plot_surface(Xg,Yg,s)
ax.plot_wireframe(Xg,Yg,s)
plt.contour(Xg,Yg,r)

# plot vectors
plt.quiver(Xg,Yg,f1,f2)
plt.streamplot(Xg,Yg,f1,f2)


def lagrangian(j: int, xp: float, xn: list[float] | np.ndarray) -> float:
    """Write a function, Lagrangian, to compute the Lagrangian polynomial
    j at a point xp, with given nodes xn.
    The function receives the values j, xp and the array of nodes xn, and
    returns the value:
    Lj(xp) = product_{k=0, k!=j}^n (xp - xk) / (xj - xk)

    Args:
        j (int): order of the polynomial
        xp (float): point to evaluate
        xn (list[float] | np.ndarray): nodes

    Returns:
        float: value of the polynomial at xp
    """
    out = 1
    for k, xk in enumerate(xn):
        if k != j:
            out *= (xp - xk) / (xn[j] - xk)

    return out

