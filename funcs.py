import numpy as np
import matplotlib.pyplot as plt

# function trapzeqd: compute numerical integration with trapezium rule with any h interval
def traprule(x,y):
    """trapezium rule for numerical integration with non-uniform intervals"""
    return np.sum((x[1:] - x[:-1]) * (y[1:] + y[:-1]) / 2)


# plot surface
ax = plt.axes(projection='3d')
ax.plot_surface(Xg,Yg,s)
ax.plot_wireframe(Xg,Yg,s)
plt.contour(Xg,Yg,r)

# plot vectors
plt.quiver(Xg,Yg,f1,f2)
plt.streamplot(Xg,Yg,f1,f2)

## Interpolation

# Lagrange interpolation
def lagrangian(j,xp,xn):
    '''takes node number j, point to be interpolated xp and known array xn and returns 
        lagrange polynomial Lj'''
    length = len(xn)
    lj = 1
    for i in range(0,length):
        if i != j and xn[j] != xn[i]:
            lj = lj * (xp-xn[i])/(xn[j]-xn[i])
    return lj

def LagrInterp(xn,yn,x):
    '''receives sets of known yn, xn, and interpolates x, returns interpolated y'''
    order = len(xn)
    length = len(x)
    y = np.empty(length)
    for i in range(length):
        tmp = 0
        for j in range(order):
            tmp += yn[j] * lagrangian(j,x[i],xn)
        y[i] = tmp
    return y

# forward Newton's divided difference
def NewtDivDiff(xn,yn):
    '''computes value of Newton's Divided difference from lists of xn and yn'''
    size=len(yn)
    xsize = len(xn)
    if size==1:
        return yn[0]
    ynn = []
    for i in range(size-1):
        ynn.append((yn[i+1]-yn[i])/(xn[i+1+xsize-size]-xn[i]))
    return NewtDivDiff(xn,ynn)


def NewtonInterp(xn,yn,x):
    '''takes known arrays xn,yn and returns interpolated array y based on x'''
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

# Cubic splines
def splines(
    xn: list[float] | np.ndarray,
    yn: list[float] | np.ndarray,
    x: list[float] | np.ndarray,
    grad_a: float,
    grad_b: float,
) -> list[float]:
    """Write a function, Splines, that receives the sets of know values,
    xn and yn, the points to be interpolated x, the clamped boundary
    conditions y'(a), y'(b), and returns the interpolated values y,
    by using cubic splines

    Args:
        xn (list[float] | np.ndarray): nodes
        yn (list[float] | np.ndarray): values
        x (list[float] | np.ndarray): points to interpolate
        grad_a (float): gradient at a
        grad_b (float): gradient at b

    Returns:
        list[float]: interpolated values
    """
    n = len(xn)
    h = np.diff(xn)
    # create matrix for gradients
    a = np.zeros((n, n))
    d = np.zeros(n)
    a[0, 0] = 1
    a[-1, -1] = 1
    d[0] = grad_a
    d[-1] = grad_b
    for j in range(1, n - 1):
        a[j, j - 1] = 1 / h[j - 1]
        a[j, j] = 2 * (1 / h[j - 1] + 1 / h[j])
        a[j, j + 1] = 1 / h[j]
        left = (yn[j] - yn[j - 1]) / (h[j - 1] ** 2)
        right = (yn[j + 1] - yn[j]) / (h[j] ** 2)
        d[j] = 3 * (left + right)
    # solve system
    v = np.linalg.solve(a, d)
    # create matrix for coefficients
    c = np.zeros((n - 1, 4))
    for j in range(n - 1):
        dy = yn[j + 1] - yn[j]
        c[j, 0] = yn[j]
        c[j, 1] = v[j]
        c[j, 2] = 3 * dy / h[j] ** 2 - (v[j + 1] + 2 * v[j]) / h[j]
        c[j, 3] = -2 * dy / h[j] ** 3 + (v[j + 1] + v[j]) / h[j] ** 2

    # interpolate
    out = []
    for xp in x:
        for j in range(n - 1):
            if xn[j] <= xp <= xn[j + 1]:
                dx = xp - xn[j]
                out.append(
                    c[j, 0] + c[j, 1] * dx + c[j, 2] * dx**2 + c[j, 3] * dx**3
                )
                break
    return out