from typing import Callable, Optional
import numpy as np
import matplotlib.pyplot as plt

# plot surface
ax = plt.axes(projection='3d')
ax.plot_surface(Xg,Yg,s)
ax.plot_wireframe(Xg,Yg,s)
plt.contour(Xg,Yg,r)

# plot vectors
plt.quiver(Xg,Yg,f1,f2)
plt.streamplot(Xg,Yg,f1,f2)

def trapz(x: np.ndarray, y: np.ndarray) -> float:
    """Trapezium integration

    Args:
        x (np.ndarray): x values
        y (np.ndarray): y values

    Returns:
    float: integral value
    """
    return sum(
        (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2 for i in range(len(x) - 1)
    )


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


def lagr_interp(
    xn: list[float] | np.ndarray,
    yn: list[float] | np.ndarray,
    x: list[float] | np.ndarray,
) -> list[float]:
    """Write a function, LagrInterp, that receives the sets of know values,
    xn and yn, the points to be interpolated x, and returns the
    interpolated values y, by using Lagrangian polynomials.

    Args:
        xn (list[float] | np.ndarray): nodes
        yn (list[float] | np.ndarray): values
        x (list[float] | np.ndarray): points to interpolate

    Returns:
        list[float]: interpolated values
    """
    out = []
    for xp in x:
        y = 0
        for j in range(len(xn)):
            y += yn[j] * lagrangian(j, xp, xn)
        out.append(y)
    return out


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


def simpson(h: float, yn: list[float] | np.ndarray) -> float:
    """Simpson integration on equally spaced nodes

    Args:
        h (float): step size
        yn (list[float] | np.ndarray): values

    Returns:
        float: integral value by Simpson method
    """
    return sum(
        h / 3 * (yn[i] + 4 * yn[i + 1] + yn[i + 2])
        for i in range(0, len(yn) - 2, 2)
    )

#adaptive simpson, takes in function f, lower and upper bounds a and b, desired absolute error e
#returns the integral value and number of nodes used
def Simpson2(f,a,b,e):
    n=2
    h = (b-a) / n
    S1 = (h / 3) * (f(a) + 4 * f((a+b)/2) + f(b))
    while True:
        n *= 2
        h = (b-a) / n
        S2 = S1
        x = np.linspace(a,b,n+1)
        S1 = (h/3)*(f(a) + 4*np.sum(f(x[1:-1:2])) + 2*np.sum(f(x[2:-1:2])) + f(b))
        if (1/15)*np.abs(S1 - S2) < e:
            break
    return S1,n+1


def factorial(n: int) -> int:
    """Factorial of an integer n

    Args:
        n (int): integer

    Returns:
        int: factorial of n
    """
    if n == 0:
        return 1
    return n * factorial(n - 1)


def choose(n: int, k: int) -> int:
    """Binomial

    Args:
        n (int): n
        k (int): k

    Returns:
        int: binomial coefficient
    """
    return factorial(n) // (factorial(k) * factorial(n - k))


def derivative_fwd(
    k: int, h: float, yn: list[float] | np.ndarray
) -> list[float] | np.ndarray:
    """k-th order derivative using forward method

    Args:
        k (int): order
        h (float): step size
        yn (list[float] | np.ndarray): values

    Returns:
        list[float] | np.ndarray: n - k derivative values
    """
    nodes = len(yn) - k
    dx = np.ndarray(nodes)
    for n in range(nodes):
        dx[n] = sum(
            (-1) ** i * choose(k, i) * yn[n + k - i] for i in range(k + 1)
        )
    return dx / h**k


def derivative_bwd(k: int, h: float, yn: list[float] | np.ndarray):
    """k-th order derivate using backward method

    Args:
        k (int): order
        h (float): step size
        yn (list[float] | np.ndarray): values

    Returns:
        list[float] | np.ndarray: n - k derivative values
    """
    nodes = len(yn) - k
    dx = np.ndarray(nodes)
    for n in range(k, len(yn)):
        dx[n] = sum((-1) ** i * choose(k, 1) * yn[n - i] for i in range(k + 1))
    return dx / h**k


def gauss_quad(func: Callable, a: float, b: float, n: int) -> float:
    """Gauss quadrature integration

    Args:
        func (Callable): function to integrate
        a (float): lower bound
        b (float): upper bound
        n (int): number of nodes

    Returns:
        float: integral value
    """
    tg = np.array(
        [
            [0],
            [1 / np.sqrt(3), -1 / np.sqrt(3)],
            [0, np.sqrt(3 / 5), -np.sqrt(3 / 5)],
            [
                np.sqrt(3 / 7 - 2 / 7 * np.sqrt(6 / 5)),
                -np.sqrt(3 / 7 - 2 / 7 * np.sqrt(6 / 5)),
                np.sqrt(3 / 7 + 2 / 7 * np.sqrt(6 / 5)),
                -np.sqrt(3 / 7 + 2 / 7 * np.sqrt(6 / 5)),
            ],
            [
                0,
                1 / 3 * np.sqrt(5 - 2 * np.sqrt(10 / 7)),
                -1 / 3 * np.sqrt(5 - 2 * np.sqrt(10 / 7)),
                1 / 3 * np.sqrt(5 + 2 * np.sqrt(10 / 7)),
                -1 / 3 * np.sqrt(5 + 2 * np.sqrt(10 / 7)),
            ],
        ]
    )

    wg = np.array(
        [
            [2],
            [1, 1],
            [8 / 9, 5 / 9, 5 / 9],
            [
                (18 + np.sqrt(30)) / 36,
                (18 + np.sqrt(30)) / 36,
                (18 - np.sqrt(30)) / 36,
                (18 - np.sqrt(30)) / 36,
            ],
            [
                128 / 225,
                (322 + 13 * np.sqrt(70)) / 900,
                (322 + 13 * np.sqrt(70)) / 900,
                (322 - 13 * np.sqrt(70)) / 900,
                (322 - 13 * np.sqrt(70)) / 900,
            ],
        ]
    )

    i: int = n - 1

    return (
        (b - a)
        / 2
        * sum(
            wg[i][j] * func((a * (1 - tg[i][j]) + b * (1 + tg[i][j])) / 2)
            for j in range(n)
        )
    )


def fwd_euler(
    func: Callable[[float, float], float],
    t0: float,
    y0: float,
    t_end: float,
    h: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Forward Euler method

    Args:
        func (Callable): function with signature func(t, y) -> float
            equal to dy/dt
        t0 (float): initial time
        y0 (float): initial value
        t_end (float): final time
        h (float): step size

    Returns:
        tuple[np.ndarray, np.ndarray]: times and values
    """
    t: np.ndarray = np.arange(t0, t_end + h, h)
    y: np.ndarray = np.zeros(len(t))
    y[0] = y0
    for i in range(len(t) - 1):
        y[i + 1] = y[i] + h * func(t[i], y[i])
    return t, y


# bwdeuler
"""rearrrange for y_n+1, proceed from there similar to fwd_euler     """


def ode_rk4(
    func: Callable, t0: float, y0: float, t_end: float, h: float
) -> tuple[np.ndarray, np.ndarray]:
    """Runge-Kutta 4th order method

    Args:
        func (Callable): function with signature func(t, y) -> float
            equal to dy/dt
        t0 (float): initial time
        y0 (float): initial value
        t_end (float): final time
        h (float): step size

    Returns:
        tuple[np.ndarray, np.ndarray]: times and values
    """
    t: np.ndarray = np.arange(t0, t_end + h, h)
    y: np.ndarray = np.zeros(len(t))
    y[0] = y0
    for i in range(len(t) - 1):
        k1 = h * func(t[i], y[i])
        k2 = h * func(t[i] + h / 2, y[i] + k1 / 2)
        k3 = h * func(t[i] + h / 2, y[i] + k2 / 2)
        k4 = h * func(t[i] + h, y[i] + k3)
        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return t, y


def fwd_euler_n(
    funcs: list[Callable],
    t0: float,
    y0: np.ndarray,
    t_end: float,
    h: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Forward Euler method for a system of ODEs.

    For a system of N ODEs, the function receives a list of functions
    [dy1/dt, dy2/dt, ..., dyN/dt]. Each dyn/dt function receives the
    current time and the values of y1, y2, ..., yN, and returns the
    derivative of the corresponding variable.

    Initial values are given in the N x 1 array y0, and the function
    returns a 1 x Nt array of times and a N x Nt array of values,
    where Nt is the number of temporal nodes computed.

    Args:
        funcs (list[Callable]): list of functions with signature
            func(t, y1, y2, ..., yN) -> float equal to dy/dt
        t0 (float): initial time
        y0 (np.ndarray): initial values
        t_end (float): final time
        h (float): step size

    Returns:
        tuple[np.ndarray, np.ndarray]: times and values
    """
    t: np.ndarray = np.arange(t0, t_end + h, h)
    y: np.ndarray = np.zeros((len(t), len(y0)))
    y[0] = y0
    for i in range(1, len(t)):
        for j, func in enumerate(funcs):
            y[i, j] = y[i - 1, j] + h * func(t[i - 1], *y[i - 1])
    return t, y



def ode_bc(
    ode: Callable,
    a: float,
    b: float,
    y_a: float,
    y_b: float,
    N: int,
    R: np.ndarray,
) -> tuple:
    """Solve a second order boundary value problem ODE with the form
    d^2y/dx^2 + f(x) dy/dx + g(x) y = p(x)
    using central differences.
    Discretised domain between a and b with N nodes and boundary
    conditions y(a) = y_a and y(b) = y_b, and a 1 x 4 array R with
    the coefficients of the boundary conditions:
    r0 dy/dx (a) + r1 y(a) = y_a
    r2 dy/dx (b) + r3 y(b) = y_b
    where r is 1 or 0 depending on the type of boundary condition.

    Args:
        ode (Callable): function with signature ode(x) -> (f, g, p)
        a (float): lower bound
        b (float): upper bound
        y_a (float): value at a
        y_b (float): value at b
        N (int): number of nodes
        R (np.ndarray): boundary conditions

    Returns:
        tuple: x and y values
    """
    x: np.ndarray
    x, _h = np.linspace(a, b, N + 1, retstep=True)
    h = float(_h)
    A = np.zeros((N + 1, N + 1))
    B = np.zeros(N + 1)
    # r0 dy/dx (a) + r1 y(a) = y_a
    # r2 dy/dx (b) + r3 y(b) = y_b
    A[0, 0] = -R[0] / h + R[1]
    A[0, 1] = R[0] / h
    A[N, N - 1] = -R[2] / h
    A[N, N] = R[2] / h + R[3]
    B[0] = y_a
    B[N] = y_b
    for i in range(1, N):
        f, g, p = ode(x[i])
        A[i, i - 1] = 1 / h**2 - f / (2 * h)
        A[i, i] = -2 / h**2 + g
        A[i, i + 1] = 1 / h**2 + f / (2 * h)
        B[i] = p
    y = np.linalg.solve(A, B)
    return x, y

def heat_eqn(
    alpha: float,
    a: float,
    b: float,
    Ta: float,
    Tb: float,
    T0: float,
    dt: float,
    dx: float,
    t_end: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve the heat equation using the backward finite difference method
    to give the form: (backwards in time, central in space)
    (T(x, t + dt) - T(x, t)) / dt
        = alpha (T(x + dx, t) - 2 T(x, t) + T(x - dx, t)) / dx^2
    with boundary conditions T(a, t) = Ta, T(b, t) = Tb, and initial
    condition T(x, 0) = T0.

    The function returns the grid points x, time points t, and the
    temperature T(x, t) at each point.

    Args:
        alpha (float): thermal diffusivity
        a (float): lower bound of x
        b (float): upper bound of x
        Ta (float): boundary condition at a
        Tb (float): boundary condition at b
        T0 (float): initial condition
        dt (float): time step
        dx (float): space step
        t_end (float): final time

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]: x, t, T(x, t)
    """
    x_range = np.arange(a, b + dx, dx)
    t_range = np.arange(0, t_end + dt, dt)
    T = np.zeros((len(x_range), len(t_range)))
    T[:, 0] = T0
    T[0, :] = Ta
    T[-1, :] = Tb
    for i in range(len(t_range) - 1):
        for j in range(1, len(x_range) - 1):
            T[j, i + 1] = (
                alpha * dt / dx**2 * (T[j + 1, i] + T[j - 1, i])
                + (1 - 2 * alpha * dt / dx**2) * T[j, i]
            )
    return x_range, t_range, T

#2D heat diffusion
def heat_conduction_2d_potato(
    alpha_air: float,
    alpha_potato: float,
    potato_xlen: float,
    potato_ylen: float,
    xa: float,
    xb: float,
    ya: float,
    yb: float,
    Txa: float,
    Txb: float,
    Tya: float,
    Tyb: float,
    T0: float,
    dt: float,
    dx: float,
    dy: float,
    t_end: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x_range = np.arange(xa, xb + dx, dx)
    y_range = np.arange(ya, yb + dy, dy)
    t_range = np.arange(0, t_end + dt, dt)
    alpha = np.zeros((len(x_range), len(y_range)))
    alpha[:, :] = alpha_air
    xap = int((xa + xb - potato_xlen) / 2 / dx)
    xbp = int((xa + xb + potato_xlen) / 2 / dx)
    yap = int((ya + yb - potato_ylen) / 2 / dy)
    ybp = int((ya + yb + potato_ylen) / 2 / dy)
    alpha[xap:xbp, yap:ybp] = alpha_potato
    T = np.zeros((len(x_range), len(y_range), len(t_range)))
    T[:, :, 0] = T0
    T[0, :, :] = Txa
    T[-1, :, :] = Txb
    T[:, 0, :] = Tya
    T[:, -1, :] = Tyb
    T[xap:xbp, yap:ybp, 0] = -15
    for i in range(1, len(t_range)):
        for j in range(1, len(x_range) - 1):
            for k in range(1, len(y_range) - 1):
                T[j, k, i] = (
                    alpha[j, k]
                    * dt
                    / dx**2
                    * (T[j + 1, k, i - 1] + T[j - 1, k, i - 1])
                    + alpha[j, k]
                    * dt
                    / dy**2
                    * (T[j, k + 1, i - 1] + T[j, k - 1, i - 1])
                    + (
                        1
                        - 2 * alpha[j, k] * dt / dx**2
                        - 2 * alpha[j, k] * dt / dy**2
                    )
                    * T[j, k, i - 1]
                )
    return x_range, y_range, t_range, T

def dft(yn: list | np.ndarray) -> np.ndarray:
    """Discrete Fourier Transform on a set of values yn

    Args:
        yn (list | np.ndarray): values

    Returns:
        np.ndarray: DFT values
    """
    N = len(yn)
    k = np.arange(N)
    FT = np.sum(
        [yn * np.exp(-1j * 2 * np.pi * k * n / N) for n in range(N)], axis=1
    )
    return FT


def dft_inv(FT: list | np.ndarray) -> np.ndarray:
    """Inverse Discrete Fourier Transform on a set of values FT"""
    N = len(FT)
    n = np.arange(N)
    yn = (
        np.sum(
            [FT * np.exp(1j * 2 * np.pi * k * n / N) for k in range(N)],
            axis=1,
        )
        / N
    )
    return yn


def barycentric(
    r1: tuple[float, float],
    r2: tuple[float, float],
    r3: tuple[float, float],
    vals: tuple[float, float, float],
    rp: tuple[float, float],
) -> float:
    """Interpolate the value of a function at a point inside a triangle using
    the barycentric coordinates method.

    Args:
        r1 (tuple[float, float]): Coordinates of the first point.
        r2 (tuple[float, float]): Coordinates of the second point.
        r3 (tuple[float, float]): Coordinates of the third point.
        vals (tuple[float, float, float]): Values at the points.
        rp (tuple[float, float]): Coordinates of the point to interpolate.

    Returns:
        float: Interpolated value of the function at the point rp.
    """
    # barycentric coordinates method
    A: np.ndarray = np.array(
        [[r1[0], r2[0], r3[0]], [r1[1], r2[1], r3[1]], [1, 1, 1]]
    )
    b: np.ndarray = np.array([rp[0], rp[1], 1])
    alpha, beta, gamma = np.linalg.solve(A, b)
    return alpha * vals[0] + beta * vals[1] + gamma * vals[2]


def bisection(f: Callable, a: float, b: float, e: float) -> float:
    """Finds the root of a function given an interval and an accuracy

    Args:
        f (function): Function to find the root of
        a (float): Left bound of the interval
        b (float): Right bound of the interval
        e (float): Desired accuracy

    Returns:
        float: Root of the function
    """
    assert b > a
    assert f(a) * f(b) < 0
    while b - a > e:
        c = (a + b) / 2
        if f(c) == 0:
            return c
        elif f(a) * f(c) < 0:
            # root is in [a, c], f(a) and f(c) have different signs
            b = c
        else:
            a = c
    return (a + b) / 2


def derivative(f: Callable, h: float = 1e-6) -> Callable:
    """Calculates the derivative of a function

    Args:
        f (function): Function to calculate the derivative of
        h (float, optional): Step size. Defaults to 1e-6.

    Returns:
        function: Derivative of the function
    """

    def df(x: float) -> float:
        return (f(x + h) - f(x)) / h

    return df


def newton_raphson(
    f: Callable,
    x0: float,
    e: float,
    df: Optional[Callable] = None,
) -> float:
    """Finds the root of a function given an initial guess and an accuracy
    using the Newton-Raphson method

    Args:
        f (function): Function to find the root of
        df (function): Derivative of the function
        x0 (float): Initial guess
        e (float): Desired accuracy

    Returns:
        float: Root of the function
    """
    if df is None:
        df = derivative(f)
    max_iter = 100
    for _ in range(max_iter):
        x1 = x0 - f(x0) / df(x0)
        if abs(x1 - x0) < e:
            return x1
        x0 = x1
    raise ValueError("Method did not converge")


def newton_raphson_system(
    f: list[Callable],
    x0: np.ndarray,
    e: float,
) -> np.ndarray:
    """Newton-Raphson method for systems of non-linear equations

    Args:
        f (list[Callable]): list of functions
        x0 (list[float]): initial guess
        e (float): tolerance

    Returns:
        list[float]: solution
    """
    max_iter: int = 100
    for _ in range(max_iter):
        J: np.ndarray = jacobian(f, x0)
        x1: np.ndarray = x0 - np.linalg.inv(J) @ np.array(
            [_f(*x0) for _f in f]
        )
        if np.linalg.norm(x1 - x0) < e:
            return x1
        x0 = x1
    raise ValueError("Method did not converge")


def jacobian(
    f: list[Callable],
    x: np.ndarray,
    h: float = 1e-5,
) -> np.ndarray:
    """Calculate the Jacobian matrix of a system of non-linear equations at
    a given list of variables

    Args:
        f (list[Callable]): list of functions of size (n x 1)
        x (list[float]): list of variables of size (n x 1)
        h (float, optional): step size. Defaults to 1e-5.

    Returns:
        np.ndarray: Jacobian matrix
    """
    n: int = len(f)
    J: np.ndarray = np.zeros((n, n))
    for i in range(n):
        x_forward: np.ndarray = x.copy()
        x_backward: np.ndarray = x.copy()
        x_forward[i] += h
        x_backward[i] -= h
        J[:, i] = [(f_i(*x_forward) - f_i(*x_backward)) / (2 * h) for f_i in f]
    return J


def golden_section_search(f, a, b, e):
    """Finds the maximum of a function using the golden section search method

    Args:
        f (function): Function to find the maximum of
        a (float): Lower bound
        b (float): Upper bound
        e (float): Desired accuracy

    Returns:
        float: Maximum of the function
    """
    ratio = (5**0.5 - 1) / 2
    x1 = a + ratio * (b - a)
    x2 = b - ratio * (b - a)
    f1 = f(x1)
    f2 = f(x2)
    while abs(b - a) > e:
        if f1 > f2:
            a = x2
            x2 = x1
            f2 = f1
            x1 = a + ratio * (b - a)
            f1 = f(x1)
        else:
            b = x1
            x1 = x2
            f1 = f2
            x2 = b - ratio * (b - a)
            f2 = f(x2)
    return (a + b) / 2


def df_dx(f, x, y, h=1e-6):
    """Calculates the derivative of a function with respect to x

    Args:
        f (function): Function to calculate the derivative of
        x (float): x value
        y (float): y value
        h (float, optional): Step size. Defaults to 1e-6.

    Returns:
        float: Derivative of the function with respect to x
    """
    return (f(x + h, y) - f(x, y)) / h


def df_dy(f, x, y, h=1e-6):
    """Calculates the derivative of a function with respect to y

    Args:
        f (function): Function to calculate the derivative of
        x (float): x value
        y (float): y value
        h (float, optional): Step size. Defaults to 1e-6.

    Returns:
        float: Derivative of the function with respect to y
    """
    return (f(x, y + h) - f(x, y)) / h


def d2f_dx2(f, x, y, h=1e-6):
    """Calculates the second derivative of a function with respect to x

    Args:
        f (function): Function to calculate the derivative of
        x (float): x value
        y (float): y value
        h (float, optional): Step size. Defaults to 1e-6.

    Returns:
        float: Second derivative of the function with respect to x
    """
    return (f(x + h, y) - 2 * f(x, y) + f(x - h, y)) / h**2


def d2f_dy2(f, x, y, h=1e-6):
    """Calculates the second derivative of a function with respect to y

    Args:
        f (function): Function to calculate the derivative of
        x (float): x value
        y (float): y value
        h (float, optional): Step size. Defaults to 1e-6.

    Returns:
        float: Second derivative of the function with respect to y
    """
    return (f(x, y + h) - 2 * f(x, y) + f(x, y - h)) / h**2


def d2f_dxdy(f, x, y, h=1e-6):
    """Calculates the second derivative of a function with respect to x and y

    Args:
        f (function): Function to calculate the derivative of
        x (float): x value
        y (float): y value
        h (float, optional): Step size. Defaults to 1e-6.

    Returns:
        float: Second derivative of the function with respect to x and y
    """
    return (f(x + h, y + h) - f(x + h, y) - f(x, y + h) + f(x, y)) / h**2


def gradient_ascent(f, x0, y0, e):
    """Finds the maximum of a function given an initial guess and an accuracy
    using the gradient ascent method

    Args:
        f (function): Function to find the maximum of
        x0 (float): Initial guess for x
        y0 (float): Initial guess for y
        e (float): Desired accuracy

    Returns:
        tuple[float, list[float]]: Maximum of the function and the point
    """
    max_iter = 100
    h = 1e-6
    for _ in range(max_iter):
        dx = df_dx(f, x0, y0, h)
        dy = df_dy(f, x0, y0, h)
        if abs(dx) < e and abs(dy) < e:
            return f(x0, y0), [x0, y0]

        def f_single(t):
            return f(x0 + t * dx, y0 + t * dy)

        t = golden_section_search(f_single, 0, 1, e)
        x0 += t * dx
        y0 += t * dy
    return f(x0, y0), [x0, y0]
