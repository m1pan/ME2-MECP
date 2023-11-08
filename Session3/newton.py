import numpy
import matplotlib.pyplot as plt
import time

def divided_difference(xn, yn):
    if len(xn) != len(yn):
        raise Exception("xn and yn do not have the same length")
    
    if len(xn) == 1:
        return yn[0]
    else:
        return (divided_difference(xn[1:], yn[1:]) - divided_difference(xn[:-1], yn[:-1]))/(xn[-1] - xn[0])

def factorial(x):
    if x == 0:
        return 1
    else:
        return x*factorial(x-1)

def binomial(r, k):
    return factorial(r)/(factorial(k)*factorial(r - k))

def newton_interp(xn, yn, x):
    out = yn[0]
    for i in range(1, len(xn)):
        total = 1
        for j in range(i):
            total *= x - xn[j]
        out += total*divided_difference(xn[:i+1], yn[:i+1])
    return out

xn_linear = numpy.linspace(1, 2, 2)
yn_linear = numpy.sin(xn_linear)
xn_quad = numpy.linspace(1, 2, 3)
yn_quad = numpy.sin(xn_quad)
xn_cube = numpy.linspace(1, 2, 4)
yn_cube = numpy.sin(xn_cube)

x_vals = numpy.arange(0, 3, 0.05)
y_sin = numpy.sin(x_vals);

t1 = time.perf_counter()

pn_linear = [newton_interp(xn_linear, yn_linear, x) for x in x_vals]
pn_quad = [newton_interp(xn_quad, yn_quad, x) for x in x_vals]
pn_cube = [newton_interp(xn_cube, yn_cube, x) for x in x_vals]

t2 = time.perf_counter()

print(f"Completed in {t2 - t1} seconds")

# Plot sin(x) interpolations

plt.plot(x_vals, pn_linear, label="linear")
plt.plot(x_vals, pn_quad, label="quadratic")
plt.plot(x_vals, pn_cube, label="cubic")
plt.plot(x_vals, y_sin, label="true")

# Runge's phenomenon

#x_vals = numpy.arange(-1, 1, 0.01)
#y_vals = numpy.empty(len(x_vals))
#for i, x in enumerate(x_vals):
#    y_vals[i] = 1/(1+25*x**2)

#plt.plot(x_vals, y_vals)

#for order in range(1, 15):
#    xn = numpy.linspace(-1, 1, order + 1)
#    yn = [1/(1+25*x**2) for x in xn]
#    for i, x in enumerate(x_vals):
#        y_vals[i] = newton_interp(xn, yn, x)
#
#    if order == 3:
#        print(f"Answer: {newton_interp(xn, yn, 0.8)}")    
#
#    plt.plot(x_vals, y_vals, label=str(order))

plt.legend()
plt.show()
