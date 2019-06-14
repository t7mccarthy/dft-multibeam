import cmath
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import time

theta = math.radians(156)
# wavelength for 5.8 GHz in cm
lmbda = 5.168835482759
k = 2 * math.pi / lmbda
d = lmbda / 2
n = 15

def get_array_factor (k, d, theta, n, alpha):
    sum = 0
    for i in range(-n, n):
        sum += cmath.exp(1j*i*(k*d*math.cos(theta)-alpha))
    return abs(sum)

def db (x):
    return 10*cmath.log10(x/2)

def main():
    t = time.time()
    global theta, d, n, k, b, alpha, plot, function

    function = lambda x: -1 * get_array_factor(k, d, theta, n, x)
    f2 = np.vectorize(function)

    res = minimize_scalar(f2)
    alpha = res.x
    opt_val = abs(get_array_factor(k, d, theta, n, alpha))
    vector = [i * alpha for i in range(0, 2 * n + 1)]
    elapsed = time.time() - t

    print (alpha)
    print (opt_val)
    print (vector)
    print ("time " + str(elapsed))


if __name__ == '__main__':
    main()
