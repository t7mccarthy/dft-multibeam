import cmath
import math
import matplotlib.pyplot as plt
import time


theta = math.radians(50)
# wavelength for 5.8 GHz in cm
lmbda = 5.168835482759
k = 2 * math.pi / lmbda
d = lmbda / 2
n = 20

def get_array_factor (k, d, theta, n, alpha):
    sum = 0
    for i in range(-n, n):
        sum += cmath.exp(1j*i*(k*d*math.cos(theta)-alpha))
    return sum

def db (x):
    return 10*cmath.log10(x/2)

def main():
    t = time.time()
    global theta, d, n, k

    analytical = k * d * math.cos(theta)
    opt_val = abs(get_array_factor(k, d, theta, n, analytical))
    result = [i * analytical for i in range(0, 2 * n + 1)]
    elapsed = time.time() - t

    print analytical
    print opt_val
    print result
    print ("time " + str(elapsed))
    #
    # rounded = 1.07448
    # rounded_val = abs(get_array_factor(k, d, theta, n, rounded))
    # print rounded_val

    x = [i for i in range(0,181)]
    y = [0] * 181
    for i in x:
        y[i] = abs(get_array_factor(k, d, math.radians(i), n, analytical))

    plt.plot(x, y)
    plt.title('uniform progressive phase factor: ' + str(analytical))
    plt.show()


if __name__ == '__main__':
    main()
