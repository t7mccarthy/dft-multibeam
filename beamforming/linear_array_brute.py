import cmath
import math
import matplotlib.pyplot as plt
import time

theta = math.radians(68)
# wavelength for 5.8 GHz in cm
lmbda = 5.168835482759
k = 2 * math.pi / lmbda
d = lmbda / 2
n = 15


def get_array_factor (k, d, theta, n, alpha):
    sum = 0
    for i in range(-n, n):
        sum += cmath.exp(1j*i*(k*d*math.cos(theta)-alpha))
    return sum

def db (x):
    # return x
    return 10*cmath.log10(x/2)

def main():
    t = time.time()
    global theta, d, n, k, b, alpha
    # print (cmath.cos(theta))
    # print (2*math.cos(.5 * (k * d * math.cos(theta) + b)))
    # print (get_array_factor(k, d, theta, n, alpha))

    max_r = float("-inf")
    max_val = float("-inf")
    rs = [i for i in range(0,361)]
    vals = [0] * 361

    for r in range(0, 361):
        temp = get_array_factor(k, d, theta, n, r)
        # print (temp.real)
        # print (temp)
        if abs(temp) > max_val:
            max_val = abs(temp)
            max_r = r

        vals[r] = temp

        # if r < 46:
        #     x = [i for i in range(0,181)]
        #     y = [0]*181
        #     for i in x:
        #         y[i] = db(get_array_factor(k, d, math.radians(i), n, r))
        #
        #     plt.plot(x, y)
        #     plt.title('uniform progressive phase factor: ' + str(r))
        #     plt.show(block=False)
        #     plt.pause(.0001)
        #     plt.close()

    result = [i * max_r for i in range(0, 2 * n + 1)]
    elapsed = time.time() - t

    print (max_r)
    print (max_val)
    print (result)
    print ("time " + str(elapsed))

    #
    # x = [i for i in range(0,181)]
    # y = [0]*181
    # for i in x:
    #     y[i] = db(get_array_factor(k, d, math.radians(i), n, max_r))
    #
    # plt.plot(x, y)
    # plt.title('uniform progressive phase factor: ' + str(max_r))
    # plt.show()
    #
    #
    # plt.plot(rs, vals)
    # plt.title('radiation intensity vs alpha' + str(max_r))
    # plt.show()


if __name__ == '__main__':
    main()
