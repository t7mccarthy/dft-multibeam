import cmath
import math
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate
import pylab as p
from matplotlib import cm
import mpl_toolkits.mplot3d.axes3d as axes3d
import time


n = 4
m = 4
N = n * 2 + 1
M = m * 2 + 1
lmbda = 5.168835482759
dx = lmbda / 2
dy = dx
wave_num = 2 * math.pi / lmbda
targets = [(30, 30)]
T = len(targets)
theta_xy = [[0] * M for _ in range(M)]
phi_xy = [[0] * M for _ in range(M)]
x_n = [0] * M
X_k = [0] * M
closest = [0] * T
show = True

# find set of angles possible to use in Fourier transform
def sample_angles():
    a = M * dx / lmbda
    b = N * dy / lmbda
    a2 = a ** 2
    b2 = b ** 2
    frac = 1 / (a * b)
    for x in range(M):
        for y in range(N):
            xm = x - m
            yn = y - n
            root = math.sqrt(b2 * (x ^ 2) + a2 * (y ^ 2))
            theta_xy[x][y] = math.degrees(math.asin(frac * root))
            if y != 0:
                phi_xy[x][y] = math.degrees(2 * (math.atan((1 / (a * y)) * (root - b * x))))
                # if phi_xy


# def get_array_factor (k, dx, dy, theta, phi, n, ax, ay):
#     sum1 = 0
#     sum2 = 0
#     # print (theta)
#     a = math.sin(theta)
#     for i in range(-n, n):
#         sum1 += cmath.exp(1j*i*(k*dx*math.sin(theta)*math.cos(phi)-ax))
#         sum2 += cmath.exp(1j*i*(k*dy*math.sin(theta)*math.sin(phi)-ay))
#     return abs(sum1 * sum2)

def main():
    sample_angles()
    print (theta_xy)
    print (phi_xy)

    pts = []

    for a in range(M):
        plt.scatter(theta_xy[a], phi_xy[a], color = 'red', marker = 'o')
        for b in range(M):
            pts.append((theta_xy[a][b], phi_xy[a][b]))

    # plt.plot(range(M), theta_xy, color = 'red', marker = 'o', label = 'Theta samples')
    # plt.plot(range(M), phi_xy, color = 'blue', marker = 'o', label = 'Phi samples')

    # plt.legend(loc='best', fontsize='small')
    plt.title('Angle Samples')
    plt.xlabel('Theta')
    plt.ylabel('Phi')
    plt.show()

    xs = []
    ys = []
    zs = []
    R = 30
    for pt in pts:
        THETA = pt[0]
        PHI = pt[1]
        z = R * math.cos(THETA)
        # if z >= 0:
        xs.append(R * math.sin(THETA) * math.cos(PHI))
        ys.append(R * math.sin(THETA) * math.sin(PHI))
        zs.append(abs(R * math.cos(THETA)))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.mplt3d([0, 0], [0,0], [-10, 10], linewidth=2, color='r')

    ax.scatter(xs, ys, zs, marker='o')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()


if __name__ == '__main__':
    main()
