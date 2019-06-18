import cmath
import math
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate
import pylab as p
from matplotlib import cm
import mpl_toolkits.mplot3d.axes3d as axes3d
import time
from sklearn.neighbors import KDTree


n = 4
m = 4
N = n * 2 + 1
M = m * 2 + 1
lmbda = 5.168835482759
dx = lmbda * 1.5
dy = dx
wave_num = 2 * math.pi / lmbda
targets = [(20, -45)]
T = len(targets)
theta_xy = [[0] * N for _ in range(M)]
phi_xy = [[0] * N for _ in range(M)]
angles = [0] * (M * N)
xy = [0] * (M * N)
f_xy = [[0] * N for _ in range(M)]
F_uv = [[0] * N for _ in range(M)]
closest = [0] * T
dists = [0] * T
show = True
ind = 0

# find set of angles possible to use in Fourier transform
def sample_angles():
    a = M * dx / lmbda
    b = N * dy / lmbda
    a2 = a ** 2
    b2 = b ** 2
    frac = 1 / (a * b)
    print (a, b, frac)
    count = 0
    for x in range(M):
        for y in range(N):
            xm = x-m
            yn = y-n
            root = math.sqrt((b2 * (xm ** 2)) + (a2 * (yn ** 2)))
            print (xm, yn)
            print (frac * root)
            theta_xy[x][y] = -math.asin(frac * root)+math.pi
            if yn != 0:
                phi_xy[x][y] = 2 * math.atan((1 / (a * yn)) * (root - (b * xm)))
            angles[count] = (theta_xy[x][y], phi_xy[x][y])
            print(xm, yn, theta_xy[x][y], phi_xy[x][y])
            count += 1


def main():
    sample_angles()
    t = theta_xy[7][8]
    p = phi_xy[7][8]
    # print(math.cos(0))
    print (M*dx*math.sin(t)*math.cos(p)/lmbda)
    print (N*dy*math.sin(t)*math.sin(p)/lmbda)
    print (theta_xy[7][8], phi_xy[7][8])


    pts = []
    for a in range(M):
        plt.scatter(theta_xy[a], phi_xy[a], color = 'red', marker = 'o')
        for b in range(M):
            pts.append((theta_xy[a][b], phi_xy[a][b]))
    plt.title('Angle Samples')
    plt.xlabel('Theta')
    plt.ylabel('Phi')
    plt.show()


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')

    xs = []
    ys = []
    zs = []
    R = 60
    for pt in angles:
        THETA = pt[0]
        PHI = pt[1]
        z = R * math.cos(THETA)
        # if z >= 0:
        xs.append(R * math.sin(THETA) * math.cos(PHI))
        ys.append(R * math.sin(THETA) * math.sin(PHI))
        zs.append(abs(R * math.cos(THETA)))

    ax.scatter(xs, ys, zs, marker='o')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()

if __name__ == '__main__':
    main()
