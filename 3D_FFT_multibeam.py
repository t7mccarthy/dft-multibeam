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
m =4
N = n * 2 + 1
M = m * 2 + 1
lmbda = 5.168835482759
dx = lmbda / 2
dy = dx
wave_num = 2 * math.pi / lmbda
targets = [(30, 30), (45, 30,), (20, 40), (10, -30)]
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
    count = 0
    for x in range(M):
        for y in range(N):
            xm = x - m
            yn = y - n
            root = math.sqrt(b2 * (x ^ 2) + a2 * (y ^ 2))
            theta_xy[x][y] = math.degrees(math.asin(frac * root))
            if y != 0:
                phi_xy[x][y] = math.degrees(2 * (math.atan((1 / (a * y)) * (root - b * x))))
            angles[count] = (theta_xy[x][y], phi_xy[x][y])
            xy[count] = (x, y)
            count += 1

def peak_approximator():
    TARGETS = np.array(targets)
    tree = KDTree(np.array(angles), leaf_size = 2)
    # print(TARGETS[:2])
    dist, ind = tree.query(TARGETS, k = 1)
    print(ind)
    print(dist)
    for i in ind:
        x, y = xy[i[0]]
        f_xy[x][y] = m * n
    dists = [d[0] for d in dist]
    closest = [angles[i[0]] for i in ind]
    print(dists)
    print(closest)
    print(f_xy)
    # print (ind[0][0])
    i = ind[0][0]
    # print (angles[i])
    # print (targets[0])

def dft():
    for u in range(M):
        for v in range(N):
            sum = 0
            for x in range(M):
                for y in range(N):
                        sum += f_xy[x][y] * cmath.exp(-2 * cmath.pi * 1j * ((u * x / M) + (v * y / N)))
            F_uv[u][v] = sum




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
    # print (angles)
    peak_approximator()
    dft()

    print (F_uv)
    pts = []

    for a in range(M):
        plt.scatter(theta_xy[a], phi_xy[a], color = 'red', marker = 'o')
        for b in range(M):
            pts.append((theta_xy[a][b], phi_xy[a][b]))

    pts.append((30, 30))


    plt.scatter([angles[15][0]], [angles[15][1]], marker='o', color = 'blue')
    plt.scatter([30], [30], marker='s', color = 'yellow')
    # print (pts == angles)
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
        THETA = math.radians(pt[0])
        PHI = math.radians(pt[1])
        z = R * math.cos(THETA)
        # if z >= 0:
        xs.append(R * math.sin(THETA) * math.cos(PHI))
        ys.append(R * math.sin(THETA) * math.sin(PHI))
        zs.append((R * math.cos(THETA)))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.mplt3d([0, 0], [0,0], [-10, 10], linewidth=2, color='r')

    ax.scatter(xs, ys, zs, marker='o')
    ax.scatter([xs[15]], [ys[15]], [zs[15]], marker='s', color = 'yellow')
    ax.scatter([xs[-1]], [ys[-1]], [zs[-1]], marker='s', color = 'red')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()


if __name__ == '__main__':
    main()
