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
dx = lmbda * 1.2
dy = dx
wave_num = 2 * math.pi / lmbda
targets = [(175, 45)]
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
    print (a, b, frac)
    for x in range(M):
        for y in range(N):
            xm = x-m
            yn = y-n
            root = math.sqrt((b2 * (xm ** 2)) + (a2 * (yn ** 2)))
            theta_xy[x][y] = math.degrees(-math.asin(frac * root)+math.pi)
            if yn != 0:
                phi_xy[x][y] = math.degrees(2 * math.atan((1 / (a * yn)) * (root - (b * xm))))
            angles[count] = (theta_xy[x][y], phi_xy[x][y])
            xy[count] = (x, y)
            count += 1

def peak_approximator():
    global closest
    # t = [(math.radians(a), math.radians(b)) for (a, b) in targets]
    # print (targets)
    # print (t)
    # print (angles)
    TARGETS = np.array(targets)
    tree = KDTree(np.array(angles), leaf_size = 2)
    # print(TARGETS[:2])
    dist, ind = tree.query(TARGETS, k = 1)
    print(ind)
    print(dist)
    for i in ind:
        x, y = xy[i[0]]
        print (x, y)
        print (theta_xy[x][y], phi_xy[x][y])
        f_xy[x][y] = (m * n / T) * M * N
    dists = [d[0] for d in dist]
    closest = [angles[i[0]] for i in ind]
    # print(dists)
    # print(closest)
    # print(f_xy)
    # # print (ind[0][0])
    i = ind[0][0]
    # print (angles[i])
    # print (targets[0])

def dft():
    print (f_xy)
    for u in range(M):
        for v in range(N):
            sum = 0
            for x in range(-m, m):
                for y in range(-n, n):
                        sum += f_xy[x+m][y+n] * cmath.exp(-2 * cmath.pi * 1j * ((u * x / M) + (v * y / N)))
            F_uv[u][v] = sum


def dft_fast():
    global F_uv
    F_uv = (np.fft.fft2(np.array(f_xy))).tolist()
    # k1 = cmath.exp(2 * cmath.pi * 1j * N / M)
    # for k in range(M):
    #     X_k[k] = X[k] * (k1 ** k)


def get_array_factor (theta, phi):
    sum = 0
    for u in range(M):
        for v in range(N):
            sum += (F_uv[u][v]/F_uv[0][0]) * cmath.exp(1j * wave_num * math.sin(theta) * (u * dx * math.cos(phi) + v * dy * math.sin(phi)))

    return abs(sum)
    # for i in range(-n, n):
    #     sum1 += cmath.exp(1j*i*(k*dx*math.sin(theta)*math.cos(phi)-ax))
    #     sum2 += cmath.exp(1j*i*(k*dy*math.sin(theta)*math.sin(phi)-ay))
    # return abs(sum1 * sum2)

def main():
    sample_angles()
    peak_approximator()
    dft()
    # print (F_uv)


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')

    xs = []
    ys = []
    zs = []
    R = 60
    for pt in angles:
        THETA = math.radians(pt[0])
        PHI = math.radians(pt[1])
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

    pts = []

    for a in range(M):
        plt.scatter(theta_xy[a], phi_xy[a], color = 'red', marker = 'o')
        for b in range(M):
            pts.append((theta_xy[a][b], phi_xy[a][b]))

    plt.title('Angle Samples')
    plt.xlabel('Theta')
    plt.ylabel('Phi')
    plt.show()

    # f2 = np.vectorize(get_array_factor)
    dim = 100
    theta, phi = np.linspace(0, np.pi, dim), np.linspace(0, 2 * np.pi, dim)
    THETA, PHI = np.meshgrid(theta, phi)
    r = [[0]*dim for _ in range(dim)]
    for i in range(dim):
        for j in range(dim):
            # print("*******************************************************")
            # print (math.degrees(THETA[i][j]), math.degrees(PHI[i][j]))
            r[i][j]=get_array_factor(THETA[i][j], PHI[i][j])
            # print (math.degrees(r[i][j]))
    R = np.array(r)
    # print (np.shape(THETA))
    # print (np.shape(PHI))
    # R = np.cos(PHI**2)
    # R = f2(THETA, PHI, k, dx, dy, n, alphax, alphay)
    # print (np.shape(R))
    # print (R)
    X = R * np.sin(THETA) * np.cos(PHI)
    Y = R * np.sin(THETA) * np.sin(PHI)
    Z = -abs(R * np.cos(THETA))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    plot = ax.plot_surface(X, Y, Z)
        # X, Y, Z, rstride=1, cstride=1, cmap=plt.get_cmap('jet'),
        # linewidth=0, antialiased=False, alpha=0.5)

    # plt.show()

    # print (F_uv)
    # pts = []
    #
    # for a in range(M):
    #     plt.scatter(theta_xy[a], phi_xy[a], color = 'red', marker = 'o')
    #     for b in range(M):
    #         pts.append((theta_xy[a][b], phi_xy[a][b]))
    #
    # pts.append((30, 30))
    #
    #
    # plt.scatter([angles[15][0]], [angles[15][1]], marker='o', color = 'blue')
    # plt.scatter([30], [30], marker='s', color = 'yellow')
    # # print (pts == angles)
    # # plt.plot(range(M), theta_xy, color = 'red', marker = 'o', label = 'Theta samples')
    # # plt.plot(range(M), phi_xy, color = 'blue', marker = 'o', label = 'Phi samples')
    #
    # # plt.legend(loc='best', fontsize='small')
    # plt.title('Angle Samples')
    # plt.xlabel('Theta')
    # plt.ylabel('Phi')
    # plt.show()
    #
    xs = []
    ys = []
    zs = []
    R = 60
    for pt in angles:
        THETA = math.radians(pt[0])
        PHI = math.radians(pt[1])
        z = R * math.cos(THETA)
        # if z >= 0:
        xs.append(R * math.sin(THETA) * math.cos(PHI))
        ys.append(R * math.sin(THETA) * math.sin(PHI))
        zs.append((R * math.cos(THETA)))

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.mplt3d([0, 0], [0,0], [-10, 10], linewidth=2, color='r')

    ax.scatter(xs, ys, zs, marker='o')
    # ax.scatter([xs[15]], [ys[15]], [zs[15]], marker='s', color = 'yellow')
    # ax.scatter([xs[-1]], [ys[-1]], [zs[-1]], marker='s', color = 'red')
    x1 = [0]
    y1 = [0]
    z1 = [0]
    c = math.radians(targets[0][0])
    d = math.radians(targets[0][1])
    x1.append(2*R * math.sin(c) * math.cos(d))
    y1.append(2*R * math.sin(c) * math.sin(d))
    z1.append((2*R * math.cos(c)))
    ax.plot(x1, y1, z1, marker='o')


    x2 = [0]
    y3 = [0]
    z4 = [0]
    print(closest)
    c = math.radians(closest[0][0])
    d = math.radians(closest[0][1])
    x2.append(2*R * math.sin(c) * math.cos(d))
    y3.append(2*R * math.sin(c) * math.sin(d))
    z4.append((2*R * math.cos(c)))
    ax.plot(x2, y3, z4, marker='s')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()


if __name__ == '__main__':
    main()
