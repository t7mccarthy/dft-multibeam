import cmath
import math
import matplotlib.pyplot as plt
import numpy as np
import pylab as p
from matplotlib import cm
# from mpl_toolkits.mplot3d import Axes3D as p3
import mpl_toolkits.mplot3d.axes3d as axes3d
import time


theta = math.radians(20)
phi = math.radians(0)
# wavelength for 5.8 GHz in cm
lmbda = 5.168835482759
k = 2 * math.pi / lmbda
dx = lmbda / 3
dy = lmbda / 3
n = 10

def get_array_factor (k, dx, dy, theta, phi, n, ax, ay):
    sum1 = 0
    sum2 = 0
    # print (theta)
    a = math.sin(theta)
    for i in range(-n, n):
        sum1 += cmath.exp(1j*i*(k*dx*math.sin(theta)*math.cos(phi)-ax))
        sum2 += cmath.exp(1j*i*(k*dy*math.sin(theta)*math.sin(phi)-ay))
    return abs(sum1 * sum2)

def db (x):
    return 10*math.log10(x/6)

def fun(x, y, k, dx, dy, n, ax, ay):
    return (get_array_factor (k, dx, dy, x, y, n, ax, ay))

def main():
    global theta, phi, dx, dy, n, k

    alphax = k * dx * math.sin(theta) * math.cos(phi)
    alphay = k * dy * math.sin(theta) * math.sin(phi)
    print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print (alphay)
    alphay = abs(k * dy * cmath.sqrt((alphax/(k * dx))**2 - math.sin(theta)**2))
    print (alphay)
    print ("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    opt_val = get_array_factor(k, dx, dy, theta, phi, n, alphax, alphay)

    resultx = [i * alphax for i in range(2 * n + 1)]
    resulty = [i * alphax for i in range(2 * n + 1)]
    result = [[0] * (2 * n + 1) for _ in range(2 * n + 1)]
    for i in range(2 * n + 1):
        for j in range(2 * n + 1):
            result[i][j] = resultx[i] + resulty[j]


    print (alphax)
    print (alphay)
    print (opt_val)
    print (result)

    f2 = np.vectorize(fun)
    dim = 100
    theta, phi = np.linspace(0, 2 * np.pi, dim), np.linspace(0, np.pi, dim)
    THETA, PHI = np.meshgrid(theta, phi)
    r = [[0]*dim for _ in range(dim)]
    for i in range(dim):
        for j in range(dim):
            # print("*******************************************************")
            # print (math.degrees(THETA[i][j]), math.degrees(PHI[i][j]))
            r[i][j]=fun(THETA[i][j], PHI[i][j], k, dx, dy, n, alphax, alphay)
            # print (math.degrees(r[i][j]))
    R = np.array(r)
    print (np.shape(THETA))
    print (np.shape(PHI))
    # R = np.cos(PHI**2)
    # R = f2(THETA, PHI, k, dx, dy, n, alphax, alphay)
    print (np.shape(R))
    print (R)
    X = R * np.sin(THETA) * np.cos(PHI)
    Y = R * np.sin(THETA) * np.sin(PHI)
    Z = abs(R * np.cos(THETA))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    plot = ax.plot_surface(X, Y, Z)
        # X, Y, Z, rstride=1, cstride=1, cmap=plt.get_cmap('jet'),
        # linewidth=0, antialiased=False, alpha=0.5)

    plt.show()

    x = [i for i in range(0,181)]
    y = [0] * 181
    for i in x:
        y[i] = fun(math.radians(i), math.radians(30), k, dx, dy, n, alphax, alphay)

    plt.plot(x, y)
    plt.title('uniform progressive phase factor: ' + str(alphax))
    plt.show()



    x = [i for i in range(0,181)]
    y = [0] * 181
    for i in x:
        y[i] = fun(math.radians(30), math.radians(i), k, dx, dy, n, alphax, alphay)

    plt.plot(x, y)
    plt.title('uniform progressive phase factor: ' + str(alphay))
    plt.show()

    # fig = p.figure()
    # ax = p3(fig)
    # f2 = np.vectorize(fun)
    #
    # theta = np.arange(0,math.pi,math.pi/100)
    #
    # for a in theta:
    #     print (a)
    #
    # phi = np.arange(0,2*math.pi,math.pi/100)
    # # r = 2 * pow(math.e, -((theta**4)/(0.25**2)))
    # r =
    # for t in theta:
    #     for p in phi:
    #         r = db(get_array_factor (k, dx, dy, theta, phi, n, alphax, alphay))
    #
    # x = r*np.outer(np.cos(phi), np.sin(theta))
    # y = r*np.outer(np.sin(phi), np.sin(theta))
    # z = r*np.outer(np.ones(phi.shape), np.cos(theta))
    #
    # print (np.shape(x), np.shape(y), np.shape(z))
    #
    # ax.plot_surface(x,y,z)
    # ax.set_xlabel("X")
    # ax.set_ylabel("Y")
    # ax.set_zlabel("Z")
    #
    # p.show()

    # print (alphax)
    # print (alphay)
    # print (opt_val)
    # print (result)
    # print ("time " + str(elapsed))
    #
    # max_rx = 0
    # max_ry = 0
    # max_val = float("-inf")
    # axs = [i for i in range(0,361)]
    # ays = [i for i in range(0,361)]
    # for u in range(0, 361):
    #     for v in range(0, 361):
    #         temp = get_array_factor(k, d, theta, n, r)
    #         if abs(temp) > max_val:
    #             max_val = abs(temp)
    #             max_r = r
    # result = [i * max_r for i in range(0, 2 * n + 1)]
    # print (max_r)
    # print (max_val)
    # print (result)

    # f2 = np.vectorize(fun)
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # x = y = np.arange(0, 360, 1)
    # R, P = np.meshgrid(x, y)
    # X, Y = R*np.cos(P), R*np.sin(P)
    # zs = np.array(f2(np.ravel(X), np.ravel(Y),k, dx, dy, n, alphax, alphay))
    # Z = zs.reshape(X.shape)

    # ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
    #                        linewidth=0, antialiased=False)
    # ax.set_xlabel('X axis')
    # ax.set_ylabel('Y axis')
    # ax.set_zlabel('Z axis')
    # plt.show()


    # r = np.linspace(0, 1.25, 50)
    # p = np.linspace(0, 2*np.pi, 50)
    # R, P = np.meshgrid(r, p)
    # Z = ((R**2 - 1)**2)
    #
    # # Express the mesh in the cartesian system.
    # X, Y = R*np.cos(P), R*np.sin(P)
    #
    # # Plot the surface.
    # ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r)
    #
    # # Tweak the limits and add latex math labels.
    # ax.set_zlim(0, 1)
    # ax.set_xlabel(r'$\phi_\mathrm{real}$')
    # ax.set_ylabel(r'$\phi_\mathrm{im}$')
    # ax.set_zlabel(r'$V(\phi)$')
    #
    # plt.show()


if __name__ == '__main__':
    main()
