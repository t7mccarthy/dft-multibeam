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


n = 2
m = 2
N = n * 2 + 1
M = m * 2 + 1
lmbda = 5.168835482759
dx = lmbda * 0.5
dy = dx
wave_num = 2 * math.pi / lmbda
targets = [(30, 70), (40, -40), (10, 0), (60, 180)]
# targets = [(45, 0), (45, 90)]

targets = np.array([(180 - x, y) for (x, y) in targets])
T = len(targets)
theta_xy = [[0] * N for _ in range(M)]
phi_xy = [[0] * N for _ in range(M)]
angles = [0] * (M * N)
xy = [0] * (M * N)
f_xy = [[0] * N for _ in range(M)]
F_uv = [[0] * N for _ in range(M)]
closest = [0] * T
dists = [0] * T

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
            xm = x-m
            yn = y-n
            root = math.sqrt((b2 * (xm ** 2)) + (a2 * (yn ** 2)))
            try:
                theta_xy[x][y] = math.degrees(-math.asin(frac * root)+math.pi)
            except:
                pass
            if yn != 0:
                phi_xy[x][y] = math.degrees(2 * math.atan((1 / (a * yn)) * (root - (b * xm))))
            angles[count] = (theta_xy[x][y], phi_xy[x][y])
            xy[count] = (x, y)
            count += 1

def peak_approximator():
    global closest, dists
    TARGETS = np.array(targets)
    tree = KDTree(np.array(angles), leaf_size = 2)
    dist, ind = tree.query(TARGETS, k = 1)
    for i in ind:
        x, y = xy[i[0]]
        f_xy[x][y] = (m * n / T) * M * N
    dists = [d[0] for d in dist]
    closest = [angles[i[0]] for i in ind]
    i = ind[0][0]

def dft():
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
    k1 = cmath.exp(2 * cmath.pi * 1j * N / M)
    for u in range(M):
        for v in range(N):
            F_uv[u][v] = F_uv[u][v] * cmath.exp(2 * cmath.pi * 1j * ((u * m / M) + (v * n / N)))


def get_array_factor (theta, phi):
    sum = 0
    for u in range(M):
        for v in range(N):
            sum += (F_uv[u][v]/F_uv[0][0]) * cmath.exp(1j * wave_num * math.sin(theta) * (u * dx * math.cos(phi) + v * dy * math.sin(phi)))

    return abs(sum)

def visualize():
    # fig = plt.figure()
    # ax = fig.add_subplot(1,1,1, projection='3d')
    # xs = []
    # ys = []
    # zs = []
    # R = 150
    # for pt in angles:
    #     THETA = math.radians(pt[0])
    #     PHI = math.radians(pt[1])
    #     z = R * math.cos(THETA)
    #     # if z >= 0:
    #     xs.append(R * math.sin(THETA) * math.cos(PHI))
    #     ys.append(R * math.sin(THETA) * math.sin(PHI))
    #     zs.append(abs(R * math.cos(THETA)))
    # ax.scatter(xs, ys, zs, marker='o')
    # ax.set_xlabel('X Label')
    # ax.set_ylabel('Y Label')
    # ax.set_zlabel('Z Label')
    #
    # plt.show()
    #
    # pts = []
    # for a in range(M):
    #     plt.scatter(theta_xy[a], phi_xy[a], color = 'red', marker = 'o')
    #     for b in range(M):
    #         pts.append((theta_xy[a][b], phi_xy[a][b]))
    # plt.title('Angle Samples')
    # plt.xlabel('Theta')
    # plt.ylabel('Phi')
    # plt.show()

    # f2 = np.vectorize(get_array_factor)
    dim = 100
    theta, phi = np.linspace(0, np.pi, dim), np.linspace(0, 2 * np.pi, dim)
    THETA, PHI = np.meshgrid(theta, phi)
    r = [[0]*dim for _ in range(dim)]
    for i in range(dim):
        for j in range(dim):
            r[i][j]=get_array_factor(THETA[i][j], PHI[i][j])
    R = np.array(r)
    X = R * np.sin(THETA) * np.cos(PHI)
    Y = R * np.sin(THETA) * np.sin(PHI)
    Z = abs(R * np.cos(THETA))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    plot = ax.plot_surface(X, Y, Z)
    xs = []
    ys = []
    zs = []
    R = 0.75 * M * N / T
    for pt in angles:
        THETA = math.radians(pt[0])
        PHI = math.radians(pt[1])
        z = -R * math.cos(THETA)
        if z >= 0:
            xs.append(R * math.sin(THETA) * math.cos(PHI))
            ys.append(R * math.sin(THETA) * math.sin(PHI))
            zs.append(z)
    # ax.scatter(xs, ys, zs, marker='o')
    for i in range(T):
        x1 = [0]
        y1 = [0]
        z1 = [0]
        c = math.radians(targets[i][0])
        d = math.radians(targets[i][1])
        x1.append(2*R * math.sin(c) * math.cos(d))
        y1.append(2*R * math.sin(c) * math.sin(d))
        z1.append((-2*R * math.cos(c)))
        ax.plot(x1, y1, z1, marker='o')
        x2 = [0]
        y3 = [0]
        z4 = [0]
        c = math.radians(closest[i][0])
        d = math.radians(closest[i][1])
        x2.append(2*R * math.sin(c) * math.cos(d))
        y3.append(2*R * math.sin(c) * math.sin(d))
        z4.append((-2*R * math.cos(c)))
        ax.plot(x2, y3, z4, marker='s')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()



def main():
    print ('\nRUNNING FFT MULTIBEAM GENERATOR...\n')

    t = time.time()
    sample_angles()
    peak_approximator()
    dft_fast()
    runtime = time.time() - t

    # Calculate phase and amplitude of each antenna based on X_k signal
    amplitudes = [0] * M * N
    phases = [0] * M * N
    count = 0
    for j in range(N):
        for i in range(M):
            x = F_uv[i][j].real
            y = F_uv[i][j].imag
            phase = math.atan2(y, x)
            phases[count] = math.degrees(phase)
            amplitudes[count] = x / math.cos(phase)
            count += 1

    # output antenna configurations
    print ('\n' + '\033[1m' + 'Required Antenna Configurations:' + '\033[0m\n')
    array_config = zip(*[range(M * N), amplitudes, phases])
    print (tabulate(array_config, headers=['Antenna #', 'Amplitude', 'Phase'], tablefmt='orgtbl') + '\n')

    # get array factor at each target angle
    target_results = []
    for target in targets:
        target_results.append(abs(get_array_factor(math.radians(target[0]), math.radians(target[1]))))
    closest_peaks = []
    for angle in closest:
        closest_peaks.append(abs(get_array_factor(math.radians(angle[0]), math.radians(angle[1]))))
    # output results
    print ('\n\033[1m Resulting Beamform:\033[0m\n')
    print ('--- Best uniform peak array factor: ' + str(M * N / T))
    print ('--- Runtime: ' + str(runtime) + ' seconds\n')
    targets_out = [(180 - x[0], x[1]) for x in targets]
    closest_out = [(180 - x[0], x[1]) for x in closest]
    results = zip(*[[x + 1 for x in range(T)], targets_out, closest_out, target_results, closest_peaks])
    print (tabulate(results, headers=['Target #', 'Target Angle', 'Closest Sample', 'Array Factor', 'AF at Closest'], tablefmt='orgtbl') + '\n')
    print ('--- Mean difference target vs. closest: ' + str(sum(dists)/T))
    print ('--- Max difference target vs. closest: ' + str(max(dists)))
    print ('--- Mean array factor at target: ' + str(sum(target_results)/T))
    print ('--- Max array factor at target: ' + str(max(target_results)))
    print ('--- Min array factor at target: ' + str(min(target_results)) + '\n')
    print ('FINISHED\n')

    visualize()




if __name__ == '__main__':
    main()
