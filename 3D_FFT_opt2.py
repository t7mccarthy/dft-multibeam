import cmath as c
from math import pi, sqrt, degrees, radians, sin, cos, asin, atan, atan2
import matplotlib.pyplot as plt
from matplotlib import cm
import mpl_toolkits.mplot3d.axes3d as axes3d
import numpy as np
import pylab as p
from scipy.optimize import minimize
from sklearn.neighbors import KDTree
from tabulate import tabulate
import time


# m = antennas in x-direction, n = antennas in y direction
# m = 2
# n = 2
# M = m * 2 + 1
# N = n * 2 + 1
M = 4
N = 4
m = (M - 1.) / 2.
n = (N - 1.) / 2.
lmbda = 5.168835482759
dx = lmbda * 0.5
dy = dx
wave_num = 2 * pi / lmbda
targets = [(30, 70), (40, -40), (10, 0), (60, 180)]
# targets = [(45, 0), (45, 90)]

targets = np.array([(180 - x, y) for (x, y) in targets])
T = len(targets)
theta_xy = np.zeros((M, N))
phi_xy = np.zeros((M, N))
angles = np.zeros(M * N, dtype=(float,2))
xy = np.zeros(M * N, dtype=(int,2))
f_xy = np.zeros((M, N))
F_uv = np.zeros((M, N))
phases = np.zeros(M * N)
targ_coeffs = [0] * T


# Find set of angles possible to use in Fourier transform
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
            root = sqrt((b2 * (xm ** 2)) + (a2 * (yn ** 2)))
            try:
                theta_xy[x][y] = degrees(-asin(frac * root)+pi)
            except:
                pass
            if yn != 0:
                phi_xy[x][y] = degrees(2 * atan((1 / (a * yn)) * (root - (b * xm))))
            angles[count] = (theta_xy[x][y], phi_xy[x][y])
            xy[count] = (x, y)
            count += 1

# Find closest sample angles to targets (to create target function for FFT)
def peak_approximator():
    tree = KDTree(angles, leaf_size = 2)
    dist, ind = tree.query(targets, k = 1)
    for i in ind:
        x, y = xy[i][0]
        f_xy[x][y] = (m * n / T) * M * N

# Run FFT to get phases/amplitudes for approximating target function
def dft_fast():
    global F_uv
    F_uv = np.fft.fft2(f_xy)
    ind = 0
    for u in range(M):
        for v in range(N):
            F_uv[u][v] = F_uv[u][v] * c.exp(2 * c.pi * 1j * ((u * m / M) + (v * n / N)))
            phases[ind] = atan2(F_uv[u][v].imag, F_uv[u][v].real)
            F_uv[u][v] = c.exp(1j * phases[ind])
            ind += 1

# Calculate AF at targets (won't recalculate coefficient vector); used during optimization
def get_array_factor_opt (theta, phi, i):
    try: targ_coeffs[i].any()
    except:
        coeffs = np.zeros((M, N), dtype = 'complex')
        def get_coeff(v, u):
            return c.exp(1j * wave_num * sin(theta) * (u * dx * cos(phi) + v * dy * sin(phi)))
        vcoeffs = np.vectorize(get_coeff)
        for u in range(N):
            coeffs[u] = vcoeffs(np.arange(N), u)
        targ_coeffs[i] = np.copy(coeffs)
    return abs(np.sum(np.multiply(F_uv / F_uv[0][0], targ_coeffs[i])))

# Maximize total AF at targets by running scipy minimization alg (BFGS)
def opt ():
    global phases, F_uv
    guess = np.copy(phases)
    def get_AFS(S):
        global phases, F_uv
        sum = 0
        phases = np.copy(S)
        F_uv = np.exp(1j * np.reshape(phases, (M, N)))
        for t in range(T):
            sum += abs(get_array_factor_opt(radians(targets[t][0]), radians(targets[t][1]), t))
        return -sum
    get_AFS(minimize(get_AFS, guess).x)

# Calculate AF at any angle (used during visualization and analysis)
def get_array_factor (theta, phi):
    sum = 0
    for u in range(M):
        for v in range(N):
            sum += (F_uv[u][v]/F_uv[0][0]) * c.exp(1j * wave_num * sin(theta) * (u * dx * cos(phi) + v * dy * sin(phi)))
    return abs(sum)

# Create 3D visualization of array factor
def visualize():
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
    R =  1.5 * M * N / T
    for i in range(T):
        x, y, z = [0], [0], [0]
        theta = radians(targets[i][0])
        phi = radians(targets[i][1])
        x.append(2 * R * sin(theta) * cos(phi))
        y.append(2 * R * sin(theta) * sin(phi))
        z.append((-2 * R * cos(theta)))
        ax.plot(x, y, z, marker='o')
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
    opt()
    runtime = time.time() - t

    # Calculate phase and amplitude of each antenna based on F_uv signal
    amplitudes = [0] * M * N
    phases = [0] * M * N
    count = 0
    for j in range(N):
        for i in range(M):
            x = F_uv[i][j].real
            y = F_uv[i][j].imag
            phase = atan2(y, x)
            phases[count] = degrees(phase)
            amplitudes[count] = x / cos(phase)
            count += 1

    # output antenna configurations
    print ('\n' + '\033[1m' + 'Required Antenna Configurations:' + '\033[0m\n')
    array_config = zip(*[range(M * N), amplitudes, phases])
    print (tabulate(array_config, headers=['Antenna #', 'Amplitude', 'Phase'], tablefmt='orgtbl') + '\n')

    # get array factor at each target angle
    target_results = []
    for target in targets:
        target_results.append(abs(get_array_factor(radians(target[0]), radians(target[1]))))

    # output results
    print ('\n\033[1m Resulting Beamform:\033[0m\n')
    print ('--- Uniform peak array factor (antennas/targets): ' + str(M * N / T))
    print ('--- Runtime: ' + str(runtime) + ' seconds\n')
    targets_out = [(180 - x[0], x[1]) for x in targets]
    results = zip(*[[x + 1 for x in range(T)], targets_out, target_results])
    print (tabulate(results, headers=['Target #', 'Target Angle', 'Array Factor'], tablefmt='orgtbl') + '\n')
    print ('--- Mean array factor at target: ' + str(sum(target_results)/T))
    print ('--- Max array factor at target: ' + str(max(target_results)))
    print ('--- Min array factor at target: ' + str(min(target_results)) + '\n')
    print ('FINISHED\n')

    visualize()


if __name__ == '__main__':
    main()
