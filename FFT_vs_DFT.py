import cmath
import math
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate
import time

N = 100000
M = N * 2 + 1
lmbda = 5.168835482759
d = lmbda / 2
wave_num = 2 * math.pi / lmbda
targets = [50, 65, 70, 100, 125]
T = len(targets)
theta_n = [0] * M
x_n = [0] * M
X_k = [0] * M
closest = [0] * T
show = True

b_bound = 4
t_bound = 32
dft_times = [0] * (t_bound - b_bound)
fft_times = [0] * (t_bound - b_bound)

# find set of angles possible to use in Fourier transform
def sample_angles():
    for n in range(M):
        theta_n[n] = math.degrees(math.acos((n-N) * lmbda / (M * d)))

# create target function where sampled angle nearest to each target is a peak
def peak_approximator():
    curr_t = T - 1
    i = 0
    last_d = float('inf')
    while curr_t >= 0:
        curr_d = targets[curr_t] - theta_n[i]
        if curr_d > 0:
            if curr_d < abs(last_d):
                x_n[i] = 2 * N / T
                closest[curr_t] = theta_n[i]
            else:
                x_n[i - 1] = 2 * N / T
                closest[curr_t] = theta_n[i-1]
            curr_t -= 1
            curr_d = targets[curr_t] - theta_n[i]
        last_d = curr_d
        i += 1

# run discrete Fourier transform to get phases/amplitudes for approximating target function
def dft():

    for k in range(M):
        sum = 0
        for n in range(M):
            sum += x_n[n] * cmath.exp(-2 * cmath.pi * 1j * n * k / M)
        X_k[k] = sum

    for k in range(M):
        X_k[k] = X_k[k] * cmath.exp(2 * cmath.pi * 1j * N * k / M)


def dft_fast():
    x = np.array(x_n)
    X = (np.fft.fft(x)).tolist()

    for k in range(M):
        X_k[k] = X[k] * cmath.exp(2 * cmath.pi * 1j * N * k / M)

# get array factor at given angle using known equation and generated phases/amplitudes
def get_array_factor (theta):
    sum = 0
    for n in range(M):
        sum += (X_k[n] / X_k[0]) * cmath.exp(1j * (n) * wave_num * d * math.cos(theta))
    return sum

# plot array factor function resulting from DFT
def visualize ():
    if show:
        x = [i * 2 + 1 for i in range(b_bound, t_bound)]
        plt.plot(x, dft_times, color = 'red', label = 'Naive')
        plt.plot(x, fft_times, color = 'blue', label = 'Fast')

        plt.legend(loc='best', fontsize='small')
        plt.title('Scalability of DFT')
        plt.xlabel('Number of Antennas')
        plt.ylabel('Runtime')
        plt.show()


def main():
    global N, M, theta_n, x_n, X_k, closest

    for n in range(b_bound, t_bound):
        # print n
        N = n
        M = N * 2 + 1
        theta_n = [0] * M
        x_n = [0] * M
        X_k = [0] * M
        closest = [0] * T
        sample_angles()
        peak_approximator()
        t = time.time()
        dft()
        dft_times[n-b_bound] = time.time() - t

        t = time.time()
        dft_fast()
        fft_times[n-b_bound] = time.time() - t


    visualize()

if __name__ == '__main__':
    main()
