import cmath
import math
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate
import time
from scipy.optimize import minimize

N = 10
M = N * 2 + 1
lmbda = 5.168835482759
d = lmbda / 2
wave_num = 2 * math.pi / lmbda
# targets = [50, 100, 125, 150, 155]
targets = [50, 100, 145]
T = len(targets)
theta_n = [0] * M
x_n = [0] * M
X_k = [0] * M
X = [0] * M
closest = [0] * T
show = True
phases = [0] * M
phasesX = [0] * M
LEN = 10001
padded_samples = [0] * LEN

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
            ind = i - (curr_d > abs(last_d))
            x_n[ind] = M / T
            closest[curr_t] = theta_n[ind]
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


# run fast Fourier transform to get phases/amplitudes for approximating target function
def dft_fast():
    global X
    X = (np.fft.fft(np.array(x_n))).tolist()
    k1 = cmath.exp(2 * cmath.pi * 1j * N / M)
    for k in range(M):
        X_k[k] = X[k] * (k1 ** k)
        phases[k] = math.atan2(X_k[k].imag, X_k[k].real)
        X_k[k] = cmath.exp(1j * phases[k])
    # for k in range(M):
    #     phasesX[k] = math.atan2(X[k].imag, X[k].real)
    #     X[k] = cmath.exp(1j * phasesX[k])

# get array factor at given angle using known equation and generated phases/amplitudes
def get_array_factor (theta):
    sum = 0
    for n in range(M):
        sum += (X_k[n] / X_k[0]) * cmath.exp(1j * (n) * wave_num * d * math.cos(theta))
    return sum

# Used for padding IFFT samples
def padded_sample_angles():
    for n in range(LEN):
        padded_samples[n] = math.degrees(math.acos((n-(LEN/2)) * lmbda / (LEN * d)))

# Works but not useful
def get_af_ifft ():
    global X
    padded_sample_angles()
    # x = M * d * math.cos(theta)
    # print (theta_n[i], x_n[i])
    X = X + ([0] * (LEN - M))
    # print (X)
    f = (np.fft.ifft(np.array(X))).tolist()
    # print (x)
    # print (f)
    # print (M * f[i])
    return f


def opt ():
    global phases, X_k
    guess = np.array(phases[:])
    def get_AFS(S):
        # r = time.time()
        sum = 0
        phases = (S.tolist())[:]
        for e in range(M):
            X_k[e] = cmath.exp(1j * phases[e])
        for t in targets:
            sum += abs(get_array_factor(math.radians(t)))
        # print (time.time() - r)
        return -sum
    get_AFS(minimize(get_AFS, guess).x)

# Significantly slower
def opt_alt ():
    global phases, X_k, X

    sample_targ_inds = [0] * T

    curr_t = T - 1
    i = 0
    last_d = float('inf')
    while curr_t >= 0:
        curr_d = targets[curr_t] - padded_samples[i]
        if curr_d > 0:
            sample_targ_inds[curr_t] = i - (curr_d > abs(last_d))
            curr_t -= 1
            curr_d = targets[curr_t] - padded_samples[i]
        last_d = curr_d
        i += 1

    guess = np.array(phases[:])
    k1 = cmath.exp(2 * cmath.pi * 1j * (LEN / 2) / LEN)
    def get_AFS(S):
        # r = time.time()
        sum = 0
        phases = (S.tolist())[:]
        for e in range(M):
            X_k[e] = cmath.exp(1j * phases[e])
            X[e] = X_k[e] / (k1 ** e)
        # print (time.time() - r)
        # r = time.time()
        f = get_af_ifft()
        for t in range(T):
            sum += abs(f[sample_targ_inds[t]])
        # print (time.time() - r)
        return -sum

    get_AFS(minimize(get_AFS, guess).x)


# plot array factor function resulting from DFT
def visualize (results):
    if show:
        plt.plot(theta_n, x_n, color = 'green', marker = 'o', label = 'Target Function')
        y = [0] * 181
        for i in range(181):
            y[i] = abs(get_array_factor(math.radians(i)))
        plt.plot(range(181), y, color = 'blue', label = 'Result', linestyle = 'dotted')
        plt.scatter(targets, results, color = 'orange', marker = '*', label = 'Target Angles', s = 100)
        plt.legend(loc='best', fontsize='small')
        plt.title('FFT-Generated Beamform')
        plt.xlabel('Angle from Array (Degrees)')
        plt.ylabel('Array Factor')
        plt.show()


def main():
    global phases
    print ('\nRUNNING FFT MULTIBEAM GENERATOR...\n')

    t = time.time()
    sample_angles()
    peak_approximator()
    dft_fast()
    opt()

    runtime = time.time() - t
    amplitudes = [1] * M

    # output antenna configurations
    print ('\n' + '\033[1m' + 'Required Antenna Configurations:' + '\033[0m\n')
    array_config = zip(*[[x + 1 for x in range(M)], amplitudes, [math.degrees(x) for x in phases]])
    print (tabulate(array_config, headers=['Antenna #', 'Amplitude', 'Phase'], tablefmt='orgtbl') + '\n')

    # get array factor at each target angle
    target_results = []
    for target in targets:
        target_results.append(abs(get_array_factor(math.radians(target))))
    closest_peaks = []
    for angle in closest:
        closest_peaks.append(abs(get_array_factor(math.radians(angle))))

    # output results
    print ('\n\033[1m Resulting Beamform:\033[0m\n')
    print ('--- Best uniform peak array factor: ' + str(M / T))
    print ('--- Runtime: ' + str(runtime) + ' seconds\n')
    results = zip(*[[x + 1 for x in range(T)], targets, closest, target_results, closest_peaks])
    print (tabulate(results, headers=['Target #', 'Target Angle', 'Closest Sample', 'Array Factor', 'AF at Closest'], tablefmt='orgtbl') + '\n')
    diffs = [abs(targets[i] - closest[i]) for i in range(T)]
    # print ('--- Mean difference target vs. closest: ' + str(sum(diffs)/T))
    # print ('--- Max difference target vs. closest: ' + str(max(diffs)))
    print ('--- Mean array factor at target: ' + str(sum(target_results)/T))
    print ('--- Max array factor at target: ' + str(max(target_results)))
    print ('--- Min array factor at target: ' + str(min(target_results)) + '\n')
    print ('FINISHED\n')

    visualize(target_results)

    # target_results = []
    # for target in targets:
    #     target_results.append(abs(get_af_ifft(math.radians(target))))

    # y = [0] * LEN
    # f = get_af_ifft()
    # for i in range(LEN):
    #     y[i] = abs(LEN * f[i])
    #     # print (i, padded_samples[i], y[i])
    # plt.plot(padded_samples, y, color = 'orange', label = 'Result', linestyle = 'dashed')
    # plt.legend(loc='best', fontsize='small')
    # plt.title('FFT-Generated Beamform')
    # plt.xlabel('Angle from Array (Degrees)')
    # plt.ylabel('Array Factor')
    # plt.show()

    # get_af_ifft(math.radians(50))
    #
    # print (get_array_factor(math.radians(50)))

if __name__ == '__main__':
    main()
