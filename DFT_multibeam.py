import cmath
import math
import matplotlib.pyplot as plt
from tabulate import tabulate
import time

N = 16
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

# find set of angles possible to use in Fourier transform
def sample_angles():
    for n in range(-N, N):
        theta_n[n + N] = math.degrees(math.acos(n * lmbda / (M * d)))

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
            sum += x_n[n] * cmath.exp(-2 * cmath.pi * 1j * (n - N) * k / M)
        X_k[k] = sum

# get array factor at given angle using known equation and generated phases/amplitudes
def get_array_factor (theta):
    sum = 0
    for n in range(-N, N):
        sum += (X_k[n + N] / X_k[0]) * cmath.exp(1j * n * wave_num * d * math.cos(theta))
    return sum

# plot array factor function resulting from DFT
def visualize ():
    if show:
        plt.plot(theta_n, x_n, color = 'red', marker = 'o', label = 'Target')

        x1 = [i for i in range(181)]
        y1 = [0] * 181
        for i in x1:
            y1[i] = abs(get_array_factor(math.radians(i)))
        plt.plot(x1, y1, color = 'blue', label = 'Result')

        for target in targets:
            plt.axvline(x=target, color = 'red')

        plt.legend(loc='best', fontsize='small')
        plt.title('DFT-Generated Beamform')
        plt.xlabel('Angle from Array (Degrees)')
        plt.ylabel('Array Factor')
        plt.show()


def main():
    print ('\nRUNNING DFT MULTIBEAM GENERATOR...')

    t = time.time()
    sample_angles()
    peak_approximator()
    dft()

    runtime = time.time() - t

    # get array factor at each target angle
    target_results = []
    for target in targets:
        target_results.append(abs(get_array_factor(math.radians(target))))

    # output results
    print ('\n' + '\033[1m' + 'Results with ' + str(M) + ' antennas and ' + str(T) + ' targets:' + '\033[0m')
    print ('--- Best uniform peak array factor: ' + str(2*N/T))
    print ('--- Runtime: ' + str(runtime) + '\n')
    results = zip(*[targets, closest, target_results])
    print (tabulate(results, headers=['Target', 'Closest Sample', 'Array Factor'], tablefmt='orgtbl') + '\n')
    diffs = [abs(targets[i] - closest[i]) for i in range(T)]
    print ('--- Mean difference target vs. closest: ' + str(sum(diffs)/T))
    print ('--- Max difference target vs. closest: ' + str(max(diffs)))
    print ('--- Mean array factor at target: ' + str(sum(target_results)/T))
    print ('--- Max array factor at target: ' + str(max(target_results)))
    print ('--- Min array factor at target: ' + str(min(target_results)) + '\n')
    print ('FINISHED\n')

    visualize()

if __name__ == '__main__':
    main()
