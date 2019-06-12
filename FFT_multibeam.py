import cmath
import math
import matplotlib.pyplot as plt
from tabulate import tabulate
import time

N = 1000
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

coeffs = [0] * M

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

def get_coeffs():
    k1 = cmath.exp(2 * cmath.pi * 1j * N / M)
    result = 1 + 0j
    for i in range(1, M):
        result *= k1
        coeffs[i] = result

# run discrete Fourier transform to get phases/amplitudes for approximating target function
def dft():
    # x_n2 = [0] * M
    # for n in range(M):
    #     x_n2[n] = x_n[n] * cmath.exp(2 * cmath.pi * 1j * N / M)
    coeffs[0] = 1 + 0j
    # get_coeffs()
    for i in range(M):
        coeffs[i] = cmath.exp(2 * cmath.pi * 1j * N * i / M)
    print coeffs

    for k in range(M):
        c = cmath.exp(2 * cmath.pi * 1j * N * k / M)
        sum = 0
        for n in range(M):
            sum += x_n[n] * cmath.exp(-2 * cmath.pi * 1j * (n) * k / M)
        X_k[k] = sum * c


# get array factor at given angle using known equation and generated phases/amplitudes
def get_array_factor (theta):
    sum = 0
    for n in range(M):
        sum += (X_k[n] / X_k[0]) * cmath.exp(1j * (n) * wave_num * d * math.cos(theta))
    return sum

# plot array factor function resulting from DFT
def visualize ():
    if show:
        plt.plot(theta_n, x_n, color = 'red', marker = 'o', label = 'Target')

        y = [0] * 181
        for i in range(181):
            y[i] = abs(get_array_factor(math.radians(i)))
        plt.plot(range(181), y, color = 'blue', label = 'Result')

        for target in targets:
            plt.axvline(x=target, color = 'red')

        plt.legend(loc='best', fontsize='small')
        plt.title('DFT-Generated Beamform')
        plt.xlabel('Angle from Array (Degrees)')
        plt.ylabel('Array Factor')
        plt.show()


def main():
    print ('\nRUNNING DFT MULTIBEAM GENERATOR...\n')

    t = time.time()
    sample_angles()
    peak_approximator()
    dft()

    runtime = time.time() - t

    # Calculate phase and amplitude of each antenna based on X_k signal
    amplitudes = [0] * M
    phases = [0] * M
    for i in range(M):
        x = X_k[i].real
        y = X_k[i].imag
        phase = math.atan2(y, x)
        phases[i] = math.degrees(phase)
        amplitudes[i] = x / math.cos(phase)

    # output antenna configurations
    print ('\n' + '\033[1m' + 'Required Antenna Configurations:' + '\033[0m\n')
    array_config = zip(*[[x + 1 for x in range(M)], amplitudes, phases])
    print (tabulate(array_config, headers=['Antenna #', 'Amplitude', 'Phase'], tablefmt='orgtbl') + '\n')

    # get array factor at each target angle
    target_results = []
    for target in targets:
        target_results.append(abs(get_array_factor(math.radians(target))))

    # output results
    print ('\n\033[1m Resulting Beamform:\033[0m\n')
    print ('--- Best uniform peak array factor: ' + str(2*N/T))
    print ('--- Runtime: ' + str(runtime) + '\n')
    results = zip(*[[x + 1 for x in range(T)], targets, closest, target_results])
    print (tabulate(results, headers=['Target #', 'Target Angle', 'Closest Sample', 'Array Factor'], tablefmt='orgtbl') + '\n')
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
