import cmath
import math
import matplotlib.pyplot as plt
import time

N = 20
M = N * 2 + 1
lmbda = 5.168835482759
d = lmbda / 2
wave_num = 2 * math.pi / lmbda
targets = [50, 65, 70, 100, 125]
T = len(targets)
theta_n = [0] * M
x_n = [0] * M
X_k = [0] * M
show = True


def sample_angles():
    for n in range(-N, N):
        theta_n[n + N] = math.degrees(math.acos(n * lmbda / (M * d)))


def peak_approximator():
    curr_t = T - 1
    i = 0
    last_d = float("inf")
    while curr_t >= 0:
        curr_d = targets[curr_t] - theta_n[i]
        if curr_d > 0:
            if curr_d < abs(last_d):
                x_n[i] = 2 * N / T
            else:
                x_n[i - 1] = 2 * N / T
            curr_t -= 1
            curr_d = targets[curr_t] - theta_n[i]
        last_d = curr_d
        i += 1


def dft():
    for k in range(M):
        sum = 0
        for n in range(M):
            sum += x_n[n] * cmath.exp(-2 * cmath.pi * 1j * (n - N) * k / M)
        X_k[k] = sum


def get_array_factor (theta):
    sum = 0
    for n in range(-N, N):
        sum += (X_k[n + N] / X_k[0])*cmath.exp(1j * n * wave_num * d * math.cos(theta))
    return sum


def visualize ():
    if show:
        plt.plot(theta_n, x_n, color = 'red', marker = "o")

        x1 = [i for i in range(181)]
        y1 = [0] * 181
        for i in x1:
            y1[i] = abs(get_array_factor(math.radians(i)))
        plt.plot(x1, y1, color = 'blue')

        for target in targets:
            plt.axvline(x=target, color = "red")
        plt.show()


def main():
    t = time.time()
    sample_angles()
    peak_approximator()
    dft()

    print (time.time() - t)

    results = []
    for target in targets:
        results.append(abs(get_array_factor(math.radians(target))))
    print(results)

    visualize()

if __name__ == '__main__':
    main()
