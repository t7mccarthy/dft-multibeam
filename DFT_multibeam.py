import cmath
import math
import matplotlib.pyplot as plt
import time

lmbda = 5.168835482759
d = lmbda / 2
M = 20
N = M * 2 + 1
k = 2 * math.pi / lmbda
targets = [50, 65, 70, 100, 125]
T = len(targets)
theta_n = [0] * N
x_n = [0] * N

def sample_angles():
    global lmbda, N, d , theta_n, M
    for n in range((-1 * M), M):
        theta_n[n+M] = math.degrees(math.acos((n*lmbda)/(N*d)))

def peak_approximator():
    global targets, T, theta_n, M, N, x_n
    curr_t = T-1
    i = 0
    last_d = float("inf")
    print(len(theta_n))
    print(N)
    while curr_t >= 0:
        curr_d = targets[curr_t] - theta_n[i]
        if curr_d > 0:
            if curr_d < abs(last_d):
                x_n[i] = 2*M/T
            else:
                x_n[i-1] = 2*M/T
            curr_t -= 1
            last_d = float("inf")
        last_d = curr_d
        i += 1

def get_array_factor (theta, I_0, I_n):
    global M, d, k
    sum = 0
    for n in range(-1 * M, M):
        sum += (I_n[n + M]/I_0)*cmath.exp(1j*n*(k*d*math.cos(theta)))
    return sum

def main():
    global theta_n, x_n, N, M

    t = time.time()
    sample_angles()
    peak_approximator()

    X_k = [0] * N
    for k in range(N):
        sum = 0
        for n in range(N):
            sum += x_n[n] * cmath.exp((-2 * cmath.pi * 1j * (n-M) * k)/N)
        X_k[k] = sum

    plt.plot(theta_n, x_n, color = 'red', marker = "o")

    x1 = [i for i in range(181)]
    y1 = [0] * 181
    for i in x1:
        y1[i] = abs(get_array_factor(math.radians(i), X_k[0], X_k))
    plt.plot(x1, y1, color = 'blue')

    t2 = time.time()
    print (t2 - t)
    results = []
    for target in targets:
        plt.axvline(x=target, color = "red")
        results.append(abs(get_array_factor(math.radians(target), X_k[0], X_k)))
    print(results)
    plt.show()



if __name__ == '__main__':
    main()
