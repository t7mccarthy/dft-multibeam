import cmath
import math
import matplotlib.pyplot as plt
import time

lmbda = 5.168835482759
# k = 2 * math.pi / lmbda
d = lmbda / 2
M =100
N = M * 2 + 1
k = 2 * math.pi / lmbda
targets = [25, 50, 65, 70, 100, 125, 160]
T = len(targets)
# theta_deg = 100
# theta = math.radians(theta_deg)
theta_n = [0] * N
x_n = [0] * N
AF = [0] * N

# # example settings: time 0.699464797974
# def func (x):
#     global targets
#     max_y = 0
#     for t in targets:
#         temp = (1/(math.sqrt(0.00018) * math.sqrt(2 * math.pi)))*math.exp((-2)*((x - t)/math.sqrt(0.5))**2)
#         max_y = max(temp, max_y)
#     return max_y

# example settings: time = 0.715882062912
def func (x):
    global targets, T
    f_i = [0] * T
    for i in range(T):
        f_i[i] = (1/(math.sqrt(0.00018) * math.sqrt(2 * math.pi)))*math.exp((-2)*((x - targets[i])/math.sqrt(0.5))**2)
    return max(f_i)

def sample_angles():
    global lmbda, N, d , theta_n, M
    for n in range((-1 * M), M):
        theta_n[n+M] = math.degrees(math.acos((n*lmbda)/(N*d)))

# def sample_curve():
#     global theta_n, x_n, M, lmbda, N, d, AF
#     for n in range((-1 * M), M):
#         theta_n[n+M] = math.degrees(math.acos((n*lmbda)/(N*d)))
#         AF[n+M] = abs(func(theta_n[n+M]))
#         x_n[n+M] = func(theta_n[n+M])

def set_output(i):
    global x_n, AF, T
    print theta_n[i]
    # AF[i-1] = max(AF[i-1], M * 0)
    AF[i] = 2*M/T
    # AF[i+1] = M/4 * 0
    # x_n[i-1] = max(x_n[i-1], M * 0)
    x_n[i] = 2*M/T
    # x_n[i+1] = M/4 * 0

def peak_approximator():
    global targets, T, theta_n, M, N, x_n
    curr_t = T-1
    i = 0
    last_d = float("inf")
    print(len(theta_n))
    print(N)
    while curr_t >= 0:
        curr_d = targets[curr_t] - theta_n[i]
        # print(targets[curr_t], theta_n[i])
        if curr_d > 0:
            if curr_d < abs(last_d):
                set_output(i)
            else:
                set_output(i-1)
            curr_t -= 1
            last_d = float("inf")
        last_d = curr_d
        i += 1

def sample_curve():
    sample_angles()
    peak_approximator()

def get_array_factor (theta, I_0, I_n):
    global M, d, k
    sum = 0
    for n in range(-1 * M, M):
        sum += (I_n[n + M]/I_0)*cmath.exp(1j*n*(k*d*math.cos(theta)))
    return sum

def main():
    global theta_n, x_n, N, M, AF

    t = time.time()
    sample_curve()
    print len(theta_n)

    X_k = [0] * N
    for k in range(N):
        sum = 0
        for n in range(N):
            sum += x_n[n] * cmath.exp((-2 * cmath.pi * 1j * (n-M) * k)/N)
        X_k[k] = sum
    # print X_k

    x = [i for i in range(0,181)]
    y = [0] * 181
    for i in x:
        y[i] = abs(func(i))

    # plt.plot(x, y, color='yellow')
    plt.plot(theta_n, AF, color = 'red', marker = "o")
    theta_n.reverse()
    # print theta_n
    x1 = [i for i in range(181)]
    y1 = [0] * 181
    for i in x1:
        y1[i] = abs(get_array_factor(math.radians(i), X_k[0], X_k))

    plt.plot(x1, y1, color = 'blue')
    t2 = time.time()
    print (t2 - t)
    plt.show()



if __name__ == '__main__':
    main()
