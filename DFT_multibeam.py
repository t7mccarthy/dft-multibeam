import cmath
import math
import matplotlib.pyplot as plt

lmbda = 5.168835482759
# k = 2 * math.pi / lmbda
d = lmbda / 2
M = 56
N = M * 2 + 1
k = 2 * math.pi / lmbda
theta_deg = 100
theta = math.radians(theta_deg)
theta_n = [0] * N
x_n = [0] * N
AF = [0] * N


def func (x):
    global theta_deg
    y1 = (1/(math.sqrt(0.00018) * math.sqrt(2 * math.pi)))*math.exp((-2)*((x - theta_deg)/math.sqrt(0.5))**2)
    # y1 = 0
    y2 = (1/(math.sqrt(0.00018) * math.sqrt(2 * math.pi)))*math.exp((-2)*((x - 160)/math.sqrt(0.5))**2)

    y3 = (1/(math.sqrt(0.00018) * math.sqrt(2 * math.pi)))*math.exp((-2)*((x - 70)/math.sqrt(0.5))**2)

    y4 = (1/(math.sqrt(0.00018) * math.sqrt(2 * math.pi)))*math.exp((-2)*((x - 50)/math.sqrt(0.5))**2)
    return max(y1, y2, y3, y4)



def sample_angles():
    global lmbda, N, d , theta_n, M
    for n in range((-1 * M), M):
        theta_n[n+M] = math.degrees(math.acos((n*lmbda)/(N*d)))

def sample_curve():
    global theta_n, x_n, M, lmbda, N, d, AF
    for n in range((-1 * M), M):
        theta_n[n+M] = math.degrees(math.acos((n*lmbda)/(N*d)))
        AF[n+M] = abs(func(theta_n[n+M]))
        x_n[n+M] = func(theta_n[n+M])/N

def get_array_factor (theta, I_0, I_n):
    global M, d, k
    sum = 0
    for n in range(-1 * M, M):
        sum += (I_n[n + M]/I_0)*cmath.exp(1j*n*(k*d*math.cos(theta)))
    return sum

def main():
    global theta_n, x_n, N, M, AF
    sample_curve()
    print len(theta_n)

    X_k = [0] * N
    for k in range(N):
        sum = 0
        for n in range(N):
            sum += x_n[n] * cmath.exp((-2 * cmath.pi * 1j * (n-M) * k)/N)
        X_k[k] = sum
    print X_k

    x = [i for i in range(0,181)]
    y = [0] * 181
    for i in x:
        y[i] = abs(func(i))

    # plt.plot(x, y, color='green')
    plt.plot(theta_n, AF, color = 'red', marker = "o")
    theta_n.reverse()
    print theta_n
    x1 = [i for i in range(181)]
    y1 = [0] * 181
    for i in x1:
        y1[i] = abs(get_array_factor(math.radians(i), X_k[0], X_k))

    plt.plot(x1, y1, color = 'blue')

    plt.show()


if __name__ == '__main__':
    main()
