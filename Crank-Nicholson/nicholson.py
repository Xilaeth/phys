import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

def gauss_seidel(A, b, k):
    size = b.shape[0]
    x = np.zeros(size, float)
    for i in range(k):
        xn = np.zeros(size, float)
        for i in range(A.shape[0]):
            a1 = np.dot(A[i, :i], xn[:i])
            a2 = np.dot(A[i, i+1 :], x[i+1 :])
            xn[i] = (b[i] - a1 - a2) / A[i, i]
            x = xn
    return x

def const_M(a1, a2, a3, size):
    v1 = np.full(size, a1)
    v2 = np.full(size-1, a2)
    v3 = np.full(size-1, a3)
    return np.diag(v2, -1) + np.diag(v1, 0) + np.diag(v3, 1)

a = 1
dx = 0.1
dt = 0.1

r = (a*dt)/(2* dx**2)

a1 = (1 + 2*r)
a2 = -r
a3 = -r
size = 10

M = const_M(a1, a2, a3, 100)

f1 = np.zeros(100)
f1[0] = 1
t = np.linspace(0, 100*dx, 100)

fig, ax = plt.subplot()
line, = ax.plot(t, f1)
ax.grid()

def gen_data(   

for i in range(5):
    plt.plot(t, f1)
    x = gauss_seidel(M, f1, 1000)
    f1 = x

plt.show()
