import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import matplotlib.animation as ani

# the matrix has to be strictly diagonally dominant, meaning a1 + a3 < a2 OR positive definite and symmetric
def gauss_seidel_tri(b, a1, a2, a3, b1, b2, b3, pot, limit):
    size = b.shape[0]
    x = np.zeros(size, complex)
    for i in range(limit):
        xn = np.ones(size, complex)
        for j in range(size):
            a = 0
            if j == 0:
                a += a3*x[j+1]
                xn[j] = ((b2*b[j] + b3*b[j+1]) - a) / a2
            elif j == size-1:
                a += a1*x[j-1]
                xn[j] = ((b1*b[j-1] + b2*b[j]) - a) / a2
            else:
                a += a1*x[j-1] + a3*x[j+1]
                xn[j] = ((b1*b[j-1] + b2*b[j] + b3*b[j+1]) - a) / a2
        # print(xn)
        if linalg.norm(x - xn) < 1e-8:
            print(i, "|", linalg.norm(x - xn))
            break
        if i == limit-1:
            print("reached limit")
            exit()
        x = xn
    return x

# limit = int(input("Itteration limit: "))
# xn = int(input("X-steps: "))
# t = float(input("Time: "))
# tn = int(input("T-steps: "))
# k = float(input("Wave number: "))

limit = 10000
tn = 1000
xn = 100
x = 1
t = 10
dt = t / tn
dx = x / xn
m = 1
hbar = 1
k = 1

fix, ax = plt.subplots()
ax.set_xlim([0, x])
ax.set_ylim([0, 1.2])

sigma = 0.01
mean = 0.3
v1 = np.array(list(map(lambda y: np.e ** (-(((y*dx)-mean)**2)/(sigma)* (np.e ** (1j * k * (dx*y)))), range(0, xn))))
t = np.linspace(0, x, xn)

a = 1j/hbar
b = (dt*(hbar**2))/(4*m*(dx**2))
a1, a2, a3 = -a*b, (1+2*a*b), -a*b
b1, b2, b3 = a*b, (1-2*a*b), a*b

def V(x):
    return 0

for i in range(50):
    if i % (tn//5) == 0:
        ax.plot(t, np.absolute(v1))
    v1 = gauss_seidel_tri(v1, a1, a2, a3, b1, b2, b3,  V, limit)
    print(i)

plt.show()
