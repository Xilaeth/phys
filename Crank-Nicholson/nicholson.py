import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import matplotlib.animation as ani

# the matrix has to be strictly diagonally dominant, meaning a1 + a3 < a2 OR positive definite and symmetric
def gauss_seidel_tri(b, a1, a2, a3, limit):
    size = b.shape[0]
    x = np.zeros(size, complex)
    for i in range(limit):
        xn = np.ones(size, complex)
        for j in range(size):
            a = 0
            if j == 0:
                a += a3*x[j+1]
            elif j == size-1:
                a += a1*x[j-1]
            else:
                a += a1*x[j-1] + a3*x[j+1]
            xn[j] = (b[j] - a) / a2
        print(xn)
        if linalg.norm(x - xn) < 1e-15:
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

limit = 100000
tn = 100
xn = 10
x = 1
t = 10
dt = t / tn
dx = x / xn
k = 0.1

# test input
# -------
# M = np.array([[-1, 0.5, 0],
              # [0.5, -1, 0.5],
              # [0, 0.5, -1]], complex)

b = np.array(list(map(lambda y: np.e ** (-(((y*dx)-0.5)**2)/(0.1)) * np.e ** (1j * (y*dx) * k), range(0, xn))))
t = np.linspace(0, x, xn)
print(np.absolute(b))
#M2 = linalg.inv(M)
# -------

v = gauss_seidel_tri(b, 0.4, 1, 0.4, limit)
# print(v)
v2 = gauss_seidel_tri(v, 0.4, 1, 0.4, limit)
# print(v2)


plt.plot(t, np.absolute(b), "r")
plt.plot(t, np.absolute(v), "b")
plt.plot(t, np.absolute(v2), "g")
plt.show()
