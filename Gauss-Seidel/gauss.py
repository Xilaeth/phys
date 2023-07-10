import numpy as np
from numpy import linalg

A = np.array([[3, 2, 3],
              [4, 5, 6],
              [7, 8, 9]], float)
b = np.array([1, 2, 3], float)

# a1 = np.dot(A[j, :j-1], xn[:j-1])
# a2 = np.dot(A[j, j:], x[j:])
# print("A: ",A[j, j+1:])
# print("x: ",x[j+1:])
# xn[j] = (b[j] - a1 - a2) / A[j, j]

def gauss_seidel(A, b, limit):
    size = b.shape[0]
    x = np.ones(size, float)
    for i in range(limit):
        xn = np.ones(size, float)
        for j in range(size):
            a1 = 0
            a2 = 0
            for k in range(j):
                a1 += A[j, k] * xn[k]
            for k in range(j+1, size):
                a2 += A[j, k] * x[k]
            xn[j] = (b[j] - a1 - a2) / A[j, j]
        if linalg.norm(x - xn) < 1e-10:
            print(i, "|", linalg.norm(x - xn))
            break            
        if i == limit-1:
            print("reached limit")
            exit()
        x = xn
    return x

v = gauss_seidel(A, b, 100)
A2 = linalg.inv(A)
#print(linalg.det(A))
print(v)
print(np.matmul(A2, b))
