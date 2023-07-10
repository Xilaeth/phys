import numpy as np

A = np.array([[2, 3], [-1, 4]])
b = np.array([1, 0])

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

v = gauss_seidel(A, b, 1000)
print(b)
print(v)
print(np.matmul(A, v))
