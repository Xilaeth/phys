import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import matplotlib.animation as ani

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

def gauss_seidel_tri(A, b, a1, a2, a3, limit):
    size = b.shape[0]
    x = np.zeros(size, float)
    for i in range(limit):
        xn = np.ones(size, float)
        for j in range(size):
            a = 0
            if j == 0:
                a += a3*x[j+1]
            elif j == size-1:
                a += a1*x[j-1]
            else:
                a += a1*x[j-1] + a3*x[j+1]
            xn[j] = (b[j] - a) / a2
        if linalg.norm(x - xn) < 1e-10:
            print(i, "|", linalg.norm(x - xn))
            break
        if i == limit-1:
            print("reached limit")
            exit()
        x = xn
    return x

tn = 100
xn = 100
x = 1
t = 10
dt = int(t / tn)
dx = int(x / xn)
