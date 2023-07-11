import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.patches as mpatches

# v is an NxN Matrix
def gauss_seidel_tri(v, l, limit):
    size = v.shape[0]
    x = np.zeros_like(v, complex)
    for rep in range(limit):
        x_new = np.zeros_like(v, complex)
        for i in range(1, size-1):
            for j in range(1, size-1):
                a = (l/2) * (x_new[i-1, j] + x_new[i+1, j] + x[i, j-1] + x[i, j+1])
                b = (l/2) * (v[i-1, j] + v[i+1, j] + v[i, j-1] + v[i, j+1]) + (1 - l) * v[i, j]
                x_new[i, j] =(b - a) / (1 + l)
        if linalg.norm(x - x_new) < 1e-8:
            print("Rep: ", rep)
            break
        if rep == limit-1:
            print("Limit reached")
            exit()
        x = x_new
    return x

limit = 1000
tn = 1000
xn = 50
yn = 50
x = 1
y = 1
t = 0.004
dt = t / tn
dx = x / xn
dy = y / yn
m = 1
hbar = 1
k = np.array([500, 500], float)

fig, ax = plt.subplots()

sigma = 0.005
mean_x = 0.5
mean_y = 0.5

lamb = 1j*(dt*hbar)/(2*m*(dx**2))

test_dat = np.random.random((16, 16))
test_dat_2 = gauss_seidel_tri(test_dat, lamb, limit)

v0 = np.array(list(map(lambda x0: list(map(lambda y0: np.e ** ((-( ((dx*x0 - mean_x)**2 + (dy*y0 - mean_y)**2)/sigma)) + 1j*np.dot(k, [dx*x0, dy*y0])), range(yn))), range(xn))))

mat = ax.matshow(np.absolute(v0))

skip = 5
interval = 20

vs = []
for i in range(tn):
    if i % skip == 0:
        vs.append(v0)
    v0 = gauss_seidel_tri(v0, lamb, limit)
    print(i, "\n----")

def animation(data):
    mat.set_data(np.absolute(vs[data]))
    return mat,

ani = anim.FuncAnimation(fig, animation, frames=(tn//skip), interval=interval, blit=True)
ani.save("quantum2d.mp4", fps=30)

# plt.show()
