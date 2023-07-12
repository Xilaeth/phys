import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import matplotlib.animation as anim

# v is an NxN Matrix
def gauss_seidel_tri_vec(v, l, limit):
    size = v.shape[0]
    vec = np.zeros_like(v, complex)
    a = np.zeros_like(v, complex)
    b = np.zeros_like(v, complex)
    for rep in range(limit):
        vec_new = np.zeros_like(v, complex)
        # calculating A-matrix
        for i in range(1, size-1):
            for j in range(1, size-1):
                if rep == 0:
                    b[i, j] = (l/2) * (v[i-1, j] + v[i+1, j] + v[i, j-1] + v[i, j+1]) + (1 - 2*l) * v[i, j]
                a[i, j] = (l/2) * (vec_new[i-1, j] + vec_new[i, j-1] + vec[i+1, j] + vec[i, j+1])
        vec_new = (b - a) / (1 + 2*l)
        if linalg.norm(vec - vec_new) < 1e-10:
            print("Rep: ", rep, " | Error: ", linalg.norm(vec-vec_new))
            break
        if rep == limit-1:
            print("Limit reached")
            exit()
        vec = vec_new
    return vec

def gauss_seidel_tri(v, l, limit):
    size = v.shape[0]
    vec = np.zeros_like(v, complex)
    b = np.zeros_like(v, complex)
    for rep in range(limit):
        vec_new = np.zeros_like(v, complex)
        for i in range(1, size-1):
            for j in range(1, size-1):
                # calculating the b-matrix
                if rep == 0:
                    b[i, j] = (l/2) * (v[i-1, j] + v[i+1, j] + v[i, j-1] + v[i, j+1]) + (1 - 2*l) * v[i, j]
                a = (l/2) * (vec_new[i-1, j] + vec_new[i, j-1] + vec[i+1, j] + vec[i, j+1])
                vec_new[i, j] = (b[i, j] - a) / (1 + 2*l)
        if linalg.norm(vec - vec_new) < 1e-10:
            print("Rep: ", rep, " | Error: ", linalg.norm(vec-vec_new))
            break
        if rep == limit-1:
            print("Limit reached")
            exit()
        vec = vec_new
    return vec

# calculating B-Matrix
def inhom(v, l):
    size = v.shape[0]
    vec = np.zeros_like(v, complex)
    for i in range(1, size-1):
        for j in range(1, size-1):
            vec[i, j] = (l/2) * (v[i-1, j] + v[i+1, j] + v[i, j-1] + v[i, j+1]) + (1 - 2*l) * v[i, j]
    return vec
    
limit = 1000
tn = 10000
xn = 200
yn = 200
x = 1
y = 1
t = 0.004
dt = complex(t / tn)
dx = complex(x / xn)
dy = complex(y / yn)
m = 1
hbar = 1
k = np.array([500, 500], complex)

fig, ax = plt.subplots()

sigma = x/100
mean_x = x/2
mean_y = x/2

lamb = 1j*((dt*hbar)/(2*m*(dx**2)))

v0 = np.array(list(map(lambda x0: list(map(lambda y0: np.e ** ((-( ((dx*x0 - mean_x)**2 + (dy*y0 - mean_y)**2)/sigma)) + 1j*np.dot(k, np.array([dx*x0, dy*y0], complex))), range(yn))), range(xn))), complex)

# plt.matshow(np.imag(v0))

mat = ax.matshow(np.absolute(v0))

skip = 5
interval = 40

vs = []
for i in range(tn):
    if i % skip == 0:
        vs.append(v0)
    v0 = gauss_seidel_tri_vec(v0, lamb, limit)
    s1 = np.array(sum(np.absolute(v0)**2))
    print(i, " | ", sum(s1), "\n----")

def animation(data):
    mat.set_data(np.absolute(vs[data])**2)
    return mat,

ani = anim.FuncAnimation(fig, animation, frames=(tn//skip), interval=interval, blit=True)
# ani.save("quantum2d-e-7.mp4", fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
