import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.patches as mpatches


# the matrix has to be strictly diagonally dominant, meaning a1 + a3 < a2 OR positive definite and symmetric
def gauss_seidel_tri(b, a1, a2, a3, b1, b2, b3, pot, limit):
    size = b.shape[0]
    x = np.zeros(size, complex)
    for i in range(limit):
        xn = np.zeros(size, complex)
        for j in range(1, size-1):
            a = 0
            a += a1*xn[j-1] + a3*x[j+1]
            xn[j] = ((b1*b[j-1] + (b2 + (dt/2)*pot(dx*j))*b[j] + b3*b[j+1]) - a) / (a2 - (dt/2)*pot(dx*j))
        if linalg.norm(x - xn) < 1e-8:
            print(i, "|", linalg.norm(x - xn))
            break
        if i == limit-1:
            print(xn)
            print("reached limit")
            exit()
        x = xn
    return x

# limit = int(input("Itteration limit: "))
# xn = int(input("X-steps: "))
# t = float(input("Time: "))
# tn = int(input("T-steps: "))
# k = float(input("Wave number: "))

limit = 1000
tn = 1000
xn = 400
x = 1
t = 0.004
dt = t / tn
dx = x / xn
m = 1
hbar = 1
k = 500

fig, ax = plt.subplots()
line, = ax.plot([], [], "k")
xdata, ydata = [], []
ax.set_xlim([0, x])
ax.set_ylim([0, 1.2])

sigma = 0.001
mean = 0.3
v = np.array(list(map(lambda y: np.e ** (-(((dx*y)-mean)**2)/sigma + 1j*k*(dx*y)), range(0, xn))))
v[0] = 0.0
v[xn-1] = 0.0
print(v)
x_range = np.linspace(0, x, xn)

a = 1j/hbar
b = (dt*(hbar**2))/(4*m*(dx**2))
a1, a2, a3 = -a*b, (1+2*a*b), -a*b
b1, b2, b3 = a*b, (1-2*a*b), a*b
print(a)

def V(x):
    if x > 0.7 and x < 0.72:
        return -50000.0
    else:
        return 0.0

skip = 5
interval = 20

vs = []
for i in range(tn):
    if i % skip == 0:
        vs.append(v)
    v = gauss_seidel_tri(v, a1, a2, a3, b1, b2, b3, V, limit)

def animation(data):
    line.set_xdata(x_range)
    line.set_ydata(np.absolute(vs[data]))
    return line,

wall = mpatches.Rectangle((0.7, 0.0), 0.03, 1.2, fill= False, color = "black", linewidth=0.5)
fig.gca().add_patch(wall)

ani = anim.FuncAnimation(fig, animation, frames=(tn//skip), interval=interval, blit=True)
# ani.save("quantum_tunneling.mp4", fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
