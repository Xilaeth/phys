import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.patches as mpatches

file = open("nicolson2d.txt", "r")
dat = file.read()
# print(dat)
file.close()

text_buf = ""
sep = []
sep2 = []
sepsep = []
for c in dat:
    if c == '>':
        break
    elif c == '\n':
        sepsep.append(np.array(sep2))
        sep2 = []
    elif c == ";":
        sep2.append(np.array(sep))
        sep = []
    elif c == '|':
        sep.append(float(text_buf))
        text_buf = ""
    else:
        text_buf += c

sepsep = np.array(sepsep)
size = len(sepsep[0])
tim = len(sepsep)
interval = 2

print(np.shape(sepsep[0]))

x = np.array(sep)
t = np.linspace(0.0, 1.0, size)

fig = plt.figure()
ax = plt.axes(projection='3d')
X = range(size)
Y = range(size)
X, Y = np.meshgrid(X, Y)
plot = ax.plot_surface(X, Y, sepsep[0], cmap=plt.cm.bwr, rstride=1, cstride=1, linewidth=0.01, color="k")
ax.grid(False)

print(sepsep[0])
print(size)
print(tim)

def animation(data, other, plot):
    ax.clear()
    plot = ax.plot_surface(X, Y, sepsep[data], cmap=plt.cm.bwr, rstride=1, cstride=1, linewidth=0.01, color="k")
    ax.grid(False)
    ax.set_axis_off()
    return plot,

ani = anim.FuncAnimation(fig, animation, frames=tim, fargs=(sepsep, plot), interval=interval, blit=False)
# ani.save("2d-cxx.mp4", fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
