import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.patches as mpatches

file = open("nicolson-k-space.txt", "r")
dat = file.read()
file.close()

text_buf = ""
sep = []
sepsep = []
for c in dat:
    if c == '>':
        break
    elif c == '\n':
        sepsep.append(np.array(sep))
        sep = []
    elif c == '|':
        sep.append(float(text_buf))
        text_buf = ""
    else:
        text_buf += c

size = len(sepsep[0])
tim = len(sepsep)
interval = 20
height = max(sepsep[0])

x = np.array(sep)
t = np.linspace(-750.0, 750.0, size)

fig, ax = plt.subplots()
line, = ax.plot([], [], "k")
xdata, ydata = [], []
ax.set_xlim([-750, 750])
ax.set_ylim([0, height])

def animation(data):
    line.set_xdata(t)
    line.set_ydata(sepsep[data])
    return line,

# wall = mpatches.Rectangle((0.7, 0.0), 0.03, 1.2, fill= False, color = "black", linewidth=0.5)
# fig.gca().add_patch(wall)

ani = anim.FuncAnimation(fig, animation, frames=tim, interval=interval, blit=True)
ani.save("k-space-1d.mp4", fps=30, extra_args=['-vcodec', 'libx264'])

# plt.show()
