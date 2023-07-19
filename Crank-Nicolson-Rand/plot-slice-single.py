import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.patches as mpatches

file = open("nicolson-rand.txt", "r")
dat = file.read()
# print(dat)
file.close()

count = 0
text_buf = ""
sep = []
sepsep = []

for c in dat:
    if c == '>':
        break
    elif c == '\n':
        sepsep.append(np.array(sep))
        sep = []
        print(count)
        count += 1
    elif c == '|':
        sep.append(float(text_buf))
        text_buf = ""
    else:
        text_buf += c

sepsep = np.array(sepsep)
size = len(sepsep[0])
tim = len(sepsep)
interval = 20
height = max(sepsep.reshape(size*tim))

t = np.linspace(0.0, 1.0, size)

# plt.plot(sepsep[0])

fig, ax = plt.subplots()
line, = ax.plot([], [], "k")
xdata, ydata = [], []
ax.set_xlim([0, 1])
ax.set_ylim([0, height])

def animation(data):
    line.set_xdata(t)
    line.set_ydata(sepsep[data])
    return line,

ani = anim.FuncAnimation(fig, animation, frames=tim, interval=interval, blit=True)
ani.save("test.mp4", fps=15, extra_args=['-vcodec', 'libx264'])

# plt.show()
