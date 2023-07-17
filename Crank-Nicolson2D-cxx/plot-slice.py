import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.patches as mpatches

file = open("nicolson2d-slice-comp-1.txt", "r")
file2 = open("nicolson2d-slice-comp-2.txt", "r")
dat = file.read()
dat2 = file2.read()
# print(dat)
file.close()
file2.close()

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

count = 0
text_buf = ""
sep2 = []
sepsep2 = []

for c in dat2:
    if c == '>':
        break
    elif c == '\n':
        sepsep2.append(np.array(sep2))
        sep2 = []
        print(count)
        count += 1
    elif c == '|':
        sep2.append(float(text_buf))
        text_buf = ""
    else:
        text_buf += c

sepsep = np.array(sepsep)
sepsep2 = np.array(sepsep2)
size = len(sepsep[0])
tim = len(sepsep)
interval = 20
height = max(sepsep2.reshape(size*tim))

t = np.linspace(0.0, 1.0, size)

fig, ax = plt.subplots()
line, = ax.plot([], [], "k")
line2, = ax.plot([], [], "r")
xdata, ydata = [], []
xdata2, ydata2 = [], []
ax.set_xlim([0, 1])
ax.set_ylim([0, height])

def animation(data):
    line.set_xdata(t)
    line.set_ydata(sepsep[data])
    line2.set_xdata(t)
    line2.set_ydata(sepsep2[data])
    return line, line2,

ani = anim.FuncAnimation(fig, animation, frames=tim, interval=interval, blit=True)
ani.save("double-slit-slice-cxx.mp4", fps=15, extra_args=['-vcodec', 'libx264'])

# plt.show()
