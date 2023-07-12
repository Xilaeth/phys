import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

file = open("data.txt", "r")
dat = file.read()
# print(dat)
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

x = np.array(sep)
t = np.linspace(0.0, 1.0, size)

fig, ax = plt.subplots()
line, = ax.plot([], [], "k")
xdata, ydata = [], []
ax.set_xlim([0, 1])
ax.set_ylim([0, 1.2])

print(sepsep[0])
print(size)

def animation(data):
    line.set_xdata(t)
    line.set_ydata(sepsep[data])
    return line,

ani = anim.FuncAnimation(fig, animation, frames=tim, interval=interval, blit=True)

plt.show()
