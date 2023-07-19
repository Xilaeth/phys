import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.patches as mpatches

file = open("gauss-seidel-comp.txt", "r")
dat = file.read()
file.close()

count = 0
text_buf = ""
sep = []
for c in dat:
    if c == '>':
        break
    elif c == '|':
        sep.append(float(text_buf))
        text_buf = ""
    else:
        text_buf += c

tim = len(sep)
interval = 20

x = np.array(sep)
t = np.linspace(50.0, 2000.0, int(2000.0/50.0))

plt.xlim([0, 2000])
plt.ylim([0, max(sep)])
plt.title("X-Steps: 500")
plt.xlabel("t-steps")
plt.ylabel("no. of computations")
plt.plot(t, x)

plt.show()
