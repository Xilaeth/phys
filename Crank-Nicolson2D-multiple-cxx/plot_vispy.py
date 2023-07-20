import numpy as np
from vispy import app, scene, io, color
import os

data_dir = "data/"

data_files = os.listdir(data_dir)

data_files.sort(key=int)
print(data_files)

count = 0
text_buf = ""
sep = []
sep2 = []
sepsep = []

for files in data_files:
    file = open(data_dir + files)
    dat = file.read()
    file.close()
    for c in dat:
        if c == '>':
            sepsep.append(np.array(sep2))
            sep2 = []
            print(count)
            count += 1
            break
        elif c == ';':
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
interval = 20
sepsep *= size

canvas = scene.SceneCanvas(keys='interactive', size=(size, size), show=True)
view = canvas.central_widget.add_view()

dx, dy = 1/size, 1/size

x = np.linspace(0, size, size)
y = np.linspace(0, size, size)
x, y = np.meshgrid(x, y)
Z = sepsep[0]

surface = scene.visuals.SurfacePlot(x=x, y=y, z=Z, shading='smooth')

view.add(surface)

view.camera = 'turntable'

def update(event):
    global count, colormap
    if count == 0:
        exit()
    new_z = sepsep[count]
    surface.set_data(z=new_z)
    frame = canvas.render()
    io.write_png(f'vispy-frames/frame-{count:03d}.png', frame)
    count -= 1
    canvas.update()

timer = app.Timer(connect=update, start=1/15)

app.run()
