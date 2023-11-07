import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure()
ax = fig.add_subplot(111)

# I like to position my colorbars this way, but you don't have to
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')

def f(x, y):
    return np.exp(x) + np.sin(y)

x = np.linspace(0, 1, 120)
y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)

# This is now a list of arrays rather than a list of artists
frames = []
for i in range(10):
    x       += 1
    curVals  = f(x, y)
    frames.append(curVals)

cv0 = frames[0]
im = ax.imshow(cv0, origin='lower') # Here make an AxesImage rather than contour
cb = fig.colorbar(im, cax=cax)
tx = ax.set_title('Frame 0')

def animate(i):
    arr = frames[i]
    vmax     = np.max(arr)
    vmin     = np.min(arr)
    im.set_data(arr)
    im.set_clim(vmin, vmax)
    tx.set_text('Frame {0}'.format(i))
    # In this version you don't have to do anything to the colorbar,
    # it updates itself when the mappable it watches (im) changes

ani = animation.FuncAnimation(fig, animate, frames=10)

plt.show()
