import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import sqrt

plt.style.use(['science'])
fig, ax = plt.subplots(figsize=(6, 4))
#fig = plt.figure()
#ax = fig.add_subplot(111)
xdata, ydata = [], []
vdata = []
vxdataR = []
vydataR = []
xdataR, ydataR = [], []

vxdataL = []
vydataL = []
xdataL, ydataL = [], []


indices = []
ln = ax.scatter(xdata, ydata, c=vdata, edgecolor = 'r', linewidth=0.3)


ax.set_xlabel(r'$x$ (m)')
ax.set_ylabel(r'$y$ (m)')
X_MAX, Y_MAX, X_MIN, Y_MIN = 0, 0, 0, 0
wXdata, wYdata = [], []
cax = fig.add_axes([0.92, 0.1, 0.03, 0.8])
cbar = fig.colorbar(ln, cax=cax)
cx = cbar.ax
cbar.mappable.set_clim(vmin=0, vmax=1)
cbar.set_label(r'$v$[m/s]' )
#tx = ax.set_title('Frame 0')
def readConfig():
    count = 0
    global X_MAX, Y_MAX, X_MIN, Y_MIN
    with open("configuration.txt", 'r') as f:
        for line in f:
            numbers = line.split(" ")
            if count == 0:
                numbers = np.asarray(numbers, dtype="float64")
                X_MIN, Y_MIN, X_MAX, Y_MAX = numbers
            else:
                numbers = np.asarray(numbers, dtype="float64")
                wXdata.append(numbers[0])
                wYdata.append(numbers[1])
                wXdata.append(numbers[2])
                wYdata.append(numbers[3])
            count += 1

def init():
    count = 0
    global X_MAX, Y_MAX, X_MIN, Y_MIN
    readConfig()
    ax.set_xlim(X_MIN - (X_MAX - X_MIN) * 0.1, X_MAX + (X_MAX - X_MIN) * 0.1)
    ax.set_ylim(Y_MIN - (Y_MAX - Y_MIN) * 0.1, Y_MAX + (Y_MAX - Y_MIN) * 0.1)
    Nw = len(wXdata)
    ax.plot([0., 7.5], [7.5, 7.5], 'r-.')
    ax.plot([0., 7.5], [-7.5, -7.5], 'r-.')
    ax.plot([7.5, 7.5], [-7.5, 7.5], 'r-.')
    for i in range(0, Nw - 1, 2):
        ax.plot(wXdata[i:i+2], wYdata[i:i+2], 'k')
    return ln,

M = ax.transData.get_matrix()
xscale = M[0, 0]
yscale = M[1, 1]
size = fig.get_size_inches() * fig.dpi
count = 0
f = open("data_pos.txt", 'r')

def update(frame):
    global count
    xdata, ydata = [], []
    vdata = []
    vxdata, vydata = [], []
    desired_data_width = []
    numbers = f.readline().split(" ")
    t = np.asarray(numbers[0], dtype='float16')
    N = np.asarray(numbers[1], dtype='int16')
    count += 1
    for i in range(0, N):
        data = f.readline().split(" ")
        data = np.asarray(data, dtype='float64')
        indices.append(data[0])
        xdata.append(data[1])
        ydata.append(data[2])
        vxdata.append(data[3])
        vydata.append(data[4])
        vdata.append(sqrt(data[3] * data[3] + data[4] * data[4]))
        desired_data_width.append(1 * xscale * ((data[-1])) ** 2)
        count += 1

    ln.set_offsets(np.column_stack((xdata, ydata)))
    vdata = np.asarray(vdata, dtype='float64')
    xdata = np.asarray(xdata)
    ydata = np.asarray(ydata)
    #qvL.set_offsets((xdata[xdata < 0], ydata[xdata < 0]))
    vx_desired = xdata[xdata < 0]#/np.sqrt(np.asarray(xdata)[xdata < 0]**2 + np.asarray(ydata)[xdata < 0]**2)
    vy_desired = ydata[xdata < 0]#/np.sqrt(np.asarray(xdata)[xdata < 0]**2 + np.asarray(ydata)[xdata < 0]**2)
    ind = np.random.randint(0, 300, size=20, dtype=int)
    qvL = ax.quiver(xdata[xdata < 0][ind], ydata[xdata < 0][ind], -xdata[xdata < 0][ind], -ydata[xdata < 0][ind], units='xy', scale = 1, color = 'r', zorder = 0, width = 0.1, headwidth = 3, ec = 'k')
    ydata_rf = 15*np.random.random_sample(20) - 7.5
    xdata_rf = 7.5*np.ones(20)
    xdata_r = np.zeros(20)
    ydata_r = (1)*np.random.random_sample(20) - 0.5
    qvR = ax.quiver(xdata_r, ydata_r, xdata_rf-xdata_r, ydata_rf, units='xy', scale = 1, color = 'b', zorder = 0, width = 0.1, headwidth = 3, ec = 'k')
    #qvR = ax.quiver(xdataR, ydataR, vxdataR, vydataL)
    ln.set_array(vdata)
    ln.set_sizes(desired_data_width)
    vmax = np.max(vdata)
    vmin = np.min(vdata)
    #ln.set_clim(vmin = vmin,vmax = vmax)
    #print(vmin, vmax)
    # Update the colorbar limits based on the new data
    #cbar.mappable.set_clim(vmin=np.min(vdata), vmax=np.max(vdata))
    #cbar.draw_all()
    #cbar.vmin = np.min(vdata)
    #cbar.vmax = np.max(vdata)
    #fig.colorbar(ln, cax=cx)
    return ln,
init()
update(0)
#ani = FuncAnimation(fig, update, frames=range(0, 1),
#                    init_func=init, blit=True, interval=200, repeat=False)
plt.show()

