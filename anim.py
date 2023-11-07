
__author__ = "Dariel Hernandez-Delfin"
__copyright__ = "Copyright 2022, Simple visulizer"
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Dariel Hernandez-Delfin"
__email__ = "dhernandez@.bcamath.org"
__status__ = "Testing"
###################################################
#input: data_pos.txt                              #
#function: visualize simply pedestrian simulations#  
###################################################

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.colors import Normalize
#from matplotlib.cm import ScalarMappable 
from matplotlib.animation import FuncAnimation
from math import sqrt
plt.style.use(['science'])
fig, ax = plt.subplots(figsize=(10,6))
xdata, ydata = [], []
vdata = []
indices = []
#ln, = plt.plot([], [], 'ro')
ln = plt.scatter(xdata, ydata, c = vdata, edgecolor = 'r', linewidth=0.3)
ax.set_xlabel(r'$x$(m)')
ax.set_ylabel(r'$y$(m)')
X_MAX, Y_MAX, X_MIN, Y_MIN = 0,0,0,0
wXdata, wYdata = [], []
cax = fig.add_axes([0.92, 0.1, 0.03, 0.8])
cbar = fig.colorbar(ln, cax=cax)
cx = cbar.ax
cbar.mappable.set_clim(vmin=0, vmax=1)
cbar.set_label(r'$v$[m/s]' )
ax.set_title(r'$t = 0.00$ s')
def readConfig():
    count = 0
    global  X_MAX, Y_MAX, X_MIN, Y_MIN
    with open("configuration.txt", 'r') as f:
        for line in f:
            numbers =  line.split(" ")
            if count==0 :
                numbers = np.asarray(numbers, dtype="float64")
                X_MIN, Y_MIN, X_MAX, Y_MAX = numbers
                #print(X_MIN, Y_MIN, X_MAX, Y_MAX)
            else:
                numbers = np.asarray(numbers, dtype="float64")
                wXdata.append(numbers[0])
                wYdata.append(numbers[1])
                wXdata.append(numbers[2])
                wYdata.append(numbers[3])
            count += 1
def init():
    count = 0
    global  X_MAX, Y_MAX, X_MIN, Y_MIN
    readConfig()
    ax.set_xlim(X_MIN - (X_MAX - X_MIN)*0.1, X_MAX + (X_MAX - X_MIN)*0.1)
    ax.set_ylim(Y_MIN - (Y_MAX - Y_MIN)*0.1, Y_MAX + (Y_MAX - Y_MIN)*0.1)
    Nw = len(wXdata)
    ax.plot([0.,7.5], [7.5, 7.5], 'r-.')
    ax.plot([0.,7.5], [-7.5, -7.5], 'r-.')
    ax.plot([7.5,7.5], [-7.5, 7.5], 'r-.')
    for i in range(0,Nw - 1,2):
        ax.plot(wXdata[i:i+2], wYdata[i:i+2], 'k')
        #print(wXdata[i:i+2], wYdata[i:i+2])
    return ln,
    
M = ax.transData.get_matrix()
xscale = M[0,0]
yscale = M[1,1]
size = fig.get_size_inches()*fig.dpi
print(xscale, yscale, size)
count = 0
f = open("data_pos.txt", 'r')
def update(frame):
    global count
    xdata, ydata = [], []
    vdata = []
    vxdata, vydata = [], []
    desired_data_width = []
    numbers = f.readline().split(" ")
    t = np.asarray(numbers[0], dtype = 'float16')
    N = np.asarray(numbers[1], dtype = 'int16')
    print(t)
    #import ipdb; ipdb.set_trace()
    count += 1
    t = frame*0.1
    for i in range(0, N):
        data = f.readline().split(" ")
        data = np.asarray(data, dtype = 'float64')
        indices.append(data[0])        
        #if data[0]==106:
        xdata.append(data[1])
        ydata.append(data[2])
        vxdata.append(data[3])
        vydata.append( data[4])
        vdata.append(sqrt(data[3]*data[3] + data[4]*data[4]))
        desired_data_width.append(1*xscale*((data[-1]))**2)
        #import ipdb; ipdb.set_trace()
        count += 1
    
    #ln.set_data(xdata, ydata)
    vdata = np.asarray(vdata, dtype = 'float64')
    
     
    ln.set_offsets(np.column_stack((xdata, ydata)))
    desired_data_width = np.asarray(desired_data_width)
    
    #Q.set_offsets(np.column_stack((xdata, ydata)))
    #Q.set_UVC(vxdata, vydata)
    vmax = np.max(vdata)
    vmin = np.min(vdata)
    """if vmax > 0:
        vdata = (vdata - vmin)/(vmax)
    else:
        vdata = 10*np.ones_like(vdata)
    """


    ln.set_array(vdata[:])
    ln.set_clim(vmin = 0.,vmax = 1.0)
    cbar.update_normal(ln)
    ax.set_title(r'$t = {:.1f}$ s'.format(t))
    
    #if update== 39:
        
    #ln.set_array(np.random.random(len(vdata))) 
    
    ln.set_sizes(desired_data_width)
    
    return ln,

ani = FuncAnimation(fig, update, frames = range(0,200),
                    init_func=init, blit=True, interval = 200, repeat = False)
ani.save('animation.gif', writer='imagemagick', fps=24)
#plt.show()
