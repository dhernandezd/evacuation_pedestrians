import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.colors import Normalize
#from matplotlib.cm import ScalarMappable 
from matplotlib.animation import FuncAnimation
from math import sqrt
#plt.style.use(['science'])
plt.style.use(['science','notebook'])
fig, ax = plt.subplots(figsize=(6,6))
xdata, ydata = [], []
ln, = ax.plot(xdata, ydata)
ax.set_xlabel(r'$t$(s)')
ax.set_ylabel(r'$N$')
ax.set_xlim(0, 20)
ax.set_ylim(0, 100)
    #plt.savefig('flux_b{}.pdf'.format(b), bbox_inches="tight" ,format = 'pdf', dpi=350)
    
f = open("flux.txt", 'r')
def update(frame):
   data = f.readline().split(" ")
   data = np.asarray(data, dtype = 'float64')
   xdata.append(data[0])
   ydata.append(data[1])
   ln.set_data(xdata, ydata)
   return ln,
ani = FuncAnimation(fig, update, frames = range(0,200), blit=True, interval = 200, repeat = False)
ani.save('animation_flux.gif', writer='imagemagick', fps=24)

