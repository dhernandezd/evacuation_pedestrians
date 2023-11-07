import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
f = np.loadtxt("kinetic_zone_wd.cg")
#plt.style.use(['science'])
vf = f[0:32*32,1]
#vf = vf.reshape(32, 32)
fig, ax = plt.subplots(figsize=(4, 4))
#ln = ax.imshow(np.rot90(vf), extent=[-15.0, 0, -7.5,7.5])
cax = fig.add_axes([0.92, 0.1, 0.03, 0.8])
cbar = fig.colorbar(ln, cax=cax)
#cbar.mappable.set_clim(vmin=0, vmax=0.9)

cbar.set_label(r'$\rho$' )
for i in range(0,1000):
    vf += f[i*32*32:(i+1)*32*32,1]
    #vf = vf.reshape(32, 32)
    #ln.set_array(np.rot90(vf))
    #return ln,
#ani = FuncAnimation(fig, animate, frames=range(0, 1000),blit=True, interval=200, repeat=False)
vf = vf.reshape(32, 32)
ln = ax.imshow(np.rot90(vf)/1000, extent=[-2.0, 0, -1,1], interpolation="bilinear")
#fig.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(4, 4))
ax.plot(vf[16,:])
plt.show()

