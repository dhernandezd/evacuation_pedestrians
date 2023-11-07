
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
from scipy import linalg
import matplotlib.pyplot as plt
#from matplotlib.colors import Normalize
#from matplotlib.cm import ScalarMappable 
from matplotlib.animation import FuncAnimation
from math import sqrt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
#ln, = plt.plot([], [], 'ro')

frames = range(0,1000)
vd = np.arange(0.5, 10.5, 0.5)
velocity = sys.argv[1]
vd = [velocity]
RHO_G = 152.788745368
D_G = 0.5

for vd_i in vd:
    rho_t = []
    vx_t = []
    vy_t = []
    trk_t = []
    trc_t = []
    trcs_t = []
    trs_t = []
    devk_t = []
    devc_t = []
    devcs_t = []
    devs_t = []
    devt_t = []
    eta_t = []
    eta2_t = []
    mu_t = []
    I_t = []
    etaC_t = []
    etaCS_t = []
    etaS_t = []
    eta2C_t = []
    eta2S_t = []
    muC_t = []
    muCS_t = []
    muS_t = []
    IC_t = []
    IS_t = []
    gamma_dot_t = []
    #root = '../clusterdisk/vd_{}/b_2.00/'.format(vd_i)
    nameK = 'kinetic_zone_wd.cg'
    nameK2 = 'kinetic_vaverage.cg'
    nameC = 'contact_zone_wd.cg'
    nameS = 'social_zone_wd.cg'
    ROWS = 1024
    COUNT = 0
    MAX_ROWS = 1024000
    NUM_TIMES= int(MAX_ROWS/ROWS)
    print(NUM_TIMES)
    vf, vx, vy, kxx, kxy, kyy= np.loadtxt(nameK, unpack = True, skiprows=COUNT, max_rows=ROWS)
    nx = 32
    ny = 32
    x = np.arange(-2, 0., 0.5/8.)
    y = np.arange(-1, 1., 0.5/8.)
    
    X,Y = np.meshgrid(x,y)
    x = X.reshape(nx,ny)
    y = Y.reshape(nx,ny)
    #zone = (x < 0) & (x > -1.0) & (y < 0.5) & (y > -0.5)
    #zone = (x < 0) & (x > -1.0) & (y < 0.5) & (y > -0.5)
    for it in range(NUM_TIMES):
        print(it)
        vf, vx, vy, kxx, kxy, kyy= np.loadtxt(nameK, unpack = True, skiprows=COUNT, max_rows=ROWS)
        cxx, cxy, cyx, cyy, csxx, csxy, csyx, csyy= np.loadtxt(nameC, unpack = True, skiprows=COUNT, max_rows=ROWS)
        sxx, sxy, syx, syy= np.loadtxt(nameS, unpack = True, skiprows=COUNT, max_rows=ROWS)
        COUNT+=ROWS
        vf = vf.reshape(nx,ny)
        vel = np.sqrt(vx**2 + vy**2)
        vel = vel.reshape(nx,ny)
        vx = vx.reshape(nx,ny)
        vy = vy.reshape(nx,ny)
        trc = (cxx + cyy)/2.
        trc = trc.reshape(nx,ny)
        trcs = (csxx + csyy)/2.
        trcs = trcs.reshape(nx,ny)
        trs = (sxx + syy)/2.
        trs = trs.reshape(nx,ny)
        trk = (kxx + kyy)/2.
        trk = trk.reshape(nx,ny)
            

        dx = 0.5/8.
        dy = 0.5/8.

        dvxdx = np.gradient(vx, dx, axis=0)
        dvxdy = np.gradient(vx, dy, axis=1)
        dvydx = np.gradient(vy, dx, axis=0)
        dvydy = np.gradient(vy, dy, axis=1)

        eta = np.zeros((nx, ny))
        eta2 = np.zeros((nx, ny))
        eta2S = np.zeros((nx, ny))
        eta2C = np.zeros((nx, ny))
        etaS = np.zeros((nx, ny))
        etaC = np.zeros((nx, ny))
        etaCS = np.zeros((nx, ny))
        mu = np.zeros((nx, ny))
        I = np.zeros((nx, ny))
        muS = np.zeros((nx, ny))
        muC = np.zeros((nx, ny))
        muCS = np.zeros((nx, ny))
        IS = np.zeros((nx, ny))
        IC = np.zeros((nx, ny))
        gamma_d = np.zeros((nx, ny))
        devt = np.zeros((nx, ny))
        devs = np.zeros((nx, ny))
        devc = np.zeros((nx, ny))
        devcs = np.zeros((nx, ny))
        devk = np.zeros((nx, ny))
        for i in range(nx)[-32:]:
            for j in range(ny)[int(ny/2) - 16: int(ny/2) + 16]:
                grad_v = [dvxdx[i,j], dvxdy[i,j], dvydx[i,j], dvydy[i,j]]
                ind = i*ny + j
                T = [RHO_G*kxx[ind] + cxx[ind] + sxx[ind], RHO_G*kxy[ind] + cxy[ind] + sxy[ind], 
                RHO_G*kxy[ind] + cyx[ind] + syx[ind], RHO_G*kyy[ind] + cyy[ind] + syy[ind]]
                T = np.asarray(T)
               	T = T.reshape(2,2)
                TC = [cxx[ind], cxy[ind], cyx[ind], cyy[ind]]
                TC = np.asarray(TC)
                TC = TC.reshape(2,2)
                TCS = [csxx[ind], csxy[ind], csyx[ind], csyy[ind]]
                TCS = np.asarray(TCS)
                TCS = TCS.reshape(2,2)
                TS = [sxx[ind], sxy[ind], syx[ind], syy[ind]]
                TS = np.asarray(TS)
                TS = TS.reshape(2,2)
                grad_v = np.asarray(grad_v)
                grad_v = grad_v.reshape(2,2)
                grad_vS = (grad_v + grad_v.T)/2.
                trace_grad = np.sum(grad_vS*grad_vS) #np.trace(np.matmul(grad_vS.T, grad_vS))
                gamma_dot = 0.5*np.sqrt((grad_v[0,1] + grad_v[1,0])**2
                + (grad_v[0,0] - grad_v[1,1])**2)
                gamma_d[i,j] = gamma_dot
                pressure = (T[0,0] + T[1,1])/2.
                pressureC = (TC[0,0] + TC[1,1])/2.
                pressureCS = (TCS[0,0] + TCS[1,1])/2.
                pressureS = (TS[0,0] + TS[1,1])/2.
                I[i,j] = D_G*gamma_dot/np.sqrt(pressure/RHO_G)
                IS[i,j] = D_G*gamma_dot/np.sqrt(pressureS/RHO_G)
                IC[i,j] = D_G*gamma_dot/np.sqrt(pressureC/RHO_G)
                sigma = np.real(linalg.eigvals(T))
                devt[i,j] = 0.5*np.sqrt((sigma[0] - sigma[1])**2)
                if pressure!=0:
                    mu[i,j] = np.sqrt((sigma[0] - sigma[1])**2)/(np.sqrt(6.)*pressure)
                else:
                    mu[i,j] = 0
                sigma = np.real(linalg.eigvals(TS))
                devs[i,j] = 0.5*np.sqrt((sigma[0] - sigma[1])**2)
                if pressureS!=0:
                    muS[i,j] = np.sqrt((sigma[0] - sigma[1])**2)/(np.sqrt(6.)*pressureS)
                else:
                    muS[i,j] = 0
                sigma = np.real(linalg.eigvals(TC))
                devc[i,j] = 0.5*np.sqrt((sigma[0] - sigma[1])**2)
                if pressureC!=0:
                    muC[i,j] = np.sqrt((sigma[0] - sigma[1])**2)/(np.sqrt(6.)*pressureC)
                else:
                    muC[i,j] = 0
                sigma = np.real(linalg.eigvals(TCS))
                devcs[i,j] = 0.5*np.sqrt((sigma[0] - sigma[1])**2)
                if pressureCS!=0:
                    muCS[i,j] = np.sqrt((sigma[0] - sigma[1])**2)/(np.sqrt(6.)*pressureCS)
                else:
                    muCS[i,j] = 0
                if trace_grad > 1.e-6:
                    #eta[i,j] = 0.5*(np.trace(np.matmul(T.T, grad_vS))/trace_grad)
                    eta[i,j] = 0.5*(np.sum(T*grad_vS)/trace_grad)
                    eta2C[i,j] = 0.5*(np.sum(TC*grad_vS)/trace_grad)
                    etaCS[i,j] = 0.5*(np.sum(TCS*grad_vS)/trace_grad)
                    eta2S[i,j] = 0.5*(np.sum(TS*grad_vS)/trace_grad)
                    eta2[i,j] = mu[i,j]*pressure/gamma_dot
                    etaS[i,j] = muS[i,j]*pressureS/gamma_dot
                    etaC[i,j] = muC[i,j]*pressureC/gamma_dot
                else:
                    eta[i,j] = 0.
                    eta2[i,j] = 0.
                    eta2S[i,j] = 0.
                    eta2C[i,j] = 0.
                    etaCS[i,j] = 0.
                    etaS[i,j] = 0.
                    etaC[i,j] = 0.
        
                
        #zone = (x < 0) & (x > -1.0) & (y < 0.5) & (y > -0.5)
        """rho_m = np.mean(vf[zone])
        vx_m = np.mean(vx[zone])
        vy_m = np.mean(vy[zone])
        trk_m = np.mean(trk[zone])
        trc_m = np.mean(trc[zone])
        trs_m = np.mean(trs[zone])
        eta_m = np.mean(eta[zone])
        eta2_m = np.mean(eta2[zone])
        mu_m = np.mean(mu[zone])
        I_m = np.mean(I[zone])
	"""
        rho_m = np.mean(vf[-32:,int(ny/2) - 16: int(ny/2) + 16])
        vx_m = np.mean(vx[-32:,int(ny/2) - 16: int(ny/2) + 16])
        vy_m = np.mean(vy[-32:,int(ny/2) - 16: int(ny/2) + 16])
        trk_m = np.mean(trk[-32:,int(ny/2) - 16: int(ny/2) + 16])
        trc_m = np.mean(trc[-32:,int(ny/2) - 16: int(ny/2) + 16])
        trs_m = np.mean(trs[-32:,int(ny/2) - 16: int(ny/2) + 16])
        trcs_m = np.mean(trcs[-32:,int(ny/2) - 16: int(ny/2) + 16])
        devt_m = np.mean(devt[-32:,int(ny/2) - 16: int(ny/2) + 16])
        devc_m = np.mean(devc[-32:,int(ny/2) - 16: int(ny/2) + 16])
        devs_m = np.mean(devs[-32:,int(ny/2) - 16: int(ny/2) + 16])
        devcs_m = np.mean(devcs[-32:,int(ny/2) - 16: int(ny/2) + 16])
        eta_m = np.mean(eta[-32:,int(ny/2) - 16: int(ny/2) + 16])
        etaS_m = np.mean(etaS[-32:,int(ny/2) - 16: int(ny/2) + 16])
        etaC_m = np.mean(etaC[-32:,int(ny/2) - 16: int(ny/2) + 16])
        etaCS_m = np.mean(etaCS[-32:,int(ny/2) - 16: int(ny/2) + 16])
        eta2_m = np.mean(eta2[-32:,int(ny/2) - 16: int(ny/2) + 16])
        eta2S_m = np.mean(eta2S[-32:,int(ny/2) - 16: int(ny/2) + 16])
        eta2C_m = np.mean(eta2C[-32:,int(ny/2) - 16: int(ny/2) + 16])
        mu_m = np.mean(mu[-32:,int(ny/2) - 16: int(ny/2) + 16])
        muS_m = np.mean(muS[-32:,int(ny/2) - 16: int(ny/2) + 16])
        muC_m = np.mean(muC[-32:,int(ny/2) - 16: int(ny/2) + 16])
        muCS_m = np.mean(muCS[-32:,int(ny/2) - 16: int(ny/2) + 16])
        I_m = np.mean(I[-32:,int(ny/2) - 16: int(ny/2) + 16])
        IS_m = np.mean(IS[-32:,int(ny/2) - 16: int(ny/2) + 16])
        IC_m = np.mean(IC[-32:,int(ny/2) - 16: int(ny/2) + 16])
        gamma_m = np.mean(gamma_d[-32:,int(ny/2) - 16: int(ny/2) + 16])
	
        rho_t.append(rho_m)
        vx_t.append(vx_m)
        vy_t.append(vy_m)
        trk_t.append(trk_m)
        trc_t.append(trc_m)
        trcs_t.append(trcs_m)
        trs_t.append(trs_m)
        devt_t.append(devt_m)
        devc_t.append(devc_m)
        devcs_t.append(devcs_m)
        devs_t.append(devs_m)
        eta_t.append(eta_m)
        etaS_t.append(etaS_m)
        etaC_t.append(etaC_m)
        etaCS_t.append(etaCS_m)
        eta2S_t.append(eta2S_m)
        eta2C_t.append(eta2C_m)
        eta2_t.append(eta2_m)
        I_t.append(I_m)
        IC_t.append(IC_m)
        IS_t.append(IS_m)
        mu_t.append(mu_m)
        muS_t.append(muS_m)
        muC_t.append(muC_m)
        muCS_t.append(muC_m)
        gamma_dot_t.append(gamma_m)
    
    np.savetxt('spatialNewwd_average_{}_D_zone_1.txt'.format(vd_i),
    np.column_stack((rho_t, vx_t, vy_t, trk_t, trc_t, trs_t, trcs_t, eta_t,
    eta2C_t, eta2S_t, etaCS_t,eta2_t, etaC_t, etaS_t, mu_t, muC_t, muS_t, muCS_t,I_t,
    IC_t, IS_t, gamma_dot_t)))
    np.savetxt('spatialNewwd_dev_{}_D_zone_1.txt'.format(vd_i),
    np.column_stack((devt_t, devc_t, devcs_t, devs_t)))


