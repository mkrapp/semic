import numpy as np
import matplotlib.pyplot as plt
import SurfacePhysics as sp

# read data
pre = 'transect'
forc = np.loadtxt('../example/data/'+pre+'_input.txt')
vali = np.loadtxt('../example/data/'+pre+'_output.txt')
var_names = ['tsurf', 'alb', 'swnet', 'smb', 'melt', 'acc', 'shf', 'lhf']
nx    = forc.shape[1]/len(var_names)
ntime = forc.shape[0]

# Load parameters from Fortran Namelist
par = sp.surface_physics.Surface_Param_Class()
sp.surface_physics.surface_physics_par_load(par,'test.nml')
par.nx = nx
sp.surface_physics.print_param(par)

# initialize boundary class
bnd = sp.surface_physics.Boundary_Opt_Class()
sp.surface_physics.surface_boundary_define(bnd,par.boundary)
sp.surface_physics.print_boundary_opt(bnd)

# initialize and allocate state class
state = sp.surface_physics.Surface_State_Class()
sp.surface_physics.surface_alloc(state,par.nx)

# set initial values for prognostic & diagnostic variables
state.tsurf = np.array([vali[0,i] for i in range(par.nx)])
state.hsnow = np.array([1.0 for i in range(par.nx)])
state.hice = np.array([0.0 for i in range(par.nx)])
state.alb = np.array([vali[0,nx+i] for i in range(par.nx)])
state.alb_snow = np.array([vali[0,nx+i] for i in range(par.nx)])
# set auxiliary variables to zero
state.qmr = np.array([0.0 for i in range(par.nx)])
state.qmr_res = np.array([0.0 for i in range(par.nx)])
state.lhf = np.array([0.0 for i in range(par.nx)])
state.shf = np.array([0.0 for i in range(par.nx)])
# set land/ice/ocean mask
state.mask = np.array([2 for i in range(par.nx)])
if nx > 1: state.mask[0:3] = 1

time = []
y1 = {}
y2 = {}
for v in var_names:
    y1[v] = np.zeros((par.nx,ntime))
    y2[v] = np.zeros((par.nx,ntime))

# populate validation data
y2['tsurf'] = vali[:,0:nx].T
y2['alb']   = vali[:,nx:2*nx].T
y2['swnet'] = vali[:,2*nx:3*nx].T
y2['smb']   = vali[:,3*nx:4*nx].T
y2['melt']  = vali[:,4*nx:5*nx].T
y2['acc']   = vali[:,5*nx:6*nx].T
y2['shf']   = vali[:,6*nx:7*nx].T
y2['lhf']   = vali[:,7*nx:].T

# loop nyears over daily time steps
nyears = 100
for y in range(nyears):
    for t in range(ntime):
        state.sf   = forc[t,0:nx]
        state.rf   = forc[t,nx:2*nx]
        state.swd  = forc[t,2*nx:3*nx]
        state.lwd  = forc[t,3*nx:4*nx]
        state.wind = forc[t,4*nx:5*nx]
        state.sp   = forc[t,5*nx:6*nx]
        state.rhoa = forc[t,6*nx:7*nx]
        state.qq   = forc[t,7*nx:8*nx]
        state.tt   = forc[t,8*nx:]

        # calculate surface energy and mass balance
        sp.surface_physics.surface_energy_and_mass_balance(state,par,bnd,t,y)

        # write output for last year
        if y == nyears-1:
            y1['tsurf'][:,t] = state.tsurf
            y1['alb'][:,t] = state.alb
            y1['swnet'][:,t] = (1.-state.alb)*y2['swnet'][:,t]
            y1['smb'][:,t] = state.smb
            y1['melt'][:,t] = state.melt
            y1['acc'][:,t] = state.acc
            y1['shf'][:,t] = state.shf
            y1['lhf'][:,t] = state.lhf
            time.append(12.*t/365.)

var = 'melt'
for n in range(nx):
    l, = plt.plot(time,y1[var][n,:],lw=2,label=n+1)
    plt.plot(time,y2[var][n,:],lw=2,c=l.get_color(),ls='--')
plt.xlim(0,12)
plt.grid()
plt.legend(loc=0)

plt.show()
