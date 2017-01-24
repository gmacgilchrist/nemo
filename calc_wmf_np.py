
# coding: utf-8

# In[1]:

# Calculate water mass formation for a range of density values from NEMO data


# In[94]:

import numpy as np
import pandas as pd
import xarray.ufuncs as xu
import scipy.io as io
from netCDF4 import Dataset
import gsw
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')


# In[91]:

# SPECIFY DATASET LOCATION
rootdir = '/home/ocean1/DRAKKAR/ORCA025.L75-GJM189-S/'
# Densities to evaluate
rrange = np.arange(1026,1028.4,0.1)
dr = 0.1
nr = rrange.shape[0]
# Timesteps to evaluate at
times = np.arange(1,5,1)
nt = times.shape[0]
# NA region
xv = np.arange(700,1300,1)
yv = np.arange(500,1021,1)
# Load grid data
grd = Dataset(rootdir+'GRID/ORCA025.L75-GJM189_mesh_hgr.nc','r')
# Grid spacings
e1t = grd.variables['e1t'][0,yv,xv]
e2t = grd.variables['e2t'][0,yv,xv]


# In[92]:

# LOAD DATA
# Load data files
F=np.empty((nr,nt,))
M=np.empty((nr,nt,))
tcount=0
for t in times:
    tstr = str(t).zfill(4)
    print tstr
    # Load dataset
    ds = Dataset(rootdir+'symbolic_links/ORCA025.L75-GJM189_'+tstr+'_gridT.nc','r')
    # Calculate sea-surface density
    r0 = gsw.rho(ds.variables['vosaline'][0,0,yv,xv],ds.variables['votemper'][0,0,yv,xv],0)
    # Calculate buoyancy flux
    Din = calc_buoyancyflux(ds,xv,yv)
    # Calculate watermass transformation and formation
    rcount = 0
    for r in rrange:
        (F[rcount,tcount],M[rcount,tcount]) = calc_wmt(Din,r0,e1t,e2t,r,dr)
        rcount+=1

    tcount+=1


# In[96]:

io.savemat(rootdir+'matfiles/wmt_r'+str(rrange.min())+'-'+str(rrange.max())+'_dr'+str(dr)+'_NA.mat',{'F':F,'M':M,'rrange':rrange})


# In[86]:

def calc_buoyancyflux(ds,xv,yv):
    import gsw
    
    # Calculate buoyancy flux
    r0 = gsw.rho(ds.variables['vosaline'][0,0,yv,xv],ds.variables['votemper'][0,0,yv,xv],0)
    alpha = gsw.alpha(ds.variables['vosaline'][0,0,yv,xv],ds.variables['votemper'][0,0,yv,xv],0)
    beta = gsw.beta(ds.variables['vosaline'][0,0,yv,xv],ds.variables['votemper'][0,0,yv,xv],0)
    D_T = -(alpha/gsw.cp0)*ds.variables['sohefldo'][0,yv,xv]
    D_S = r0*beta*ds.variables['vosaline'][0,0,yv,xv]*ds.variables['sowaflup'][0,yv,xv]/1000
    
    return D_T + D_S

def calc_wmt(Din,r0,e1t,e2t,r,dr):
    import numpy as np
    
    # Multiply by area of grid cell
    Din_area = Din*e1t*e2t
    
    # Calculate summed buoyancy flux for surface densities greater than r
    # Generate new array for values above and below the target density
    ind = r0>r-dr
    Dbelow = Din_area[ind].sum()
    ind = r0>r
    Dtarget = Din_area[ind].sum()
    ind = r0>r+dr
    Dabove = Din_area[ind].sum()
    # Perform differentiation on these arrays
    
    # TRANSFORMATION
    # First order central difference of D wrt rho
    F = -(1e-6)*(Dabove-Dbelow)/(2*dr)
    # FORMATION
    # Second order central differences of D wrt rho
    M = (1e-6)*(Dabove-2*Dtarget+Dbelow)/(dr**2)
    
    return (F,M)

