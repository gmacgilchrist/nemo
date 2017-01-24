
# coding: utf-8

# In[ ]:

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

