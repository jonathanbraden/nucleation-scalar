#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from field_analysis import *

###### Data Listings #####


def energy_density_plot(rho,dt,dx,a=None):
    xv = dx*np.arange(rho.shape[0])
    tv = dt*np.arange(rho.shape[1])

    if a==None:
        f,a = plt.subplots()
    cs = a.contourf(xv,yv,rho,51,cmap='cividis',zorder=-10)
    a.set_rasterization_zorder(-5)
    #for c_ in cs.collections:
    #    c_.set_rasterized(True)

    a.set_xlabel(r'$\sqrt{V_0}\phi_0^{-1}x$')
    a.set_ylabel(r'$\sqrt{V_0}\phi_0^{-1}t$')

    cb = f.colorbar(cs,pad=0.01,frac=0.03,ax=a)
    cb.solids.set_rasterized(True)
    
    return a.get_figure(),a
    
if __name__=="__main__":
    pass
