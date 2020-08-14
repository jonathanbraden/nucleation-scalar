import numpy as np
import matplotlib.pyplot as plt

def cutout_bubble(dat,fv,tind,xmin,xmax,w):
    win = planck_taper(dat.shape[0],xmin,xmax,w)
    dd = np.empty(dat.shape)
    dd[:,0] = win*(dat[:,0]-fv)+fv; dd[:,1] = win*dat[:,1]
    return dd

def planck_taper(sz,xmin,xmax,w):
    ii = np.arange(1,w); jj = w*(1./ii-1./(w-ii))
    win = np.zeros(sz)
    win[xmin:xmax] = 1.
    win[xmin-w+1:xmin] = 1./(1.+np.exp(jj))
    win[xmax:xmax+w-1] = 1./(1.+np.exp(jj[-1::-1]))
    return win

if __name__=="__main__":
    dat = np.genfromtxt('field-base.dat').reshape((-1,1024,6))
