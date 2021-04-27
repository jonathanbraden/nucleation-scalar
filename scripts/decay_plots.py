#!/usr/bin/env python

# A collection of useful Python functions for making various plots to
# explore the decays of false vacuum trajectories
#
# This file only contains subroutines for plotting.
# For file reading and data manipulation, please see
# decay_plots.py

import numpy as np
import matplotlib.pyplot as plt
from decay_rates import *

def _decay_labels(ax):
    """
    Add time label and decay probability labels
    """
    ax.set_xlabel(r'$\bar{t}$'); ax.set_ylabel(r'$P_{\rm decay}$')
    return

def _survive_labels(ax):
    """
    Add time label and survival probability label
    """
    ax.set_xlabel(r'$\bar{t}$')
    ax.set_ylabel(r'$P_{\rm survive}$')
    return

def decay_comparison_scatter(t_coarse,t_fine):
    i_dd = np.where( (t_coarse>=0.) & (t_fine>=0.) )
    i_du = np.where( (t_coarse>=0.) & (t_fine<0.) )
    i_ud = np.where( (t_coarse<0.) & (t_fine>=0.) )
    i_uu = np.where( (t_coarse<0.) & (t_fine<0.) )
    
    f,a = plt.subplots()
    a.scatter(t1[i_dd],t2[i_dd],alpha=0.2)
    a.scatter(t1[i_du],t2[i_du],'rx',alpha=0.2)
    a.scatter(t1[i_ud],t2[i_ud],'kx',alpha=0.2)
    return f,a

if __name__=="__main__":
    pass
