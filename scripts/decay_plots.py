#!/usr/bin/env python

# A collection of useful Python functions for making various plots to
# explore the decays of false vacuum trajectories
#
# This file only contains subroutines for plotting.
# For file reading and data manipulation, please see
# decay_plots.py

import numpy as np
import matplotlib.pyplot as plt

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
