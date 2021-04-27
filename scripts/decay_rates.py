#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

# Convenience routine for fitting analysis
def log_p(tc,ax):
    print(tc[1])
    ax.plot(np.sort(tc[0]),np.log(survive_prob(tc[0],tc[1])))
    return

# Combined returns to combine various steps in the pipeline
def create_times(files,sig_cut=-5):
    d = [read_mean_field(f) for f in files]
    mu = []; sig = []
    for dc in d:
        mu.append(np.mean(dc[:,0,1]))
        sig.append(np.std(dc[:,0,1]))
        myThresh = mu[i]+sig_cut*sig[i]
    return [extract_decay_times(tc[:,:,1],thresh=myThresh,cut=myThresh,interp=True,dt=tc[0,0,0]) for i,tc in enumerate(d) ]

def get_means_and_sigma(d):
    mu = []; sig = []
    for dc in d:
        mu.append(np.mean(dc[:,0,1]))
        sig.append(np.std(dc[:,0,1]))
    return np.array(mu), np.array(sig)

def decays_from_dat(dat,sigCut=-10,full=False):
    mu, sig = get_means_and_sigma(dat)
    th = mu + sigCut*sig
    if full:
        return [extract_decay_times_full(dc[:,:,1],dt=dc[0,0,0],interp=True,thresh=th,cut=th) for dc in d]
    return [extract_decay_times(dc[:,:,1],dt=dc[0,0,0],interp=True,thresh=th,cut=th) for dc in d]
        
# Move this into plotting file
def plot_survival(times,a=None):
    if a == None:
        f,a = plt.subplots()
    a.plot(np.sort(times[0]),survive_prob(times[0],times[1]))
    return a.get_figure(),a

def lin_fit_times(times,tmin,tmax,o=1):
    """
    Given a collection of decay times, do a linear fit to
    the logarithmic survival probability between given times
    
    Input
      times : Times object (first index is array of decay times, 2nd is original number of samples
      tmin  : minimum time to fit inside
      tmax  : maximum time to fit inside
    """
    t = np.sort(times[0])
    p = np.log( survive_prob(times[0],times[1]) )
    ii = np.where( (t>tmin) & (t<tmax) )
    return np.polyfit(t[ii],p[ii],o)

def lin_fit_probs(times,pmin,pmax):
    """
    """
    t = np.sort(times[0])
    p = np.survive_prob(times[0],times[1])
    ii = np.where( (p<pmax) & (p>pmin) )
    return np.polyfit( t[ii],np.log(p[ii]) )

# File reading and manipulation
def read_mean_field(fName):
    """
    Read in time-stream data of mean field from ASCII file.

    Assumed format of the file is a file with 4 columns:
    Time   <phi>  rho_spec  rho_fd

    With trajectories separated by a carriage return.
    
    This function assumes each sample has the same number of time-steps, which is inferred from the time outputs in the file.  If there is a partial trajectory at the end, it is removed.
    """
    tind = 0
    d = np.loadtxt(fName)
    dt = np.diff(d[:,tind])
    nTstep = np.where(dt<0)[0][0]+1
    return d[:d.shape[0]//nTstep*nTstep,:].reshape((-1,nTstep,d.shape[1]))

# start debugging this thing
# Issues:
#   Assumes each files has trajectories with the same number of time-steps
def get_trajectories(files,num_tsteps,ax=1):
    """
    Extract the trajectores from a collection of files.

    Input:
     files      - A list of files to read the trajectories from
     num_tsteps - Number of time steps in each trajectory
     ax         - Axis storing the trajectories
    Output:
     d     - A list of numpy arrays containing the decay trajectories of length (num files).  Each numpy array has shape (num_traj,num_tstep)
    """
    d = []
    for f in files:
        a = np.loadtxt(f)
        d.append( a[:a.shape[0]//num_tsteps*num_tsteps,ax].reshape((-1,num_tsteps)) )
    return d

# Need to debug the interpolation subroutine
def extract_decay_times(time_streams,thresh=0.5,cut=0.5,interp=True,up_cross=False,dt=1.,**kwargs):
    """
    Extract the first crossing of a time-stream past some threshold.

    time_streams - The time streams to get the upcrossings from.  
      Has shape (num_samples,num_tsteps)

    thresh - Crossing threshold

    cut - During the time-stream, only time-streams that exceed cut are considered to have decayed.

    up_cross   - If True look for first upcrossing
               - If False look for first downcrossing

    interp - If True will linearly interpolate to get the time
           - If False will just get the lower bin on the time

    dt     - Time step between outputs

    Returns:
        times    - numpy array of decay times (unsorted)
        num_traj - total number of trajectories 

        Note: Because not all trajectories need decay, times.size may be less than the number of trajectories
    """
    if up_cross:
        td = np.where(np.max(time_streams[:,:],axis=-1) > cut)
        ti = np.argmax(time_streams[td,:] > thresh,axis=-1)[0]  # Assumes increading
    else:
        td = np.where(np.min(time_streams[:,:],axis=-1) < cut)
        ti = np.argmax(time_streams[td,:] < thresh,axis=-1)[0]

    # This needs to be fixed to work when ti = 0, or fucky slope
    if interp:
        t = ti + (thresh-time_streams[td,ti]) / (time_streams[td,ti]-time_streams[td,ti-1])
        t = t[0]
    else:
        t = ti
    return dt*t, time_streams.shape[0]

def extract_decay_times_full(time_streams,thresh=0.5,cut=0.5,interp=True,up_cross=False,dt=1.,**kwargs):
    """
    Extract the decay times, including -1 for undecayed trajectories.
    This allows a direct comparison between different resolutions using
    identical initial conditions.
    """
    t = -np.ones(time_streams.shape[0])/dt  # For future normalisation
    if up_cross:
        td = np.where(np.max(time_streams[:,:],axis=-1) > cut)
        ti = np.argmax(time_streams[td,:] > thresh,axis=-1)[0]
    else:
        td = np.where(np.min(time_streams[:,:],axis=-1) < cut)
        ti = np.argmax(time_streams[td,:] < thresh,axis=-1)[0]

    #### Need to add interpolation
    if interp:
        print("Interpolation not yet fully tested")
        t[td] = ti + (thresh-time_streams[td,ti]) / (time_streams[td,ti]-time_streams[td,ti-1])
    else:
        t[td] = ti
    return dt*t

def decay_indices(t1,t2):
    """
    Given two streams of input decay times, return the indices where both decay, the first decays but the second doesn't, the second decays but the first doesn't, and where neither decays

    Input:
     t1 - First set of decay times
     t2 - Second set of decay times

    Output:
     i_dd - Indices where both decay
     i_du - Indices where first decays but second doesn't
     i_ud - Indices where second decays but first doesn't
     i_dd - Indices where neither decays
    """
    i_dd = np.where( (t1>=0.) & (t2>=0.) )
    i_du = np.where( (t1>=0.) & (t2<0.) )
    i_ud = np.where( (t1<0.) & (t2>=0.) )
    i_uu = np.where( (t1<0.) & (t2<0.) )        
    return i_dd, i_du, i_ud, i_uu

def ind_count(t0,t1,ind,thresh):
    """
    Input: 
      ind    - The indices where both ensembles decay, one does and the other doesn't, and neither do.  In the form of the output of decay_indices
      thresh - Threshold dt for which two decaying trajectories are determined to not match

    Returns:
      n - Counts of the 4 indices
      n_error - Number of misatched trajectories that decayed in both sims
    """
    n = np.array([ic[0].size for ic in ind])
    n_error = np.where(np.abs(t0[ind[0]]-t1[ind[0]]) > thresh)[0].size
    return n, n_error

def error_prob(t_bad,t_res,thresh=0.5):
    """
    Return the fraction of 'bad' decay times given an unresolved and compared resolved simulation

    Input:
      t_bad : Decay times (including undecayed) times for unresolved sims
      t_res : Decay times for resolved simulations
      thresh : Error in dt to be considered a bad decay time

    Output:
      Fraction of bad decay times defined as:
        (number_wrong_decays + number_decays_in_unresolved + number_decays_in_resolved) / (number_decays + number_decays_in_unresolved + number_decays_in_resolved)
    """
    ind = decay_indices(t_bad,t_res)
    n,n_error = ind_count(t_bad,t_res,ind,thresh)
    return (n_error+n[1]+n[2])*1./(n[0]+n[1]+n[2])

# Return the distribution of decay trajectories
def decay_frac_dist(times):
    uv_ax = 1; ir_ax = 0
    nIR = times.shape[ir_ax]; nUV = times.shape[uv_ax]
    nDecay = np.sum(times >= 0.,axis=uv_ax) # Check the axes are correct

    decayDist = np.array( [np.sum(nDecay==i) for i in range(nUV+1)] )
    return nDecay, decayDist

def vary_uv_stats(times):
    # First, extract only the trajectories that all decay
    nDecay = np.sum(times >= 0., axis=1)
    ii_decay = np.where(nDecay==nUV)
    ii_survive = np.where(nDecay==0)
    # Now get variance, etc. for the decayed trajectories
    s = np.std(times[ii_decay],axis=-1)
    d = np.max(times[ii_decay],axis=-1)-np.min(times[ii_decay],axis=-1)
    return s,d

def bin_decay_times():
    """
    When this is written, I will slice and dice the decay times for identical
    choices of simulations
    """
    return

# To do: Debug more to ensure all offsets are correct.
# I've done a first go through and I think they're ok
def survive_prob(t_decay,num_samp):
    """
    Return the survival probability as a function of time.

    Input:
      t_decay  : Decay times of trajectories
      num_samp : Total number of samples in Monte Carlo

    Note: Since some trajectories may not decay, t_decay.size isn't always num_sampe

    Output:
      prob     : Survival Probability

    These can be plotted as plt.plot(t_sort,prob) to get a survival probability graph.
    """
    frac_remain = float(num_samp-t_decay.size)/float(num_samp)
    prob = 1.-np.linspace(1./num_samp,1.-frac_remain,t_decay.size,endpoint=True)
    return prob
    
if __name__=="__main__":
    # Temporary to reduce typing
#    fName = 'files-vary-cut.txt'
#    fName = 'files-hertzberg-varyL.txt'
#    fName = 'files-cut-pairs.txt'
#    fName = 'files-cut-base-to-converged.txt'
#    fName = 'files-check-converged-pairs.txt'
    fName = 'files-vary-kir-fix-nrat.txt'
    
    with open(fName) as f:
        files = f.read().splitlines()
        
    d = [read_mean_field(f) for f in files]
    mu,sig = get_means_and_sigma(d)
    nsig = -10; th = mu+nsig*sig

    tf = [ extract_decay_times_full(dc[:,:,1],dt=dc[0,0,0],thresh=th[i],cut=th[i]) for i,dc in enumerate(d) ]
#    t = [ extract_decay_times(dc[:,:,1],dt=dc[0,0,0],thresh=th[i],cut=th[i]) for i,dc in enumerate(d) ]
#    for i in range(3):
#        print( error_prob(tf[2*i],tf[2*i+1]) )
    
#    fName = 'files-cut-base-to-converged.txt'
#    with open(fName) as f:
#        files = f.read().splitlines()
        
#    d = [read_mean_field(f) for f in files]
#    mu,sig = get_means_and_sigma(d)
#    for i in range(4):
#        mu[2*i] = mu[2*i+1]
#        sig[2*i] = sig[2*i+1]
#    nsig = -10; th = mu+nsig*sig
        
#    tf = [ extract_decay_times_full(dc[:,:,1],dt=dc[0,1,0]-dc[0,0,0],thresh=th[i],cut=th[i]) for i,dc in enumerate(d) ]
#    t = [ extract_decay_times(dc[:,:,1],dt=dc[0,1,0]-dc[0,0,0],thresh=th[i],cut=th[i]) for i,dc in enumerate(d) ]
#    for i in range(4):
#        print( error_prob(tf[2*i],tf[2*i+1]) )
    
