#!/usr/bin/env python
import numpy as np

# Combined returns to combine various steps in the pipeline
def create_times(files):
    d = [read_mean_field(f) for f in files]
    return [extract_decay_times(tc[:,:,1],thresh=0.5,cut=0.5,interp=True,dt=tc[0,0,0]) for tc in d]

# Move this into plotting file
def plot_survival(times,f=None,a=None):
    if f == None:
        f,a = plt.subplots()
    a.plot(np.sort(times[0]),survive_prob(times[0],times[1]))
    return f,a

def lin_fit_times(times,tmin,tmax):
    """
    Given a collection of decay times, do a linear fit to
    the logarithmic survival probability between given times
    
    Input
      times : Times object (first index is array of decay times, 2nd is original number of samples
      tmin  : minimum time to fit inside
      tmax  : maximum time to fit inside
    """
    t = np.sort(times[0])
    p = survive_prob(times[0],times[1])
    ii = np.where( (t>tmin) & (t<tmax) )
    return np.polyfit(t[ii],p[ii])

def lin_fit_probs(times,pmin,pmax):
    """
    """
    t = np.sort(times[0])
    p = survive_prob(times[0],times[1])
    ii = np.where( (p<pmax) & (p>pmin) )
    return np.polyfit(t[ii],p[ii])

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
    return d[:d.shape[0]//nTstep*nTstep,:].reshape((-1,nTstep,4))

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
    pass
