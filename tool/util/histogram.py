import numpy as np
import scipy.stats as st

"""
    Returns 2 arrays (length : nbin) containing bin centers 
    and bin counts (depends on 'norm' flag)
"""
def histogram_1d(x, nbin = 1, min = None, max = None, norm = False) :
    # Re-assign range if required
    if not min :
        min = np.amin(x)
    if not max :
        max = np.amax(x)

    bin_count, bin_edge = np.histogram(x, bins = nbin, 
                  range = (min, max), density = norm)

    bin_center = 0.5*(bin_edge[:-1] + bin_edge[1:])

    return bin_center.flatten(), bin_count.flatten()

"""
    Returns 3 arrays (length : nbin) containing bin centers,
    mean and errors for each bin
"""
def histogram_profile(x, y, nbin = 1, min = None, max = None) :
    # Re-assign range if required
    if not min :
        min = np.amin(x)
    if not max :
        max = np.amax(x)

    y2 = [i**2 for i in y]
    bin_stat = st.binned_statistic(x, [y, y2], bins = nbin,
                    range = (min, max), statistic = 'mean')
    
    bin_mean, bin_rms = bin_stat.statistic
    bin_center, bin_count = histogram_1d(x, nbin)

    # Calculate standard deviation for each bin
    bin_err = np.sqrt((bin_rms - bin_mean**2)/bin_count)
    
    return bin_center, bin_mean.flatten(), bin_err.flatten()