import numpy as np
from scipy import stats
from scipy import optimize

__author__ = "Steven Gough-Kelly"
__copyright__ = "Copyright 2022, UCLan Galaxy Dynamics"
__credits__ = ["Steven Gough-Kelly"]
__license__ = "MIT"
__version__ = "1.3.1"
__maintainer__ = "Steven Gough-Kelly"
__email__ = "sgoughkelly@gmail.com"
__status__ = "Production"

def sn_bins(var, n, w=None, cent='avg', order='asc', leftover='join', rtn_xerr=False):
    """
    Defines bin edges with equal 'n' in each bin (i.e. where you can assume 
    noise = sqrt(n)). Useful for none uniform distribution of data. Order 
    controls the direction bin definition as either the last or first bin 
    will have less than the target n. Leftover controls the handling of 
    these remaining elements. Function also returns bin centers.
    This allows for quick easy integration into binned_statistic and plotting.

    Parameters
    ----------
    var : list or numpy arrary of int or float
        Variable to be binned.
    n : int
        Number of elements to put in a bin
    w : list or numpy arrary of int or float
        Weight of each element in var. Used only when cent=='avg'
    cent : str
        Takes values of 'avg' or 'mid' to chose definition of bin centers.
        Either weighted average of each bin or midpoint between bins.
    order : str
        Takes values of 'asc' or 'des' to set direction in which you define
        bins. In 'asc' mode the final bin will have n_i <= n. The first bin in
        'des' mode will have n_i <= n.
        If 'ValueError: The smallest edge difference is numerically 0.' is
        raised in binned_statistic try changing order.
    leftover : str
        Takes value of 'join' to join final bin which has len(elements) < n
        to the penultimate bin to create a combined bin with len(elements) > n
        or 'drop' to exclude the bin with len(elements) < n. Any other value
        will default to leaving the final bin with len(elements) < n in 
        the output.
    rtn_xerr : bool
        Flag to return the uncertainty on the bin_cents when using cent='avg'.
    Returns
    -------
    bin_edges : numpy array
        the bin edges
    bin_cents : numpy array
        the bin centers
    bin_xerr : numpy array
        (optional) the uncertaintiy on the bin center average
    """

    if order not in ['asc','des']:
        order = 'asc'
        raise Exception("Use only 'asc' or 'des' values for order. Assuming 'asc'.")

    # Create bin arrays

    # First Catch Case - perfectly divisible
    # defaults order as would be equivient
    # sets leftover to None
    if len(var)%n == 0:
        order = 'asc'
        leftover = None

        sorted = np.argsort(var)
        bin_edges = var[sorted[::n]]
        bin_edges = np.append(bin_edges,var[sorted[-1]])
    else:
        # Second Catch Case - edge difference 0
        # Can't have two bin edges that are the same value
        # Has to default to 'join'
        if len(var)%n == 1:
            leftover = 'join'

        #find length of remainder
        l = len(var)%n

        if order=='asc':
            sorted = np.argsort(var)
            tmp = sorted[:-l]
            bin_edges = var[tmp[::n]]
            if leftover=='join':
                bin_edges = np.append(bin_edges,var[sorted[-1]])
            elif leftover=='drop':
                bin_edges = np.append(bin_edges,var[sorted[-l-1]])
            else:
                bin_edges = var[sorted[::n]]
                bin_edges = np.append(bin_edges,var[sorted[-1]])

        elif order=='des':
            sorted = np.argsort(var)
            tmp = sorted[l:]
            bin_edges = var[tmp[::n]]
            if leftover=='join':
                bin_edges = np.append(var[sorted[0]],bin_edges[1:])
                bin_edges = np.append(bin_edges, var[sorted[-1]])
            elif leftover=='drop':
                bin_edges = np.append(bin_edges, var[sorted[-1]])
            else:
                bin_edges = np.append(var[sorted[0]],bin_edges)
                bin_edges = np.append(bin_edges, var[sorted[-1]])
    
    # Calculate Bin Centers

    if cent not in ['avg','mid']:
        cent = 'avg'
        raise Exception("Use only 'avg' or 'mid' values for bin centers. Assuming 'avg'.")

    if (cent=='avg'):
        if w is None:
            w = np.ones(len(var))
        bin_cents, _, _ = stats.binned_statistic(var,w*var,statistic='sum',bins=bin_edges)
        counts, _, _ = stats.binned_statistic(var,w,statistic='sum',bins=bin_edges)
        bin_cents = bin_cents/counts
        
        if rtn_xerr:
            bin_xsig, _, _ = stats.binned_statistic(var,var,statistic='std',bins=bin_edges)
            nbin, _, _ = stats.binned_statistic(var,None,statistic='count',bins=bin_edges)
            bin_xerr = bin_xsig/np.sqrt(nbin)
            return bin_edges, bin_cents, bin_xerr

    elif (cent=='mid'):
        bin_cents = bin_edges[:-1] + np.diff(bin_edges)/2

    return bin_edges, bin_cents

def lin_sn_bins(var,nmin,range=None,mode='min'):
    """
    Defines linear bins with a minimum 'nmin' elements in each bin using optimize.
    Two modes allow for high resolution (max number of bins) or minimum number of
    bins possible (default).

    Parameters
    ----------
    var : list or numpy arrary of int or float
        Variable to be binned.
    nmin : int
        Minimum number of elements required for each bin
    range : tuple or list of two elements.
        The Min and Max values to determine bins over.
    mode : str
        Takes values of min or max to allow for the minimum number of bins
        that fits the nmin criteria or the maximum possible to still ensure
        nmin in each bin.
    Returns
    -------
    nbins : min/max number of bins possible to ensure each bin has nmin elements.
    """
    if range is None:
        range = (np.min(var),np.max(var))
        
    def sn_min(nbins,var,nmin,range):
        count, edges = np.histogram(var,range=range,bins=int(nbins))
        if np.nanmin(count)>nmin:
            return (1/(np.diff(edges)[0]))
        elif np.nanmin(count)==nmin:
            return -1e-10
        else:
            return -(1/(np.diff(edges)[0]))
    
    def sn_max(nbins,var,nmin,range):
        count, edges = np.histogram(var,range=range,bins=int(nbins))
        if np.nanmin(count)<nmin:
            return (1/(np.diff(edges)[0]))
        elif np.nanmin(count)==nmin:
            return -1e-10
        else:
            return -(1/(np.diff(edges)[0]))
    
    if mode=='min':
        sol = optimize.root_scalar(sn_min, args=(var,nmin,range),
                                    bracket=[1, len(var)], method='brentq')
    if mode=='max':
        sol = optimize.root_scalar(sn_max, args=(var,nmin,range),
                                    bracket=[1, len(var)], method='brentq')
        
    return int(sol.root)