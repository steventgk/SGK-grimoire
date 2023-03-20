import numpy as np
from scipy import stats
from scipy import optimize

__author__ = "Steven Gough-Kelly"
__copyright__ = "Copyright 2022, UCLan Galaxy Dynamics"
__credits__ = ["Steven Gough-Kelly"]
__license__ = "MIT"
__version__ = "1.2.1"
__maintainer__ = "Steven Gough-Kelly"
__email__ = "sgoughkelly@gmail.com"
__status__ = "Production"

def sn_bins(var, n, w=None, cent='avg', order='asc', leftover='join'):
    """
    Defines bins with equal 'n' in each bin. Useful for sparse data. Order
    allows for the control of direction of definition of the bins as either the
    the last or first bin will have less than the target n. If leftover is join
    then this bin gets joined to the previous creating a bin with > n.

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
        Takes values of 'asc' or 'dec' to set direction in which you define
        bins. In 'asc' mode the final bin will have n_i <= n. The first bin in
        'dec' mode will have n_i <= n.
        If 'ValueError: The smallest edge difference is numerically 0.' is
        raised in binned_statistic try changing order.
    leftover: str
        Takes value of 'join' to join final bin which has len(elements) < n.
    Returns
    -------
    bins : the bin edges
    bin_cents : the bin centers
    """
    if order=='asc':
        sorted = np.argsort(var)
        bins = var[sorted[::n]]
        if leftover=='join':
            bins = np.append(bins[:-1],var[sorted[-1]])
        else:
            bins = np.append(bins,var[sorted[-1]])
    elif order=='dec':
        sorted = np.argsort(var)[::-1]
        bins = var[sorted[::n]]
        if leftover=='join':
            bins = np.append(bins[:-1],var[sorted[-1]])[::-1]
        else:
            bins = np.append(bins,var[sorted[-1]])[::-1]
        
    if cent=='avg':
        if w is None:
            w = np.ones(len(var))
        bin_cents, _, _ = stats.binned_statistic(var,w*var,statistic='sum',bins=bins)
        counts, _, _ = stats.binned_statistic(var,w,statistic='sum',bins=bins)
        bin_cents = bin_cents/counts
    if cent=='mid':
        bin_cents = bins[:-1] + np.diff(bins)/2

    return bins, bin_cents

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
        count, edges = np.histogram(var,range=arange,bins=int(nbins))
        if np.nanmin(count)>nmin:
            return (1/(np.diff(edges)[0]))
        elif np.nanmin(count)==nmin:
            return -1e-10
        else:
            return -(1/(np.diff(edges)[0]))
    
    def sn_max(nbins,var,nmin,range):
        count, edges = np.histogram(var,range=arange,bins=int(nbins))
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