import numpy as np
from scipy import stats

__author__ = "Steven Gough-Kelly"
__copyright__ = "Copyright 2022, UCLan Galaxy Dynamics"
__credits__ = ["Steven Gough-Kelly"]
__license__ = "MIT"
__version__ = "1.1.1"
__maintainer__ = "Steven Gough-Kelly"
__email__ = "sgoughkelly@gmail.com"
__status__ = "Production"

def sn_bins(var, n, w=None, cent='avg', order='asc', leftover='join'):
    """
    Defines bins with equal 'n' in each bin. Useful for sparse data. Order
    allows for the control of direction of definition of the bins as either the
    the last or first bin will have less than the target n.

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
