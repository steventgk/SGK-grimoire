import numpy as np

def sn_bins(var, n, order='asc'):
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
    order : str
        Takes values of 'asc' or 'dec' to set direction in which you define
        bins. In 'asc' mode the final bin will have n_i <= n. The first bin in
        'dec' mode will have n_i <= n
    Returns
    -------
    None
    """
    if order=='asc':
        sorted = np.argsort(var)
        bins=np.append(var[sorted[0::n]],var[sorted[-1]])
    elif order=='dec':
        sorted = np.argsort(var)#[::-1]
        bins=np.append(var[sorted[-1]],var[sorted[-n::-n]])
        bins=np.append(bins,var[sorted[0]])[::-1]
    bin_cents = bins[:-1] + np.diff(bins)/2

    return bins, bin_cents
