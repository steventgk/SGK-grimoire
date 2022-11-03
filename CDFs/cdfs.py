import numpy as np

__author__ = "Steven Gough-Kelly"
__copyright__ = "Copyright 2022, UCLan Galaxy Dynamics"
__credits__ = ["Steven Gough-Kelly"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Steven Gough-Kelly"
__email__ = "sgoughkelly@gmail.com"
__status__ = "Production"

def CDF(arr,bins=10_000,rng=None, norm=True):
    """
    Produces cumulative distribution functions of input array.
    Large number of bins produce curves which step up for each element where
    len(arr) < bins. For len(arr) > bins, default number of bins is still large
    enough to produce smooth lines.

    Parameters
    ----------
    arr : list or numpy arrary of int or float
        Variable to be binned.
    bins : int
        Number of bins
    rng : tuple or list of two elements
        upper and lower limits over which to calculate the CDF
    norm : bool
        Flag to determine if CDF ranges from 0-max(cumsum) or 0-1
    Returns
    -------
    x : the bin centres
    cumulative : the cumulative sum values
    """

    if type(rng)==type(None):
        rng = (np.nanmin(arr),np.nanmax(arr))

    values, base = np.histogram(arr, bins=10000,range=rng)
    cumulative = np.cumsum(values)
    if norm:
        cumulative = cumulative/max(cumulative)
    x = base[:-1] + np.diff(base)[0]

    return x, cumulative
