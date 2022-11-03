import numpy as np
from scipy import signal
from scipy.interpolate import CubicSpline

__author__ = "Steven Gough-Kelly"
__copyright__ = "Copyright 2022, UCLan Galaxy Dynamics"
__credits__ = ["Steven Gough-Kelly", "Stuart R. Anderson"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Steven Gough-Kelly"
__email__ = "sgoughkelly@gmail.com"
__status__ = "Production"

def buttersmooth(x, y, order=2, crit=None, interp=True, cs_np=1000,log=True):

    """
    A lowpass butterwooth frequency filter and optional CubicSpline
    interpolation to smooth a 1D arr input.

    Parameters
    ----------
    x : list or numpy arrary of int or float
        x coordinates of arr
    y : list or numpy arrary of int or float
        y coordinates of arr
    order : int
        Order of butterworth filter which controls the amplitude of suppresion
        of higher frequencies than crit frequency.
    crit : float or None
        If None the function calculates the critical frequency estimate by
        calculating a window size of len(x)/10. Can be used as a 'good first
        guess'.
    interp : bool
        Option to implement CubicSpline interpolation.
    cs_np : int
        The number of data points to interpolate the arr to. Does nothing if
        interp=False
    log : bool
        Option to print calculated critical frequency if crit=None.
    Returns
    -------
    xout : smoothed x coordinates with len(xout)==len(x) if interp=False,
           len(xout)==len(cs_np) if interp=True.
    yout : smoothed y coordinates with len(yout)==len(x) if interp=False,
           len(yout)==len(cs_np) if interp=True.
    """

    if type(crit)==type(None):
        crit = np.sum(np.diff(x)[:int(len(x)/10)])
        if log:
            print('Using Critical frequency of {}'.format(crit))

    b, a = signal.butter(order, crit, 'lowpass')

    ys = signal.filtfilt(b, a, y)

    if interp:
        x_cs = np.linspace(np.min(x), np.max(x), cs_np)
        cs = CubicSpline(x, ys)
        y_cs = cs(x_cs)

        xout, yout = x_cs, y_cs
    else:
        xout, yout = x, ys

    return xout, yout
