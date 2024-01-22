from fast_histogram import histogram2d
import matplotlib.pyplot as plt
import matplotlib.colors as colors

__author__ = "Steven Gough-Kelly"
__copyright__ = "Copyright 2024, UCLan Galaxy Dynamics"
__credits__ = ["Steven Gough-Kelly"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Steven Gough-Kelly"
__email__ = "sgoughkelly@gmail.com"
__status__ = "Production"

# Solution from Paul Gavrikov with minor edits to make it more flexible
# https://towardsdatascience.com/how-to-create-fast-and-accurate-scatter-plots-with-lots-of-data-in-python-a1d3f578e551
# Implementation only has square bins for speed
# uses https://pypi.org/project/fast-histogram/


def quick_hist(x, y, nbin=1000, cmap='magma', xlim=None, ylim=None, vmin=1, vmax=None):
    """
    Produces a simple 2D histogram with colourbar to quickly sutdy the distribution of two variables.

    Slightly faster than an numpy implementation and significantly faster than Scipy.

    Parameters
    ----------
    x : list or array
        x variable

    y : list or array
        y variable

    nbin : int
        The number of bins in BOTH x and y.

    cmap : string
        Matplotlib colormap to display.
    
    xlim : array like shape (2,1)
        The x limits of the plot. Default None will show full range of x.

    ylim : array like shape (2,1)
        The y limits of the plot. Default None will show full range of y.

    vmin : float
        The limit of minium intensity of the plot. Default will impose a lower limit of 1.

    vmax : float
        The limit of maximum intensity of the plot. Default None will be
        largest bin count of the histogram.


    Returns
    -------
    None

    """

    bounds = [[x.min(), x.max()], [y.min(), y.max()]]
    extent = [x.min(), x.max(), y.min(), y.max()]
    h = histogram2d(x, y, range=bounds, bins=nbin)
    
    if vmax is None:
        vmax = h.max()
    
    plt.imshow(h.T, extent=extent, origin='lower',
               norm=colors.LogNorm(vmin=vmin, vmax=vmax),
               cmap=cmap)
    
    plt.xlabel('x')
    plt.ylabel('y')
    
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    
    plt.colorbar(label='N')
    plt.show()