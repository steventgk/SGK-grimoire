import numpy as np
import pynbody as pb
import scipy
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from scipy import signal
from scipy import stats
from scipy.interpolate import CubicSpline, griddata
from scipy.ndimage import gaussian_filter
from fast_histogram import histogram2d
from scipy import optimize
import matplotlib.colors as colors

__author__ = "Steven Gough-Kelly"
__copyright__ = "Copyright 2022, UCLan Galaxy Dynamics"
__credits__ = ["Steven Gough-Kelly", "Victor P. Debattista", "Stuart R. Anderson"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Steven Gough-Kelly"
__email__ = "sgoughkelly@gmail.com"
__status__ = "Production"

def bar_align(galaxy, rbar, barfrac = 0.5, zlim=0.5, log=False):
    """
    Aligns the bar of pynbody galaxy simulation with the x-axis assuming the
    galaxy disc is already aligned to the XY plane using the inertial tensor.

    Function does not return any values/object. Pynbody functions effect the
    global variable which stores 'galaxy' so rotations within the functions
    are applied to input variable 'galaxy'.

    Parameters
    ----------
    galaxy : pynbody simulation object
        Galaxy object in the XY plane to be aligned.

    rbar : float
        Bar radius in simulation units e.g. kpc.

    barfrac : float
        Fraction of bar length to calculate the inertial tensor within in
        simulation units e.g. kpc.

    zlim : float
        Vertical limit to calculate intertial tensor within in simulation units
        e.g. kpc. Useful in galaxies with thick discs and weak bars.

    log : Bool
        Flag to output print statements.

    Returns
    -------
    None

    """
    if np.isnan(rbar):
        if log:
            print('* Bar undefined, using 1 kpc *')
        rbar = 1.0
    elif rbar*barfrac < 1.:
        rbar = 1
        if log:
            print('* Short Bar, using 1 kpc *')
    else:
        rbar = rbar*barfrac
        if log:
            print('* Bar defined, aligning to {} kpc *'.format(rbar))

    if log:
        print('* Realigning bar using |z| < {} *'.format(zlim))

    zfilt = pb.filt.LowPass('z',zlim)&pb.filt.HighPass('z',-zlim)
    rfilt = pb.filt.LowPass('rxy',rbar)

    x = np.array(galaxy[zfilt&rfilt].star['pos'].in_units('kpc'))[:,0]
    y = np.array(galaxy[zfilt&rfilt].star['pos'].in_units('kpc'))[:,1]
    m = np.array(galaxy.star[zfilt&rfilt]['mass'])

    #Calculate the inertia tensor
    I_yy, I_xx, I_xy = np.sum(m*y**2),np.sum(m*x**2),np.sum(m*x*y)
    I = np.array([[I_yy, -I_xy], [-I_xy, I_xx]])

    #Calculate the eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(I)
    lowest = eigenvalues.argmin()
    maj_axis = eigenvectors[:, lowest]

    #Get the angle we need to rotate by
    r_angle = np.degrees(np.arctan2(maj_axis[1], maj_axis[0]))

    galaxy.rotate_z(-r_angle)

    if log:
        print('* Bar realigned by {} degrees*'.format(r_angle))

    return None

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

    values, base = np.histogram(arr, bins=bins,range=rng)
    cumulative = np.cumsum(values)
    if norm:
        cumulative = cumulative/max(cumulative)
    x = base[:-1] + np.diff(base)[0]

    return x, cumulative

def distmod(app_mag=None,abs_mag=None,dist=None,extinc=None,output=None):
    """
    Conversions through distance modulus equation for apparent, absolute
    and distance parameters with optional reddening.

    Parameters
    ----------
    app_mag : float or ndarray of float
        apparent magnitude(s) of source(s) in mag.

    abs_mag : float or ndarray of float
        absolute magnitude(s) of source(s) in mag.

    dist : float or ndarray of float
        distances(s) of source(s) in kpc.

    extinc : float or ndarray of float
        extinction value(s) in relevant band of source(s) in mag.

    output : str
        changes appropriate output. Also changes required inputs.
        can take values:
            * 'abs'    : calculate the absolute magnitude.
                         Requires app_mag and dist.
            * 'app'    : calculate the apparent magnitude.
                         Requires abs_mag and dist.
            * 'dist'   : calculate the distance.
                         Requires app_mag and dist.
            * 'extinc' : calculate the extinction.
                         Requires app_mag, abs_mag and dist.

    Returns
    -------
    calculation: float or ndarray of float
        dependent on the output parameter.
    """

    def app2abs(app_mag=None,dist=None,extinc=None):
        if type(extinc) is type(None):
            extinc = np.zeros(len(app_mag))
        return app_mag - extinc - 5*(np.log10(dist)+2)

    def abs2app(abs_mag=None,dist=None,extinc=None):
        if type(extinc) is type(None):
            extinc = np.zeros(len(abs_mag))
        return abs_mag + extinc + 5*(np.log10(dist)+2)

    def mag2dist(app_mag=None,abs_mag=None,extinc=None):
        if type(extinc) is type(None):
            extinc = np.zeros(len(app_mag))
        return 10**(((app_mag - abs_mag - extinc)/5)-2)

    def magdist2extinc(app_mag=None,abs_mag=None,dist=None):
        return app_mag - 5*(np.log10(dist)+2) - abs_mag

    if output == 'abs':
            return app2abs(app_mag=app_mag,dist=dist,extinc=extinc)
    elif output=='app':
            return abs2app(abs_mag=abs_mag,dist=dist,extinc=extinc)
    elif output=='dist':
            return mag2dist(app_mag=app_mag,abs_mag=abs_mag,extinc=extinc)
    elif output=='extinc':
            return magdist2extinc(app_mag=app_mag,abs_mag=abs_mag,dist=dist)

    else:
        return print('output type not found: '+ \
                   'Options "abs","app","dist" and "extinc"')

def Gauss_Hermite(w, n):
    """
    Return the normalised Gauss Hermite function of order n, weights w
    Gerhard MNRAS (1993) 265, 213-230
    Equations 3.1 - 3.7

    Parameters
    ----------
    w : list or numpy arrary of int or float
        Weights of velocities/positions.
    n : int
        nth order of the Gauss Hermite polynomial.
    Returns
    -------
    numpy array of Gauss Hermite of order n for input weights.
    """
    w = np.array(w)
    p = scipy.special.hermite(n, monic=False) #hermite poly1d obj
    norm = np.sqrt((2**(n+1))*np.pi*np.math.factorial(n)) # N_n Eqn 3.1

    return (p(w)/norm) * np.exp( -0.5 * w * w )

def GaussHermiteMoment(v, n):
    """
    Calculate the Gauss Hermite moment of order n for input distribution

    Parameters
    ----------
    v : list or numpy arrary of int or float
        input distribution of velocities/positions.
    n : int
        nth order of the Gauss Hermite polynomial.
    Returns
    -------
    float Gauss Hermite moment of order n for input distribution.
    """
    v = v[np.isfinite(v)] # remove nans&inf
    if len(v) <= 1: # Added SL speed+error catch when used in binned_statistic
        return np.nan

    v_dash = (v - np.mean(v))/np.std(v) # center on 0, norm width to 1sig
    hn = np.sum(Gauss_Hermite(v_dash, n))
    return np.sqrt(4*np.pi) * hn / len(v)

def linint(x,y,intercept,axis,pout=False):
    """
    Finds the intercept point between arbitrary line and a horizontal or vertical slice.
    Note this does not work when a slice intersects the line multiple times.
    Most commonly used with cumulative distribution functions.

    Parameters
    ----------
    x : list or array
        x coordinates of line.
    y : list or array
       y coordinates of line.
    intercept : float or int
        x or y axis position of the slice to calculate the intercept for.
    axis : str
        Takes values of 'x' or 'y' to define which axis the slice is taken from.
        Pass 'x' for vertical slice or 'y' for horizontal slice.
    pout : bool
        Set to true if you want to print the (x,y) coordinate of the intercept point.
    Returns
    -------
    shapely point object of the intercept. Use intercept.x or intercept.y to access the coordinates.
    """
    line_1 = LineString(np.column_stack((x, y)))

    if axis=='x':
        line_2 = LineString(np.column_stack(([intercept,intercept],[min(y),max(y)])))
    if axis=='y':
        line_2 = LineString(np.column_stack(([min(x),max(x)],[intercept,intercept])))

    if pout:
        print((line_1.intersection(line_2).x,line_1.intersection(line_2).y))

    return line_1.intersection(line_2)

def plot_linint(x,y,intercept,axis,plot_axis,pout=False,**plt_kwargs):
    """
    Plot the intercept of an arbitrary line and a slice as defined using linint on a given plot axis.

    Parameters
    ----------
    x : list or array
        x coordinates of line.
    y : list or array
       y coordinates of line.
    intercept : float or int
        x or y axis position of the slice to calculate the intercept for.
    axis : str
        Takes values of 'x' or 'y' to define which axis the slice is taken from.
        Pass 'x' for vertical slice or 'y' for horizontal slice.
    plot_axis : axis object
        The axis object of a matplotlib plot on which you wish to plot the intercept lines.
    pout : bool
        Set to true if you want to print the (x,y) coordinate of the intercept point.

    Returns
    -------
    None
    """
    intercept = linint(x,y,intercept,axis,pout)
    plot_axis.plot([intercept.x,intercept.x],[min(y),intercept.y],**plt_kwargs)
    plot_axis.plot([min(x),intercept.x],[intercept.y,intercept.y],**plt_kwargs)

    return None

def sn_bins(var, n, w=None, cent='avg', order='asc', leftover='join', rtn_xerr=False):
    """
    Defines bins with equal 'n' in each bin. Useful for sparse data. Order
    allows for the control of direction of definition of the bins as either the
    the last or first bin will have less than the target n. If leftover is join
    (default) then this bin gets joined to the previous creating a bin with > n.

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

def interpolate_missing_pixels(
        image: np.ndarray,
        mask: np.ndarray,
        method: str = 'nearest',
        fill_value: float = np.nan
):
    """
    interpolate missing pixels in an image to allow for smoothing 
    and unsharp mask. Adapted from 
    https://stackoverflow.com/questions/37662180/interpolate-missing-values-2d-python
    
    :param image: a 2D image
    :param mask: a 2D boolean image, True indicates missing values
    :param method: interpolation method, one of
        'nearest', 'linear', 'cubic'.
    :param fill_value: which value to use for filling up data outside the
        convex hull of known pixel values.
        Default is 0, Has no effect for 'nearest'.
    :return: the image with missing values interpolated
    """

    h, w = image.shape[:2]
    xx, yy = np.meshgrid(np.arange(w), np.arange(h))

    known_x = xx[~mask]
    known_y = yy[~mask]
    known_v = image[~mask]
    missing_x = xx[mask]
    missing_y = yy[mask]

    interp_values = griddata(
        (known_x, known_y), known_v, (missing_x, missing_y),
        method=method, fill_value=fill_value
    )

    interp_image = image.copy()
    interp_image[missing_y, missing_x] = interp_values

    return interp_image

def fill_holes(
    image: np.ndarray,
    interpolate: bool = True,
    fill_value: float = 0.0
):
    """
    fill nan/inf values in an image using either interpolation or fill value
    
    :param image: a 2D image
    :param mask: a 2D boolean image, True indicates missing values
    :param interpolate: flag to use either interpolation or single value fill
    :param fill_value: (optional) which value to use for filling missing data
                        has no effect when interpolate==True
    :return: the image with missing values filled
    """
    fullimg = np.copy(image)
    
    if interpolate:
        fullimg = interpolate_missing_pixels(fullimg, ~np.isfinite(image))
    elif type(fill_value)==float:
        fullimg[~np.isfinite(image)] = fill_value
    else:
        fullimg = image
    return fullimg

def unsharp(
    image: np.ndarray,
    sigma: int = 0,
    rtn_smooth: bool = False,
    log: bool = True
):
    """
    Produce astronomy unsharp image as fractional difference
    from smoothed image
    
    :param image: a 2D image
    :param sigma: (optional) pixel sigma to smooth over, default is 1/10 
                    of the (w+h)/2
    :param rtn_smooth: flag to also return smoothed image
    :param log: flag to control print statements
    :return: the unsharp (and optional smoothed) images from input
    """  
    if sigma==0:
        sigma = int(np.mean(image.shape)/10)
        if log:
            print('Using sigma of:',sigma)
        
    smoothed = gaussian_filter(image, sigma=sigma)
    unsharp = (image - smoothed)/smoothed
    if rtn_smooth:
        return smoothed, unsharp
    else:
        return unsharp


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