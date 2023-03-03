import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata

__author__ = "Steven Gough-Kelly"
__copyright__ = "Copyright 2022, UCLan Galaxy Dynamics"
__credits__ = ["Steven Gough-Kelly"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Steven Gough-Kelly"
__email__ = "sgoughkelly@gmail.com"
__status__ = "Production"

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

# Example
# unsharp(fill_holes(image))