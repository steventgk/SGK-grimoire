import numpy as np

__author__ = "Steven Gough-Kelly"
__copyright__ = "Copyright 2022, UCLan Galaxy Dynamics"
__credits__ = ["Steven Gough-Kelly"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Steven Gough-Kelly"
__email__ = "sgoughkelly@gmail.com"
__status__ = "Production"

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
