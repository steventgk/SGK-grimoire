import numpy as np
import scipy

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
