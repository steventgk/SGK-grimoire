import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import LineString

__author__ = "Steven Gough-Kelly"
__copyright__ = "Copyright 2022, UCLan Galaxy Dynamics"
__credits__ = ["Steven Gough-Kelly"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Steven Gough-Kelly"
__email__ = "sgoughkelly@gmail.com"
__status__ = "Production"

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