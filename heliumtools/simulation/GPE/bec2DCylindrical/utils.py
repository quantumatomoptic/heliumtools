import numpy as np
from scipy.interpolate import make_interp_spline


def derivative_rui(func,x,order,extrapolate,k):
    """ Computes the 5-point derivative of the function 
    ---------------------------------------------------------
    Parameters
    func : np.array
        function to compute derivative
    x : np.array
        points at which function is evaluated
    order : int
        order of derivative to compute. Currently one 1st and 2nd orders are implemented
    extrapolate : bool
        if True, it extrapolates the points at the edge, by interpolating the derivative. 
        If false, it computes the near-edge points using the 3-point derivative and the edge points using a second order approximation.
    k : int
        order of interpolation for when extrapolate is set to True.
    -----------------------------------------------------------------
    Return
        The derivative of func at order order.
    """

    order = int(order)
    # compute step
    step = x[1] - x[0]
    # pad array
    Npad = 2
    pad = np.pad(func,(Npad,Npad), mode='constant' , constant_values = (func[0],func[-1]))
    N = len(func)
    # compute derivative 1st order
    if order == 1:
        # 5-point derivative
        first = (pad[Npad-2:Npad+N-2] - 8*pad[Npad-1:Npad+N-1] + 8*pad[Npad+1:Npad+N+1] - pad[Npad+2:Npad+N+2])/12
        # if true, we extrapolate the edges
        if extrapolate:
            first = first/step
            first = extrapolate_derivative(x,first,k)
            return first
        else:
            # compute edges using second order approximation
            first[0] = -3/2*func[0] + 2*func[1] - 0.5*func[2]
            first[-1] = 3/2*func[-1] - 2*func[-2] + 0.5*func[-3]
            # compute near-edge using 3-point derivative
            first[1] = (func[2] - func[0])/2
            first[-2] = (func[-1] - func[-3])/2
            return first/step
    # compute second order derivative
    elif order == 2:
        # 5-point derivative
        second = (-pad[Npad-2:Npad+N-2] + 16*pad[Npad-1:Npad+N-1] - 30*pad[Npad:Npad+N] + 16*pad[Npad+1:Npad+N+1] - pad[Npad+2:Npad+N+2])/12
        # if true, we extrapolate the edges
        if extrapolate:
            second = second/step**2
            second = extrapolate_derivative(x,second,k)
            return second
        else:
            # compute edges using 2-order approximation
            second[0] = 0.5*func[2] - func[1] + 0.5*func[0]
            second[-1] = 0.5*func[-1] - func[-2] + 0.5*func[-3]
            # compute near-edge using 3-point derivative
            second[1] = func[2] - 2*func[1] + func[0]
            second[-2] = func[-1] - 2*func[-2] + func[-3]
            return second/step**2
    else:
        print("Invalid order!")
        return np.zeros(len(x),dtype=float)   

def extrapolate_derivative(x,func,k):
    """ Given the array func, it extrapolates the edges, by interpolating the rest of the array 
    ------------------------------------------------------------------------------
    Parameters
    x : np.array
        point coordinates
    func : np.array
        array with function values at x
    k : int
        order of interpolation
    ---------------------------------------------------------------------------------
    Return
        func with edges extrapolated
    """
    N = len(func)
    # interpolate func points, except the edges
    spline = make_interp_spline(x[2:N-2], func[2:N-2],k = k)
    # extrapolate the edges and update func
    func[0:2] = spline(x[0:2])
    func[N-2:N] = spline(x[N-2:N])
    return func

def laplacian(func,x):
    """ Computes the radial part of the laplacian. It uses the 5-point derivative
    and it extrapolates the values at the edge using a 9th order interpolator.
    ----------------------------------------------------------------------------
    Parameters
    func : np.array
        array with function values
    x : np.array
        values of coordinates at which func is evaluated
    ----------------------------------------------------------------
    Return 
        radial part of the laplacian of func
    """

    """ compute first derivative """
    first = derivative_rui(func,x,1,False,9)

    """ compute second derivative """
    second = derivative_rui(func,x,2,False,9)

    """ compute laplacian """
    lap = second + first/x

    """ extrapolate edges """
    lap = extrapolate_derivative(x,lap,9)

    return lap

