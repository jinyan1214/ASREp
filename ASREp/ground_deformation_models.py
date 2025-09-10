import numpy as np
from scipy.stats import norm

def ground_disp_Mair_1993(x, z0, K, D, vl):
    """
    Calculate the 2D Gaussian based ground deformation. 
    The vertical displacement comes from Mair, R. J., Taylor, R. N., & Bracegirdle,
    A. (1993). Subsurface settlement profiles above tunnels in clays. Geotechnique,
    43(2), 315-320.
    The horizontal displacement comes from O'reilly, M. P., & New, B. M. (1982).
    Settlements above tunnels in the United Kingdom-their magnitude and prediction
    (No. Monograph).

    Parameters
    ----------
    x : float or np.array(dtype=float)
        Distance from tunnel center. Unit: m
    z0 : float or np.array(dtype=float)
        Tunnel center depth. Unit: m
    K : float or np.array(dtype=float)
        Trough width parameter. Unit: Unitless
    D : float
        Tunnel diameter. Unit: m
    vl : float or np.array(dtype=float)
        Tunnel volume loss. Unit: the actual quantity, not percentages

    Returns
    -------
    tuple
        A tuple containing:
        - The vertical displacement (float  or np.array(dtype=float)) in the unit
          of m and upward as positive value.
        - The horizontal displacement (float  or np.array(dtype=float)) in the unit
          of m and positive direction is the same with x.
    """
    # Eq. 3 of Mair
    i = K*z0
    # Eq. 6a of Mair
    s_max = 0.313 * vl * D**2 / i
    # EQ. 1 of Mair
    dv = - s_max * np.exp(-np.power(x, 2)/2/np.power(i,2))
    # Eq. 6 of O'Reilly
    dh = dv * x / z0
    # Return the vertical and horizontal disp
    return dv, dh


def ground_disp_Zhao_2023(x, y, z, vl, d, z0, ys, yf, k, delta):
    """
    Calculate the 3D Gaussian based ground deformation used in Zhao and DeJong (2023). 
    This is a combination of Peck 1969, Attewell 1982, Mair 1993 and Camos et al. 2014,2015,2016.

    Parameters
    ----------
    x : float or np.array(dtype=float)
        Coordinate in the x direction. Unit: m
    y : float or np.array(dtype=float)
        Coordinate in the y direction. Unit: m
    z : float or np.array(dtype=float)
        Coordinate in the z direction. Unit: m
    vl : float or np.array(dtype=float)
        Tunnel volume loss. Unit: the actual quantity, not percentages
    d : float
        Tunnel diameter. Unit: m
    z0 : float
        Tunnel center depth. Unit: m
    ys : float
        Horizontal distance from the tunnel face to the origin.
        Unit: m
    yf : float
        Horizontal distance from the tunnel start to the origin.
        Unit: m
    k : float
        Trough width parameter. Unit: Unitless
    delta : float
        Vertical displacement factor. Unit: Unitless

    Returns
    -------
    tuple
        A tuple containing:
        - The vertical displacement (float or np.array(dtype=float)) in the unit of m and upward as positive value.
        - The horizontal displacement in the x direction (float or np.array(dtype=float)) in the unit of m and positive direction is the same with x.
        - The horizontal displacement in the y direction (float or np.array(dtype=float)) in the unit of m and positive direction is the same with y.
    """

    # Calculate max vertical displacement.
    uz_max = - (np.pi * vl * d**2) / (4 * np.sqrt(2 * np.pi) * k * (z0 - np.min(z)))

    y0 = -norm.ppf(delta) * k * z0

    u_z = uz_max * np.exp(- x**2 / (2 * k**2 * (z0 - z)**2)) * (
        norm.cdf((y - (ys + y0)) / (k * (z0 - z))) - 
        norm.cdf((y - yf) / (k * (z0 - z)))
    )
        
    u_x = (x / (z0 - z)) * u_z

    u_y = vl * d**2 / (8 * (z0 - z)) * (
        np.exp(-((y - (ys + y0))**2 + x**2) / (2 * k**2 * (z0 - z)**2)) - 
        np.exp(-((y - yf)**2 + x**2) / (2 * k**2 * (z0 - z)**2))
    )
    
    return u_x, u_y, u_z