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

# ---- 2D greenfield displacement: retaining wall cavity superposition ----
def u_GON(z, x, depth, radius, epsilon, soil_ovalization, volumetric):
    """
    Zheng et al. analytical displacement increment from a single equivalent cavity.
    Used internally by wall_deflection to superpose contributions along the wall depth.

    Parameters
    ----------
    z : float or np.ndarray
        Depth of observation point below ground surface. Unit: m
    x : float or np.ndarray
        Horizontal distance from the retaining wall. Unit: m
    depth : float
        Depth of the equivalent cavity centre. Unit: m
    radius : float
        Equivalent cavity radius. Unit: m
    epsilon : float
        Sign of wall deflection at this depth (±1). Unit: unitless
    soil_ovalization : float
        Ovalization parameter (soil_ovalization * epsilon from caller). Unit: unitless
    volumetric : float
        Volumetric parameter (1 = no volume change). Unit: unitless

    Returns
    -------
    tuple
        A tuple containing:
        - ux (float or np.ndarray): Horizontal displacement increment. Unit: m, positive away from wall.
        - uz (float or np.ndarray): Vertical displacement increment. Unit: m, positive downward.
    """
    ovalization_ratio = soil_ovalization / epsilon
    normalized_x = x / depth
    depth_below = z - depth
    depth_above = z + depth
    normalized_depth_below = depth_below / depth
    normalized_depth_above = depth_above / depth
    normalized_depth = z / depth
    radius_below = np.sqrt(x ** 2 + depth_below ** 2) / depth
    radius_above = np.sqrt(x ** 2 + depth_above ** 2) / depth

    ux = 2 * epsilon * radius * (radius / depth) ** (2 * volumetric - 1) * (
        -((normalized_x * (1 - (
            ovalization_ratio * (normalized_x ** 2 - normalized_depth_below ** 2)) / radius_below ** 2)) / (
                2 * radius_below ** (2 * volumetric))) -
        (normalized_x * (1 - (
            ovalization_ratio * (normalized_x ** 2 - normalized_depth_above ** 2)) / radius_above ** 2)) / (
                2 * radius_above ** (2 * volumetric)) +
        (4 * normalized_x * normalized_depth * (normalized_depth_above / radius_above ** 2 - (ovalization_ratio * (
            normalized_x ** 2 - 3 * normalized_depth_above ** 2)) / radius_above ** 4)) / (
                2 * radius_above ** (2 * volumetric)))

    uz = 2 * epsilon * radius * (radius / depth) ** (2 * volumetric - 1) * (
        -((normalized_depth_below * (1 - (
            ovalization_ratio * (normalized_x ** 2 - normalized_depth_below ** 2)) / radius_below ** 2)) / (
                2 * radius_below ** (2 * volumetric))) +
        (normalized_depth_above * (1 + (
            ovalization_ratio * (normalized_x ** 2 - normalized_depth_above ** 2)) / radius_above ** 2)) / (
                2 * radius_above ** (2 * volumetric)) -
        ((2 * (normalized_depth + ovalization_ratio) * (
            normalized_x ** 2 - normalized_depth_above ** 2)) / radius_above ** 2 +
         (4 * ovalization_ratio * normalized_depth * normalized_depth_above * (
             3 * normalized_x ** 2 - normalized_depth_above ** 2)) / radius_above ** 4) / (
                2 * radius_above ** (2 * volumetric)))

    return ux, uz


def wall_deflection(retaining_wall_depth, avg_wall_displacement, wall_deflection_shape,
                    soil_ovalization, volumetric, x_coords=None, z_coords=None,
                    C1_val=None, C2_val=None):
    """
    Compute 2D greenfield displacements behind a retaining wall using Zheng et al.'s
    analytical solution. The wall deflection profile is idealised as one of five
    parametric shapes (M0–M5); each depth slice is converted to an equivalent
    circular cavity whose contribution is superposed via u_GON.

    Reference: Zheng et al. (see Wall_deflection_v38 in MastersThesis repo).

    Parameters
    ----------
    retaining_wall_depth : float
        Total embedded depth of the retaining wall. Unit: m
    avg_wall_displacement : float
        Average wall displacement normalised by the wall height. Unit: unitless
    wall_deflection_shape : int
        Deflection mode: 0=Uniform, 1=Cantilever, 2=Parabola,
        3=Composite (35% cantilever + 65% parabola), 4=Kick-in, 5=Custom (C1/C2).
        Unit: unitless
    soil_ovalization : float
        Soil ovalization parameter. Unit: unitless
    volumetric : float
        Volumetric parameter (1 = elastic, no contraction/expansion). Unit: unitless
    x_coords : np.ndarray, optional
        Horizontal distances from wall at which to evaluate displacements. Unit: m
        Defaults to 100 points from 0.01*Hwall to 2*Hwall.
    z_coords : np.ndarray, optional
        Depths below ground surface at which to evaluate displacements. Unit: m
        Defaults to ground surface only (z=0).
    C1_val : float, optional
        Weight of cantilever component for shape 5. Must satisfy C1+C2+C3=1.
    C2_val : float, optional
        Weight of parabolic component for shape 5. Must satisfy C1+C2+C3=1.

    Returns
    -------
    tuple
        A tuple containing:
        - horizontal_displacement (np.ndarray): Shape (len(z_coords), len(x_coords)). Unit: m
        - vertical_displacement (np.ndarray): Shape (len(z_coords), len(x_coords)). Unit: m
        - wall_deflection_all (np.ndarray): Wall deflection profile over full depth. Unit: m
        - volume_loss_wall (float): Integrated wall volume loss. Unit: m²
        - max_wall_deflection (float): Peak wall deflection. Unit: m
        - check_ratio (float): Ratio of numerical to analytical volume loss (should be ~1). Unit: unitless
        - volume_loss_surface (float): Integrated vertical displacement at z_coords[0]. Unit: m²
        - ratio_volume (float): Ratio of surface to wall volume loss. Unit: unitless
    """
    deltaH = retaining_wall_depth / 100

    if x_coords is None:
        x_coords = np.arange(retaining_wall_depth * 0.01, 2.0 * retaining_wall_depth,
                             retaining_wall_depth * 0.01)

    if z_coords is None:
        z_coords = np.array([0.0])

    def wall_uniform_deflection(depth):
        return 1 * avg_wall_displacement * retaining_wall_depth + 0 * depth

    def wall_cantilever_deflection(depth):
        return 2 * avg_wall_displacement * retaining_wall_depth - 2 * avg_wall_displacement * depth

    def wall_parabola_deflection(depth):
        return -6 * avg_wall_displacement / retaining_wall_depth * depth ** 2 + 6 * avg_wall_displacement * depth

    def wall_composite_deflection(depth):
        return (0.35 * (2 * avg_wall_displacement * retaining_wall_depth - 2 * avg_wall_displacement * depth) +
                0.65 * (-6 * avg_wall_displacement / retaining_wall_depth * depth ** 2 + 6 * avg_wall_displacement * depth))

    def wall_kickin_deflection(depth):
        return -4 * avg_wall_displacement / retaining_wall_depth * depth ** 2 + 4.7 * avg_wall_displacement * depth

    def wall_custom_deflection(depth):
        C3_val = 1 - C1_val - C2_val
        if not (0.99999 < C1_val + C2_val + C3_val < 1.00001):
            raise ValueError(f'C1 + C2 + C3 must equal 1, got {C1_val + C2_val + C3_val}')
        return (C1_val * (2 * avg_wall_displacement * retaining_wall_depth - 2 * avg_wall_displacement * depth) +
                C2_val * (-6 * avg_wall_displacement / retaining_wall_depth * depth ** 2 + 6 * avg_wall_displacement * depth) +
                C3_val * (2 * avg_wall_displacement * depth))

    shape_map = {
        0: wall_uniform_deflection,
        1: wall_cantilever_deflection,
        2: wall_parabola_deflection,
        3: wall_composite_deflection,
        4: wall_kickin_deflection,
        5: wall_custom_deflection,
    }
    if wall_deflection_shape not in shape_map:
        raise ValueError(f'wall_deflection_shape must be 0–5, got {wall_deflection_shape}')
    wall_deflection_func = shape_map[wall_deflection_shape]

    depth_vec = np.arange(0, retaining_wall_depth + deltaH, deltaH)
    cavity_depth_vec = np.arange(deltaH / 2, retaining_wall_depth - deltaH / 2 + deltaH, deltaH)
    num_cavities = len(cavity_depth_vec)

    x_section, z_section = np.meshgrid(x_coords, z_coords)
    horizontal_displacement = np.zeros_like(x_section)
    vertical_displacement = np.zeros_like(x_section)

    deflection_vec = wall_deflection_func(cavity_depth_vec)
    sign_vector = np.sign(deflection_vec)
    radius_vector = np.sqrt(np.abs(deflection_vec * deltaH / np.pi))

    for i in range(num_cavities):
        epsilon = sign_vector[i]
        radius = radius_vector[i]
        cavity_depth = cavity_depth_vec[i]
        u_xinc, u_zinc = u_GON(z_section, x_section, cavity_depth, radius, epsilon,
                               soil_ovalization * epsilon, volumetric)
        horizontal_displacement += u_xinc
        vertical_displacement += u_zinc

    wall_deflection_all = wall_deflection_func(depth_vec)
    volume_loss_wall = np.trapz(wall_deflection_all, depth_vec)
    max_wall_deflection = np.max(wall_deflection_all)
    volume_loss_wall_v2 = avg_wall_displacement * retaining_wall_depth ** 2
    check_ratio = volume_loss_wall / volume_loss_wall_v2

    volume_loss_surface = np.trapz(vertical_displacement[0, :], x_section[0, :])
    ratio_volume = volume_loss_surface / volume_loss_wall

    return (horizontal_displacement, vertical_displacement, wall_deflection_all,
            volume_loss_wall, max_wall_deflection, check_ratio, volume_loss_surface, ratio_volume)
