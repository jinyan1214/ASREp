import numpy as np

# %% WALL PARAMETERS
retaining_wall_depth    = 20.7          # [m]  Embedded depth of the retaining wall
avg_wall_displacement   = 0.027 / 100   # [-]  Average wall displacement / wall height (0.027%)
wall_deflection_shape   = 5             # [-]  0=Uniform, 1=Cantilever, 2=Parabola,
                                        #      3=Composite, 4=Kick-in, 5=Custom (C1/C2)
C1_val                  = -0.36         # [-]  Cantilever weight  (shape 5)
C2_val                  =  1.29         # [-]  Parabola weight    (shape 5)
                                        #      C3 = 1 - C1 - C2 = 0.07 (kick-in)

# %% SOIL PARAMETERS
soil_ovalization        = 0.0           # [-]  Ovalization parameter
volumetric              = 1.0           # [-]  1 = elastic (no volume change)
soil_poisson            = 0.5           # [-]  Undrained Poisson's ratio

# %% BUILDING PARAMETERS
building_length         = 12.0          # [m]  Length of building
building_height         = 36.0          # [m]  Equivalent building height
foundation_depth        = 2.2           # [m]  Foundation depth below ground surface
building_offset         = 11.0          # [m]  Horizontal distance from wall to building start
building_width          = 1.0           # [m]  Foundation width (per running metre)
num_nodes               = 101           # [-]  Nodes along building length
num_elements            = num_nodes - 1

neutral_line            = 18.0          # [m]  Neutral axis from bottom fibre (pure bending)

# %% STRUCTURAL PROPERTIES (Timoshenko equivalent beam)
# Stiffness values as per MastersThesis uncertainty parameters
EA      = 25.9          # [GN]    Axial stiffness
EI      = 2777.0        # [GN·m²] Bending stiffness
GAs     = 5.95          # [GN]    Shear stiffness

EA      = EA   * 1e9    # [N]
EI      = EI   * 1e9    # [N·m²]
GAs     = GAs  * 1e9    # [N]

poissons_ratio  = 0.3               # [-]  Building material Poisson's ratio

Ab      = 2.76                      # [m²] Cross-sectional area of equivalent beam
Eb      = EA  / Ab                  # [Pa] Equivalent Young's modulus
Ib      = EI  / Eb                  # [m⁴] Second moment of area
Gb      = GAs / 1.65                # [Pa] Shear modulus (1.65 = effective shear area As [m²])
EoverG  = Eb  / Gb                  # [-]  = 2*(1 + poissons_ratio) ≈ 2.6

d_NA    = neutral_line              # [m]  Distance from neutral axis to beam-soil interface
q_foot  = 3.2 * 10 * 1000          # [N/m] Uniform foundation load

# %% SOIL–STRUCTURE INTERFACE
Es_isotropic    = 175 * 1e6         # [Pa]  Isotropic soil stiffness (Es = 175 MN/m²)
mu_int          = 0.0               # [-]   Friction coefficient (0 for elastic solver)

# %% GRID AND COORDINATES
Hwall = retaining_wall_depth

# 2D grid for contour plots — x horizontal from wall, z depth below surface
x_grid = np.linspace(0.05 * Hwall, 5.0 * Hwall, 150)   # [m]
z_grid = np.linspace(0.0,          5.0 * Hwall, 150)    # [m]

# Building line at foundation depth
x_coords_building = building_offset + np.linspace(0.0, building_length, num_nodes)  # [m]
z_coords_building = np.array([foundation_depth])                                     # [m]
