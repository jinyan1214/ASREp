import os
import sys
import matplotlib.pyplot as plt
import numpy as np
# import ASRE models
file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(file_dir)

cur_dir = os.getcwd()
ASREp_dir = os.path.join(os.path.dirname(os.path.dirname(cur_dir)))
sys.path.append(ASREp_dir)
import ASREp.ground_deformation_models as gdm
import ASREp
import pandas as pd
from scipy import interpolate
from scipy.io import loadmat
np.set_printoptions(threshold=np.inf)
import pickle
import time



# Read the mesh file
input_dir = "./input"
output_dir = "./output"
mesh_file = os.path.join(input_dir, "burd_mat_20250901T162943.mat")
mesh_data = loadmat(mesh_file)

mesh_data['wholeElem2n'] = mesh_data['wholeElem2n'] - 1  # convert to zero-based index
mesh_data['elem2nInter'] = mesh_data['elem2nInter'] - 1  # convert to zero-based index
# Define building properties
Eb = 3e9  # Young's modulus of the building (Pa)
nu = 0.2  # Poisson's ratio of the building
rho = 23.75*10**3  # Density of the building (N/m^3) burd et al. 2022 table 1
# Define soil properties
EsNominal = 14.5e6 # Nominal Young's modulus of the soil (Pa)
nus = 0.49 # Poisson's ratio of the soil
# Define interface properties
mu_int = 0.4 # friction coefficient of the interface
lim_t_int = 0 # tensile strength of the interface (Pa)
lim_c_int = -348000 # compressive strength of the interface (Pa). It must


wholeNodesXYZ = mesh_data['wholeNodesXYZ']  # (N, 3) array
min_z = np.min(wholeNodesXYZ[:, 2])

groundNodeInd = np.where(wholeNodesXYZ[:, 2] == min_z)[0]  # indices of ground nodes
x = wholeNodesXYZ[groundNodeInd, 0]
y = wholeNodesXYZ[groundNodeInd, 1]
# z = wholeNodesXYZ[groundNodeInd, 2]
z = np.zeros_like(x)  # set z to zero for ground surface
vl = 1.65/100.0
d = 11.0
z0 = 23.0
ys = 200.0 # assume the tunnel has fully passed the building
yf = -200.0 # assume the tunnel start point is far away from the building
k = 0.57
delta = 0.3
u_x, u_y, u_z = gdm.ground_disp_Zhao_2023(x, y, z, vl, d, z0, ys, yf, k, delta)

np.savetxt(os.path.join(output_dir, 'greenfield_ground_disp.txt'), np.column_stack((u_x, u_y, u_z)))
print("Greenfield ground displacements saved to 'greenfield_ground_disp.txt'.")

# Prepare inputs for ASRElib3DSolid
model_el = ASREp.ASRE_3D_solid_model(
    mesh_data['wholeNodesXYZ'],
    mesh_data['wholeElem2n'],
    mesh_data['interNodesXYZ'],
    mesh_data['elem2nInter'],
    solver = 'elastic'
)
model_el.set_building_properties(
    Eb = Eb,
    nu = nu,
    rho = rho,
)
model_el.set_soil_properties(
    EsNominal = EsNominal,
    nus = nus
)
model_el.set_interface_properties(
    mu_int = mu_int,
    lim_t_int = lim_t_int,
    lim_c_int = lim_c_int

)

start_time = time.time()
print("Running elastic model...")
model_el.run_model(u_x, u_y, u_z)
end_time = time.time()
print(f"Elastic model completed in {end_time - start_time:.2f} seconds.")

# Pickle the model_el for post-processing
with open(os.path.join(output_dir, 'model_el.pkl'), 'wb') as f:
    pickle.dump(model_el, f)
print("Elastic model results saved to 'model_el.pkl'.")

# Repeat for elastoplastic model
model_ep = ASREp.ASRE_3D_solid_model(
    mesh_data['wholeNodesXYZ'],
    mesh_data['wholeElem2n'],
    mesh_data['interNodesXYZ'],
    mesh_data['elem2nInter'],
    solver = 'elasto-plastic'
)
model_ep.set_building_properties(
    Eb = Eb,
    nu = nu,
    rho = rho,
)
model_ep.set_soil_properties(
    EsNominal = EsNominal,
    nus = nus
)
model_ep.set_interface_properties(
    mu_int = mu_int,
    lim_t_int = lim_t_int,
    lim_c_int = lim_c_int

)

start_time = time.time()
print("Running elastoplastic model...")
model_ep.run_model(u_x, u_y, u_z)
end_time = time.time()
print(f"Elastoplastic model completed in {end_time - start_time:.2f} seconds.")

# Pickle the model_ep for post-processing
with open(os.path.join(output_dir, 'model_ep.pkl'), 'wb') as f:
    pickle.dump(model_ep, f)
print("Elastoplastic model results saved to 'model_ep.pkl'.")