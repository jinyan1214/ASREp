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
Gs = 6e6  # Shear modulus of the soil (Pa)
nus = 0.2 # Poisson's ratio of the soil
EsNominal = 2*Gs*(1+nus) # Nominal Young's modulus of the soil (Pa)
# Define interface properties
mu_int = 0.3 # friction coefficient of the interface
lim_t_int = 0 # tensile strength of the interface (Pa)
lim_c_int = -347.8261*1000 # compressive strength of the interface (Pa). It must

# The current soil flexibility model only works for ground surface loading
mesh_data['wholeNodesXYZ'][:, 2] = mesh_data['wholeNodesXYZ'][:, 2] - np.min(mesh_data['wholeNodesXYZ'][:, 2])
mesh_data['interNodesXYZ'][:, 2] = 0.0
wholeNodesXYZ = mesh_data['wholeNodesXYZ']  # (N, 3) array
min_z = np.min(wholeNodesXYZ[:, 2])
print(f"Minimum z value (ground surface): {min_z} m")
groundNodeInd = np.where(wholeNodesXYZ[:, 2] == min_z)[0]  # indices of ground nodes
x = wholeNodesXYZ[groundNodeInd, 0]
y = wholeNodesXYZ[groundNodeInd, 1]
z = wholeNodesXYZ[groundNodeInd, 2]
vl = 1.65/100.0
d = 11.0
z0 = 23.0 - 1.0 # The 1.0 is to model the embedding of building foundation
ys = -50.0 # assume the tunnel has fully passed the building
yf = 500.0 # assume the tunnel start point is far away from the building
k = 0.57
delta = 0.3
# If use empirical GF model
# u_x, u_y, u_z = gdm.ground_disp_Zhao_2023(x, y, z, vl, d, z0, ys, yf, k, delta)
# If use Burd et al. 2022 GF (based on FEM)
u_z_yiu_raw = loadmat(os.path.join(input_dir, "Yiu_GF.mat"))['Yiu_GF']
f = interpolate.interp1d(u_z_yiu_raw[:,0], u_z_yiu_raw[:,1])
u_z = -f(x)/1000.0
u_x = x / (z0 - z) * u_z
u_y = np.zeros_like(u_z)

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
print("Calculating principal strains...")
model_el.calculate_principal_strain()
end_time = time.time()
print(f"Strain calculation completed in {end_time - start_time:.2f} seconds.")

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
model_ep.run_model(u_x, u_y, u_z, print_iteration=1)
end_time = time.time()
print(f"Elastoplastic model completed in {end_time - start_time:.2f} seconds.")
print("Calculating principal strains...")
model_ep.calculate_principal_strain()
end_time = time.time()
print(f"Strain calculation completed in {end_time - start_time:.2f} seconds.")

# Pickle the model_ep for post-processing
with open(os.path.join(output_dir, 'model_ep.pkl'), 'wb') as f:
    pickle.dump(model_ep, f)
print("Elastoplastic model results saved to 'model_ep.pkl'.")

# Save principal strains to text files
np.savetxt(os.path.join(output_dir, 'principal_strain_elastic.txt'), 
           np.column_stack((model_el.result_tensile_ptr, model_el.result_compressive_ptr)))
np.savetxt(os.path.join(output_dir, 'principal_strain_elastoplastic.txt'), 
           np.column_stack((model_ep.result_tensile_ptr, model_ep.result_compressive_ptr)))
print("Principal strains saved to 'principal_strain_elastic.txt' and 'principal_strain_elastoplastic.txt'.")

# Save nodal displacements to text files
np.savetxt(os.path.join(output_dir, 'nodal_displacement_elastic.txt'),
           model_el.result_array_ptr)
np.savetxt(os.path.join(output_dir, 'nodal_displacement_elastoplastic.txt'),
           model_ep.result_array_ptr)
print("Nodal displacements saved to 'nodal_displacement_elastic.txt' and 'nodal_displacement_elastoplastic.txt'.")

print("99% quantile of elastoplastic principal tensile strain:", 
      np.quantile(model_ep.result_tensile_ptr, 0.99))