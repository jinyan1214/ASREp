import openseespy.opensees as ops
import opsvis as opsv
import matplotlib.pyplot as plt
# Clear existing model
ops.wipe()

# Define Model Builder
ops.model('basic', '-ndm', 3, '-ndf', 6)  # 3D model with 6 DOF per node

# Define parameters
n = 3        # Number of bays
l = 5.0      # Length of each bay
h = 3.0      # Height of each story
numStories = 2  # Number of stories

# Material properties
bfoot = 10
dfoot = 1
E = 200e9  # Young's Modulus in Pa
A = bfoot * dfoot   # Area of the element in m^2
EoverG = 0.001
G = E/EoverG
Iy = bfoot * dfoot**3/12
Iz = dfoot**3 * bfoot/12  # Moment of Inertia in m^4
J = Iy + Iz

# Create nodes
for i in range(numStories + 1):
    for j in range(n + 1):
        nodeTag = i * (n + 1) + j + 1
        x = j * l
        z = i * h
        ops.node(nodeTag, x, 0, z)

# Define geometric transformation
horizontal_gTTag = 1
vertical_gTTag = 2
ops.geomTransf('Linear', horizontal_gTTag, 0, 1, 0)
ops.geomTransf('Linear', vertical_gTTag, 0, 1, 0)
opsv.plot_model()


ops.element('elasticBeamColumn', 1, 1, 2, A, E, G, J, Iy, Iz, 1)