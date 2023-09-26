"""
This example shows how to initialize pyTACS directly from a pyNASTRAN BDF object.
In this example we create an annulus out of shell elements using pyNASTRAN's BDF class.
We then use this class object to initialize pyTACS rather than reading in a file.
We analyze two load cases: A static case where we apply a load at the spider-webbed RBE2 centered inside the annulus
and modal case where we evaluate the mode shapes of the structure.
"""

import numpy as np
from pyNastran.bdf.bdf import BDF

from tacs import pyTACS, functions


# Geometric properties
Ri = 0.5
Ro = 1.0
z = 0.0
t = 0.05
# Material properties
E = 70e9
nu = 0.3
ys = 270e6
rho = 2700.0

# number of elements in r/theta
nr = 10
ntheta = 20
# Center node ID
center_nid = 1000

# Initialize our BDF object
bdfInfo = BDF()
# Create a cylindrical coord system
bdfInfo.add_cord2c(1, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0])
# Discretize in r and theta
r_array = np.linspace(Ri, Ro, nr)
theta_array = np.linspace(0.0, 360.0, ntheta, endpoint=False)


# Loop through and add each node
curID = 1
nodeIDs = {}
for r in r_array:
    for theta in theta_array:
        bdfInfo.add_grid(curID, [r, theta, z], 1)
        nodeIDs[r, theta] = curID
        curID += 1

# Loop through and add each element
curID = 1
for i in range(nr - 1):
    for j in range(ntheta):
        r0 = r_array[i]
        t0 = theta_array[j - 1]
        r1 = r_array[i + 1]
        t1 = theta_array[j]
        conn = [nodeIDs[r0, t0], nodeIDs[r0, t1], nodeIDs[r1, t1], nodeIDs[r1, t0]]
        bdfInfo.add_cquad4(curID, 1, conn)
        curID += 1

# Simply support the outer edge
for theta in theta_array:
    bdfInfo.add_spc(0, nodeIDs[Ro, theta], "123", 0.0)

# Add node at center for applying load
bdfInfo.add_grid(center_nid, [0.0, 0.0, 0.0], 0)

# Connect center node with RBE to inner edge
rbe_dep_nodes = []
for theta in theta_array:
    rbe_dep_nodes.append(nodeIDs[Ri, theta])
bdfInfo.add_rbe2(curID, center_nid, "123456", rbe_dep_nodes)

# Define material card for Aluminum
bdfInfo.add_mat1(1, E, None, nu, rho, St=ys)
# Define shell thickness
bdfInfo.add_pshell(1, 1, t)

# Initialize pyTACS from our BDF object
FEAAssembler = pyTACS(bdfInfo)
FEAAssembler.initialize()

# Apply a static load case with a x moment at the center
prob = FEAAssembler.createStaticProblem("center_moment")
prob.addLoadToNodes(center_nid, [0.0, 0.0, 0.0, 1e7, 0.0, 0.0], nastranOrdering=True)
prob.addFunction(
    "Izz", functions.MomentOfInertia, direction1=[0, 0, 1], direction2=[0, 0, 1]
)
prob.solve()
f = {}
prob.evalFunctions(f)
print(f)
prob.writeSolution()

# Perform modal analysis
modal = FEAAssembler.createModalProblem("modes", 1.0, 10)
modal.solve()
modal.writeSolution()
