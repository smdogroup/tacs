from pyNastran.bdf.bdf import read_bdf
import numpy as np

"""
Using pynastran open up original bdf file loop through all nodes on the edge of the plate 
and add a temperature boundary condition such that:
 T(theta) = T0 + dT * sin(2*theta)
Finally write out the modified bdf file so we can run it in the analysis.py script
"""
DEG2RAD = np.pi / 180
# Temperature constants for boundary function
T0 = 70
dT = 30
# Radius of plate
R = 1.0
# tolerance for nodes on edge
rtol = 1e-3
# read in original bdf
bdfInfo = read_bdf("circ-plate-no-bcs.bdf")
# Loop through every node
for nodeID in bdfInfo.node_ids:
    node = bdfInfo.nodes[nodeID]
    r = node.xyz[0]
    theta = node.xyz[1] * DEG2RAD
    # If node is on edge, add a new spc boundary condition
    if abs(r - R) < rtol:
        T = T0 + dT * np.sin(2 * theta)
        bdfInfo.add_spc(0, [nodeID], ["1"], [T])
# Write out new bdf with bcs
bdfInfo.write_bdf("circ-plate-dirichlet-bcs.bdf")
