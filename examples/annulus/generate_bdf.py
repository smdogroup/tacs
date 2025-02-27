"""
Create a .bdf file for a plane stress problem for a circular annulus

This only creates one quarter of problem.
"""

import argparse
import numpy as np

# Note that for performance results we use either:
# nyelems = 750
# nyelems = 375


def writeElementDefinition(fp, elemType, elemID, compID, nodeIDs):
    elementDef = [elemType, elemID, compID] + nodeIDs
    entries = 0
    defString = ""
    for entry in elementDef:
        defString += f"{entry:<8}"
        entries += 1
        if entries % 9 == 0:
            defString += "\n" + 8 * " "
            entries += 1
    if defString[-1] != "\n":
        defString += "\n"
    fp.write(defString)


# Set up the argument parser object
parser = argparse.ArgumentParser(
    description="Generate a .bdf file for a circular annulus"
)
parser.add_argument(
    "--nyelems", type=int, default=100, help="number of elements in the y-direction"
)
parser.add_argument("--order", type=int, default=2, help="element order")
args = parser.parse_args()

# Extract the arguments
nyelems = args.nyelems
order = args.order

# Compute the number of element circumferentially
nxelems = int(1.8 * nyelems)

# Compute the number of nodes
nx = (order - 1) * nxelems + 1
ny = (order - 1) * nyelems + 1

# Generate a bdf file
fp = open("circular_annulus%dx%d_order%d.bdf" % (nx, ny, order), "w")
fp.write("$ Input file for a plane stress problem\n")
fp.write("SOL 103\nCEND\nBEGIN BULK\n")

# Set the outer and inner radii
Router = 10.0
Rinner = 4.0

# Set up the first mapped section
nodes = np.zeros((nx, ny), dtype=int)
x = np.zeros((nx, ny))
y = np.zeros((nx, ny))

# Compute the node numbers
nodes[:] = np.arange(1, nx * ny + 1).reshape((nx, ny))

# Set the locations of theta/u
theta = np.linspace(np.pi / 2.0, 0.0, nx)
u = np.linspace(Rinner, Router, ny)

# Set the nodal coordinates the coordinates for the center section
for j in range(ny):
    for i in range(nx):
        r = u[j]
        x[i, j] = r * np.cos(theta[i])
        y[i, j] = r * np.sin(theta[i])

# Write the grid points to a file
for j in range(ny):
    for i in range(nx):
        # Write the nodal data
        spc = " "
        coord_disp = 0
        coord_id = 0
        seid = 0

        fp.write(
            "%-8s%16d%16d%16.9e%16.9e*       \n"
            % ("GRID*", nodes[i, j], coord_id, x[i, j], y[i, j])
        )
        fp.write("*       %16.9e%16d%16s%16d        \n" % (0.0, coord_disp, spc, seid))

elem = 1
part_id = 1
if order == 2:
    # Ouput 2nd order elements
    for j in range(nodes.shape[1] - 1):
        for i in range(nodes.shape[0] - 1):
            # Write the connectivity data
            fp.write(
                "%-8s%8d%8d%8d%8d%8d%8d\n"
                % (
                    "CQUAD4",
                    elem,
                    1,
                    nodes[i, j],
                    nodes[i + 1, j],
                    nodes[i + 1, j + 1],
                    nodes[i, j + 1],
                )
            )
elif order == 3:
    # Output 3rd order elements
    for j in range(0, nodes.shape[1] - 1, order - 1):
        for i in range(0, nodes.shape[0] - 1, order - 1):
            nodeOrdering = [
                [0, 0],
                [2, 0],
                [2, 2],
                [0, 2],
                [1, 0],
                [2, 1],
                [1, 2],
                [0, 1],
                [1, 1],
            ]
            elementDef = ["CQUAD9", elem, elem]
            nodeIDs = [nodes[i + node[0], j + node[1]] for node in nodeOrdering]
            writeElementDefinition(fp, "CQUAD9", elem, 1, nodeIDs)
            elem += 1
elif order == 4:
    # Output 3rd order elements
    for j in range(0, nodes.shape[1] - 1, order - 1):
        for i in range(0, nodes.shape[0] - 1, order - 1):
            # Write the connectivity data
            # CQUAD16 elem id n1 n2 n3 n4 n5 n6
            #         n7   n8 n9 n10 n11 n12 n13
            #         n14  n15 n16
            nodeOrdering = [
                [0, 0],
                [3, 0],
                [3, 3],
                [0, 3],
                [1, 0],
                [2, 0],
                [3, 1],
                [3, 2],
                [2, 3],
                [1, 3],
                [0, 2],
                [0, 1],
                [1, 1],
                [2, 1],
                [2, 2],
                [1, 2],
            ]
            elementDef = ["CQUAD16", elem, elem]
            nodeIDs = [nodes[i + node[0], j + node[1]] for node in nodeOrdering]
            writeElementDefinition(fp, "CQUAD16", elem, 1, nodeIDs)
            elem += 1
elif order == 5:
    for j in range(0, nodes.shape[1] - 1, order - 1):
        for i in range(0, nodes.shape[0] - 1, order - 1):
            nodeOrdering = [
                [0, 0],
                [4, 0],
                [4, 4],
                [0, 4],
                [1, 0],
                [2, 0],
                [3, 0],
                [4, 1],
                [4, 2],
                [4, 3],
                [3, 4],
                [2, 4],
                [1, 4],
                [0, 3],
                [0, 2],
                [0, 1],
                [1, 1],
                [3, 1],
                [3, 3],
                [1, 3],
                [2, 1],
                [3, 2],
                [2, 3],
                [1, 2],
                [2, 2],
            ]
            elementDef = ["CQUAD25", elem, elem]
            nodeIDs = [nodes[i + node[0], j + node[1]] for node in nodeOrdering]
            writeElementDefinition(fp, "CQUAD25", elem, 1, nodeIDs)
            elem += 1


# Set up the plate so that it is fully clamped
for k in range(ny):
    # Set the y = const edges
    # Clamp w, roty, rotz
    spc = "2"
    fp.write("%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[0, k], spc, 0.0))

for k in range(ny):
    spc = "1"
    fp.write("%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[-1, k], spc, 0.0))

fp.write("END BULK")
fp.close()
