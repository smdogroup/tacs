# Generate a plate mesh
import numpy as np

nx = 12
ny = 12
nodes = np.arange(1, (2 * nx + 1) * (2 * ny + 1) + 1, dtype=np.int32).reshape(
    2 * nx + 1, 2 * ny + 1
)

Lx = 1.0
Ly = 0.7

x = Lx * np.linspace(-0.5, 0.5, 2 * nx + 1)
y = Ly * np.linspace(-0.5, 0.5, 2 * ny + 1)

fp = open("plate.bdf", "w")
fp.write("$ Input file for a square shear-disp BC plate\n")
fp.write("SOL 103\nCEND\nBEGIN BULK\n")

# Write the grid points to a file
for j in range(2 * ny + 1):
    for i in range(2 * nx + 1):
        # Write the nodal data
        spc = " "
        coord_disp = 0
        coord_id = 0
        seid = 0

        fp.write(
            "%-8s%16d%16d%16.9e%16.9e*       \n"
            % ("GRID*", nodes[i, j], coord_id, x[i], y[j])
        )
        fp.write("*       %16.9e%16d%16s%16d        \n" % (0.0, coord_disp, spc, seid))

# Output 3rd order elements
elem = 1
part_id = 1
for j in range(0, nodes.shape[1] - 1, 2):
    for i in range(0, nodes.shape[0] - 1, 2):
        # Write the connectivity data
        # CQUAD9 elem id n1 n2 n3 n4 n5 n6
        #        n7   n8 n9
        fp.write(
            "%-8s%8d%8d%8d%8d%8d%8d%8d%8d\n"
            % (
                "CQUAD",
                elem,
                part_id,
                nodes[i, j],
                nodes[i + 2, j],
                nodes[i + 2, j + 2],
                nodes[i, j + 2],
                nodes[i + 1, j],
                nodes[i + 2, j + 1],
            )
        )
        fp.write(
            "        %8d%8d%8d\n"
            % (nodes[i + 1, j + 2], nodes[i, j + 1], nodes[i + 1, j + 1])
        )
        elem += 1

# Set up the plate BCs so that it has u = uhat, for shear disp control
# u = eps * y, v = eps * x, w = 0
eps = 1e-3
for j in range(2 * ny + 1):
    for i in range(2 * nx + 1):
        # check on boundary
        on_bndry = False
        if i == 0 or i == 2 * nx:
            on_bndry = True
        if j == 0 or j == 2 * ny:
            on_bndry = True
        if on_bndry:
            fp.write(
                "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "1", eps * y[j])
            )  # u = eps * y
            fp.write(
                "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "2", eps * x[i])
            )  # v = eps * x
            fp.write(
                "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "3", 0.0)
            )  # w = 0

            # x-rot, y-rot = 0 from paper 6
            # "An Adaptive Model reduction strategy for post-buckling analysis of stiffened structures"
            fp.write("%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "4", 0.0)) # w = 0
            fp.write("%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "5", 0.0)) # w = 0
        
        # on all nodes eliminate drill DOF
        fp.write("%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "6", 0.0)) # w = 0

fp.write("ENDDATA")
fp.close()
