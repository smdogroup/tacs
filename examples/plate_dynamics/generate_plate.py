# Generate a plate mesh
import numpy as np

nx = 3
ny = 3
nodes = np.arange(1, (2 * nx + 1) * (2 * ny + 1) + 1, dtype=np.int).reshape(
    2 * nx + 1, 2 * ny + 1
)

x = np.linspace(0, 1, 2 * nx + 1)
y = np.linspace(0, 1, 2 * ny + 1)

fp = open("plate.bdf", "w")
fp.write("$ Input file for a square clamped plate\n")
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

# Set up the plate so that it is fully clamped
for i in range(2 * nx + 1):
    # Set the y = const edges
    spc = "123"
    fp.write("%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, 0], spc, 0.0))

fp.write("END BULK")
fp.close()
