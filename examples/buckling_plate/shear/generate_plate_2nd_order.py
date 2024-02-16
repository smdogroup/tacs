# Generate a plate mesh with CQUAD4 elements
import numpy as np
import matplotlib.pyplot as plt

nx = 12
ny = 12
nodes = np.arange(1, (nx + 1) * (ny + 1) + 1, dtype=np.int32).reshape(
    nx + 1, ny + 1
)

Lx = 1.0
Ly = 0.7

x = Lx * np.linspace(-0.5, 0.5, nx + 1)
y = Ly * np.linspace(-0.5, 0.5, ny + 1)

fp = open("plate.bdf", "w")
fp.write("$ Input file for a square shear-disp BC plate\n")
fp.write("SOL 103\nCEND\nBEGIN BULK\n")

# Write the grid points to a file
for j in range(ny + 1):
    for i in range(nx + 1):
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

# Output 2nd order elements
elem = 1
part_id = 1
for j in range(ny):
    for i in range(nx):
        # Write the connectivity data
        # CQUAD4 elem id n1 n2 n3 n4
        fp.write(
            "%-8s%8d%8d%8d%8d%8d%8d\n"
            % (
                "CQUAD4",
                elem,
                part_id,
                nodes[i, j],
                nodes[i + 1, j],
                nodes[i + 1, j + 1],
                nodes[i, j+1],
            )
        )
        elem += 1

# Set up the plate BCs so that it has u = uhat, for shear disp control
# u = eps * y, v = eps * x, w = 0
eps = 1e-3
for j in range(ny + 1):
    for i in range(nx + 1):
        # check on boundary
        on_bndry = False
        if i == 0 or i == nx:
            on_bndry = True
        if j == 0 or j == ny:
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
        #fp.write("%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "6", 0.0)) # w = 0

# plot the mesh to make sure it makes sense
X,Y = np.meshgrid(x,y)
W = X * 0.0
for j in range(ny+1):
    for i in range(nx+1):
        if i == 0 or i == nx or j == 0 or j == ny:
            W[i,j] = 1.0

#plt.scatter(X,Y)
#plt.contour(X,Y,W, corner_mask=True, antialiased=True)
#plt.show()

fp.write("ENDDATA")
fp.close()
