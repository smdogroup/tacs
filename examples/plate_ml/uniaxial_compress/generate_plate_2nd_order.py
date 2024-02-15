# Generate a plate mesh with CQUAD4 elements
import numpy as np
import matplotlib.pyplot as plt

def generate_plate(Lx=1.0,Ly=0.7,nx=12,ny=12, displacement_control=False):

    nodes = np.arange(1, (nx + 1) * (ny + 1) + 1, dtype=np.int32).reshape(
        nx + 1, ny + 1
    )

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
            #uhat = 0.001
            if i == 0 and j == 0: # node 1
                fp.write(
                    "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "123456", 0.0)
                ) # u = v = w = theta_x = theta_y = theta_z = 0 on right edge
            elif j == 0:
                fp.write(
                    "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "236", 0.0)
                )  # v = 0 on bottom edge
            elif i == 0:
                fp.write(
                    "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "36", 0.0)
                )  # w = theta_z = 0 on all edges
        
            if i == nx and j != 0:
                fp.write(
                    "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "36", 0.0)
                ) # nothing on right edge, was u = 0 on right edge

            if j == ny:
                if displacement_control:
                    fp.write(
                        "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "2", -0.001)
                    ) # v = vhat on top edge
                if i != 0  and i != nx:
                    fp.write(
                        "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "36", 0.0)
                    ) # w = theta_z = 0

    # apply LOAD cards on upper and lower edges..
    # since v = 0 on lower edge, first try just on upper edge
    if not displacement_control:
        for i in range(nx+1):
            if i == 0 or i == nx:
                Fy = 1000.0
            else:
                Fy = 2000.0
            # j = 0
            # don't need top edge since it's constrained there
            # fp.write(
            #     "%-8s%8d%8d%8d%8.1f%8.1f%8.2f%8.2f\n" % ("FORCE", 1, nodes[i,j], 0, 1.0, 0, Fy, 0)
            # ) # bottom edge

            j = ny
            fp.write(
                "%-8s%8d%8d%8d%8.1f%8.1f%8.2f%8.2f\n" % ("FORCE", 1, nodes[i,j], 0, 1.0, 0, -Fy, 0)
            ) # top edge

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


if __name__=="__main__":
    generate_plate(Lx=1.0, Ly=0.7, nx=12, ny=12, displacement_control=False)