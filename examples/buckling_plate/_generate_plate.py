# Generate a plate mesh with CQUAD4 elements
import numpy as np
import matplotlib.pyplot as plt


def generate_plate(
    Lx=1.0, Ly=1.0, nx=30, ny=30, exx=0.0, eyy=0.0, exy=0.0, clamped=True
):
    """
    create mixed compression and shear load problem with
    """

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
            fp.write(
                "*       %16.9e%16d%16s%16d        \n" % (0.0, coord_disp, spc, seid)
            )

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
                    nodes[i, j + 1],
                )
            )
            elem += 1

    # Set up the plate BCs so that it has u = uhat, for shear disp control
    # u = eps * y, v = eps * x, w = 0
    for j in range(ny + 1):
        for i in range(nx + 1):
            u = exy * y[j]
            v = exy * x[i]

            if i == nx or exy != 0:
                u -= exx * x[i]
            elif j == ny:
                v -= eyy * Ly
            elif i == 0 or j == 0:
                pass

            # check on boundary
            if i == 0 or j == 0 or i == nx or j == ny:
                if clamped or (i == 0 and j == 0):
                    fp.write(
                        "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "3456", 0.0)
                    )  # w = theta_x = theta_y
                else:
                    fp.write(
                        "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "36", 0.0)
                    )  # w = theta_x = theta_y
                if exy != 0 or i == 0 or i == nx:
                    fp.write(
                        "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "1", u)
                    )  # u = eps_xy * y
                if exy != 0.0 or j == 0:
                    fp.write(
                        "%-8s%8d%8d%8s%8.6f\n" % ("SPC", 1, nodes[i, j], "2", v)
                    )  # v = eps_xy * x

    # # plot the mesh to make sure it makes sense
    # X, Y = np.meshgrid(x, y)
    # W = X * 0.0
    # for j in range(ny + 1):
    #     for i in range(nx + 1):
    #         if i == 0 or i == nx or j == 0 or j == ny:
    #             W[i, j] = 1.0

    # plt.scatter(X,Y)
    # plt.contour(X,Y,W, corner_mask=True, antialiased=True)
    # plt.show()

    fp.write("ENDDATA")
    fp.close()


if __name__ == "__main__":
    generate_plate(
        Lx=1.0, Ly=0.7, nx=30, ny=20, exx=0.0, eyy=0.0, exy=0.001, clamped=False
    )
