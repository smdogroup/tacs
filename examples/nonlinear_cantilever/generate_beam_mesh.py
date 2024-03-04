import numpy
import argparse


def write_80(line):
    newline = "{:80s}\n".format(line.strip("\n"))
    fout.write(newline)


parser = argparse.ArgumentParser()
parser.add_argument("--lx", type=float, default=10.0, help="Beam length")
parser.add_argument("--ly", type=float, default=1.0, help="Beam width")
parser.add_argument("--ncx", type=int, default=1, help="# components along length")
parser.add_argument("--ncy", type=int, default=1, help="# components along width")
parser.add_argument("--nex", type=int, default=16, help="# elements along length")
parser.add_argument("--ney", type=int, default=1, help="# elements along width")
parser.add_argument("--name", type=str, default="Beam")
args = parser.parse_args()

# Overall plate dimensions
lx = args.lx
ly = args.ly

# Number of components in each direction
ncx = args.ncx
ncy = args.ncy
ncomp = ncx * ncy

# Number of elements along each edge of a single panel
nex = args.nex
ney = args.ney

# Nodes
nx = ncx * nex + 1
ny = ncy * ney + 1
xtmp = numpy.linspace(0, lx, nx)
ytmp = numpy.linspace(0, ly, ny)
X, Y = numpy.meshgrid(xtmp, ytmp)
Z = numpy.zeros_like(X)
nodes = numpy.stack((X, Y, Z), axis=2)
nmat = nodes.reshape((nx * ny, 3))

# Node numbering
nid = numpy.zeros((nx, ny), dtype="intc")
bcnodes = []
count = 1
for i in range(ny):
    for j in range(nx):
        nid[j, i] = count
        if j == 0:
            bcnodes.append(count)
        count += 1

# Connectivity
nex = nx - 1
ney = ny - 1
ne = nex * ney
conn = {i + 1: [] for i in range(ncomp)}
ie = 1
for i in range(ney):
    for j in range(nex):
        compID = i // ney * ncx + j // nex + 1
        conn[compID].append(
            [ie, nid[j, i], nid[j + 1, i], nid[j + 1, i + 1], nid[j, i + 1]]
        )
        ie += 1


# Write BDF
output_file = args.name + ".bdf"

with open(output_file, "w") as fout:
    write_80("SOL 103")
    write_80("CEND")
    write_80("BEGIN BULK")

    # Make component names
    compNames = {}
    compID = 1
    for i in range(ncy):
        for j in range(ncx):
            compNames[compID] = "PLATE.{:03d}/SEG.{:02d}".format(i, j)
            compID += 1

    def write_bulk_line(key, items, format="small"):
        if format == "small":
            width = 8
            writekey = key
        elif format == "large":
            width = 16
            writekey = key + "*"
        line = "{:8s}".format(writekey)
        for item in items:
            if type(item) in [int, numpy.int64, numpy.int32]:
                line += "{:{width}d}".format(item, width=width)[:width]
            elif type(item) in [float, numpy.float64]:
                line += "{: {width}f}".format(item, width=width)[:width]
            elif type(item) is str:
                line += "{:{width}s}".format(item, width=width)[:width]
            else:
                print(type(item), item)
            if len(line) == 72:
                write_80(line)
                line = " " * 8
        if len(line) > 8:
            write_80(line)

    # Write nodes
    for i in range(nx * ny):
        write_bulk_line("GRID", [i + 1, 0, nmat[i, 0], nmat[i, 1], nmat[i, 2], 0, 0, 0])

    # Write elements
    compID = 1
    for key in conn:
        famPrefix = "$       Shell element data for family    "
        famString = "{}{:39s}".format(famPrefix, compNames[compID])
        write_80(famString)
        compID += 1
        for element in conn[key]:
            element.insert(1, key)
            write_bulk_line("CQUAD4", element)

    # Write boundary conditions
    for node in bcnodes:
        write_bulk_line("SPC", [1, node, "123456", 0.0])

    write_80("ENDDATA")
