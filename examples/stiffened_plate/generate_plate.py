from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pprint import pprint
from scipy.spatial.transform import Rotation

parser = argparse.ArgumentParser()
parser.add_argument("--nex", type=int, default=20, help="# elements from top to bottom")
parser.add_argument(
    "--ney", type=int, default=40, help="# elements around circumference"
)
parser.add_argument("--pinEdges", action="store_true", help="Pin edges of plate")
args = parser.parse_args()

# Overall dimensions
width = 0.75
length = 1.5

# Number of elements along each edge of a single panel
nex = args.nex
ney = args.ney

np1 = 1
np2 = 1

# Nodes
n1 = nex + 1
n2 = ney + 1
x = np.linspace(0, width, n1)
y = np.linspace(0, length, n2)
X, Y = np.meshgrid(x, y)
Z = np.zeros_like(X)

# Node numbering
nid = np.zeros((n2, n1), dtype="intc")
bcnodes = []
count = 1
for i in range(n2):
    isYmin = i == 0
    isYmax = i == n2 - 1
    for j in range(n1):
        isXmin = j == 0
        isXmax = j == n1 - 1
        nid[i, j] = count
        # Fix bottom left corner in x,y, z
        if isYmin and isXmin:
            bcnodes.append({"nodenum": count, "fixedDOF": "123"})
        # fix bottom right corner in y
        elif isYmin and isXmax:
            bcnodes.append({"nodenum": count, "fixedDOF": "23"})
        # Fix top left corner in z
        elif isYmax and isXmin:
            bcnodes.append({"nodenum": count, "fixedDOF": "3"})
        # Fix all other edge nodes in z
        elif args.pinEdges and (isYmax or isYmin or isXmin or isXmax):
            bcnodes.append({"nodenum": count, "fixedDOF": "3"})
        count += 1
nodes = np.stack((X, Y, Z), axis=2)
nmat = nodes.reshape((n1 * n2, 3))

# Connectivity
nex = n1 - 1
ney = n2 - 1
ne = nex * ney
ncomp = 1
conn = {i + 1: [] for i in range(ncomp)}
ie = 1
for i in range(ney):
    for j in range(nex):
        compID = i // ney * np1 + j // nex + 1
        conn[compID].append(
            [ie, nid[i, j], nid[i + 1, j], nid[i + 1, j + 1], nid[i, j + 1]]
        )
        ie += 1


# Write BDF
output_file = "plate.bdf" if not args.pinEdges else "plate_pinned_edges.bdf"
fout = open(output_file, "w")


def write_80(line):
    newline = "{:80s}\n".format(line.strip("\n"))
    fout.write(newline)


write_80("SOL 103")
write_80("CEND")
write_80("BEGIN BULK")

# Make component names
compNames = {}
compID = 1
for i in range(np2):
    for j in range(np1):
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
        if type(item) in [int, np.int64, np.int32]:
            line += "{:{width}d}".format(item, width=width)[:width]
        elif type(item) in [float, np.float64]:
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
for i in range(n1 * n2):
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
    write_bulk_line("SPC", [1, node["nodenum"], node["fixedDOF"], 0.0])

write_80("ENDDATA")

fout.close()
