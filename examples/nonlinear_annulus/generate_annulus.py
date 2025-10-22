from __future__ import print_function
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--ne1", type=int, default=10, help="# elements in radial direction"
)
parser.add_argument(
    "--ne2", type=int, default=80, help="# elements incircumferential direction"
)
args = parser.parse_args()

# Overall dimensions
Ri = 6.0
Ro = 10.0


# Number of elements along each edge of a single panel
ne1 = args.ne1
ne2 = args.ne2

np1 = 1
np2 = 1

# Nodes
n1 = ne1 + 1
n2 = ne2 + 1
r = np.linspace(Ri, Ro, n1)
theta = np.linspace(0, 2 * np.pi, n2)
R, Theta = np.meshgrid(r, theta)
X = np.zeros_like(R)
Y = np.zeros_like(R)
Z = np.zeros_like(R)

# Node numbering
nid = np.zeros((n2, n1), dtype="intc")
bcnodes = []
count = 1
for i in range(n2):
    for j in range(n1):
        nid[i, j] = count
        # Compute node location on annulus
        R = r[j]
        t = theta[i]
        X[i, j] = R * np.cos(t)
        Y[i, j] = -R * np.sin(t)
        Z[i, j] = 0.0
        if i == 0:
            bcnodes.append({"nodenum": count, "fixedDOF": "123456"})
        count += 1
nodes = np.stack((X, Y, Z), axis=2)
nmat = nodes.reshape((n1 * n2, 3))

# Connectivity
ne1 = n1 - 1
ne2 = n2 - 1
ne = ne1 * ne2
ncomp = 1
conn = {i + 1: [] for i in range(ncomp)}
ie = 1
for i in range(ne2):
    for j in range(ne1):
        compID = i // ne2 * np1 + j // ne1 + 1
        conn[compID].append(
            [ie, nid[i, j], nid[i + 1, j], nid[i + 1, j + 1], nid[i, j + 1]]
        )
        ie += 1


# Write BDF
output_file = "annulus.bdf"
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
