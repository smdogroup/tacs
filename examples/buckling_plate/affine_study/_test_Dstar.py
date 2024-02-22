from tacs import buckling_surrogate
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD

"""
test the range of Dstar values I can get by rotating these plies
"""

_materials = buckling_surrogate.FlatPlateAnalysis.get_materials()
print(f"materials = {_materials}")
for material in _materials:
    for angle in np.linspace(0.0, 90.0, 10):
        comp_plate = material(
            comm=comm,
            bdf_file="plate.bdf",
            a=1.0,
            b=1.0,
            h=0.1,
            ply_angle=angle,
        )

        print(
            f"mat = {comp_plate.material_name}, angle = {angle}, Dstar = {comp_plate.Dstar}"
        )
