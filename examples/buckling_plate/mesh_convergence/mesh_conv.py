from tacs import buckling_surrogate
import numpy as np, os, time
from mpi4py import MPI
import pandas as pd

comm = MPI.COMM_WORLD

nelem_list = []
eig_value_list = []

# remove the previous version of the csv file
_clear_csv = False
if _clear_csv and comm.rank == 0:
    if os.path.exists("mesh_convergence.csv"):
        os.remove("mesh_convergence.csv")

ct = 0
for AR in [5.0]:  # [1.0, 3.0, 5.0]
    if AR == 1.0:
        max_scale = 1.3
    elif AR == 3.0:
        max_scale = 1.15
    elif AR == 5.0:
        max_scale = 1.05

    # max # elems about 150x150 grid or 22,500 elements

    for log10_scale in np.linspace(0.0, max_scale, 10):
        scale = 10**log10_scale
        _start_time = time.time()

        flat_plate = buckling_surrogate.FlatPlateAnalysis(
            comm=comm,
            bdf_file="plate.bdf",
            a=AR,
            b=1.0,
            h=0.01,
            E11=82.14e9,
            nu12=0.1487,
        )

        flat_plate.generate_bdf(
            nx=int(5 * AR * scale),
            ny=int(5 * scale),
            exx=flat_plate.affine_exx,
            eyy=0.0,
            exy=0.0,
            clamped=False,
        )

        # avg_stresses = flat_plate.run_static_analysis(write_soln=True)

        tacs_eigvals, errors = flat_plate.run_buckling_analysis(
            sigma=5.0, num_eig=12, write_soln=False
        )

        nelem_list += [flat_plate.num_elements]
        eigvals = [np.real(_) for _ in tacs_eigvals[:5]]
        eig_value_list += [eigvals]

        _dt = time.time() - _start_time

        data_dict = {
            "AR": [AR],
            "nx": [flat_plate._nx],
            "ny": [flat_plate._ny],
            "nelem": [flat_plate.num_elements],
            "dt[s]": [_dt],
        }

        del flat_plate

        # write to a csv file
        if comm.rank == 0 and abs(errors[0]) < 1e-9:
            for i in range(5):
                data_dict[f"lam{i}"] = [eigvals[i]]
            df = pd.DataFrame(data_dict)
            if ct == 0 and _clear_csv:
                df.to_csv("mesh_convergence.csv", index=False)
            else:
                df.to_csv("mesh_convergence.csv", mode="a", index=False, header=False)

            ct += 1
