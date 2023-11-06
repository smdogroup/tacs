from mpi4py import MPI

# make sure your system has procs 1 - N using mpirun command (not all rank zero)
# ended up uninstalling pip mpi4py, conda uninstalling mpi4py and then installing mpi4py through conda-forge channel
# this comes with its own mpirun binary located in your conda-env/bin folder (then I alias that mpirun as conda-mpirun)
# so that my OS mpirun doesn't get messed up. -Sean Engelstad
# command is then conda-mpirun -n 4 python _parallel_test.py (using conda-mpirun alias)
comm = MPI.COMM_WORLD
print(f"rank = {comm.rank}")
