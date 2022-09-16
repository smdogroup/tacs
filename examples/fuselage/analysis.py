"""
Fuselage section example provided in AFEM library.
The fuselage is subjected to a cabin pressure of 1 atm and 1G gravity load.
User can run this example using aluminum or composite material properties.
"""
# ==============================================================================
# Standard Python modules
# ==============================================================================
from __future__ import print_function
import os
import argparse

# ==============================================================================
# External Python modules
# ==============================================================================
from mpi4py import MPI

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import functions, pyTACS

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--composite", action="store_true", default=False)
    args = parser.parse_args()
    use_composite = args.composite
else:
    use_composite = False

comm = MPI.COMM_WORLD

# Instantiate FEAAssembler
structOptions = {
    "printtiming": True,
}

bdfFile = os.path.join(os.path.dirname(__file__), "fuselage.bdf")
FEAAssembler = pyTACS(bdfFile, options=structOptions, comm=comm)

# Load a composite or isotropic elemCallBack depending on user input
if use_composite:
    from composite_setup import elemCallBack
else:
    from aluminum_setup import elemCallBack

# Set up elements and TACS assembler
FEAAssembler.initialize(elemCallBack)

# Setup static problem
problem = FEAAssembler.createStaticProblem("pressurization")

# Add loads to problem
# Apply cabin pressure to fuselage skins
compIDs = FEAAssembler.selectCompIDs("SKIN")
p = 101.3e3  # Pa
problem.addPressureToComponents(compIDs, p)
# Apply gravity load to fuselage
problem.addInertialLoad([0.0, -9.81, 0.0])

# Add eval functions to problem
problem.addFunction("mass", functions.StructuralMass)
problem.addFunction("ks_failure", functions.KSFailure, ksWeight=100.0, safetyFactor=1.5)

# Solve problem
problem.solve()

# Write solution
problem.writeSolution(outputDir=os.path.dirname(__file__))

# Evaluate functions
funcs = {}
problem.evalFunctions(funcs)

if MPI.COMM_WORLD.rank == 0:
    print(funcs)
