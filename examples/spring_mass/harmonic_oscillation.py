"""
Same spring-mass case as analysis.py.
This time we first apply an oscilating force on the mass as a baseline case.
The chosen frequency ends up being far away from the systems natural frequency.
We then create a modal analysis problem to find the system's natural frequencies.
We then create a second transient problem and use the lowest natural frequency
and eigenmode found from the modal problem to drive the new transient problem at
the model's first natural frequency. We then plot the respone of both
the baseline and resonating transient problem.
"""
# ==============================================================================
# Standard Python modules
# ==============================================================================
from __future__ import print_function
import os
import matplotlib.pyplot as plt

# ==============================================================================
# External Python modules
# ==============================================================================
from pprint import pprint
from mpi4py import MPI
import numpy as np

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import pyTACS, functions, constitutive, elements

# Beginning and end time of transient problem
tStart = 0.0  # s
tEnd = 4 * np.pi  # s
# Number of steps for transient problem
nSteps = 100

# Force vector to apply to mass
Fmax = 1.0  # N
freqBaseline = 1.5  # rad/s

comm = MPI.COMM_WORLD

# Instantiate FEAAssembler
structOptions = {
    "printtiming": True,
}

bdfFile = os.path.join(os.path.dirname(__file__), "mass_spring.bdf")
FEAAssembler = pyTACS(bdfFile, comm, options=structOptions)


def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    # Setup spring element
    # NOTE: Mass elements are always setup automatically for users outside of elemCallBack
    k = [1, 4, 9, 16, 25, 36]
    con = constitutive.DOFSpringConstitutive(k=k)
    # Stiffness components are in global coordinate system
    # We don't need a transform
    transform = None
    elem = elements.SpringElement(transform, con)
    return elem


# Set up TACS Assembler
FEAAssembler.initialize(elemCallBack)

# create baseline transient problems
baselineProb = FEAAssembler.createTransientProblem("baseline", tStart, tEnd, nSteps)

# Apply sinusoidal force at each time step
timeSteps = baselineProb.getTimeSteps()
for step_i, time in enumerate(timeSteps):
    force = FEAAssembler.createVec()
    force[:] = Fmax * np.sin(freqBaseline * time)
    baselineProb.addLoadToRHS(step_i, force)

# Solve problems
baselineProb.solve()

# Pick off state variable values at each time step
nnodes = FEAAssembler.getNumOwnedNodes()
vpn = FEAAssembler.getVarsPerNode()
stateHistory = np.zeros([len(timeSteps), nnodes, vpn])
for step_i in range(len(timeSteps)):
    baselineProb.getVariables(step_i, states=stateHistory[step_i, :, :].reshape(-1))

# Create a modal problem to find natural frequencies of system
modalProb = FEAAssembler.createModalProblem("freq_analysis", sigma=0.8, numEigs=6)
# Solve the modal problem
modalProb.solve()
# Print out each found eigenfrequency
numFreqs = modalProb.getNumEigs()
for i in range(numFreqs):
    eigVal, _ = modalProb.getVariables(i)
    # Frequency is the sqrt of the eigenvalue
    freq = np.sqrt(eigVal)
    print(f"Mode {i+1}:")
    print(f"Frequency: {freq} (rad/s)")
    print(" ")

# Pull out eigen vector information for first mode
eigVal0, eigVec0 = modalProb.getVariables(0)
freq0 = np.sqrt(eigVal0)

# Create a new transient problem where we force the model at resonance
resonanceProb = FEAAssembler.createTransientProblem("resonance", tStart, tEnd, nSteps)
# Apply sinusoidal force at each time step
timeSteps = resonanceProb.getTimeSteps()
for step_i, time in enumerate(timeSteps):
    force = Fmax * np.sin(freq0 * time) * eigVec0
    resonanceProb.addLoadToRHS(step_i, force)

# Solve the problem
resonanceProb.solve()

# Pick off state variable values at each time step
newStateHistory = np.zeros([len(timeSteps), nnodes, vpn])
for step_i in range(len(timeSteps)):
    resonanceProb.getVariables(step_i, states=newStateHistory[step_i, :, :].reshape(-1))

# Plot results for first 3 dofs
plt.plot(
    timeSteps,
    stateHistory[:, 1, 0],
    "r-",
    timeSteps,
    stateHistory[:, 1, 1],
    "g-",
    timeSteps,
    stateHistory[:, 1, 2],
    "b-",
    timeSteps,
    newStateHistory[:, 1, 0],
    "r--",
    timeSteps,
    newStateHistory[:, 1, 1],
    "g--",
    timeSteps,
    newStateHistory[:, 1, 2],
    "b--",
)
plt.legend(
    [
        "dof 1 (baseline)",
        "dof 2 (baseline)",
        "dof 3 (baseline)",
        "dof 1 (resonance)",
        "dof 2 (resonance)",
        "dof 3 (resonance)",
    ]
)
plt.ylabel("displacement (m)")
plt.xlabel("time (s)")
# Show transient plot
plt.show()
