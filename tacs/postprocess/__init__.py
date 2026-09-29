"""
Post-processing utilities for TACS ``.f5`` solution files.

Currently provides shear / bending-moment / torque (VMT) diagrams for
wing-like structures, for single load cases and as an envelope over several
cases, see :mod:`tacs.postprocess.vmt`. matplotlib is only imported when a
plot is requested.
"""

from .vmt import (
    F5Data,
    Silhouette,
    VMTEnvelope,
    VMTResult,
    buildSilhouette,
    computeEnvelope,
    computeVMT,
    computeVMTCases,
    computeVMTFromF5,
    listComponents,
    loadF5,
    nodalCoordinates,
    nodalLoads,
    plotVMT,
    plotVMTEnvelope,
    plotVMTReport,
    scaleVMTResult,
    selectComponents,
    writeVMTCsv,
    writeVMTEnvelopeCsv,
)

__all__ = [
    "F5Data",
    "Silhouette",
    "VMTEnvelope",
    "VMTResult",
    "buildSilhouette",
    "computeEnvelope",
    "computeVMT",
    "computeVMTCases",
    "computeVMTFromF5",
    "listComponents",
    "loadF5",
    "nodalCoordinates",
    "nodalLoads",
    "plotVMT",
    "plotVMTEnvelope",
    "plotVMTReport",
    "scaleVMTResult",
    "selectComponents",
    "writeVMTCsv",
    "writeVMTEnvelopeCsv",
]
