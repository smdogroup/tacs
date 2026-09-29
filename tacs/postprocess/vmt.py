r"""
Shear, bending-moment and torque (VMT) diagrams from TACS ``.f5`` files.

This module reduces the nodal load field stored in a TACS ``.f5`` solution
file to the classical beam quantities used to sanity-check a wing-like
structure: shear force ``V``, bending moment ``M`` and torque ``T`` along a
user-defined, piecewise-linear beam axis.

Method
------
Every node is projected onto the beam axis polyline, which gives it an
arc-length coordinate ``s``. At each spanwise station the applied nodal loads
(``fx, fy, fz, mx, my, mz``) and, when present, the support reactions
(``rfx ... rmz``) of all nodes *outboard* of the station are summed:

.. math::

    \mathbf{F}(s) = \sum_{s_i > s} \mathbf{f}_i, \qquad
    \mathbf{Mom}(s) = \sum_{s_i > s} \left[\mathbf{m}_i +
    (\mathbf{x}_i - \mathbf{x}_s) \times \mathbf{f}_i\right]

where :math:`\mathbf{x}_s` is the station point on the axis. The resultants
are then resolved in a local frame made of the axis tangent :math:`\hat{t}`
(root to tip), the shear direction :math:`\hat{s}` (the user-supplied shear
direction with its component along :math:`\hat{t}` removed) and the bending
axis :math:`\hat{b} = \hat{s} \times \hat{t}`:

.. math::

    V = \mathbf{F} \cdot \hat{s}, \qquad
    M = \mathbf{Mom} \cdot \hat{b}, \qquad
    T = \mathbf{Mom} \cdot \hat{t}

With this convention :math:`(\hat{t}, \hat{b}, \hat{s})` is right-handed and
``dM/ds = V`` holds exactly. Example: a cantilever along ``+x`` with a tip load
``P`` in ``+z`` and the default shear direction ``(0, 0, -1)`` gives
``V = -P`` and ``M = +P (L - s)``.

Notes
-----
* The loads written to the f5 file (pyTACS option ``writeLoads``) are the
  total external nodal loads the solution is in equilibrium with, including
  the consistent nodal loads of pressure and inertial loads. They are only
  meaningful after the problem has been solved.
* Reactions (``writeReactions``) are genuine external forces acting on the
  structure, so summing loads and reactions over the whole model gives zero.
  Nodes at or inboard of the first axis point project to ``s = 0`` exactly and
  are therefore excluded from the root station, which reports the total
  applied load rather than zero.
* A station that coincides with a ring of nodes treats that ring as inboard
  (the reported value is the outboard limit). Use ``stationMode="midpoint"``
  or explicit ``stations`` if the inboard limit is wanted.
* The f5 writer zeroes the load rows of constrained degrees of freedom, and
  the reactions balance only the remaining loads. External load applied
  directly at constrained nodes (for example the inertial load lumped on the
  clamped root nodes) is therefore not recoverable from the file. For a
  cantilever root this only affects the load share of the root node ring.
* Coordinates and loads are stored as ``float32`` in the f5 file; all sums are
  performed in ``float64``.

Examples
--------
Python::

    from tacs.postprocess import computeVMTFromF5, plotVMT, loadF5

    data = loadF5("pullup_000.f5")
    result = computeVMTFromF5(
        data, axisPts=[[4.0, 0.0, 0.0], [4.0, 0.0, 13.8]], shearDir=[0, -1, 0]
    )
    plotVMT(result, data, fileName="pullup_vmt.pdf")

Several load cases with an envelope, written as a multi-page PDF::

    from tacs.postprocess import computeVMTCases, plotVMTReport

    results = computeVMTCases(
        ["pullup_000.f5", "pushover_000.f5", "gust_000.f5"],
        axisPts=[[4.0, 0.0, 0.0], [4.0, 0.0, 13.8]],
        shearDir=[0, -1, 0],
    )
    plotVMTReport(results, "vmt_report.pdf", data=data)

Command line::

    python -m tacs.postprocess.vmt pullup_000.f5 pushover_000.f5 --axis 4 0 0  4 0 13.8 \
        --shear-dir 0 -1 0 --num-stations 30 --output vmt_report.pdf

Label the plot axes with units, or convert the plotted numbers first (here
from N and m to kN and mm)::

    python -m tacs.postprocess.vmt pullup_000.f5 --axis 4 0 0  4 0 13.8 --shear-dir 0 -1 0 \
        --output vmt.pdf --force-unit N --length-unit m
    python -m tacs.postprocess.vmt pullup_000.f5 --axis 4 0 0  4 0 13.8 --shear-dir 0 -1 0 \
        --output vmt.pdf --force-scale 1e-3 --length-scale 1e3 --force-unit kN --length-unit mm
"""

import argparse
import fnmatch
import os
import struct
import sys
import warnings
from dataclasses import dataclass, replace

import numpy as np

# ---------------------------------------------------------------------------
# Element layout constants (mirror of the ElementLayout enum, see
# src/elements/TACSElementTypes.h). Kept here so the numpy kernel and the
# silhouette builder do not need the compiled extension at import time.
# ---------------------------------------------------------------------------
LAYOUT_NONE = 0
POINT_ELEMENT = 1
LINE_ELEMENT = 2
LINE_QUADRATIC_ELEMENT = 3
LINE_CUBIC_ELEMENT = 4
TRI_ELEMENT = 5
TRI_QUADRATIC_ELEMENT = 6
TRI_CUBIC_ELEMENT = 7
QUAD_ELEMENT = 8
QUAD_QUADRATIC_ELEMENT = 9
QUAD_CUBIC_ELEMENT = 10
QUAD_QUARTIC_ELEMENT = 11
QUAD_QUINTIC_ELEMENT = 12
RBE2_ELEMENT = 24
RBE3_ELEMENT = 25

_LINE_LAYOUTS = (LINE_ELEMENT, LINE_QUADRATIC_ELEMENT, LINE_CUBIC_ELEMENT)
_TRI_LAYOUTS = (TRI_ELEMENT, TRI_QUADRATIC_ELEMENT, TRI_CUBIC_ELEMENT)
_QUAD_LAYOUTS = (
    QUAD_ELEMENT,
    QUAD_QUADRATIC_ELEMENT,
    QUAD_CUBIC_ELEMENT,
    QUAD_QUARTIC_ELEMENT,
    QUAD_QUINTIC_ELEMENT,
)
_SOLID_LAYOUTS = tuple(range(13, 24))

# f5 binary format constants (src/io/TACSFH5.h / TACSFH5.cpp)
_F5_DTYPE_ITEMSIZE = {0: 4, 1: 8, 2: 4}  # FH5_INT, FH5_DOUBLE, FH5_FLOAT
_F5_MAX_NAME_LEN = 1 << 20
_F5_REQUIRED_ZONES = ("components", "ltypes", "ptr", "connectivity")
_F5_CONTINUOUS_PREFIX = "continuous data"

_COORD_NAMES = ("X", "Y", "Z")
_FORCE_NAMES = ("fx", "fy", "fz")
_MOMENT_NAMES = ("mx", "my", "mz")
_RFORCE_NAMES = ("rfx", "rfy", "rfz")
_RMOMENT_NAMES = ("rmx", "rmy", "rmz")


# ---------------------------------------------------------------------------
# Layer 1: f5 file -> numpy
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class F5Zone:
    """Header information of one zone in an f5 file."""

    name: str
    varNames: str
    dtype: int
    dim1: int
    dim2: int


@dataclass(frozen=True)
class F5Data:
    """
    Nodal data and connectivity read from a TACS f5 file.

    Attributes
    ----------
    fileName : str
        Path of the file the data was read from.
    varNames : tuple of str
        Names of the continuous (nodal) variables, one per column of ``data``.
    data : numpy.ndarray
        Continuous data, shape ``(numNodes, numVars)``, ``float64``.
    comps : numpy.ndarray
        Component id of every element, shape ``(numElements,)``.
    ltypes : numpy.ndarray
        Element layout id (``ElementLayout`` enum) of every element.
    ptr : numpy.ndarray
        CSR offsets into ``conn``, shape ``(numElements + 1,)``.
    conn : numpy.ndarray
        Flattened element connectivity, indexing rows of ``data``.
    componentNames : tuple of str
        Name of every component, indexed by the values in ``comps``.
    """

    fileName: str
    varNames: tuple
    data: np.ndarray
    comps: np.ndarray
    ltypes: np.ndarray
    ptr: np.ndarray
    conn: np.ndarray
    componentNames: tuple

    @property
    def numNodes(self):
        """Number of nodes (rows of ``data``)."""
        return self.data.shape[0]

    @property
    def numElements(self):
        """Number of elements."""
        return self.comps.shape[0]

    def has(self, *names):
        """Return True if every named variable is present."""
        return all(name in self.varNames for name in names)

    def column(self, name):
        """Return one nodal variable as a 1-D array, or None if absent."""
        if name not in self.varNames:
            return None
        return self.data[:, self.varNames.index(name)]

    def columns(self, names):
        """Return several nodal variables as an ``(numNodes, len(names))`` array."""
        idx = [self.varNames.index(name) for name in names]
        return self.data[:, idx]


def _scanF5Header(fileName):
    """
    Walk the header of an f5 file in pure Python.

    The compiled loader does not report read failures and dereferences NULL
    pointers afterwards, so this scan is used to reject invalid or truncated
    files before the loader is touched.

    Parameters
    ----------
    fileName : str
        Path to the f5 file.

    Returns
    -------
    componentNames : list of str
    zones : list of F5Zone

    Raises
    ------
    ValueError
        If the file is not a structurally valid TACS f5 file.
    """
    fileSize = os.path.getsize(fileName)

    def bad(reason):
        return ValueError(f"'{fileName}' is not a valid TACS .f5 file ({reason})")

    with open(fileName, "rb") as fp:

        def readInts(count):
            raw = fp.read(4 * count)
            if len(raw) != 4 * count:
                raise bad("unexpected end of file in header")
            return struct.unpack(f"{count}i", raw)

        def readName(length):
            if not 0 < length <= _F5_MAX_NAME_LEN:
                raise bad(f"invalid name length {length}")
            raw = fp.read(length)
            if len(raw) != length:
                raise bad("unexpected end of file in name")
            return raw.split(b"\0", 1)[0].decode("utf-8", errors="replace")

        (numComp,) = readInts(1)
        if not 0 <= numComp <= _F5_MAX_NAME_LEN:
            raise bad(f"invalid component count {numComp}")
        componentNames = []
        for _ in range(numComp):
            (slen,) = readInts(1)
            componentNames.append(readName(slen))

        zones = []
        while fp.tell() < fileSize:
            if fileSize - fp.tell() < 20:
                raise bad("truncated zone header")
            dtype, dim1, dim2, lenZone, lenVars = readInts(5)
            if dtype not in _F5_DTYPE_ITEMSIZE:
                raise bad(f"unknown data type {dtype}")
            if dim1 < 0 or dim2 < 0:
                raise bad(f"negative zone dimensions {dim1} x {dim2}")
            zoneName = readName(lenZone)
            varNames = readName(lenVars)
            nbytes = dim1 * dim2 * _F5_DTYPE_ITEMSIZE[dtype]
            if fp.tell() + nbytes > fileSize:
                raise bad(f"zone '{zoneName}' extends past the end of the file")
            zones.append(F5Zone(zoneName, varNames, dtype, dim1, dim2))
            fp.seek(nbytes, os.SEEK_CUR)

    return componentNames, zones


def loadF5(fileName):
    """
    Read the nodal data and connectivity of a TACS f5 file.

    Parameters
    ----------
    fileName : str
        Path to the ``.f5`` file written by :class:`tacs.TACS.ToFH5` or
        pyTACS ``writeSolution``.

    Returns
    -------
    F5Data
        Connectivity and continuous nodal data as ``float64`` copies.

    Raises
    ------
    FileNotFoundError
        If ``fileName`` does not exist.
    ValueError
        If the file is not a valid f5 file, or was written without
        connectivity or continuous nodal data.
    """
    if not os.path.isfile(fileName):
        raise FileNotFoundError(f"f5 file '{fileName}' does not exist")

    componentNames, zones = _scanF5Header(fileName)
    zoneNames = [zone.name for zone in zones]
    for required in _F5_REQUIRED_ZONES:
        if required not in zoneNames:
            raise ValueError(
                f"'{fileName}' has no '{required}' zone; re-write the file with "
                "the pyTACS option writeConnectivity=True "
                "(TACS_OUTPUT_CONNECTIVITY)"
            )
    if not any(name.startswith(_F5_CONTINUOUS_PREFIX) for name in zoneNames):
        raise ValueError(
            f"'{fileName}' contains no continuous nodal data; re-write the file "
            "with the pyTACS options writeNodes=True and writeLoads=True"
        )

    # Only import the compiled extension once the file has been validated
    from tacs import TACS

    loader = TACS.FH5Loader()
    loader.loadData(fileName)
    comps, ltypes, ptr, conn = loader.getConnectivity()
    varString, fdata = loader.getContinuousData()

    # Copy out of the loader-owned float32 buffer
    data = np.array(fdata, dtype=np.float64)
    varNames = tuple(name.strip() for name in varString.split(","))
    numComp = loader.getNumComponents()
    names = tuple(loader.getComponentName(i) for i in range(numComp))

    return F5Data(
        fileName=fileName,
        varNames=varNames,
        data=data,
        comps=np.asarray(comps, dtype=np.intp),
        ltypes=np.asarray(ltypes, dtype=np.intp),
        ptr=np.asarray(ptr, dtype=np.intp),
        conn=np.asarray(conn, dtype=np.intp),
        componentNames=names,
    )


def nodalCoordinates(data):
    """
    Return the nodal coordinates ``(numNodes, 3)`` stored in an f5 file.

    Raises
    ------
    ValueError
        If the file was written without nodal coordinates.
    """
    if not data.has(*_COORD_NAMES):
        raise ValueError(
            f"'{data.fileName}' has no nodal coordinates; re-write the file with "
            "the pyTACS option writeNodes=True (TACS_OUTPUT_NODES)"
        )
    return np.array(data.columns(_COORD_NAMES))


def nodalLoads(data, includeReactions=None):
    """
    Return the external nodal forces and moments stored in an f5 file.

    Parameters
    ----------
    data : F5Data
        Data read with :func:`loadF5`.
    includeReactions : bool or None
        ``True`` adds the support reactions (``rfx ... rmz``) to the applied
        loads and raises if they are absent; ``False`` never adds them;
        ``None`` (default) adds them when present and warns otherwise.

    Returns
    -------
    F : numpy.ndarray
        Nodal forces, shape ``(numNodes, 3)``.
    M : numpy.ndarray
        Nodal moments, shape ``(numNodes, 3)`` (zeros if the element type
        has no rotational degrees of freedom).
    usedReactions : bool
        Whether reactions were added.

    Raises
    ------
    ValueError
        If the applied loads are missing, or reactions were requested
        explicitly but are missing.
    """
    if not data.has(*_FORCE_NAMES):
        raise ValueError(
            f"'{data.fileName}' has no nodal loads (fx, fy, fz); re-write the "
            "file with the pyTACS option writeLoads=True (TACS_OUTPUT_LOADS)"
        )
    F = np.array(data.columns(_FORCE_NAMES))
    if data.has(*_MOMENT_NAMES):
        M = np.array(data.columns(_MOMENT_NAMES))
    else:
        warnings.warn(
            f"'{data.fileName}' has no nodal moments (mx, my, mz); assuming zero",
            stacklevel=2,
        )
        M = np.zeros_like(F)

    hasReactions = data.has(*_RFORCE_NAMES)
    if includeReactions is None:
        usedReactions = hasReactions
        if not hasReactions:
            warnings.warn(
                f"'{data.fileName}' has no reactions (rfx, rfy, rfz); continuing "
                "with applied loads only. Write the file with writeReactions=True "
                "to include them",
                stacklevel=2,
            )
    elif includeReactions and not hasReactions:
        raise ValueError(
            f"'{data.fileName}' has no reactions (rfx, rfy, rfz); re-write the "
            "file with the pyTACS option writeReactions=True "
            "(TACS_OUTPUT_REACTIONS) or pass includeReactions=False"
        )
    else:
        usedReactions = bool(includeReactions)

    if usedReactions:
        F += data.columns(_RFORCE_NAMES)
        if data.has(*_RMOMENT_NAMES):
            M += data.columns(_RMOMENT_NAMES)

    if not (np.all(np.isfinite(F)) and np.all(np.isfinite(M))):
        raise ValueError(
            f"'{data.fileName}' contains non-finite nodal loads; the solution "
            "may not have converged"
        )
    return F, M, usedReactions


# ---------------------------------------------------------------------------
# Layer 2: component and node selection
# ---------------------------------------------------------------------------
def _elementNodeCounts(data):
    """Return the number of connectivity slots of every element."""
    return data.ptr[1:] - data.ptr[:-1]


def _loadCarryingSlotCounts(data):
    """
    Return, per element, the number of leading connectivity slots that hold
    physical nodes. RBE2/RBE3 elements append Lagrange-multiplier nodes that
    carry constraint residuals instead of loads; these are dropped.
    """
    counts = _elementNodeCounts(data).copy()
    isRBE2 = data.ltypes == RBE2_ELEMENT
    isRBE3 = data.ltypes == RBE3_ELEMENT
    counts[isRBE2] -= (counts[isRBE2] - 1) // 2
    counts[isRBE3] -= 1
    return counts


def _slotsOfElements(data, elemMask, slotCounts):
    """Return the flattened connectivity slot indices of the masked elements."""
    elems = np.flatnonzero(elemMask)
    starts = data.ptr[elems]
    lens = slotCounts[elems]
    total = int(lens.sum())
    if total == 0:
        return np.zeros(0, dtype=np.intp)
    offsets = np.arange(total) - np.repeat(np.cumsum(lens) - lens, lens)
    return np.repeat(starts, lens) + offsets


def listComponents(data):
    """
    Return ``[(componentName, numElements), ...]`` for an f5 data set.

    Parameters
    ----------
    data : F5Data or str
        Loaded data or the path of an f5 file.
    """
    if isinstance(data, str):
        data = loadF5(data)
    counts = np.bincount(data.comps, minlength=len(data.componentNames))
    return [(name, int(counts[i])) for i, name in enumerate(data.componentNames)]


def selectComponents(data, components=None):
    """
    Build element and node masks for a subset of components.

    Parameters
    ----------
    data : F5Data
        Data read with :func:`loadF5`.
    components : None, str, int or sequence of str/int
        Components to keep. ``None`` keeps everything. Strings are matched
        case-insensitively, first exactly and then as ``fnmatch`` glob
        patterns (``"WING_U_SKIN*"``). Integers are component ids.

    Returns
    -------
    elemMask : numpy.ndarray
        Boolean mask over elements.
    nodeMask : numpy.ndarray
        Boolean mask over nodes. A node is selected if it is a load-carrying
        node of any selected element, whatever the element type. The
        Lagrange-multiplier slots of RBE2/RBE3 elements are excluded.

    Raises
    ------
    ValueError
        If a component name or id matches nothing.
    """
    names = data.componentNames
    if components is None:
        elemMask = np.ones(data.numElements, dtype=bool)
    else:
        if isinstance(components, (str, int, np.integer)):
            components = [components]
        lowered = [name.lower() for name in names]
        selected = set()
        for item in components:
            if isinstance(item, (int, np.integer)):
                if not 0 <= int(item) < len(names):
                    raise ValueError(
                        f"component id {item} out of range [0, {len(names) - 1}]"
                    )
                selected.add(int(item))
                continue
            pattern = str(item).lower()
            matches = [i for i, name in enumerate(lowered) if name == pattern]
            if not matches:
                matches = [
                    i
                    for i, name in enumerate(lowered)
                    if fnmatch.fnmatchcase(name, pattern)
                ]
            if not matches:
                preview = ", ".join(names[:20])
                more = ", ..." if len(names) > 20 else ""
                raise ValueError(
                    f"component '{item}' matched none of the {len(names)} "
                    f"components in '{data.fileName}': {preview}{more}"
                )
            selected.update(matches)
        elemMask = np.isin(data.comps, np.fromiter(selected, dtype=np.intp))

    slots = _slotsOfElements(data, elemMask, _loadCarryingSlotCounts(data))
    nodeMask = np.zeros(data.numNodes, dtype=bool)
    nodeMask[data.conn[slots]] = True
    if not nodeMask.any():
        raise ValueError("no nodes selected")
    return elemMask, nodeMask


# ---------------------------------------------------------------------------
# Layer 3: pure-numpy VMT kernel
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class VMTResult:
    """
    Shear, bending moment and torque at spanwise stations.

    Attributes
    ----------
    s : numpy.ndarray
        Arc length of every station along the beam axis, shape ``(N,)``.
    V : numpy.ndarray
        Shear force ``F . shearAxis`` at every station.
    M : numpy.ndarray
        Bending moment ``Mom . bendingAxis`` at every station.
    T : numpy.ndarray
        Torque ``Mom . tangent`` at every station.
    force : numpy.ndarray
        Resultant force of all nodes outboard of the station, ``(N, 3)``.
    moment : numpy.ndarray
        Resultant moment about the station point, ``(N, 3)``.
    stationPoints, tangent, bendingAxis, shearAxis : numpy.ndarray
        Station point and local unit frame, each ``(N, 3)``.
    nodeS : numpy.ndarray
        Arc-length coordinate of every node passed to the kernel.
    axisPoints : numpy.ndarray
        Beam axis vertices after removing zero-length segments.
    axisLength : float
        Total arc length of the beam axis.
    includeReactions : bool or None
        Whether reactions were part of the summed loads (``None`` when the
        kernel was called directly with user arrays).
    """

    s: np.ndarray
    V: np.ndarray
    M: np.ndarray
    T: np.ndarray
    force: np.ndarray
    moment: np.ndarray
    stationPoints: np.ndarray
    tangent: np.ndarray
    bendingAxis: np.ndarray
    shearAxis: np.ndarray
    nodeS: np.ndarray
    axisPoints: np.ndarray
    axisLength: float
    includeReactions: object = None

    def asDict(self):
        """Return the station arrays as a dict keyed by name."""
        return {
            "s": self.s,
            "V": self.V,
            "M": self.M,
            "T": self.T,
            "Fx": self.force[:, 0],
            "Fy": self.force[:, 1],
            "Fz": self.force[:, 2],
            "Mx": self.moment[:, 0],
            "My": self.moment[:, 1],
            "Mz": self.moment[:, 2],
            "x": self.stationPoints[:, 0],
            "y": self.stationPoints[:, 1],
            "z": self.stationPoints[:, 2],
        }


def _prepareAxis(axisPts):
    """Validate the axis polyline and return its segment data."""
    pts = np.asarray(axisPts, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3 or pts.shape[0] < 2:
        raise ValueError(
            "axisPts must be an (npts, 3) array with at least two points, got "
            f"shape {pts.shape}"
        )
    if not np.all(np.isfinite(pts)):
        raise ValueError("axisPts contains non-finite values")
    d = np.diff(pts, axis=0)
    segLen = np.linalg.norm(d, axis=1)
    if segLen.max() <= 0.0:
        raise ValueError("all axis points coincide; the beam axis has zero length")
    keep = segLen > 1e-12 * segLen.max()
    starts = pts[:-1][keep]
    d = d[keep]
    segLen = segLen[keep]
    cleanPts = np.vstack([starts, pts[-1:]])
    tangents = d / segLen[:, None]
    sAxis = np.concatenate([[0.0], np.cumsum(segLen)])
    return starts, d, segLen, tangents, sAxis, cleanPts


def _projectToPolyline(X, starts, d, segLen, sAxis):
    """Return the arc-length coordinate of the closest axis point per node."""
    bestD2 = np.full(X.shape[0], np.inf)
    bestS = np.zeros(X.shape[0])
    for j in range(len(segLen)):
        w = X - starts[j]
        tpar = np.clip((w @ d[j]) / (segLen[j] * segLen[j]), 0.0, 1.0)
        diff = w - tpar[:, None] * d[j]
        d2 = np.einsum("ij,ij->i", diff, diff)
        update = d2 < bestD2
        bestD2 = np.where(update, d2, bestD2)
        bestS = np.where(update, sAxis[j] + tpar * segLen[j], bestS)
    return bestS


def _makeStations(axisLength, numStations, stations, stationMode):
    """Return the sorted station arc lengths."""
    if stations is not None:
        s = np.asarray(stations, dtype=np.float64).ravel()
        if s.size == 0 or not np.all(np.isfinite(s)):
            raise ValueError("stations must be a non-empty array of finite values")
        if s.min() < 0.0 or s.max() > axisLength:
            raise ValueError(
                f"stations must lie within [0, {axisLength}] (the axis arc length)"
            )
        return np.sort(s)
    numStations = int(numStations)
    if numStations < 1:
        raise ValueError("numStations must be at least 1")
    if stationMode == "endpoint":
        if numStations == 1:
            return np.zeros(1)
        return np.linspace(0.0, axisLength, numStations)
    if stationMode == "midpoint":
        return (np.arange(numStations) + 0.5) * axisLength / numStations
    raise ValueError(
        f"stationMode must be 'endpoint' or 'midpoint', got '{stationMode}'"
    )


def _localFrame(shearDir, tangent, segmentIndex):
    """Return the unit shear and bending axes for one station."""
    sh = shearDir - np.dot(shearDir, tangent) * tangent
    norm = np.linalg.norm(sh)
    if norm <= 1e-8 * np.linalg.norm(shearDir):
        raise ValueError(
            f"shearDir {shearDir.tolist()} is parallel to axis segment "
            f"{segmentIndex}; choose a shear direction perpendicular to the axis"
        )
    sh = sh / norm
    b = np.cross(sh, tangent)
    return sh, b


def computeVMT(
    X,
    F,
    axisPts,
    M=None,
    shearDir=(0.0, 0.0, -1.0),
    numStations=20,
    stations=None,
    nodeMask=None,
    stationMode="endpoint",
    stationTol=1e-6,
):
    """
    Compute shear, bending moment and torque along a piecewise-linear axis.

    Parameters
    ----------
    X : array_like
        Nodal coordinates, shape ``(n, 3)``.
    F : array_like
        External nodal forces, shape ``(n, 3)``.
    axisPts : array_like
        Beam axis vertices from root to tip, shape ``(npts, 3)``, ``npts >= 2``.
    M : array_like, optional
        External nodal moments, shape ``(n, 3)``. Defaults to zero.
    shearDir : array_like
        Shear direction. Its component along the local axis tangent is
        removed at every station. Default ``(0, 0, -1)``.
    numStations : int
        Number of stations when ``stations`` is not given. Default 20.
    stations : array_like, optional
        Explicit station arc lengths in ``[0, axisLength]``; overrides
        ``numStations`` and ``stationMode``.
    nodeMask : array_like of bool, optional
        Only nodes where the mask is True are summed.
    stationMode : {"endpoint", "midpoint"}
        ``"endpoint"`` spaces stations uniformly including the root and tip;
        ``"midpoint"`` places them at the centres of ``numStations`` equal bins.
    stationTol : float
        A node contributes to a station when ``sNode > s + stationTol * scale``
        with ``scale = max(axisLength, max|X|)``. Absorbs the ``float32``
        round-off of coordinates stored in f5 files.

    Returns
    -------
    VMTResult

    Raises
    ------
    ValueError
        On inconsistent shapes, a degenerate axis, a shear direction parallel
        to an axis segment, or non-finite input.
    """
    X = np.asarray(X, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)
    M = np.zeros_like(F) if M is None else np.asarray(M, dtype=np.float64)
    if X.ndim != 2 or X.shape[1] != 3:
        raise ValueError(f"X must have shape (n, 3), got {X.shape}")
    if F.shape != X.shape or M.shape != X.shape:
        raise ValueError(
            f"F and M must have the same shape as X {X.shape}, got {F.shape} and {M.shape}"
        )
    if nodeMask is not None:
        nodeMask = np.asarray(nodeMask, dtype=bool)
        if nodeMask.shape != (X.shape[0],):
            raise ValueError(
                f"nodeMask must have shape ({X.shape[0]},), got {nodeMask.shape}"
            )
        X, F, M = X[nodeMask], F[nodeMask], M[nodeMask]
    if X.shape[0] == 0:
        raise ValueError("no nodes selected")
    if not (
        np.all(np.isfinite(X)) and np.all(np.isfinite(F)) and np.all(np.isfinite(M))
    ):
        raise ValueError("X, F and M must be finite")

    shearDir = np.asarray(shearDir, dtype=np.float64).ravel()
    if shearDir.shape != (3,) or not np.all(np.isfinite(shearDir)):
        raise ValueError("shearDir must be a finite 3-vector")
    if np.linalg.norm(shearDir) == 0.0:
        raise ValueError("shearDir must be non-zero")

    starts, d, segLen, tangents, sAxis, cleanPts = _prepareAxis(axisPts)
    axisLength = float(sAxis[-1])
    nodeS = _projectToPolyline(X, starts, d, segLen, sAxis)
    if X.shape[0] > 1 and np.ptp(nodeS) < 1e-6 * axisLength:
        warnings.warn(
            "all nodes project to the same point of the beam axis; check that the "
            "axis runs along the structure",
            stacklevel=2,
        )
    stationS = _makeStations(axisLength, numStations, stations, stationMode)
    tol = stationTol * max(axisLength, float(np.abs(X).max()))

    numSeg = len(segLen)
    N = stationS.shape[0]
    V = np.zeros(N)
    Mb = np.zeros(N)
    T = np.zeros(N)
    force = np.zeros((N, 3))
    moment = np.zeros((N, 3))
    stationPoints = np.zeros((N, 3))
    tangent = np.zeros((N, 3))
    bendingAxis = np.zeros((N, 3))
    shearAxis = np.zeros((N, 3))

    for k, s in enumerate(stationS):
        j = int(np.clip(np.searchsorted(sAxis, s, side="right") - 1, 0, numSeg - 1))
        t = tangents[j]
        xs = starts[j] + (s - sAxis[j]) * t
        sh, b = _localFrame(shearDir, t, j)

        outboard = nodeS > s + tol
        Fk = F[outboard].sum(axis=0)
        Mk = M[outboard].sum(axis=0) + np.cross(X[outboard] - xs, F[outboard]).sum(
            axis=0
        )

        V[k] = Fk @ sh
        Mb[k] = Mk @ b
        T[k] = Mk @ t
        force[k] = Fk
        moment[k] = Mk
        stationPoints[k] = xs
        tangent[k] = t
        bendingAxis[k] = b
        shearAxis[k] = sh

    return VMTResult(
        s=stationS,
        V=V,
        M=Mb,
        T=T,
        force=force,
        moment=moment,
        stationPoints=stationPoints,
        tangent=tangent,
        bendingAxis=bendingAxis,
        shearAxis=shearAxis,
        nodeS=nodeS,
        axisPoints=cleanPts,
        axisLength=axisLength,
    )


def computeVMTFromF5(
    source,
    axisPts,
    shearDir=(0.0, 0.0, -1.0),
    numStations=20,
    components=None,
    includeReactions=None,
    **kwargs,
):
    """
    Compute shear, bending moment and torque from a TACS f5 file.

    Parameters
    ----------
    source : str or F5Data
        Path of an f5 file, or data already read with :func:`loadF5`.
    axisPts : array_like
        Beam axis vertices from root to tip, shape ``(npts, 3)``.
    shearDir : array_like
        Shear direction, default ``(0, 0, -1)``.
    numStations : int
        Number of stations, default 20.
    components : None, str, int or sequence
        Components whose nodes are summed, see :func:`selectComponents`.
    includeReactions : bool or None
        Whether to add support reactions, see :func:`nodalLoads`.
    **kwargs
        Passed to :func:`computeVMT` (``stations``, ``stationMode``,
        ``stationTol``).

    Returns
    -------
    VMTResult
    """
    data = loadF5(source) if isinstance(source, str) else source
    X = nodalCoordinates(data)
    F, M, usedReactions = nodalLoads(data, includeReactions)
    _, nodeMask = selectComponents(data, components)
    result = computeVMT(
        X,
        F,
        axisPts,
        M=M,
        shearDir=shearDir,
        numStations=numStations,
        nodeMask=nodeMask,
        **kwargs,
    )
    return VMTResult(**{**result.__dict__, "includeReactions": usedReactions})


def _checkScale(name, value):
    """Return ``value`` as a float after checking it is a positive finite number."""
    try:
        scale = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a positive number, got {value!r}") from exc
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError(f"{name} must be a positive number, got {value!r}")
    return scale


def scaleVMTResult(result, forceScale=1.0, lengthScale=1.0):
    """
    Return a copy of ``result`` with all forces and lengths multiplied by constants.

    This is a plain unit conversion of the numbers: forces (``V``, ``force``)
    are multiplied by ``forceScale``, lengths (``s``, ``nodeS``,
    ``stationPoints``, ``axisPoints``, ``axisLength``) by ``lengthScale`` and
    moments (``M``, ``T``, ``moment``) by their product. For example a model
    in N and m is converted to kN and mm with ``forceScale=1e-3`` and
    ``lengthScale=1e3``.

    Parameters
    ----------
    result : VMTResult
    forceScale : float
        Multiplier applied to every force, must be positive.
    lengthScale : float
        Multiplier applied to every length, must be positive.

    Returns
    -------
    VMTResult
        The same object when both factors are ``1``, otherwise a scaled copy.

    Raises
    ------
    ValueError
        If a scale factor is not a positive finite number.
    """
    forceScale = _checkScale("forceScale", forceScale)
    lengthScale = _checkScale("lengthScale", lengthScale)
    if forceScale == 1.0 and lengthScale == 1.0:
        return result
    momentScale = forceScale * lengthScale
    return replace(
        result,
        s=result.s * lengthScale,
        V=result.V * forceScale,
        M=result.M * momentScale,
        T=result.T * momentScale,
        force=result.force * forceScale,
        moment=result.moment * momentScale,
        stationPoints=result.stationPoints * lengthScale,
        nodeS=result.nodeS * lengthScale,
        axisPoints=result.axisPoints * lengthScale,
        axisLength=result.axisLength * lengthScale,
    )


# ---------------------------------------------------------------------------
# Layer 4: planform silhouette and plotting
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Silhouette:
    """
    Planform geometry projected along the shear direction.

    Attributes
    ----------
    quads : numpy.ndarray
        Corner polygons of quadrilateral elements, ``(m, 4, 2)``.
    tris : numpy.ndarray
        Corner polygons of triangular elements, ``(m, 3, 2)``.
    lines : numpy.ndarray
        Segments of line (beam) elements, ``(k, 2, 2)``.
    axisXY : numpy.ndarray
        Projected beam axis vertices, ``(npts, 2)``.
    stationXY : numpy.ndarray
        Projected station points, ``(N, 2)``.
    stationDirXY : numpy.ndarray
        Projected station tangents (not normalised), ``(N, 2)``.
    e1, e2, origin : numpy.ndarray
        In-plane basis vectors and origin used for the projection.
    """

    quads: np.ndarray
    tris: np.ndarray
    lines: np.ndarray
    axisXY: np.ndarray
    stationXY: np.ndarray
    stationDirXY: np.ndarray
    e1: np.ndarray
    e2: np.ndarray
    origin: np.ndarray


def _projectionBasis(axisPts, shearDir):
    """
    Return the in-plane basis ``(e1, e2)`` of the plane normal to ``shearDir``.

    ``e1`` follows the root-to-tip direction of the axis and ``e2`` is the
    screen-up direction for a viewer looking along ``shearDir``, so the
    default ``(0, 0, -1)`` gives the usual top view with ``x`` right, ``y`` up.
    """
    sHat = np.asarray(shearDir, dtype=np.float64)
    sHat = sHat / np.linalg.norm(sHat)
    span = axisPts[-1] - axisPts[0]
    e1 = span - np.dot(span, sHat) * sHat
    if np.linalg.norm(e1) <= 1e-8 * max(np.linalg.norm(span), 1.0):
        # Axis parallel to the shear direction: use any perpendicular
        trial = np.zeros(3)
        trial[np.argmin(np.abs(sHat))] = 1.0
        e1 = np.cross(sHat, trial)
    e1 = e1 / np.linalg.norm(e1)
    # Screen-up direction for a viewer looking along sHat (e1 x e2 points
    # back towards the viewer)
    e2 = np.cross(e1, sHat)
    return e1, e2


def _cornerOffsets(layout, numNodes):
    """Return the connectivity slots of an element's corner nodes in polygon order."""
    if layout in _TRI_LAYOUTS:
        return np.array([0, 1, 2], dtype=np.intp)
    if layout in _QUAD_LAYOUTS:
        p = int(round(np.sqrt(numNodes)))
        return np.array([0, p - 1, p * p - 1, p * (p - 1)], dtype=np.intp)
    raise ValueError(f"layout {layout} has no polygon representation")


def buildSilhouette(data, result, shearDir=(0.0, 0.0, -1.0), elemMask=None):
    """
    Project the selected elements onto the plane normal to the shear direction.

    Parameters
    ----------
    data : F5Data
        Data read with :func:`loadF5`.
    result : VMTResult
        Result whose axis and stations are projected alongside the elements.
    shearDir : array_like
        Projection direction, default ``(0, 0, -1)``.
    elemMask : array_like of bool, optional
        Elements to draw; defaults to all. Point, RBE and solid layouts are
        never drawn.

    Returns
    -------
    Silhouette
    """
    X = nodalCoordinates(data)
    axisPts = result.axisPoints
    e1, e2 = _projectionBasis(axisPts, shearDir)
    origin = axisPts[0]
    basis = np.column_stack([e1, e2])
    xy = (X - origin) @ basis

    if elemMask is None:
        elemMask = np.ones(data.numElements, dtype=bool)
    elemMask = np.asarray(elemMask, dtype=bool)
    counts = _elementNodeCounts(data)

    quads = []
    tris = []
    lines = []
    solidSeen = False
    for layout in np.unique(data.ltypes[elemMask]):
        layout = int(layout)
        elems = np.flatnonzero(elemMask & (data.ltypes == layout))
        if elems.size == 0:
            continue
        if layout in _SOLID_LAYOUTS:
            solidSeen = True
            continue
        if layout not in _TRI_LAYOUTS + _QUAD_LAYOUTS + _LINE_LAYOUTS:
            continue
        numNodes = int(counts[elems[0]])
        if layout in _LINE_LAYOUTS:
            slots = data.ptr[elems][:, None] + np.arange(numNodes)[None, :]
            pts = xy[data.conn[slots]]  # (m, numNodes, 2)
            segs = np.stack([pts[:, :-1, :], pts[:, 1:, :]], axis=2)
            lines.append(segs.reshape(-1, 2, 2))
            continue
        offsets = _cornerOffsets(layout, numNodes)
        slots = data.ptr[elems][:, None] + offsets[None, :]
        polys = xy[data.conn[slots]]
        if layout in _QUAD_LAYOUTS:
            quads.append(polys)
        else:
            tris.append(polys)
    if solidSeen:
        warnings.warn(
            "solid elements are not drawn in the planform silhouette", stacklevel=2
        )

    def stack(parts, k):
        return np.concatenate(parts, axis=0) if parts else np.zeros((0, k, 2))

    return Silhouette(
        quads=stack(quads, 4),
        tris=stack(tris, 3),
        lines=stack(lines, 2),
        axisXY=(axisPts - origin) @ basis,
        stationXY=(result.stationPoints - origin) @ basis,
        stationDirXY=result.tangent @ basis,
        e1=e1,
        e2=e2,
        origin=origin,
    )


_DIAGRAM_LABELS = (
    ("V", r"Shear $V = \mathbf{F}\cdot\hat{s}$"),
    ("M", r"Bending $M = \mathbf{Mom}\cdot\hat{b}$"),
    ("T", r"Torque $T = \mathbf{Mom}\cdot\hat{t}$"),
)


def _unitSuffix(unit, sep=" "):
    """Return ``"[unit]"`` preceded by ``sep`` for an axis label, or ``""`` without a unit."""
    return f"{sep}[{unit}]" if unit else ""


def _momentUnit(forceUnit, lengthUnit):
    """Return the moment unit, e.g. ``"N·m"``, or ``None`` unless both are given."""
    if forceUnit and lengthUnit:
        return f"{forceUnit}·{lengthUnit}"
    return None


def _diagramLabels(forceUnit=None, lengthUnit=None):
    """Return ``_DIAGRAM_LABELS`` with the unit of every quantity appended."""
    momentUnit = _momentUnit(forceUnit, lengthUnit)
    units = (forceUnit, momentUnit, momentUnit)
    # The unit goes on a second line so the long y labels fit the panel height
    return tuple(
        (key, label + _unitSuffix(unit, sep="\n"))
        for (key, label), unit in zip(_DIAGRAM_LABELS, units, strict=True)
    )


def _scaleSilhouette(silhouette, lengthScale):
    """Return the silhouette with every in-plane coordinate multiplied by ``lengthScale``."""
    lengthScale = _checkScale("lengthScale", lengthScale)
    if lengthScale == 1.0:
        return silhouette
    return replace(
        silhouette,
        quads=silhouette.quads * lengthScale,
        tris=silhouette.tris * lengthScale,
        lines=silhouette.lines * lengthScale,
        axisXY=silhouette.axisXY * lengthScale,
        stationXY=silhouette.stationXY * lengthScale,
        origin=silhouette.origin * lengthScale,
    )


def _importPyplot():
    """Import matplotlib lazily with an installation hint on failure."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError(
            "plotting requires matplotlib; install it with 'pip install tacs[postprocess]'"
        ) from exc
    return plt


def _arcColormap(cmap, cmapRange, axisLength):
    """Return the (colormap, norm) pair used for the arc length in every panel."""
    from matplotlib.colors import ListedColormap, Normalize

    plt = _importPyplot()
    baseCmap = plt.get_cmap(cmap)
    low, high = cmapRange
    arcCmap = ListedColormap(baseCmap(np.linspace(low, high, 256)), name="arcLength")
    arcNorm = Normalize(vmin=0.0, vmax=axisLength)
    return arcCmap, arcNorm


def _resolveSilhouette(result, data, silhouette, shearDir, elemMask):
    """Return the silhouette to draw, building an axis-only one if needed."""
    if silhouette is not None:
        return silhouette
    if data is not None:
        return buildSilhouette(data, result, shearDir, elemMask)
    axisPts = result.axisPoints
    e1, e2 = _projectionBasis(axisPts, shearDir)
    basis = np.column_stack([e1, e2])
    return Silhouette(
        quads=np.zeros((0, 4, 2)),
        tris=np.zeros((0, 3, 2)),
        lines=np.zeros((0, 2, 2)),
        axisXY=(axisPts - axisPts[0]) @ basis,
        stationXY=(result.stationPoints - axisPts[0]) @ basis,
        stationDirXY=result.tangent @ basis,
        e1=e1,
        e2=e2,
        origin=axisPts[0],
    )


def _newFigure(figsize, forceUnit=None, lengthUnit=None):
    """Create the figure with the planform panel above three shared-x diagrams."""
    from matplotlib.gridspec import GridSpec

    plt = _importPyplot()
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(4, 1, height_ratios=[2.2, 1.0, 1.0, 1.0], hspace=0.35, figure=fig)
    axPlan = fig.add_subplot(gs[0])
    axV = fig.add_subplot(gs[1])
    axM = fig.add_subplot(gs[2], sharex=axV)
    axT = fig.add_subplot(gs[3], sharex=axV)
    labels = _diagramLabels(forceUnit, lengthUnit)
    for ax, (_, label) in zip((axV, axM, axT), labels, strict=True):
        ax.axhline(0.0, color="0.6", linewidth=0.8)
        ax.set_ylabel(label)
        ax.grid(True, alpha=0.3)
    axT.set_xlabel("arc length along beam axis, s" + _unitSuffix(lengthUnit))
    return fig, (axPlan, axV, axM, axT)


def _drawPlanform(
    fig,
    axPlan,
    silhouette,
    result,
    arcCmap,
    arcNorm,
    faceColor,
    lineColor,
    title,
    lengthUnit=None,
):
    """Draw the silhouette, the arc-length coloured axis, station ticks and colourbar."""
    from matplotlib.cm import ScalarMappable
    from matplotlib.collections import LineCollection, PolyCollection

    for polys in (silhouette.quads, silhouette.tris):
        if polys.shape[0] > 0:
            axPlan.add_collection(
                PolyCollection(polys, facecolors=faceColor, edgecolors="none", zorder=1)
            )
    if silhouette.lines.shape[0] > 0:
        axPlan.add_collection(
            LineCollection(silhouette.lines, colors=lineColor, linewidths=0.8, zorder=2)
        )

    # Beam axis drawn as short segments coloured by arc length
    axisPts = result.axisPoints
    sAxis = np.concatenate(
        [[0.0], np.cumsum(np.linalg.norm(np.diff(axisPts, axis=0), axis=1))]
    )
    sDense = np.union1d(np.linspace(0.0, result.axisLength, 256), sAxis)
    dense3d = np.column_stack(
        [np.interp(sDense, sAxis, axisPts[:, k]) for k in range(3)]
    )
    denseXY = (dense3d - silhouette.origin) @ np.column_stack(
        [silhouette.e1, silhouette.e2]
    )
    axisLine = LineCollection(
        np.stack([denseXY[:-1], denseXY[1:]], axis=1),
        cmap=arcCmap,
        norm=arcNorm,
        linewidths=3.0,
        capstyle="round",
        zorder=5,
    )
    axisLine.set_array(0.5 * (sDense[:-1] + sDense[1:]))
    axPlan.add_collection(axisLine)
    axPlan.plot(
        silhouette.axisXY[:, 0],
        silhouette.axisXY[:, 1],
        "o",
        markersize=5,
        markerfacecolor="white",
        markeredgecolor="0.25",
        markeredgewidth=0.8,
        zorder=6,
    )

    # Station tick marks perpendicular to the projected axis tangent
    extent = np.ptp(silhouette.axisXY, axis=0).max()
    if silhouette.quads.shape[0] + silhouette.tris.shape[0] + silhouette.lines.shape[0]:
        allPts = np.concatenate(
            [
                silhouette.quads.reshape(-1, 2),
                silhouette.tris.reshape(-1, 2),
                silhouette.lines.reshape(-1, 2),
            ]
        )
        extent = max(extent, np.ptp(allPts, axis=0).max())
    tickLen = 0.02 * extent if extent > 0 else 1.0
    ticks = []
    for xy, dxy in zip(silhouette.stationXY, silhouette.stationDirXY, strict=True):
        norm = np.linalg.norm(dxy)
        perp = np.array([-dxy[1], dxy[0]]) / norm if norm > 0 else np.array([0.0, 1.0])
        ticks.append([xy - tickLen * perp, xy + tickLen * perp])
    if ticks:
        axPlan.add_collection(
            LineCollection(
                np.array(ticks),
                colors=arcCmap(arcNorm(result.s)),
                linewidths=1.2,
                zorder=4,
            )
        )
    axPlan.autoscale_view()
    axPlan.set_aspect("equal", adjustable="datalim")
    axPlan.set_xlabel("in-plane coordinate along axis" + _unitSuffix(lengthUnit))
    axPlan.set_ylabel("in-plane coordinate" + _unitSuffix(lengthUnit))
    axPlan.set_title(
        title if title is not None else "Planform projected along the shear direction"
    )
    fig.colorbar(
        ScalarMappable(norm=arcNorm, cmap=arcCmap),
        ax=axPlan,
        fraction=0.04,
        pad=0.02,
        label="beam axis arc length, s" + _unitSuffix(lengthUnit),
    )


def plotVMT(
    result,
    data=None,
    silhouette=None,
    shearDir=(0.0, 0.0, -1.0),
    elemMask=None,
    fileName=None,
    show=False,
    figsize=(8.0, 10.0),
    title=None,
    faceColor="0.85",
    lineColor="0.45",
    cmap="viridis",
    cmapRange=(0.0, 1.0),
    forceUnit=None,
    lengthUnit=None,
    forceScale=1.0,
    lengthScale=1.0,
):
    """
    Plot the planform silhouette with the beam axis, and V, M, T below it.

    The beam axis in the planform panel and the station markers in the three
    diagrams share one sequential colour ramp over the arc length ``s``, so a
    point in a diagram can be traced back to its physical location on the
    axis.

    Parameters
    ----------
    result : VMTResult
        Result of :func:`computeVMT` or :func:`computeVMTFromF5`.
    data : F5Data, optional
        Used to build the silhouette when ``silhouette`` is not given. If
        both are omitted only the beam axis and stations are drawn.
    silhouette : Silhouette, optional
        Pre-built silhouette from :func:`buildSilhouette`.
    shearDir : array_like
        Projection direction used when building the silhouette.
    elemMask : array_like of bool, optional
        Elements to draw, see :func:`buildSilhouette`.
    fileName : str, optional
        If given, the figure is saved to this path (any format matplotlib
        supports, e.g. ``.pdf`` or ``.png``).
    show : bool
        Call ``matplotlib.pyplot.show()`` before returning.
    figsize : tuple of float
        Figure size in inches.
    title : str, optional
        Figure title.
    faceColor, lineColor : color
        Colours of the shell fill and of beam element lines.
    cmap : str or matplotlib.colors.Colormap
        Sequential colour map used for the arc length along the beam axis.
    cmapRange : tuple of float
        Fraction of ``cmap`` that is used, ``(low, high)`` in ``[0, 1]``. Useful
        to skip the near-white end of single-hue maps such as ``"Blues"``.
    forceUnit, lengthUnit : str, optional
        Unit labels appended to the axis labels, e.g. ``"N"`` and ``"m"`` or
        ``"lbf"`` and ``"in"``. The moment and torque axes are labelled
        ``forceUnit·lengthUnit`` and only get a unit when both are given.
        The labels are cosmetic; combine them with the scale factors below
        when the plotted numbers have to be converted first.
    forceScale, lengthScale : float
        Multipliers applied to the plotted forces and lengths (moments and
        the planform geometry scale accordingly), see
        :func:`scaleVMTResult`. ``result``, ``data`` and ``silhouette`` are
        always given in model units. For example a model in N and m is
        plotted in kN and mm with ``forceScale=1e-3, lengthScale=1e3,
        forceUnit="kN", lengthUnit="mm"``.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : tuple of matplotlib.axes.Axes
        ``(planform, shear, moment, torque)`` axes.

    Raises
    ------
    ImportError
        If matplotlib is not installed.
    ValueError
        If a scale factor is not a positive finite number.
    """
    plt = _importPyplot()
    silhouette = _resolveSilhouette(result, data, silhouette, shearDir, elemMask)
    result = scaleVMTResult(result, forceScale, lengthScale)
    silhouette = _scaleSilhouette(silhouette, lengthScale)
    arcCmap, arcNorm = _arcColormap(cmap, cmapRange, result.axisLength)
    fig, axes = _newFigure(figsize, forceUnit, lengthUnit)
    axPlan, axV, axM, axT = axes
    _drawPlanform(
        fig,
        axPlan,
        silhouette,
        result,
        arcCmap,
        arcNorm,
        faceColor,
        lineColor,
        title,
        lengthUnit,
    )

    for ax, (key, _) in zip((axV, axM, axT), _DIAGRAM_LABELS, strict=True):
        values = getattr(result, key)
        ax.plot(result.s, values, "-", color="0.35", linewidth=1.2, zorder=2)
        ax.scatter(
            result.s,
            values,
            c=result.s,
            cmap=arcCmap,
            norm=arcNorm,
            s=36,
            edgecolors="white",
            linewidths=0.6,
            zorder=3,
        )
    axV.set_xlim(0.0, result.axisLength)

    if fileName is not None:
        fig.savefig(fileName, dpi=150, bbox_inches="tight")
    if show:  # pragma: no cover - interactive
        plt.show()
    return fig, axes


# ---------------------------------------------------------------------------
# Layer 4b: several load cases and their envelope
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class VMTEnvelope:
    """
    Station-wise minimum and maximum of V, M and T over several load cases.

    Attributes
    ----------
    s : numpy.ndarray
        Station arc lengths shared by all cases, shape ``(N,)``.
    caseNames : tuple of str
        Names of the load cases in the order they were given.
    Vmin, Vmax, Mmin, Mmax, Tmin, Tmax : numpy.ndarray
        Envelope values at every station.
    VminCase, VmaxCase, MminCase, MmaxCase, TminCase, TmaxCase : numpy.ndarray
        Name of the case that governs each envelope value.
    axisPoints : numpy.ndarray
        Beam axis vertices of the first case.
    axisLength : float
        Arc length of the beam axis.
    """

    s: np.ndarray
    caseNames: tuple
    Vmin: np.ndarray
    Vmax: np.ndarray
    Mmin: np.ndarray
    Mmax: np.ndarray
    Tmin: np.ndarray
    Tmax: np.ndarray
    VminCase: np.ndarray
    VmaxCase: np.ndarray
    MminCase: np.ndarray
    MmaxCase: np.ndarray
    TminCase: np.ndarray
    TmaxCase: np.ndarray
    axisPoints: np.ndarray
    axisLength: float

    def asDict(self):
        """Return the station arrays as a dict keyed by name."""
        return {
            "s": self.s,
            "Vmin": self.Vmin,
            "Vmax": self.Vmax,
            "Mmin": self.Mmin,
            "Mmax": self.Mmax,
            "Tmin": self.Tmin,
            "Tmax": self.Tmax,
            "VminCase": self.VminCase,
            "VmaxCase": self.VmaxCase,
            "MminCase": self.MminCase,
            "MmaxCase": self.MmaxCase,
            "TminCase": self.TminCase,
            "TmaxCase": self.TmaxCase,
        }


def _caseNames(sources, names=None):
    """Return unique case names, defaulting to the f5 file names without extension."""
    if names is None:
        names = []
        for src in sources:
            fileName = src if isinstance(src, str) else src.fileName
            names.append(os.path.splitext(os.path.basename(fileName))[0])
    names = [str(name) for name in names]
    if len(names) != len(sources):
        raise ValueError(
            f"{len(names)} case names were given for {len(sources)} sources"
        )
    unique = []
    for name in names:
        candidate = name
        k = 2
        while candidate in unique:
            candidate = f"{name}_{k}"
            k += 1
        unique.append(candidate)
    return unique


def computeVMTCases(
    sources,
    axisPts,
    names=None,
    shearDir=(0.0, 0.0, -1.0),
    numStations=20,
    components=None,
    includeReactions=None,
    **kwargs,
):
    """
    Compute V, M and T for several load cases on the same beam axis.

    Parameters
    ----------
    sources : sequence of str or F5Data
        One f5 file (or loaded data set) per load case.
    axisPts : array_like
        Beam axis vertices from root to tip, shape ``(npts, 3)``.
    names : sequence of str, optional
        Case names; default to the file names without extension. Duplicates
        get a numeric suffix.
    shearDir, numStations, components, includeReactions, **kwargs
        Passed to :func:`computeVMTFromF5`.

    Returns
    -------
    dict
        ``{caseName: VMTResult}`` in the order of ``sources``.
    """
    sources = list(sources)
    if not sources:
        raise ValueError("at least one f5 source is required")
    names = _caseNames(sources, names)
    return {
        name: computeVMTFromF5(
            source,
            axisPts,
            shearDir=shearDir,
            numStations=numStations,
            components=components,
            includeReactions=includeReactions,
            **kwargs,
        )
        for name, source in zip(names, sources, strict=True)
    }


def _asCaseDict(results):
    """Accept a dict or a sequence of results and return an ordered dict."""
    if isinstance(results, dict):
        cases = dict(results)
    else:
        results = list(results)
        cases = {f"case {i + 1}": res for i, res in enumerate(results)}
    if not cases:
        raise ValueError("at least one VMTResult is required")
    return cases


def computeEnvelope(results):
    """
    Compute the station-wise min/max envelope of several results.

    Parameters
    ----------
    results : dict or sequence of VMTResult
        Results computed on the same axis with the same stations, e.g. from
        :func:`computeVMTCases`.

    Returns
    -------
    VMTEnvelope

    Raises
    ------
    ValueError
        If the results do not share the same stations.
    """
    cases = _asCaseDict(results)
    names = tuple(cases.keys())
    first = next(iter(cases.values()))
    for name, res in cases.items():
        if res.s.shape != first.s.shape or not np.allclose(res.s, first.s):
            raise ValueError(
                f"case '{name}' has different stations from case '{names[0]}'; "
                "compute all cases with the same axis and station settings"
            )
    nameArray = np.array(names, dtype=object)
    fields = {}
    for key, _ in _DIAGRAM_LABELS:
        stack = np.vstack([getattr(res, key) for res in cases.values()])  # (ncase, N)
        imin = np.argmin(stack, axis=0)
        imax = np.argmax(stack, axis=0)
        fields[f"{key}min"] = stack.min(axis=0)
        fields[f"{key}max"] = stack.max(axis=0)
        fields[f"{key}minCase"] = nameArray[imin]
        fields[f"{key}maxCase"] = nameArray[imax]
    return VMTEnvelope(
        s=first.s.copy(),
        caseNames=names,
        axisPoints=first.axisPoints,
        axisLength=first.axisLength,
        **fields,
    )


def plotVMTEnvelope(
    results,
    data=None,
    silhouette=None,
    shearDir=(0.0, 0.0, -1.0),
    elemMask=None,
    fileName=None,
    show=False,
    figsize=(8.0, 10.0),
    title=None,
    faceColor="0.85",
    lineColor="0.45",
    cmap="viridis",
    cmapRange=(0.0, 1.0),
    bandColor="0.55",
    maxCaseLines=10,
    forceUnit=None,
    lengthUnit=None,
    forceScale=1.0,
    lengthScale=1.0,
):
    """
    Plot the min/max envelope of V, M and T over several load cases.

    The planform panel is identical to :func:`plotVMT`. Each diagram shows
    the envelope as a shaded band with its min and max curves, and the
    individual cases as thin lines with a legend (up to ``maxCaseLines``
    cases; beyond that the cases are drawn in grey without a legend).

    Parameters
    ----------
    results : dict or sequence of VMTResult
        Results of :func:`computeVMTCases` (``{name: result}``).
    data, silhouette, shearDir, elemMask, fileName, show, figsize, title, faceColor, lineColor, cmap, cmapRange
        As in :func:`plotVMT`.
    bandColor : color
        Colour of the envelope band and of its bounding curves.
    maxCaseLines : int
        Largest number of cases that get individual colours and a legend.
    forceUnit, lengthUnit, forceScale, lengthScale
        Unit labels and unit-conversion factors, as in :func:`plotVMT`. The
        envelope is computed from the scaled cases.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : tuple of matplotlib.axes.Axes
        ``(planform, shear, moment, torque)`` axes.
    """
    plt = _importPyplot()
    cases = _asCaseDict(results)
    first = next(iter(cases.values()))
    silhouette = _resolveSilhouette(first, data, silhouette, shearDir, elemMask)
    silhouette = _scaleSilhouette(silhouette, lengthScale)
    cases = {
        name: scaleVMTResult(res, forceScale, lengthScale)
        for name, res in cases.items()
    }
    first = next(iter(cases.values()))
    envelope = computeEnvelope(cases)
    arcCmap, arcNorm = _arcColormap(cmap, cmapRange, first.axisLength)
    fig, axes = _newFigure(figsize, forceUnit, lengthUnit)
    axPlan, axV, axM, axT = axes
    _drawPlanform(
        fig,
        axPlan,
        silhouette,
        first,
        arcCmap,
        arcNorm,
        faceColor,
        lineColor,
        title,
        lengthUnit,
    )

    numCases = len(cases)
    if numCases <= maxCaseLines:
        caseColors = plt.get_cmap("tab10")(np.arange(numCases) % 10)
    else:
        caseColors = [(0.5, 0.5, 0.5, 0.6)] * numCases

    for ax, (key, _) in zip((axV, axM, axT), _DIAGRAM_LABELS, strict=True):
        lo = getattr(envelope, f"{key}min")
        hi = getattr(envelope, f"{key}max")
        ax.fill_between(
            envelope.s,
            lo,
            hi,
            color=bandColor,
            alpha=0.25,
            linewidth=0,
            zorder=1,
            label="envelope (min/max)" if key == "V" else None,
        )
        ax.plot(envelope.s, lo, "-", color=bandColor, linewidth=1.6, zorder=2)
        ax.plot(envelope.s, hi, "-", color=bandColor, linewidth=1.6, zorder=2)
        for (name, res), color in zip(cases.items(), caseColors, strict=True):
            ax.plot(
                res.s,
                getattr(res, key),
                "-",
                color=color,
                linewidth=1.0,
                zorder=3,
                label=name if key == "V" and numCases <= maxCaseLines else None,
            )
    axV.set_xlim(0.0, first.axisLength)
    axV.legend(loc="best", fontsize="x-small", ncol=2 if numCases > 4 else 1)

    if fileName is not None:
        fig.savefig(fileName, dpi=150, bbox_inches="tight")
    if show:  # pragma: no cover - interactive
        plt.show()
    return fig, axes


def plotVMTReport(
    results,
    fileName,
    data=None,
    silhouette=None,
    shearDir=(0.0, 0.0, -1.0),
    elemMask=None,
    figsize=(8.0, 10.0),
    title=None,
    **kwargs,
):
    """
    Write a multi-page PDF with the envelope followed by one page per case.

    Parameters
    ----------
    results : dict or sequence of VMTResult
        Results of :func:`computeVMTCases`.
    fileName : str
        Output ``.pdf`` path.
    data, silhouette, shearDir, elemMask, figsize
        As in :func:`plotVMT`. The silhouette of the first case is used on
        every page.
    title : str, optional
        Prefix for the page titles.
    **kwargs
        Passed to :func:`plotVMT` and :func:`plotVMTEnvelope`: colours, unit
        labels (``forceUnit``, ``lengthUnit``) and unit-conversion factors
        (``forceScale``, ``lengthScale``).

    Returns
    -------
    int
        Number of pages written (``len(results) + 1`` when there are several
        cases, ``1`` for a single case).
    """
    from matplotlib.backends.backend_pdf import PdfPages

    plt = _importPyplot()
    cases = _asCaseDict(results)
    first = next(iter(cases.values()))
    silhouette = _resolveSilhouette(first, data, silhouette, shearDir, elemMask)
    prefix = f"{title}: " if title else ""
    pages = 0
    with PdfPages(fileName) as pdf:
        if len(cases) > 1:
            fig, _ = plotVMTEnvelope(
                cases,
                silhouette=silhouette,
                shearDir=shearDir,
                figsize=figsize,
                title=f"{prefix}envelope of {len(cases)} load cases",
                **kwargs,
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            pages += 1
        for name, res in cases.items():
            fig, _ = plotVMT(
                res,
                silhouette=silhouette,
                shearDir=shearDir,
                figsize=figsize,
                title=f"{prefix}{name}",
                **kwargs,
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            pages += 1
    return pages


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------
def _writeCsv(fileName, table, meta):
    """Write a dict of equal-length columns as CSV with ``#`` metadata lines."""
    header = list(table.keys())
    numRows = len(next(iter(table.values())))
    with open(fileName, "w") as fp:
        fp.write("# tacs.postprocess.vmt\n")
        for key, value in meta.items():
            fp.write(f"# {key}={value}\n")
        fp.write(",".join(header) + "\n")
        for i in range(numRows):
            cells = []
            for key in header:
                value = table[key][i]
                cells.append(
                    f"{value:.12e}"
                    if isinstance(value, (float, np.floating))
                    else str(value)
                )
            fp.write(",".join(cells) + "\n")


def writeVMTCsv(result, fileName, metadata=None):
    """
    Write the station results to a CSV file.

    The file starts with ``#``-prefixed metadata lines, followed by a header
    row ``s,V,M,T,Fx,Fy,Fz,Mx,My,Mz,x,y,z`` and one row per station. It can be
    read back with ``numpy.loadtxt(fileName, delimiter=",", comments="#")``.

    Parameters
    ----------
    result : VMTResult
    fileName : str
    metadata : dict, optional
        Extra ``key=value`` pairs written to the metadata lines.
    """
    meta = {
        "axisPoints": result.axisPoints.tolist(),
        "axisLength": result.axisLength,
        "includeReactions": result.includeReactions,
    }
    if metadata:
        meta.update(metadata)
    _writeCsv(fileName, result.asDict(), meta)


def writeVMTEnvelopeCsv(envelope, fileName, metadata=None):
    """
    Write an envelope to a CSV file.

    Columns are ``s, Vmin, Vmax, Mmin, Mmax, Tmin, Tmax`` followed by the
    governing case name of each envelope value. Because of the text columns,
    read it back with pandas (``comment="#"``) or with
    ``numpy.genfromtxt(..., delimiter=",", names=True, dtype=None,
    skip_header=<number of leading # lines>, encoding=None)``.

    Parameters
    ----------
    envelope : VMTEnvelope
    fileName : str
    metadata : dict, optional
        Extra ``key=value`` pairs written to the metadata lines.
    """
    meta = {
        "cases": list(envelope.caseNames),
        "axisPoints": envelope.axisPoints.tolist(),
        "axisLength": envelope.axisLength,
    }
    if metadata:
        meta.update(metadata)
    _writeCsv(fileName, envelope.asDict(), meta)


# ---------------------------------------------------------------------------
# Layer 5: command-line interface
# ---------------------------------------------------------------------------
def buildParser():
    """Return the argparse parser of the command-line interface."""
    parser = argparse.ArgumentParser(
        prog="python -m tacs.postprocess.vmt",
        description=(
            "Compute shear (V), bending moment (M) and torque (T) diagrams along "
            "a piecewise-linear beam axis from the nodal loads in TACS f5 files. "
            "Several f5 files are treated as load cases on the same model and "
            "also produce a min/max envelope."
        ),
    )
    parser.add_argument(
        "f5Files", nargs="+", metavar="F5FILE", help="TACS .f5 solution file(s)"
    )
    parser.add_argument(
        "--axis",
        type=float,
        nargs="+",
        metavar="COORD",
        help=(
            "beam axis vertices from root to tip as a flat list of x y z triplets "
            "(at least two points), e.g. --axis 0 0 0  0 10 0"
        ),
    )
    parser.add_argument(
        "--shear-dir",
        type=float,
        nargs=3,
        default=[0.0, 0.0, -1.0],
        metavar=("X", "Y", "Z"),
        help="shear direction (default: 0 0 -1)",
    )
    parser.add_argument(
        "--num-stations", type=int, default=20, help="number of stations (default 20)"
    )
    parser.add_argument(
        "--station-mode",
        choices=("endpoint", "midpoint"),
        default="endpoint",
        help="station placement rule (default endpoint)",
    )
    parser.add_argument(
        "--components",
        nargs="+",
        metavar="NAME",
        help="component names, glob patterns or ids to include (default: all)",
    )
    parser.add_argument(
        "--no-reactions",
        action="store_true",
        help="sum applied loads only, ignore support reactions",
    )
    parser.add_argument(
        "--output",
        metavar="FILE",
        help=(
            "save the figure(s) to this file. A .pdf becomes a multi-page report "
            "(envelope first, then one page per case); any other extension writes "
            "one file per case plus an '_envelope' file when there are several"
        ),
    )
    parser.add_argument(
        "--csv",
        metavar="CSV",
        help=(
            "write station results to this file; with several cases one file per "
            "case plus an '_envelope' file are written"
        ),
    )
    parser.add_argument("--show", action="store_true", help="display the figure")
    parser.add_argument(
        "--force-unit",
        metavar="UNIT",
        help=(
            "force unit written in the plot axis labels, e.g. N or lbf; moments "
            "are labelled UNIT·LENGTH when --length-unit is also given"
        ),
    )
    parser.add_argument(
        "--length-unit",
        metavar="UNIT",
        help="length unit written in the plot axis labels, e.g. m or in",
    )
    parser.add_argument(
        "--force-scale",
        type=_positiveFloat,
        default=1.0,
        metavar="FACTOR",
        help=(
            "multiply the plotted forces by FACTOR (moments scale with it) to "
            "convert units, e.g. 0.001 for N to kN; CSV and table output stay "
            "in model units (default 1)"
        ),
    )
    parser.add_argument(
        "--length-scale",
        type=_positiveFloat,
        default=1.0,
        metavar="FACTOR",
        help=(
            "multiply the plotted lengths by FACTOR (moments and the planform "
            "scale with it) to convert units, e.g. 39.3701 for m to in; CSV and "
            "table output stay in model units (default 1)"
        ),
    )
    parser.add_argument(
        "--list-components",
        action="store_true",
        help="list the components in the first file and exit",
    )
    return parser


def _positiveFloat(text):
    """Parse an argparse value as a strictly positive, finite float."""
    try:
        return _checkScale("scale factor", text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def _parseComponents(tokens):
    if tokens is None:
        return None
    return [int(tok) if tok.lstrip("-").isdigit() else tok for tok in tokens]


def _printTable(result, stream):
    stream.write(f"{'s':>14s} {'V':>16s} {'M':>16s} {'T':>16s}\n")
    for s, V, M, T in zip(result.s, result.V, result.M, result.T, strict=True):
        stream.write(f"{s:14.6e} {V:16.8e} {M:16.8e} {T:16.8e}\n")


def _printEnvelopeTable(envelope, stream):
    cols = ("Vmin", "Vmax", "Mmin", "Mmax", "Tmin", "Tmax")
    stream.write(f"{'s':>14s}" + "".join(f" {c:>16s}" for c in cols) + "\n")
    for i, s in enumerate(envelope.s):
        row = "".join(f" {getattr(envelope, c)[i]:16.8e}" for c in cols)
        stream.write(f"{s:14.6e}{row}\n")


def _withSuffix(fileName, suffix):
    stem, ext = os.path.splitext(fileName)
    return f"{stem}_{suffix}{ext}"


def main(argv=None):
    """
    Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Arguments without the program name; defaults to ``sys.argv[1:]``.

    Returns
    -------
    int
        ``0`` on success, ``2`` on a user or data error.
    """
    parser = buildParser()
    args = parser.parse_args(argv)
    try:
        datas = [loadF5(fileName) for fileName in args.f5Files]
        if args.list_components:
            for name, count in listComponents(datas[0]):
                sys.stdout.write(f"{name}\t{count} elements\n")
            return 0
        if args.axis is None:
            raise ValueError("--axis is required")
        if len(args.axis) % 3 != 0 or len(args.axis) < 6:
            raise ValueError(
                "--axis expects a flat list of x y z triplets with at least two points"
            )
        axisPts = np.asarray(args.axis, dtype=np.float64).reshape(-1, 3)
        components = _parseComponents(args.components)
        includeReactions = False if args.no_reactions else None
        meta = {"shearDir": list(args.shear_dir), "components": components}

        results = computeVMTCases(
            datas,
            axisPts,
            names=_caseNames(args.f5Files),
            shearDir=args.shear_dir,
            numStations=args.num_stations,
            components=components,
            includeReactions=includeReactions,
            stationMode=args.station_mode,
        )
        multi = len(results) > 1
        envelope = computeEnvelope(results) if multi else None

        if args.csv:
            if multi:
                for name, res in results.items():
                    writeVMTCsv(
                        res,
                        _withSuffix(args.csv, name),
                        metadata={"case": name, **meta},
                    )
                writeVMTEnvelopeCsv(
                    envelope, _withSuffix(args.csv, "envelope"), metadata=meta
                )
            else:
                (res,) = results.values()
                writeVMTCsv(res, args.csv, metadata={"file": args.f5Files[0], **meta})

        if args.output or args.show:
            if args.output and not args.show:
                import matplotlib

                matplotlib.use("Agg")
            plt = _importPyplot()
            elemMask, _ = selectComponents(datas[0], components)
            first = next(iter(results.values()))
            silhouette = buildSilhouette(datas[0], first, args.shear_dir, elemMask)
            plotKwargs = {
                "silhouette": silhouette,
                "shearDir": args.shear_dir,
                "forceUnit": args.force_unit,
                "lengthUnit": args.length_unit,
                "forceScale": args.force_scale,
                "lengthScale": args.length_scale,
            }
            isPdf = bool(args.output) and args.output.lower().endswith(".pdf")
            if args.output and isPdf:
                plotVMTReport(results, args.output, **plotKwargs)
            elif args.output and multi:
                for name, res in results.items():
                    fig, _ = plotVMT(
                        res,
                        fileName=_withSuffix(args.output, name),
                        title=name,
                        **plotKwargs,
                    )
                    plt.close(fig)
                fig, _ = plotVMTEnvelope(
                    results,
                    fileName=_withSuffix(args.output, "envelope"),
                    title=f"envelope of {len(results)} load cases",
                    **plotKwargs,
                )
                plt.close(fig)
            if args.show or (args.output and not isPdf and not multi):
                (name, res), *_ = results.items()
                if multi:
                    fig, _ = plotVMTEnvelope(
                        results,
                        show=args.show,
                        title=f"envelope of {len(results)} load cases",
                        **plotKwargs,
                    )
                else:
                    fig, _ = plotVMT(
                        res,
                        fileName=None if isPdf else args.output,
                        show=args.show,
                        title=name,
                        **plotKwargs,
                    )
                plt.close(fig)

        if not (args.csv or args.output or args.show):
            for name, res in results.items():
                if multi:
                    sys.stdout.write(f"case: {name}\n")
                _printTable(res, sys.stdout)
            if multi:
                sys.stdout.write("envelope\n")
                _printEnvelopeTable(envelope, sys.stdout)
        return 0
    except (OSError, ValueError, ImportError) as exc:
        sys.stderr.write(f"error: {exc}\n")
        return 2


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
