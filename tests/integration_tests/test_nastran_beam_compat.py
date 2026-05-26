import os
import unittest

import numpy as np
from mpi4py import MPI

from tacs import pyTACS

"""
Nastran beam compatibility regression suite.

Compares TACS results against Nastran reference data for every supported beam
element/property card combination and every supported field on those cards.
Each (configuration, analysis) pair is its own unit test so a single failure
identifies exactly which card, feature combination, and analysis is broken.

The 252 tests (36 configurations x 7 analyses) self-skip until the Nastran
reference CSVs have been generated and committed. See the directory README and
DESIGN.md for the maintainer workflow that produces the references.
"""

_BASE_DIR = os.path.dirname(os.path.abspath(__file__))
_INPUT_DIR = os.path.join(_BASE_DIR, "input_files", "nastran_beam_compat")

_LADDERS = {
    ("CBAR", "PBAR", None): ["rotation", "wa", "nsm"],
    ("CBAR", "PBARL", "BAR"): ["rotation", "wa", "nsm"],
    ("CBAR", "PBARL", "TUBE"): ["rotation", "nsm"],
    ("CBAR", "PBARL", "T"): ["rotation", "wa", "nsm"],
    ("CBEAM", "PBEAM", None): ["taper", "rotation", "wa", "n", "nsm", "m"],
    ("CBEAM", "PBEAML", "BAR"): ["taper", "rotation", "wa", "nsm"],
    ("CBEAM", "PBEAML", "TUBE"): ["taper", "rotation", "nsm"],
    ("CBEAM", "PBEAML", "T"): ["taper", "rotation", "wa", "nsm"],
}


def _max_normalize(shape):
    """Divide by the element with largest |value| so max(|shape|) = 1."""
    flat = shape.flatten()
    idx = np.argmax(np.abs(flat))
    return shape / flat[idx]


def _align_mode_sign(reference, candidate):
    """Flip candidate's sign if its inner product with reference is negative."""
    if np.dot(reference.flatten(), candidate.flatten()) < 0:
        return -candidate
    return candidate


class NastranCompatBeamBase(unittest.TestCase):
    SOL101_BDF = None
    SOL103_BDF = None
    REF_DIR = None
    STATIC_RTOL = 1e-3
    STATIC_ATOL = 1e-4
    MODAL_FREQ_RTOL = 0.05
    MODAL_SHAPE_RTOL = 1e-3
    MODAL_SHAPE_ATOL = 1e-4

    @classmethod
    def setUpClass(cls):
        if cls.REF_DIR is None or not os.path.isdir(cls.REF_DIR):
            # testflo runs each test method as its own case and does not always
            # propagate a setUpClass skip, so the per-test setUp also guards on
            # REF_DIR. Bail out of the (expensive) pyTACS build here regardless.
            return
        comm = MPI.COMM_WORLD
        cls.fea_static = pyTACS(cls.SOL101_BDF, comm)
        cls.fea_static.initialize()
        cls.static_problems = cls.fea_static.createTACSProbsFromBDF()
        for prob in cls.static_problems.values():
            prob.setOption("L2Convergence", 1e-20)
            prob.setOption("L2ConvergenceRel", 1e-20)
        cls.fea_modal = pyTACS(cls.SOL103_BDF, comm)
        cls.fea_modal.initialize()
        cls.modal_problems = cls.fea_modal.createTACSProbsFromBDF()

    def setUp(self):
        if self.REF_DIR is None or not os.path.isdir(self.REF_DIR):
            raise unittest.SkipTest(
                f"Reference CSVs not found at {self.REF_DIR}. "
                "Run generate_inputs.py -> run_nastran.sh -> extract_nastran_refs.py first."
            )

    def _static_by_name(self, name):
        for prob in self.static_problems.values():
            if prob.name == name:
                return prob
        raise KeyError(f"No static problem named {name!r}")

    def _solve_and_get_disps(self, name):
        prob = self._static_by_name(name)
        prob.solve()
        u_global = self.fea_static.localToGlobalArray(prob.getVariables())
        return u_global.reshape(-1, 6)

    def _assert_static(self, name):
        actual = self._solve_and_get_disps(name)
        ref_path = os.path.join(self.REF_DIR, f"{name}.csv")
        expected = np.loadtxt(ref_path, delimiter=",", comments="#")
        np.testing.assert_allclose(
            actual,
            expected,
            rtol=self.STATIC_RTOL,
            atol=self.STATIC_ATOL,
        )

    def test_static_fx(self):
        self._assert_static("static_fx")

    def test_static_fy(self):
        self._assert_static("static_fy")

    def test_static_fz(self):
        self._assert_static("static_fz")

    def test_static_mx(self):
        self._assert_static("static_mx")

    def test_static_my(self):
        self._assert_static("static_my")

    def test_static_mz(self):
        self._assert_static("static_mz")

    def test_modal(self):
        prob = next(iter(self.modal_problems.values()))
        prob.solve()
        freqs_actual = np.array(
            [
                np.sqrt(prob.getVariables(ii)[0]) / (2 * np.pi)
                for ii in range(prob.getNumEigs())
            ]
        )
        freq_path = os.path.join(self.REF_DIR, "modal_frequencies.csv")
        freqs_expected = np.loadtxt(freq_path, delimiter=",", comments="#").flatten()
        np.testing.assert_allclose(
            freqs_actual,
            freqs_expected,
            rtol=self.MODAL_FREQ_RTOL,
        )
        for ii in range(prob.getNumEigs()):
            _, shape_local = prob.getVariables(ii)
            shape_actual = self.fea_modal.localToGlobalArray(shape_local).reshape(-1, 6)
            shape_actual = _max_normalize(shape_actual)
            mode_path = os.path.join(self.REF_DIR, f"modal_mode{ii + 1:02d}.csv")
            shape_expected = np.loadtxt(mode_path, delimiter=",", comments="#")
            shape_actual = _align_mode_sign(shape_expected, shape_actual)
            np.testing.assert_allclose(
                shape_actual,
                shape_expected,
                rtol=self.MODAL_SHAPE_RTOL,
                atol=self.MODAL_SHAPE_ATOL,
                err_msg=f"Mode {ii + 1} shape mismatch",
            )


for (_element, _prop, _section), _ladder in _LADDERS.items():
    _parts = [_element.lower(), _prop.lower()]
    if _section is not None:
        _parts.append(_section.lower())
    _stem_prefix = "_".join(_parts)

    for _rung_len in range(len(_ladder) + 1):
        _active = _ladder[:_rung_len]
        _rung_name = "baseline" if not _active else "_".join(_active)
        _stem = f"{_stem_prefix}_{_rung_name}"

        _cls = type(
            f"Test_{_stem}",
            (NastranCompatBeamBase,),
            {
                "SOL101_BDF": os.path.join(_INPUT_DIR, f"{_stem}_sol101.bdf"),
                "SOL103_BDF": os.path.join(_INPUT_DIR, f"{_stem}_sol103.bdf"),
                "REF_DIR": os.path.join(_INPUT_DIR, _stem),
                "ACTIVE_FEATURES": tuple(_active),
            },
        )
        globals()[_cls.__name__] = _cls

# Remove the abstract base from the module namespace so test runners (testflo)
# don't collect its method stubs as a standalone test case. The generated
# subclasses retain their own reference to it via __bases__.
del NastranCompatBeamBase
