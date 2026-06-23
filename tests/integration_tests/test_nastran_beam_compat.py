import contextlib
import io
import os
import sys
import tempfile
import unittest

import numpy as np
from mpi4py import MPI
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.properties.bars import PBARL
from pyNastran.bdf.cards.properties.beam import PBEAML

from tacs import pyTACS
from tacs.utilities import Error

"""
Nastran beam compatibility regression suite.

Compares TACS results against Nastran reference data for every supported beam
element/property card combination and every supported field on those cards.
Each (configuration, analysis) pair is its own unit test so a single failure
identifies exactly which card, feature combination, and analysis is broken.

The 308 positive tests (44 configurations x 7 analyses) self-skip until the
Nastran reference CSVs have been generated and committed. See the directory
README and DESIGN.md for the maintainer workflow that produces the references.
A separate set of negative tests asserts that every unsupportable card/section
combination raises a clear error; these need no Nastran reference and always
run.

Set DEBUG = True below to save (a) 6-panel TACS-vs-Nastran comparison plots for
every static load case, (b) 6-panel per-mode comparison plots for every Nastran
reference mode, and (c) a MAC heatmap summarising the cross-correlation between
TACS and Nastran mode shapes. Plots are written to ``<REF_DIR>/debug_plots/``.
"""

# Set to True to save comparison plots for each test into <REF_DIR>/debug_plots/.
DEBUG = False

_BASE_DIR = os.path.dirname(os.path.abspath(__file__))
_INPUT_DIR = os.path.join(_BASE_DIR, "input_files", "nastran_beam_compat")
_BDF_DIR = os.path.join(_INPUT_DIR, "nastran_input_files")

sys.path.insert(0, _INPUT_DIR)
from cases import iterCases, SUPPORTED_SECTIONS  # noqa: E402

_DOF_LABELS = ("u", "v", "w", "rotx", "roty", "rotz")

if DEBUG:
    import matplotlib.pyplot as plt
    import niceplots

    plt.style.use(niceplots.get_style())


def _plot_beam_disps(x, actual, expected, title, output_path):
    """Save a 6-panel figure comparing TACS vs Nastran displacement components along the beam."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    for ax, dof_label, act_comp, exp_comp in zip(
        axes.flatten(), _DOF_LABELS, actual.T, expected.T
    ):
        ax.plot(x, act_comp, label="TACS", lw=1.5, clip_on=False)
        ax.plot(x, exp_comp, "--", label="Nastran", lw=1.5, clip_on=False)
        ax.set_xlabel("x (m)")
        ax.set_ylabel(dof_label)
        ax.legend(fontsize=7)
        niceplots.adjust_spines(ax)
    fig.suptitle(title, fontsize=13)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def _plot_beam_mode_shape(x, actual, expected, mode_num, output_path):
    """Save a 6-panel figure comparing TACS vs Nastran mode shape components along the beam."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    for ax, dof_label, act_comp, exp_comp in zip(
        axes.flatten(), _DOF_LABELS, actual.T, expected.T
    ):
        ax.plot(x, act_comp, label="TACS", lw=1.5, clip_on=False)
        ax.plot(x, exp_comp, "--", label="Nastran", lw=1.5, clip_on=False)
        ax.set_xlabel("x (m)")
        ax.set_ylabel(dof_label)
        ax.legend(fontsize=7)
        niceplots.adjust_spines(ax)
    fig.suptitle(f"Mode {mode_num}", fontsize=13)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


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


def _compute_mac_matrix(shapes_a, shapes_b):
    """Compute the full n_a x n_b Modal Assurance Criterion matrix.

    Each row/col represents one mode; cell (i, j) is the MAC between
    shapes_a[i] and shapes_b[j], a scalar in [0, 1] that is 1 for
    parallel (same-shape) modes and 0 for orthogonal ones.

    Parameters
    ----------
    shapes_a, shapes_b : ndarray, shape (n_modes, n_nodes, 6)

    Returns
    -------
    mac : ndarray, shape (n_a, n_b)
    """
    n_a = shapes_a.shape[0]
    n_b = shapes_b.shape[0]
    mac = np.zeros((n_a, n_b))
    for ii in range(n_a):
        vec_a = shapes_a[ii].flatten()
        for jj in range(n_b):
            vec_b = shapes_b[jj].flatten()
            mac[ii, jj] = np.dot(vec_a, vec_b) ** 2 / (
                np.dot(vec_a, vec_a) * np.dot(vec_b, vec_b)
            )
    return mac


def _plot_mac_heatmap(mac, title, output_path):
    """Save an annotated viridis heatmap of an n x n MAC matrix."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    n_rows, n_cols = mac.shape
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.matshow(mac, cmap="viridis", vmin=0, vmax=1)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_xlabel("Nastran mode")
    ax.set_ylabel("TACS mode")
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([str(jj + 1) for jj in range(n_cols)], fontsize=8)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([str(ii + 1) for ii in range(n_rows)], fontsize=8)
    ax.xaxis.set_label_position("bottom")
    ax.xaxis.tick_bottom()
    ax.set_title(title, fontsize=11)
    for ii in range(n_rows):
        for jj in range(n_cols):
            ax.text(
                jj,
                ii,
                f"{mac[ii, jj]:.2f}",
                ha="center",
                va="center",
                fontsize=7,
                color="white" if mac[ii, jj] < 0.5 else "black",
            )
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


class NastranCompatBeamBase(unittest.TestCase):
    SOL101_BDF = None
    SOL103_BDF = None
    REF_DIR = None
    STATIC_RTOL = 1e-3
    STATIC_ATOL = 1e-6
    MODAL_FREQ_RTOL = 0.05
    MAC_THRESHOLD = 0.99

    @classmethod
    def setUpClass(cls):
        if cls.REF_DIR is None or not os.path.isdir(cls.REF_DIR):
            # testflo runs each test method as its own case and does not always
            # propagate a setUpClass skip, so the per-test setUp also guards on
            # REF_DIR. Bail out of the (expensive) pyTACS build here regardless.
            return
        if not os.path.isfile(cls.SOL101_BDF) or not os.path.isfile(cls.SOL103_BDF):
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
        if not os.path.isfile(self.SOL101_BDF) or not os.path.isfile(self.SOL103_BDF):
            raise unittest.SkipTest(
                f"BDF inputs not found in {_BDF_DIR}. Run generate_inputs.py first."
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
        ref_path = os.path.join(self.REF_DIR, f"{name}.csv")
        if not os.path.isfile(ref_path):
            raise unittest.SkipTest(f"Reference not found: {ref_path}")
        actual = self._solve_and_get_disps(name)
        expected = np.loadtxt(ref_path, delimiter=",", comments="#")
        if DEBUG:
            plot_dir = os.path.join(self.REF_DIR, "debug_plots")
            x = np.linspace(0, 1, actual.shape[0])
            stem = os.path.basename(self.REF_DIR)
            _plot_beam_disps(
                x,
                actual,
                expected,
                f"{stem} — {name}",
                os.path.join(plot_dir, f"{name}.png"),
            )
            # Also dump the TACS solution (.f5) next to the plots for inspection.
            os.makedirs(plot_dir, exist_ok=True)
            self._static_by_name(name).writeSolution(outputDir=plot_dir)
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
        freq_path = os.path.join(self.REF_DIR, "modal_frequencies.csv")
        if not os.path.isfile(freq_path):
            raise unittest.SkipTest(f"Reference not found: {freq_path}")
        prob = next(iter(self.modal_problems.values()))
        prob.solve()
        n_modes = prob.getNumEigs()

        # Collect frequencies and mode shapes before any assertions so that
        # debug plots are always generated when DEBUG=True, regardless of
        # which check fails.
        freqs_actual = np.array(
            [np.sqrt(prob.getVariables(ii)[0]) / (2 * np.pi) for ii in range(n_modes)]
        )
        freqs_expected = np.loadtxt(freq_path, delimiter=",", comments="#").flatten()

        shapes_actual = []
        shapes_expected = []
        for ii in range(n_modes):
            _, shape_local = prob.getVariables(ii)
            shape_actual = self.fea_modal.localToGlobalArray(shape_local).reshape(-1, 6)
            shapes_actual.append(shape_actual)
            mode_path = os.path.join(self.REF_DIR, f"modal_mode{ii + 1:02d}.csv")
            shapes_expected.append(np.loadtxt(mode_path, delimiter=",", comments="#"))
        shapes_actual = np.stack(shapes_actual)
        shapes_expected = np.stack(shapes_expected)

        mac = _compute_mac_matrix(shapes_actual, shapes_expected)

        if DEBUG:
            plot_dir = os.path.join(self.REF_DIR, "debug_plots")
            stem = os.path.basename(self.REF_DIR)
            _plot_mac_heatmap(
                mac,
                f"{stem} — MAC (TACS rows, Nastran cols)",
                os.path.join(plot_dir, "mac.png"),
            )
            x = np.linspace(0, 1, shapes_actual.shape[1])
            for ii in range(n_modes):
                shape_a = _max_normalize(shapes_actual[ii])
                shape_e = shapes_expected[ii]
                shape_a = _align_mode_sign(shape_e, shape_a)
                _plot_beam_mode_shape(
                    x,
                    shape_a,
                    shape_e,
                    ii + 1,
                    os.path.join(plot_dir, f"modal_mode{ii + 1:02d}.png"),
                )
            # Also dump the TACS modal solution (.f5) next to the plots.
            os.makedirs(plot_dir, exist_ok=True)
            prob.writeSolution(outputDir=plot_dir)

        # Frequencies: loose 5% rtol is intentional — known TACS/Nastran beam
        # element formulation differences cause up to 3-4% frequency error even
        # when mode shapes agree very well.
        np.testing.assert_allclose(
            freqs_actual,
            freqs_expected,
            rtol=self.MODAL_FREQ_RTOL,
        )

        # MAC: full n x n matrix. Assert each TACS mode's best Nastran match is
        # its index-paired mode AND the diagonal MAC clears the threshold.
        for ii in range(n_modes):
            best_match = int(np.argmax(mac[ii]))
            if best_match != ii:
                self.fail(
                    f"Mode {ii + 1}: best Nastran match is mode {best_match + 1} "
                    f"(MAC = {mac[ii, best_match]:.3f}), not mode {ii + 1} "
                    f"(MAC = {mac[ii, ii]:.3f}). "
                    "Modes appear to be swapped — see mac.png in debug_plots."
                )
            if mac[ii, ii] < self.MAC_THRESHOLD:
                self.fail(
                    f"Mode {ii + 1}: MAC = {mac[ii, ii]:.4f} < {self.MAC_THRESHOLD:.4f}. "
                    "Modes are paired correctly but the shapes disagree — "
                    "see mac.png and modal_mode*.png in debug_plots."
                )


for _element, _prop, _section, _stem, _features in iterCases():
    _cls = type(
        f"Test_{_stem}",
        (NastranCompatBeamBase,),
        {
            "SOL101_BDF": os.path.join(_BDF_DIR, f"{_stem}_sol101.bdf"),
            "SOL103_BDF": os.path.join(_BDF_DIR, f"{_stem}_sol103.bdf"),
            "REF_DIR": os.path.join(_INPUT_DIR, "nastran_ref_results", _stem),
            "ACTIVE_FEATURES": _features,
        },
    )
    globals()[_cls.__name__] = _cls

# Remove the abstract base from the module namespace so testflo does not
# collect its method stubs as a standalone test case.
del NastranCompatBeamBase


# ==============================================================================
# Negative tests: unsupportable card/section combinations must raise.
# ==============================================================================
# These guard against a routing regression that silently starts "supporting" a
# section for which pyTACS cannot obtain a correct torsion constant J. Each test
# builds a minimal single-element BDF in memory with pyNastran and asserts that
# pyTACS.initialize() raises tacs.utilities.Error. The section list is derived
# from pyNastran's valid_types minus SUPPORTED_SECTIONS, so it tracks pyNastran.
# See docs/adr/0001-nastran-beam-section-support-bounded-by-torsion-constant.md.

# Minimal beam material; values are irrelevant since the error is raised during
# constitutive translation, before any material property is used.
_NEG_E, _NEG_NU, _NEG_RHO = 71.7e9, 0.33, 2810.0

# Valid example dimensions for each unsupported section. PBEAML validates section
# geometry in its constructor, so generic dims are rejected for many types; these
# were vetted against that validation (seeded from pyNastran's own test suite).
# PBARL does not validate geometry, so any correctly-sized dims work there.
_UNSUPPORTED_SECTION_DIMS = {
    "BOX": [1.0, 2.0, 0.1, 0.1],
    "BOX1": [2.0, 2.0, 0.1, 0.1, 0.1, 0.1],
    "CHAN": [2.0, 2.0, 0.1, 0.1],
    "CHAN1": [3.0, 3.0, 2.0, 3.0],
    "CHAN2": [3.0, 3.0, 3.0, 3.0],
    "CROSS": [2.0, 2.0, 0.1, 0.1],
    "DBOX": [3.0, 2.0, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
    "H": [2.0, 2.0, 0.1, 0.1],
    "HAT": [3.0, 1.0, 3.0, 3.0],
    "HAT1": [2.0, 2.0, 0.1, 0.1, 0.1],
    "HEXA": [1.0, 3.0, 3.0],
    "I": [2.0, 2.0, 0.1, 0.1, 0.1, 0.1],
    "I1": [2.0, 2.0, 0.1, 0.1],
    "L": [2.0, 2.0, 0.1, 0.1],
    "T": [1.0, 2.0, 0.1, 0.2],
    "T1": [1.0, 2.0, 0.2, 0.1],
    "T2": [1.0, 2.0, 0.2, 0.1],
    "Z": [2.0, 2.0, 0.1, 0.1],
    # Circular sections; used only for the PBEAML-with-offset cases.
    "ROD": [2.0],
    "TUBE": [2.0, 1.0],
    "TUBE2": [2.0, 0.5],
}

# A shear-centre (WA/WB) offset with non-zero local y and z components, so it
# projects onto the section axes and trips the offset-handling path.
_NEG_OFFSET_VEC = [0.0, 0.1, 0.1]


def _sectionDims(prop, section):
    """Return valid dims for an unsupported section, falling back to a generic
    correctly-sized list for any section not in the vetted table (works for
    PBARL; a new PBEAML section would surface as a visible test error)."""
    if section in _UNSUPPORTED_SECTION_DIMS:
        return list(_UNSUPPORTED_SECTION_DIMS[section])
    validTypes = PBARL.valid_types if prop == "PBARL" else PBEAML.valid_types
    numDims = validTypes[section]
    return [round(0.12 - 0.01 * ii, 4) for ii in range(numDims)]


def _writeSingleElementBeamBdf(path, element, prop, section, dims, withOffset):
    """Write a minimal single-element CBAR/CBEAM + PBARL/PBEAML BDF to ``path``."""
    model = BDF(debug=False)
    model.add_grid(1, [0.0, 0.0, 0.0])
    model.add_grid(2, [1.0, 0.0, 0.0])
    model.add_mat1(1, _NEG_E, None, _NEG_NU, _NEG_RHO)
    orientationVec = [0.0, 1.0, 0.0]
    wa = _NEG_OFFSET_VEC if withOffset else None
    wb = _NEG_OFFSET_VEC if withOffset else None
    if prop == "PBARL":
        model.add_pbarl(10, 1, section, dims)
        model.add_cbar(100, 10, [1, 2], orientationVec, None, wa=wa, wb=wb)
    else:
        model.add_pbeaml(10, 1, section, xxb=[0.0, 1.0], dims=[dims, dims])
        model.add_cbeam(100, 10, [1, 2], orientationVec, None, wa=wa, wb=wb)
    model.write_bdf(path)


class NastranCompatUnsupportedBase(unittest.TestCase):
    # Runs on a single process; the in-memory BDF has one element.
    ELEMENT = None
    PROP = None
    SECTION = None
    WITH_OFFSET = False

    def test_raises(self):
        dims = _sectionDims(self.PROP, self.SECTION)
        fd, path = tempfile.mkstemp(suffix=".bdf")
        os.close(fd)
        try:
            _writeSingleElementBeamBdf(
                path, self.ELEMENT, self.PROP, self.SECTION, dims, self.WITH_OFFSET
            )
            # tacs.utilities.Error prints its message rather than storing it on
            # the exception, so capture stdout to confirm we hit *our* routing
            # branch (the message names the offending section type).
            captured = io.StringIO()
            with self.assertRaises(Error):
                with contextlib.redirect_stdout(captured):
                    fea = pyTACS(path, MPI.COMM_WORLD)
                    fea.initialize()
            self.assertIn(
                self.SECTION,
                captured.getvalue(),
                msg=f"{self.PROP} {self.SECTION} raised, but the message did not "
                "name the section — it may be a different error than the "
                "unsupported-section routing error.",
            )
        finally:
            if os.path.isfile(path):
                os.remove(path)


def _iterUnsupportedCases():
    """Yield (element, prop, section, withOffset) for every unsupportable combo."""
    for prop, element, validTypes in (
        ("PBARL", "CBAR", PBARL.valid_types),
        ("PBEAML", "CBEAM", PBEAML.valid_types),
    ):
        for section in sorted(set(validTypes) - set(SUPPORTED_SECTIONS)):
            yield element, prop, section, False
    # Circular sections with a shear-centre offset: supported for PBARL (PBARL.J()
    # is a correct circular J) but unsupported for PBEAML (PBEAML.J() is None).
    for section in ("ROD", "TUBE", "TUBE2"):
        yield "CBEAM", "PBEAML", section, True


for _element, _prop, _section, _withOffset in _iterUnsupportedCases():
    _suffix = f"{_prop.lower()}_{_section.lower()}" + ("_offset" if _withOffset else "")
    _negName = f"Test_unsupported_{_suffix}"
    _negCls = type(
        _negName,
        (NastranCompatUnsupportedBase,),
        {
            "ELEMENT": _element,
            "PROP": _prop,
            "SECTION": _section,
            "WITH_OFFSET": _withOffset,
        },
    )
    globals()[_negName] = _negCls

del NastranCompatUnsupportedBase

if __name__ == "__main__":
    unittest.main()
