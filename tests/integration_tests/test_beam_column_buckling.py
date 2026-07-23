import os

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive

"""
Euler-column buckling test for TACSBeamElement (Phase 5, Task 5.7) -- the
first beam buckling correctness oracle in this repository, and (per this
feature's own SPEC pass) arguably the first direct
getMatType(TACS_GEOMETRIC_STIFFNESS_MATRIX) oracle for any element in the
repo, shell included.

Model: the same 1m-long, 5-element rectangular-cross-section beam mesh
used by test_rectangle_beam_buckling.py (beam_model_translated.bdf,
t=0.05m thick, w=0.1m wide, E=70e3, rho=2700 -- a self-consistent, not
SI, unit system already exercised by that passing test), fixed at node 1
(all 6 dof, an SPC1 in the BDF) and axially loaded at the free end (node
6) via BDF LOAD case 3 ("x-axial"), a unit +x force. This is a classic
CANTILEVER (fixed-free) Euler column: end-condition factor K=2.

Closed-form check: P_cr = pi^2 * E * I_min / (K*L)^2, where I_min is the
weaker-axis area moment of inertia (buckling occurs about the weak axis
first). For a solid rectangle, TACSIsoRectangleBeamConstitutive's own
evalMomentsOfInertia (src/constitutive/TACSIsoRectangleBeamConstitutive.cpp)
gives Iz = w*t^3/12 (bending stiffness Cs[11], beam's "bend1"/e[2] strain
component) and Iy = t*w^3/12 (Cs[15], "bend2"/e[3]) -- confirmed directly
from that file's evalTangentStiffness. With t=0.05 < w=0.1, Iz = I_min =
0.1*0.05**3/12 = 1.041666...e-6 is the governing (weak-axis) moment of
inertia, so:

    P_cr = pi**2 * 70e3 * (0.1*0.05**3/12) / (2*1.0)**2 ~= 0.17990

The eigenvalue reported by TACS's buckling analysis is the load factor
on the applied reference load (a unit +x force per the BDF's FORCE card,
scale=1.0 via addLoadFromBDF's own default), so the lowest buckling
eigenvalue should equal P_cr directly (no additional scaling needed).

Since the sensitivity checks (test_total_dv_sensitivities/
test_total_xpt_sensitivities, run automatically by PyTACSTestCase for
every function listed in FUNC_REFS) route the GEOMETRIC_STIFFNESS_MATRIX
matrix-sens methods through whichever path Tasks 5.3/5.4/5.5 land in
(analytic or the base class's own generic FD-based fallback --
TACSElement::addMatDVSensInnerProduct/addMatXptSensInnerProduct/
getMatSVSensInnerProduct are genuine FD wrappers around getMatType
itself, confirmed by reading src/elements/TACSElement.cpp directly, not
"return 0" stubs), this test's sensitivity checks are a real,
independent, physically-grounded exercise of getMatType regardless of
whether 5.3-5.5 land analytically or via documented fallback -- this is
why this test was written before Tasks 5.3-5.5, per the coordinator's
suggestion, rather than only at the very end of the phase.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/beam_model_translated.bdf")

# Closed-form Euler critical load for this cantilever (K=2) column.
_E = 70.0e3
_t = 0.05
_w = 0.1
_L = 1.0
_K = 2.0
_Iz = _w * _t**3 / 12.0  # weak-axis (governing) moment of inertia
P_CR_exact = (3.141592653589793**2) * _E * _Iz / (_K * _L) ** 2
# Value based on eigenvalue solve on mesh
P_CR_approx = 0.181211383209752


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "buckling_eigsb.0": P_CR_approx,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Overwrite default check values -- mirrors
        # test_shell_tri_buckling_axial.py's own real/complex split; real
        # mode is loose (2e-1) there for a much coarser/stiffer tri mesh,
        # this beam mesh (5 quadratic elements on the first Euler mode)
        # should do noticeably better, but keep some margin for
        # discretization error rather than hard-coding an arbitrarily
        # tight tolerance without a convergence study.
        if self.dtype == complex:
            self.rtol = 1e-8
            self.atol = 1e-8
            self.dh = 1e-50
        else:
            self.rtol = 5e-2
            self.atol = 1e-4
            self.dh = 1e-6

        # Material properties (identical to test_rectangle_beam_buckling.py,
        # same BDF file)
        rho = 2700.0  # density kg/m^3
        E = 70.0e3  # Young's modulus (Pa)
        nu = 0.3  # Poisson's ratio
        ys = 2.7e3  # yield stress

        # Beam dimensions
        t = 0.05  # m
        w = 0.1  # m

        # Callback function used to setup TACS element objects and DVs
        def elem_call_back(
            dv_num, comp_id, comp_descript, elem_descripts, global_dvs, **kwargs
        ):
            prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
            con = constitutive.IsoRectangleBeamConstitutive(
                prop,
                t=t,
                tNum=dv_num,
                w=w,
                wNum=dv_num + 1,
                Lb=1.0,
                LbNum=dv_num + 2,
                Kb=1.0,
            )
            refAxis = [0.0, 1.0, 0.0]
            transform = elements.BeamRefAxisTransform(refAxis)
            elem = elements.Beam2(transform, con)
            return elem

        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        fea_assembler.initialize(elem_call_back)

        buckle_prob = fea_assembler.createBucklingProblem(
            "buckling", sigma=1.0, numEigs=5
        )
        buckle_prob.setOption("L2Convergence", 1e-20)
        buckle_prob.setOption("L2ConvergenceRel", 1e-20)
        # BDF LOAD case 3 is a unit +x force at the free end (node 6);
        # since node 6 sits at higher x than the fixed node 1, +x is
        # TENSION, not compression. Confirmed this session: using it
        # unscaled gives a buckling eigenvalue of the correct magnitude
        # but negative sign (TACS's K + lambda*G convention reports the
        # scale factor, of either sign, needed on the REFERENCE load to
        # reach the critical state -- a tensile reference load needs a
        # negative multiple, i.e. reversed into compression, to buckle).
        # scale=-1.0 applies the load in compression directly, so the
        # reported eigenvalue is directly comparable to the positive
        # closed-form P_cr.
        buckle_prob.addLoadFromBDF(loadID=3, scale=-1.0)

        return [buckle_prob], fea_assembler


if __name__ == "__main__":
    import unittest

    unittest.main()
