Composite failure criteria
**************************

:class:`~tacs.constitutive.OrthotropicPly` supports five failure criteria, selected with ``setFailureCriterion``.
In every case the returned value ``fail`` is interpreted as

.. math::

    \text{fail} < 1 \Rightarrow \text{no failure}, \qquad \text{fail} \geq 1 \Rightarrow \text{failure}

The ply strains are first rotated into the ply coordinate system, and for the stress-based criteria the ply stresses :math:`\mathbf{s} = [s_1, s_2, s_{12}]` follow from the ply stiffness.

Tsai-Wu
=======

The Tsai-Wu criterion is a quadratic tensor polynomial in the ply stresses,

.. math::

    F(\mathbf{s}) = F_1 s_1 + F_2 s_2 + F_{11} s_1^2 + F_{22} s_2^2 + 2 F_{12} s_1 s_2 + F_{66} s_{12}^2 \leq 1

Select it with ``setFailureCriterion(CompositeFailureCriterion.TSAI_WU)``.

The coefficients are derived from the ply strengths, where :math:`X_t, X_c` are the tensile and compressive strengths along the fibre direction, :math:`Y_t, Y_c` the transverse strengths, and :math:`S_{12}` the in-plane shear strength:

.. math::

    F_1 = \frac{X_c - X_t}{X_t X_c}, \quad
    F_2 = \frac{Y_c - Y_t}{Y_t Y_c}, \quad
    F_{11} = \frac{1}{X_t X_c}, \quad
    F_{22} = \frac{1}{Y_t Y_c}, \quad
    F_{66} = \frac{1}{S_{12}^2}

The interaction coefficient :math:`F_{12}` cannot be derived from the uniaxial strengths alone.
For orthotropic materials TACS sets :math:`F_{12} = 0`, matching the default of the Nastran ``MAT8`` card.
For isotropic materials it is set to the value that recovers the von Mises criterion, as derived below.

Any value of :math:`F_{12}` must satisfy the stability criterion

.. math::

    F_{12}^2 < F_{11} F_{22}

which ensures the failure surface is a closed ellipsoid rather than an open hyperboloid.

Modified Tsai-Wu (strength ratio)
=================================

This is the default criterion, selected with ``setFailureCriterion(CompositeFailureCriterion.TSAI_WU_MODIFIED)``.
Writing :math:`b` for the linear part and :math:`a` for the quadratic part of the Tsai-Wu polynomial, it returns

.. math::

    \text{fail} = \frac{1}{2} \left( b + \sqrt{b^2 + 4a} \right)

This is equivalent to the standard form at the failure boundary, but scales linearly when the stress state is scaled uniformly.
It can therefore be used directly to compute safety factors, which is not true of the quadratic form.

Equivalence with von Mises for isotropic materials
--------------------------------------------------

For an isotropic material with yield stress :math:`\sigma_y`, TACS derives the ply strengths as

.. math::

    X_t = X_c = Y_t = Y_c = \sigma_y, \qquad S_{12} = \frac{\sigma_y}{\sqrt{3}}

The shear strength follows from von Mises: in pure shear :math:`\sigma_{vm} = \sqrt{3}\,\tau`, so yielding occurs at :math:`\tau = \sigma_y/\sqrt{3}`.
The Tsai-Wu coefficients then become

.. math::

    F_1 = F_2 = 0, \qquad
    F_{11} = F_{22} = \frac{1}{\sigma_y^2}, \qquad
    F_{66} = \frac{3}{\sigma_y^2}

and the interaction coefficient is set to

.. math::

    F_{12} = -\frac{1}{2}\sqrt{F_{11} F_{22}} = -\frac{1}{2\sigma_y^2}

which is the normalised value :math:`F_{12}^* = F_{12}/\sqrt{F_{11}F_{22}} = -1/2`.
It satisfies the stability criterion by construction, since :math:`F_{12}^2 = \tfrac{1}{4} F_{11} F_{22} < F_{11} F_{22}`.

Because :math:`F_1 = F_2 = 0` the linear part :math:`b` vanishes, so the modified Tsai-Wu value reduces to the square root of the quadratic part:

.. math::

    \text{fail} = \sqrt{a}
                = \frac{1}{\sigma_y}\sqrt{s_1^2 + s_2^2 - s_1 s_2 + 3 s_{12}^2}

which is exactly the plane-stress von Mises failure index.

Note that the equivalence is exact in *value* only for the modified form.
The standard ``TSAI_WU`` criterion returns :math:`a`, the square of the von Mises index — the same failure surface at :math:`\text{fail} = 1`, but a different value elsewhere.

Cuntze
======

The Cuntze Failure Mode Concept computes a global material stressing effort :math:`\text{Eff}` by interacting several independent failure modes through an exponent :math:`m`.
Fibre failure modes are labelled ``FF`` and inter-fibre failure modes ``IFF``.
If any mode reaches :math:`\text{Eff} \geq 1` the material fails in that mode; the global effort reflects the interaction of simultaneously active modes.

The concept is defined for general 3D stress states; the plane-stress specialisations implemented in TACS are given below.

Unidirectional plies
--------------------

Select with ``setFailureCriterion(CompositeFailureCriterion.CUNTZE_UD)``.

.. math::

    \text{Eff} = \left( \text{Eff}_{FF1}^m + \text{Eff}_{FF2}^m + \text{Eff}_{IFF1}^m
                 + \text{Eff}_{IFF2}^m + \text{Eff}_{IFF3}^m \right)^{1/m} \leq 1

with

.. math::

    \text{Eff}_{FF1} &= \frac{e_1 E_1}{X_t} \quad (e_1 \geq 0) \\
    \text{Eff}_{FF2} &= \frac{-e_1 E_1}{X_c} \quad (e_1 < 0) \\
    \text{Eff}_{IFF1} &= \frac{s_2}{Y_t} \quad (s_2 \geq 0) \\
    \text{Eff}_{IFF2} &= \frac{-s_2}{Y_c} \quad (s_2 < 0) \\
    \text{Eff}_{IFF3} &= \sqrt{\frac{s_{12}^2 \left( b_{tl} s_2
                         + \sqrt{b_{tl}^2 s_2^2 + S_{12}^2} \right)}{S_{12}^3}}

The plane-stress equations are derived from the 3D formulation presented in
E. Petersen, R. G. Cuntze, and C. Huehne, "Experimental determination of material parameters in
Cuntze's Failure-Mode-Concept-based UD strength failure conditions,"
*Compos. Sci. Technol.*, vol. 134, pp. 12-25, Oct. 2016,
doi: `10.1016/j.compscitech.2016.08.006 <https://doi.org/10.1016/j.compscitech.2016.08.006>`_.

Woven plies
-----------

Select with ``setFailureCriterion(CompositeFailureCriterion.CUNTZE_WOVEN)``.

.. math::

    \text{Eff} = \Big( \text{Eff}_{FF1}^m + \text{Eff}_{FF2}^m + \text{Eff}_{FF3}^m
                 + \text{Eff}_{FF4}^m + \text{Eff}_{IFF1}^m + \text{Eff}_{IFF2}^m \\
                 + \text{Eff}_{IFF3}^m + \text{Eff}_{IFF4}^m + \text{Eff}_{IFF5}^m \Big)^{1/m} \leq 1

with

.. math::

    \text{Eff}_{FF1} &= \frac{e_1 E_1}{X_t} \quad (e_1 \geq 0) \\
    \text{Eff}_{FF2} &= \frac{-e_1 E_1}{X_c} \quad (e_1 < 0) \\
    \text{Eff}_{FF3} &= \frac{e_2 E_2}{Y_t} \quad (e_2 \geq 0) \\
    \text{Eff}_{FF4} &= \frac{-e_2 E_2}{Y_c} \quad (e_2 < 0) \\
    \text{Eff}_{IFF3} &= \frac{|s_{12}|}{S_{12} - \mu_{WF}(s_1 + s_2)}

For a plane stress state :math:`\text{Eff}_{IFF1} = \text{Eff}_{IFF2} = \text{Eff}_{IFF4} = \text{Eff}_{IFF5} = 0`.

Adapted from J. Bold, "Vergleich des Impaktverhaltens von monolithischer und hybrider CFK-Platte
unter Verwendung eines neuen Werkstoffmodells," Dr.-Ing. dissertation,
Technische Universitaet Braunschweig, Braunschweig, Germany, 2018.

Maximum strain
==============

Select with ``setFailureCriterion(CompositeFailureCriterion.MAX_STRAIN)``.

.. math::

    \text{KS}\left( \frac{e_i}{e_{i,\max}^{\pm}}, \rho \right) \leq 1

where :math:`e_{i,\max}^{\pm}` are the positive and negative strains at failure for component :math:`i`, and :math:`\rho` is the KS aggregation weight set by ``setKSWeight``.
The KS aggregation gives a smooth, conservative approximation to the maximum of the six ratios, which is desirable for gradient-based optimisation.
