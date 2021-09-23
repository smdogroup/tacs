Beam and shell elements in TACS
*******************************

The beam and shell elements in TACS are designed to provide both linear and geometrically nonlinear analysis for static and transient analysis.
A key feature of these elements is their flexible parameterization.
In the case of beam elements, the rotational parametrization captures the change in orientation of the normal sections of the beam.
For shell elements, the rotational parametrization captures the change in orientation of the shell normal.
Small, moderate and finite rotations can be modeled.

Director parametrization
------------------------

A director is a vector that defines the orientation and length of the reference shell normal in the deformed configuration.
In the beam and shell elements implemented in TACS, the director itself is is not required, but can be evaluated via :math:`\mathbf{d} + \mathbf{n}`.
The vector field :math:`\mathbf{d}` defines the rate of change of displacement in a reference direction :math:`\mathbf{t}` induced by an exact or approximate rotation.
The relationship between :math:`\mathbf{d}` and the reference direction is

.. math::

    \mathbf{d} = \left(\mathbf{C}(\mathbf{q})^{T} - \mathbf{1}\right) \mathbf{t} = \mathbf{Q}(\mathbf{q}) \mathbf{t}

where :math:`\mathbf{C}(\mathbf{q})` is a rotation matrix parametrized by the variables :math:`\mathbf{q}`.
The reference direction :math:`\mathbf{t}` is generated from the initial, undeformed geometry of the beam or shell.

The spatial and temporal derivatives of :math:`\mathbf{d}` at each node are used to determine the potential and kinetic energies, respectively, of the beam and shell elements.
Here we focus on the kinetic energy term, since the spatial derivatives are a result of the interpolation between nodes.
The first and second time derivatives are required for computing the kinetic energy and governing equations.
When the matrix :math:`\mathbf{C}(\mathbf{q})` is an exact rotation matrix, the relationship :math:`\dot{\mathbf{C}} = -\omega^{\times} \mathbf{C}` holds.
Therefore :math:`\dot{\mathbf{C}}^{T} = \mathbf{C}^{T}\omega^{\times}` and the director satisfies

.. math::

    \begin{aligned}
    \dot{\mathbf{d}} & = \mathbf{C}^{T} \omega^{\times} \mathbf{t} \\
    \ddot{\mathbf{d}} & = \mathbf{C}^{T} \left( \dot{\omega}^{\times} + \omega^{\times}\omega^{\times} \right) \mathbf{t} \\
    \end{aligned}

There are several parametrizations of the field :math:`\mathbf{d}` that are implemented to enable linear and geometrically nonlinear analysis that are described below.

The equations of motion in TACS are derived from Lagrange's equations.
To provide a concrete explanation of the implementation, consider the kinetic energy expression

.. math::

    T(\dot{\mathbf{u}}_{0}, \dot{\mathbf{d}}) = \frac{1}{2} \left( m_{0} \dot{\mathbf{u}}_{0}^{T} \dot{\mathbf{u}}_{0} +
    2 m_{1} \dot{\mathbf{u}}_{0}^{T} \dot{\mathbf{d}} +
    m_{2} \dot{\mathbf{d}}^{T} \dot{\mathbf{d}} \right)

where :math:`\mathbf{u}_{0} \in \mathbb{R}^{3}` are the displacements from a reference point and :math:`m_{1}, m_{2}, m_{3}` are mass and inertia coefficients.
In this example, the potential energy expression takes the form

.. math::

    U = U(\mathbf{u}_{0}, \mathbf{u}_{0,\xi}, \mathbf{d}, \mathbf{d}_{,\xi})

where :math:`\mathbf{u}_{0}, \mathbf{u}_{0,\xi}, \mathbf{d}, \mathbf{d}_{,\xi}` represent the displacement values and their derivatives with respect to the parametrization of the beam or shell.

Assuming that the Lagrangian can be written as :math:`L = T(\dot{\mathbf{q}}, \dot{\mathbf{u}}_{0}, \mathbf{q}, \mathbf{u}_{0}) - U(\mathbf{q})`, the contributions to the governing equations are derived from Lagrange's equations, giving

.. math::

    \begin{aligned}
    \mathbf{R}_{u} =& \dfrac{d}{dt} \left( \dfrac{\partial T}{\partial \dot{\mathbf{u}}_{0}}\right) - \dfrac{\partial L}{\partial \mathbf{u}_{0}} \\
    %
    \mathbf{R}_{q} =& \dfrac{d}{dt} \left( \dfrac{\partial T}{\partial \dot{\mathbf{q}}}\right) - \dfrac{\partial L}{\partial \mathbf{q}} \\
    %
    =& \dfrac{d}{dt}\left( \dfrac{\partial T}{\partial \dot{\mathbf{d}}} \right) \dfrac{\partial \dot{\mathbf{d}}}{\partial \dot{\mathbf{q}}} +
    \dfrac{\partial T}{\partial \dot{\mathbf{d}}} \left( \dfrac{d}{dt} \left( \dfrac{\partial \dot{\mathbf{d}}}{\partial \dot{\mathbf{q}}} \right) -
    \dfrac{\partial \dot{\mathbf{d}}}{\partial \mathbf{q}} \right)
    - \dfrac{\partial L}{\partial \mathbf{d}} \dfrac{\partial \mathbf{d}}{\partial \mathbf{q}}\\
    %
    =& \dfrac{d}{dt}\left( \dfrac{\partial T}{\partial \dot{\mathbf{d}}} \right) \dfrac{\partial \dot{\mathbf{d}}}{\partial \dot{\mathbf{q}}} -
    \dfrac{\partial L}{\partial \mathbf{d}} \dfrac{\partial \mathbf{d}}{\partial \mathbf{q}}\\
    \end{aligned}


The full expression for the kinetic engery is integrated over the quadrature points in the element via the relationship

.. math::

    T = \sum_{k} w_{k} T_{k}( \sum_{i} N_{i}(\xi_{k}) \dot{\mathbf{u}}_{0i}, \sum_{i} N_{i}(\xi_{k}) \dot{\mathbf{d}}_{i})

A linear parametrization is based on a 3 parameter parameter model, where :math:`\mathbf{q} \in \mathbb{R}^{3}` that takes the form

.. math::

    \mathbf{Q}(\mathbf{q}) = {\mathbf{q}^{\times}}

.. math::

    \mathbf{\omega} = \dot{\mathbf{q}}

In this case, the parameters represent small rotations about the global coordinate axes.
A quadratic approximation of the rotation matrix can also be used, again with :math:`\mathbf{q} \in \mathbb{R}^{3}` where

.. math::

    \mathbf{Q}(\mathbf{q}) = \mathbf{q}^{\times} + \frac{1}{2} \mathbf{q}^{\times} \mathbf{q}^{\times}

This model enables a more accurate prediction of moderate rotations of the shell and beam elements.

An exact rotational parametrization is also available using a quaternion parametrization.
In this case, the parametrization utilizes a state vector :math:`\mathbf{q} \in \mathbb{R}^{5}` that contains both the quaternion variables and a Lagrange multiplier variable used to enforce the constraint.

.. math::

    \mathbf{Q}(\mathbf{q}) = 2 \epsilon^{\times} \epsilon^{\times} + 2 \eta \epsilon^{\times}

where :math:`\mathbf{q} = (\epsilon, \eta, \lambda)` and satisifies the constraint

.. math::

    \epsilon^{T}\epsilon + \eta^2 = 1

For the quaternion parametrization, this constraint is imposed at each node in the finite-element mesh.

Beam volume parametrization
---------------------------

The beam volume is parametrized by a single parameter which defines the reference line of the beam :math:`\xi \in \mathbb{R}`.
The beam reference line location is given as

.. math::

    \mathbf{X}_{0}(\xi) = \sum_{i} N_{i}(\xi) \mathbf{X}_{i}

where :math:`\mathbf{X}_{i}` are the element node locations and :math:`N_{i}(\xi)` are the element shape functions.

At each point in the beam element, a local reference frame is constructed by utilizing a reference direction, :math:`\mathbf{e}_{ref}`.
This reference direction must have a component perpendicular to the centerline of the beam.
The total


.. Kinetic energy for the shell
.. ----------------------------

.. The kinetic energy in the shell element is computed as

.. .. math::

..     T = \frac{1}{2} \int_{\Omega} \rho \dot{\mathbf{u}}^{T} \dot{\mathbf{u}} \mathbf{u} d\Omega

.. .. math::

..     T = \frac{1}{2} \int_{A} \rho_{0} \dot{\mathbf{u}}^{T}  \dot{\mathbf{u}} + \rho_{1} \omega^{T} \mathbf{J} \omega^{T} \; dA

.. here :math:`\mathbf{J} = \mathbf{1} - \mathbf{n}\mathbf{n}^{T}`.

.. We assume that the midsurface of the shell is co-located with the center of mass of the shell.

.. .. math::

..     \dot{\mathbf{C}} = - \omega^{\times} \mathbf{C}


Shell volume parametrization
----------------------------

The shell is parameterized by two coordinates which define in the mid-surface of the shell :math:`\xi = (\xi_{1}, \xi_{2})`.
The mid-surface of the shell is computed based on the element node locations and the element shape functions

.. math::

    \mathbf{X}_{0}(\xi) = \sum_{i} N_{i}(\xi) \mathbf{X}_{i}

where :math:`\mathbf{X}_{i}` are the element node locations and :math:`N_{i}(\xi)` are the shape functions.

The shell normal is computed based on the mid-surface tangents

.. math::

    \hat{\mathbf{n}} = \dfrac{\mathbf{X}_{0,\xi_{1}} \times \mathbf{X}_{0,\xi_{2}}}{||\mathbf{X}_{0,\xi_{1}} \times \mathbf{X}_{0,\xi_{2}}||_{2}}

The through-thickness volume of the shell is parametrized by interpolating the normal between points.
This interpolation enables an exact preservation of the rigid body rotations.
To form this interpolation, the surface normals are computed at the nodes of the finite element given by :math:`\hat{\mathbf{n}}_{i}`.
With these normal directions defined, the full parametrized volume is given as

.. math::

    \mathbf{X}(\eta) = \mathbf{X}_{0}(\xi) + \zeta \mathbf{n}(\xi) = \sum_{i} N_{i}(\xi)(\mathbf{X}_{i} + \zeta \hat{\mathbf{n}}_{i})

Here :math:`\zeta` is the through-thickness direction for the shell.
The mid-surface parameters and through thickness parameter are conveniently collected in the vector :math:`\eta = (\xi_{1}, \xi_{2}, \zeta)`.

The derivative of the position with respect to the volume parameterization :math:`\eta` is

.. math::

    \mathbf{X}_{,\eta} = \sum_{i} \begin{bmatrix} N_{i,\xi_1} (\mathbf{X}_{i} + \zeta \hat{\mathbf{n}}_{i}) &
    N_{i,\xi_2} (\mathbf{X}_{i} + \zeta \hat{\mathbf{n}}_{i}) &
    N_{i} \hat{\mathbf{n}}_{i} \end{bmatrix}

Note that this varies through the thickness of the shell.

The goal in the analysis of shell behavior is to reduce the response to data on the shell mid-surface.
The Jacobian transformation from derivatives with respect to the shell volume transformation at the mid-surface to the global coordinates is

.. math::

    \eta_{\mathbf{X}}^{0} = \left. \mathbf{X}_{,\eta}^{-1} \right|_{\zeta = 0} =
    \left[ \sum_{i} \begin{bmatrix}
    N_{i,\xi_1} \mathbf{X}_{i} &
    N_{i,\xi_2} \mathbf{X}_{i} &
    N_{i} \hat{\mathbf{n}}_{i} \end{bmatrix} \right]^{-1}

The Jacobian transformation varies through the thickness of the shell.
It is often required ot consider the rate of change of the Jacobian transformation through the thickness of the shell at the mid-surface

.. math::

    \eta_{\mathbf{X}\zeta}^{0} =  \left. \dfrac{\partial \mathbf{X}_{,\eta}^{-1}}{\partial \zeta} \right|_{\zeta = 0} =
    - \eta_{\mathbf{X}}^{0}
    \left[ \sum_{i}
    \begin{bmatrix}
    N_{i,\xi_1} \hat{\mathbf{n}}_{i} &
    N_{i,\xi_2} \hat{\mathbf{n}}_{i} & 0 \end{bmatrix} \right]
    \eta_{\mathbf{X}}^{0}

These quantities express the derivatives of the parameters with respect to the global coordinates.
Later, a transformation will be introduced to a local shell-oriented coordinate systen.

Displacement parametrization
----------------------------

The displacement field in the shell is parameterized using a combination of the mid-plane deflections and the through-thickness rate of deformation parameterized by :math:`\mathbf{d}`.
At each node in the shell element, :math:`\mathbf{d}` is parametrized based on the rotational variables at each node, :math:`\mathbf{q}_{i}`.
The field :math:`\mathbf{d}` gives the rate of change of the displacement in the through-thickness direction and is computed at each node :math:`i`

.. math::

    \mathbf{d}_{i}(\mathbf{q}_{i}) = (\mathbf{Q}^{T}(\mathbf{q}_{i}) - \mathbf{I})\mathbf{n}_{i}

The matrix :math:`\mathbf{Q}(\mathbf{q}_{i})` is either an exact or approximate rotation matrix.
Note that this matrix is only ever evaluated at the nodes and is never interpolated directly, only the :math:`\mathbf{d}` field itself is interpolated.

The displacement field is a combination of the mid-surface displacements at each node :math:`\mathbf{u}_{0i}` and the :math:`\mathbf{d}_{i}` values at each node

.. math::

    \mathbf{u}(\eta) = \sum_{i} N_{i}(\xi) \left( \mathbf{u}_{0i} + \zeta \mathbf{d}_{i}(\mathbf{q}_{i}) \right)

The gradient of the displacement field with respect to the parameters :math:`\eta` is required to compute the strain.
This gradient involves a combination of the in-plane and through-thickness parameters as follow

.. math::

    \mathbf{u}_{,\eta} = \sum_{i} \begin{bmatrix} N_{i,\xi_{1}} \left( \mathbf{u}_{0i} + \zeta \mathbf{d}_{i}(\mathbf{q}_{i}) \right) &
    N_{i,\xi_{2}} \left( \mathbf{u}_{0i} + \zeta \mathbf{d}_{i}(\mathbf{q}_{i}) \right) &
    N_{i} \mathbf{d}_{i}(\mathbf{q}_{i}) \end{bmatrix}

It will be important to consider the rate of change at the mid-surface of the shell as

.. math::

    \mathbf{u}^{0}_{,\eta} = \sum_{i} \begin{bmatrix} N_{i,\xi_{1}} \mathbf{u}_{0i} &
    N_{i,\xi_{2}} \mathbf{u}_{0i} &
    N_{i} \mathbf{d}_{i}(\mathbf{q}_{i}) \end{bmatrix}

The derivative depends on both

.. math::

    \mathbf{u}_{,\eta\zeta} = \sum_{i} \begin{bmatrix} N_{i,\xi_{1}} \mathbf{d}_{i}(\mathbf{q}_{i}) &
    N_{i,\xi_{2}} \mathbf{d}_{i}(\mathbf{q}_{i}) & 0 \end{bmatrix}

The derivative of the displacement with respect to the global coordinate system is

.. math::

    \mathbf{u}_{,\mathbf{X}} = \mathbf{u}_{,\eta} \mathbf{X}_{,\eta}^{-1}

This nonlinear expression is approximated using a linearization through the thickness as follows

.. math::

    \mathbf{u}_{,\mathbf{X}} \approx
    \mathbf{u}^{0}_{\mathbf{X}} + \zeta \mathbf{u}^{1}_{\mathbf{X}} =
    \mathbf{u}_{,\eta}^{0} \eta_{\mathbf{X}} + \zeta\left( \mathbf{u}_{,\eta\zeta}\eta_{\mathbf{X}}^{0} + \mathbf{u}_{,\eta}^{0} \eta_{\mathbf{X}\zeta}^{0} \right)

The zeroth and first order displacement gradient expressions are used to construct the shell-aligned strain expressions.

Transformation to local shell-attached frame
--------------------------------------------

At each point on the shell, we construct a transformation :math:`\mathbf{T}` that transforms the displacements from a local shell-attached reference frame to the global reference frame.
This transformation preserves the normal direction such that

.. math::

    \mathbf{T} \mathbf{e}_{3} = \hat{\mathbf{n}}

The transformation is computed at quadrature points and to visualize the results.

There are two methods that are implemented to compute the local shell transformation described below.
Both methods assemble the transformation by finding tangent directions :math:`\mathbf{t}_{1}` and :math:`\mathbf{t}_{2}`.
After these vectors have been computed, the full transformation matrix is

.. math::

    \mathbf{T} = \begin{bmatrix} \mathbf{t}_{1} & \mathbf{t}_{2} & \hat{\mathbf{n}} \end{bmatrix}

Reference axis projection transform
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first transformation utilizes a unit reference direction, denoted :math:`\mathbf{e}_{ref}`.
Note: *the reference direction cannot be normal to the shell surface*.
The reference direction is projected onto the surface of the shell to construct the local 1-direction.
This local reference direction is then combined with the normal to create the 2-direction.
The projection of the reference direction onto the shell surface takes the form

.. math::

    \mathbf{t}_{1} = \dfrac{\mathbf{e}_{ref} - \hat{\mathbf{n}}^{T}\mathbf{e}_{ref} \hat{\mathbf{n}}}{||\mathbf{e}_{ref} - \hat{\mathbf{n}}^{T}\mathbf{e}_{ref} \hat{\mathbf{n}}||_{2}}

The 2-direction is then computed by combining the reference direction with the surface normal to give the 2-direction :math:`\mathbf{t}_{2}`

.. math::

    \mathbf{t}_{2} = \hat{\mathbf{n}} \times \mathbf{t}_{1}.

Natural transform
^^^^^^^^^^^^^^^^^

The second transformation method utilizes the first tangent direction captured via the parametrization of the shell.
This tangent direction is always well defined and computed as

.. math::

    \mathbf{t}_{1} = \dfrac{\mathbf{X}_{0,\xi_{1}}}{|| \mathbf{X}_{0,\xi_{1}} ||_{2}}

The 2-direction is taken as the vector that completes the orthogonal basis

.. math::

    \mathbf{t}_{2} = \hat{\mathbf{n}} \times \mathbf{t}_{1}

Caution should be used when utilizing this transformation, since it will vary between shell elements depending on their orientation.
When the shell material is orthotropic the reference direction method should be used.

Strain expressions
------------------

The displacement gradient is transformed into the local reference frame as

.. math::

    \mathbf{u}_{,x} =
    \mathbf{u}_{,x}^{0} + \mathbf{u}_{,x}^{1} =
    \mathbf{T}^{T} \mathbf{u}^{0}_{\mathbf{X}} \mathbf{T} +
    \zeta \mathbf{T}^{T} \mathbf{u}^{1}_{\mathbf{X}} \mathbf{T}

The strain distribution throughout the shell is

.. math::

    \epsilon = \frac{1}{2} \left[ \mathbf{u}_{,x}^{0} + {\mathbf{u}_{,x}^{0}}^{T} + {\mathbf{u}_{,x}^{0}}^{T} \mathbf{u}_{,x}^{0} +
    \zeta \left( \mathbf{u}_{,x}^{1} + {\mathbf{u}_{,x}^{1}}^{T} +
    {\mathbf{u}_{,x}^{1}}^{T}\mathbf{u}_{,x}^{0} +
    {\mathbf{u}_{,x}^{0}}^{T}\mathbf{u}_{,x}^{1}\right) \right] +
    \mathcal{O}(\zeta^{2})

For analysis, the strain is split into the zeroth order and bending components.
The zeroth order strain terms consist of both in-plane normal and shear strains and out-of-plane shear strains

.. math::

    \epsilon^{0} = \frac{1}{2} \left[ \mathbf{u}_{,x}^{0} + {\mathbf{u}_{,x}^{0}}^{T} + {\mathbf{u}_{,x}^{0}}^{T} \mathbf{u}_{,x}^{0} \right]

For linear analysis, the zeroth order strains are

.. math::

    \epsilon^{0} = \frac{1}{2} \left[ \mathbf{u}_{,x}^{0} + {\mathbf{u}_{,x}^{0}}^{T} \right]

The bending strains consist of the normal and twisting bending components and are computed as

.. math::

    \kappa = \frac{1}{2} \left[ \mathbf{u}_{,x}^{1} + {\mathbf{u}_{,x}^{1}}^{T} +
    {\mathbf{u}_{,x}^{1}}^{T}\mathbf{u}_{,x}^{0} +
    {\mathbf{u}_{,x}^{0}}^{T}\mathbf{u}_{,x}^{1} \right]

For linear analysis, the bending strain components are

.. math::

    \kappa = \frac{1}{2} \left[ \mathbf{u}_{,x}^{1} + {\mathbf{u}_{,x}^{1}}^{T} \right]

Thermal strain formulation
--------------------------

For thermoelastic analysis, we add a scalar temperature variable at each node, :math:`\theta`.
The scalar field is invariant under the local shell frame transformation.
However, the derivatives of the temperature with respect to the local coordinates depend on this transofrmation.
The temperature field gradient is

.. math::

    \theta_{,x} = \begin{bmatrix} \theta_{,\xi_{1}} & \theta_{,\xi_{1}} & 0 \end{bmatrix} \mathbf{X}_{,\eta}^{-1} \mathbf{T}

The change in temperature causes a strain due to thermal expansion

.. math::

    \epsilon_{T} = \theta \begin{bmatrix}
    \alpha_{xx}(\zeta) & \alpha_{xy}(\zeta) & 0 \\
    \alpha_{xy}(\zeta) & \alpha_{yy}(\zeta) & 0 \\
    0 & 0 & 0 \\
    \end{bmatrix}

where :math:`\alpha_{xx}`, :math:`\alpha_{xy}` and :math:`\alpha_{yy}` are thermal coefficients of expansion from lamination theory.
Note that the term :math:`\alpha_{xy}` arises due to the transformation between material reference frame and the shell-aligned local reference frame.
In addition, there are coupling terms that arise due to the dependence of the thermal strain through-thickness :math:`\zeta`.

Under thermoelastic analysis, the in-plane mechanical strain for linear elements is

.. math::

    \epsilon^{0} = \frac{1}{2} \left[ \mathbf{u}_{,x}^{0} + {\mathbf{u}_{,x}^{0}}^{T} \right] - \epsilon_{T}^{0}

The bending components of the mechanical strain for the linear elements is

.. math::

    \kappa = \frac{1}{2} \left[ \mathbf{u}_{,x}^{1} + {\mathbf{u}_{,x}^{1}}^{T} \right] - \kappa_{T}

Here :math:`\epsilon_{T}^{0}` and :math:`\kappa_{T}` are the components of the strain due to thermal expansion from lamination theory.

Drilling rotation
-----------------

The rotation of the shell about the shell normal is called the drill rotation.
In this formulation, we add a penalization between the rotation normal to the shell and the rotation computed from the in-plane rotation of the displacement.
This penalization adds stiffness to the shell.
The value of the penalization is taken from the shell constitutive object.

Given the rotation matrix :math:`\mathbf{C}(q_{i})`, at each node, the rotation penalty term is computed as

.. math::

    \epsilon_{t} =
    \mathbf{e}_{2}^{T} \mathbf{T}^{T} \left[ \sum_{i} N_{i} \mathbf{C}(\mathbf{q}_{i}) \right] \mathbf{T} \left(\mathbf{e}_{1} + \mathbf{u}_{,x}^{0} \mathbf{e}_{1} \right) -
    \mathbf{e}_{1}^{T} \mathbf{T}^{T} \left[ \sum_{i} N_{i} \mathbf{C}(\mathbf{q}_{i}) \right] \mathbf{T} \left(\mathbf{e}_{2} + \mathbf{u}_{,x}^{0} \mathbf{e}_{2} \right)

here :math:`\mathbf{e}_{1}` and :math:`\mathbf{e}_{2}` denote the cartesian basis, and :math:`\mathbf{u}_{,x}^{0}` is the derivative of the mid surface displacements in the locally attached reference frame.
This deviation is treated by adding a strain energy penalty term to the total potential energy of the element :math:`\frac{1}{2} k_{t} \epsilon_{t}^2`.

In the case of the linear rotation matrix :math:`\mathbf{C}(\mathbf{q}) = \mathbf{1} - \mathbf{q}^{\times}` for :math:`\mathbf{q} \in \mathbb{R}^{3}`.
Linearizing the expression for :math:`\epsilon_{t}` gives

.. math::

    \begin{aligned}
    \epsilon_{t} &=
    \mathbf{t}_{2}^{T} (\mathbf{1} -  \mathbf{q}(\xi)^{\times}) \mathbf{t}_{1} + \mathbf{t}_{2}^{T} \mathbf{T} \mathbf{u}_{,x}^{0} \mathbf{e}_{1}
    - \mathbf{t}_{1}^{T} (\mathbf{1} -  \mathbf{q}(\xi)^{\times}) \mathbf{t}_{2} - \mathbf{t}_{1}^{T} \mathbf{T} \mathbf{u}_{,x}^{0} \mathbf{e}_{2} \\
    &= \mathbf{e}_{2} \mathbf{u}_{,x}^{0} \mathbf{e}_{1} - \mathbf{e}_{1} \mathbf{u}_{,x}^{0} \mathbf{e}_{2} - 2 \hat{\mathbf{n}}^{T} \mathbf{q}(\xi) \\
    \end{aligned}

Note that for a plate in the :math:`x-y` plane this simplifies to the relationship :math:`\epsilon_{t} = v_{,x} - u_{,y} - 2 q_{z}`.

Mixed Interpolation of Tensorial Components
-------------------------------------------

Shell and beam elements can suffer from locking behavior where the predictive capability of the shell or beam elements suffers.
This locking phenomena is due to an inability of some elements to capture pure bending behavior without producing shear artificially.
To alleviate shear and in-plane locking behavior, the shell and beam elements in TACS utilize an mixed interpolation of tensorial components (MITC) formulation.
This formulation naturally extends to higher-order element implementations.

The MITC approach works by evaluating the displacement-based expressions for the strain at tying points within the element.
These strain components are then interpolated across the element with an assumed strain distribution.
When selected appropriately, the modified element exhibits locking-free behavior.

The tensorial components of the strain are interpolated within the element.
In this context, the interpolated tensorial components are given by the zeroth order strain terms and are

.. math::

    \tilde{\epsilon} = \frac{1}{2}\left( \mathbf{X}_{,\eta}^{T}\mathbf{u}_{,\eta} + \mathbf{u}_{,\eta}^{T} \mathbf{X}_{,\eta}
    + \mathbf{u}_{,\eta}^{T}\mathbf{u}_{,\eta} \right)

The tensorial components of the strain can be transformed to the Green strain in the global coordinate systen using

.. math::

    \epsilon = \mathbf{X}_{,\eta}^{-T} \tilde{\epsilon} \mathbf{X}_{,\eta}^{-1}

The tying points are given by the parametric points :math:`\eta_{t}`, and the interpolation for the strain is given by the basis :math:`N^{as}(\xi)`.
With these definitions, the zeroth order strain components can be computed from the tying strain values as

.. math::

    \epsilon_{as}^{0} = \mathbf{T}^{T} {\eta_{\mathbf{X}}^{0}}^{T} \left[ \sum_{t} N^{as}_{t}(\xi) \tilde{\epsilon}(\xi_{t}) \right] \eta_{\mathbf{X}}^{0} \mathbf{T}

Constitutive relationships for the shell element
------------------------------------------------

For the shell element, the constitutive relationship is


Equations of motion
-------------------


  The equations of motion for this element are derived using
  Lagrange's equations with constraints. Each node imposes a
  constraint that its own quaternions satisfy the required unit norm
  constraint.

  The equations of motion can be divided into two parts: (1) the
  motion of the deformed surface and (2) the motion of the normals of
  the surface (or directors since they are no longer normal to the
  deformed surface during deformation). The linear motion takes the
  form:

  M*ddot{u} + dU/dx - fg = 0

  where M is a mass matrix. The rotational degrees of freedom satisfy
  the following equations of motion:

  S^{T}*J*d{omega} + 2*dot{S}^{T}*J*omega + dU/dq - A^{T}*lamb = 0

  where J = (I - n*n^{T}) is a rotational inertia term and the
  constraints produce the term A^{T}*lamb.




Director implementation
***********************



Beam element implementation
***************************




Shell element implementation
****************************

The shell element implementation consists of the following

1. A shell element basis that defines the shape functions and the mixed interpolation functions required for the strain interpolation.
2. A director parametrization that computes the director field as a function of the element variables.
3. A transformation that computes the local shell element coordinates
4. A constitutive object that computes the stress resultants as a function of the strains.

Shell element basis
-------------------

The shell element basis handles the element parametrization and quadrature points.
It inherits from the
computes the local tangents at each node in the mesh

.. code-block:: cpp

    // Given the nodal coordinates, Xpts, compute the shell coordinate frame at each
    // node in the element and store each one in Xf.
    virtual void computeNodalFrames( const TacsScalar Xpts[], TacsScalar Xf[] ) = 0;

    // Given a parametric point in the element, typically a quadrature point,
    // evaluate the local frame Xd
    void interpolateFrame( const int n, const double pt[],
                           const TacsScalar Xpts[], const TacsScalar Xf[],
                           TacsScalar Xd[], TacsScalar Xdz[] );

    // Compute the parametric derivatives of the displacement field and director
    // fields
    virtual void computeDeriv( const int npts, const double pts[],
                               )

Director field parametrization
------------------------------

.. code-block:: cpp

    virtual void computeDirectors( const int vars_per_node, const int offset,
                                   const int num_nodes, const TacsScalar Xf[],
                                   const TacsScalar vars[], TacsScalar dirs[] );



Transformation
--------------

The following computes the transformations at each of the quadrature points in the element

.. code-block:: cpp

    virtual void computeTransform( const TacsScalar Xd[], TacsScalar T[] ) = 0;


Strain computation
------------------

.. code-block:: cpp

    // Compute the natural curvilinear reference frame at each node
    TacsScalar Xf[9*nnodes];
    computeNodalFrames(Xpts, Xf);

    // Interpolate the frame to the parametric point
    TacsScalar Xd[9], Xdz[9];
    interpolateFrame(n, pt, Xpts, Xf, Xd, Xdz);

    // Compute the transformation at the node
    TacsScalar T[9];
    computeTransform(Xd, T);

    // Compute the inverse of the 3x3 transformation
    TacsScalar Xdinv[9];
    TacsScalar detXd = inv3x3(Xd, Xdinv);

    //
    TacsScalar zXdinv[9], tmp[9]
    // zXdinv = -Xdinv*Xdz*Xdinv


    // Compute the transformation ux0 = T*ueta*Xdinv*T^{T}
    // u0x = T*u0d*Xdinv*T^{T}
    TacsScalar u0x[9];
    matMatMult(u0d, Xdinv, u0x);
    matMatMult(u0d, T, tmp);
    matTransMatMult(T, tmp, u0x);

    // u1x = T*(u0d*zXdinv + u1d*Xdinv)*T^{T}
    TacsScalar u1x[9];
    matMatMult(u0d, zXdinv, u1x);
    matMatMultAdd(u1d, Xdinv, u1x);
    matMatMult(u1x, T, tmp);
    matTransMatMult(T, tmp, u1x);


