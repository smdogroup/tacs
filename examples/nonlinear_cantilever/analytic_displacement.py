"""
==============================================================================

==============================================================================
@Description : Functions for computing the geometrically exact tip displacement
of a cantilever beam with an applied tip force
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================

# ==============================================================================
# Extension modules
# ==============================================================================


def analyticCantileverDisplacement(alpha):
    if len(alpha) > 1:
        thetaTip = np.zeros_like(alpha)
        XDispNorm = np.zeros_like(alpha)
        ZDispNorm = np.zeros_like(alpha)
    else:
        alpha = np.array(alpha)
    for a in range(len(alpha)):
        if alpha[a] != 0.0:
            # Use the linear cantilever tip rotation as a good starting guess for the nonlinear rotation
            ResidualFun = lambda theta: tipRotationResidual(theta, alpha[a])
            sol = root_scalar(
                ResidualFun,
                x0=1e-4,
                x1=min(alpha[a] / 2.0, 0.5),
                xtol=1e-12,
                rtol=1e-14,
            )
            theta = sol.root

            mu = (1.0 + np.sin(theta)) / 2.0
            phi = np.arcsin(1.0 / (np.sqrt(2.0 * mu)))
            XDisp = 1.0 - np.sqrt((2.0 * np.sin(theta)) / alpha[a])
            ZDisp = 1.0 - 2.0 / np.sqrt(alpha[a]) * (ellipe(mu) - ellipeinc(phi, mu))
        else:
            theta = 0.0
            XDisp = 0.0
            ZDisp = 0.0

        if len(alpha) > 1:
            thetaTip[a] = theta
            XDispNorm[a] = XDisp
            ZDispNorm[a] = ZDisp
        else:
            thetaTip = theta
            XDispNorm = XDisp
            ZDispNorm = ZDisp

    return thetaTip, XDispNorm, ZDispNorm


# This function computes the residual of the transcendental equation that is solved to compute the tip rotation angle theta
def tipRotationResidual(theta, alpha):
    mu = (1.0 + np.sin(theta)) / 2.0
    phi = np.arcsin(1.0 / (np.sqrt(2.0 * mu)))
    R = ellipk(mu) - ellipkinc(phi, mu) - np.sqrt(alpha)
    print(
        "Alpha = %f, theta = %f, mu = %f, phi = %f, R = %f" % (alpha, theta, mu, phi, R)
    )
    return R
