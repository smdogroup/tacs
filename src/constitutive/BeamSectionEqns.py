import sympy as sp

# IMPORTANT: strain sign convention.
# The TACS beam element (TACSBeamElementModel.h::evalStrain) computes the
# bending strains as the axial derivatives of the directors:
#     e[2] = d1x[0] = (kappa x e_y)_x = -kappa_z   (director d1 ~ local y)
#     e[3] = d2x[0] = (kappa x e_z)_x = +kappa_y   (director d2 ~ local z)
# so the e[2] slot carries -kappa_z, NOT the textbook +kappa_z. Here xi3 is the
# symbol occupying that e[2] slot, hence xi3 = -kappa_z and xi2 = +kappa_y. Every
# expression below is written in this element convention so the resulting C
# matches populateMats directly (no ad-hoc xc2 negation required downstream).
e11, xi1, xi2, xi3, y12, y13 = sp.symbols("e11, xi1, xi2, xi3, y12, y13")
elementStrains = sp.Array([e11, xi1, xi3, xi2, y12, y13])

xo2, xo3 = sp.symbols("xo2 xo3")
offset = sp.Array([xo2, xo3])

# Beam strains at offset point. Physical axial strain at (y=xo2, z=xo3) is
#   eps = e11 + z*kappa_y - y*kappa_z = e11 + xo3*xi2 - xo2*(-xi3),
# and with xi3 = -kappa_z this becomes e11 + xo3*xi2 + xo2*xi3.
e11c = e11 + xo3 * xi2 + xo2 * xi3
xi1c = xi1
xi2c = xi2
xi3c = xi3
y12c = y12 - xo3 * xi1
y13c = y13 + xo2 * xi1
offsetStrains = sp.Array([e11c, xi1c, xi3c, xi2c, y12c, y13c])

# --- Express the strain transformation as a matrix ---
# sympy uses the wrong convention for multidimensional derivatives, so the strain transformation is the transpose of the sympy derivative of the beam strains w.r.t the shell strains
strainTransform = sp.permutedims(sp.diff(offsetStrains, elementStrains), (1, 0))
print("\nstrainTransform = ")
for row in strainTransform:
    print(row)
strainTransform = strainTransform.tomatrix()

xc2, xc3 = sp.symbols("xc2 xc3")
xk2, xk3 = sp.symbols("xk2 xk3")
# Compute transformation matrices for centroid and shear center
centroidTransform = strainTransform.subs({xo2: xc2, xo3: xc3})
shearCenterTransform = strainTransform.subs({xo2: xk2, xo3: xk3})

# Stiffness terms
EA, EI22, EI33, EI23, GJ, kGA22, kGA33, kGA23 = sp.symbols(
    "EA, EI22, EI33, EI23, GJ, kGA22, kGA33, kGA23"
)

# Axial and bending stiffness matrix at centroid.
# Note that the EI23 terms are positive (rather than negative as most texts define) because of the same xi3 = -kappa_z
# sign convention mentioned above.
axialBendingMat = sp.Matrix(
    [
        [EA, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, EI33, EI23, 0, 0],
        [0, 0, EI23, EI22, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
    ]
)

# Shear and torsion stiffness matrix at shear center
shearTorsionMat = sp.Matrix(
    [
        [0, 0, 0, 0, 0, 0],
        [0, GJ, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, kGA22, -kGA23],
        [0, 0, 0, 0, -kGA23, kGA33],
    ]
)

totalStiffness = sp.simplify(
    centroidTransform.T @ axialBendingMat @ centroidTransform
    + shearCenterTransform.T @ shearTorsionMat @ shearCenterTransform
)

flatInd = 0
for ii in range(totalStiffness.shape[0]):
    for jj in range(ii, totalStiffness.shape[1]):
        if totalStiffness[ii, jj] != 0:
            print(f"C[{flatInd}] = {totalStiffness[ii, jj]}")
        flatInd += 1
