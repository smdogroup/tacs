"""
The following script generates the coefficients for the linear and nonlinear elasticity
element implementations. This uses sympy for symbolic manipulation.

This makes extensive use of the trace. Note that by the transpose property of the
trace (since the matrix and its transpose have the same trace)

tr(A^{T}*B) = tr(B^{T}*A)

And by the cyclic permutation property of the trace

tr(A^{T}*B) = tr(B*A^{T}) = tr(A*B^{T})

The integrand in the principle of virtual work can be written as

tr(S*dE) - dot(dU,f)

where S is the stress tensor, dE is the first variation of the strain tensor,
dU is the first variation of the displacement and dU is the variation of the
displacements.

Our goal is to find the coefficients of dU and dUx (the variation of U and
the variation of the derivative of U w.r.t. x). This can be achieved with the
trace by finding the coefficients such that

tr(coef*dUx^{T}) = tr(coef^{T}*dUx) = coef_{ij}*dUx_{ij}

Examining the first term, we find that

tr(S*dE)
= 0.5*tr(S*d(Ux + Ux^{T} + Ux^{T}*Ux))
= 0.5*tr(S*(dUx + dUx^{T} + dUx^{T}*Ux + Ux^{T}*dUx))
= tr(S*dUx) + 0.5*tr(S*dUx^{T}*Ux) + 0.5*(S*Ux^{T}*dUx)
= tr(S*dUx) + 0.5*tr(Ux*S*dUx^{T}) + (S*Ux^{T}*dUx)
= tr(S*dUx) + 0.5*tr(S*Ux^{T}*dUx)
= tr(S*(I + Ux^{T})*dUx)
= tr(coef^{T}*dUx)

Therefore:

coef = (S*(I + Ux^{T}))^{T} = (I + Ux)*S
"""

import sympy as sym

def get_const(dim=3):
    """Create a symmetric constitutive matrix with the specified dimensions"""
    index = 0
    C = []
    for i in range(dim):
        C.append([])
        for j in range(0, i):
            C[i].append(C[j][i])
        for j in range(i, dim):
            C[i].append(sym.Symbol('C[%d]'%(index)))
            index += 1
    return sym.Matrix(C)

def str_coef(coef):
    c = str(coef).replace('1.0*', '')
    c = c.replace('+ 1', '+ 1.0')
    c = c.replace('+ 2', '+ 2.0')
    return c

def print_coef(coef):
    print(str_coef(coef) + ';')

def get_3d_stress(s='s'):
    """Create a symbolic stress matrix for 3D"""
    s0 = sym.Symbol('%s[0]'%(s))
    s1 = sym.Symbol('%s[1]'%(s))
    s2 = sym.Symbol('%s[2]'%(s))
    s3 = sym.Symbol('%s[3]'%(s))
    s4 = sym.Symbol('%s[4]'%(s))
    s5 = sym.Symbol('%s[5]'%(s))
    S = sym.Matrix([[s0, s5, s4],
                    [s5, s1, s3],
                    [s4, s3, s2]])
    return S

def get_3d_disp_grad():
    """Create a symbolic displacement gradient for 3D"""
    index = 0
    Ux = []
    for i in range(3):
        row = []
        for j in range(3):
            row.append(sym.Symbol('Ux[%d]'%(index)))
            index += 1
        Ux.append(row)
    return sym.Matrix(Ux)

def get_3d_adjoint_grad(phi=False):
    """Create a symbolic adjoint gradient for 3D"""
    index = 0
    Psix = []
    for i in range(3):
        row = []
        for j in range(3):
            if phi:
                row.append(sym.Symbol('Phix[%d]'%(index)))
            else:
                row.append(sym.Symbol('Psix[%d]'%(index)))
            index += 1
        Psix.append(row)
    return sym.Matrix(Psix)

def get_3d_strains(Ux=None):
    """Get the 3D strains"""
    if Ux is None:
        Ux = get_3d_disp_grad()
    Ex = 0.5*(Ux + Ux.T + Ux.T*Ux)
    return sym.Matrix([Ex[0,0], Ex[1,1], Ex[2,2], 2*Ex[1,2], 2*Ex[0,2], 2*Ex[0,1]])

def print_3d_linear_coef():
    """Write out the coefficients and Jacobian matrix for the 3D nonlinear equations of elasticity"""
    voigt = ((0,0), (1,1), (2,2), (1,2), (0,2), (0,1))
    Ux = get_3d_disp_grad()
    Ex = 0.5*(Ux + Ux.T)
    e = sym.Matrix([Ex[0,0], Ex[1,1], Ex[2,2], 2*Ex[1,2], 2*Ex[0,2], 2*Ex[0,1]])
    C = get_const(dim=6)
    s = C*e

    # Compute the expression for the coefficients Coef = S + 0.5*(S*U,x^{T} + U,x*S)
    coef = sym.Matrix([[s[0], s[5], s[4]], [s[5], s[1], s[3]], [s[4], s[3], s[2]]])

    # Compute the derivative of the coefficients w.r.t. Ux
    index = 0
    for arg in ((0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2)):
        for i in range(3):
            for j in range(3):
                expr = sym.diff(coef[arg], Ux[i,j])
                print('Jac[%d] = '%(index) + str_coef(expr)  + ';')
                index += 1
    return

def print_3d_coef():
    """Write out the coefficients and Jacobian matrix for the 3D nonlinear equations of elasticity"""
    all = ((0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2))
    voigt = ((0,0), (1,1), (2,2), (1,2), (0,2), (0,1))
    S = get_3d_stress()
    Ux = get_3d_disp_grad()

    coef = S + Ux*S

    print('3D expressions:')
    for k, args in enumerate(all):
        coef[args] = sym.collect(coef[args], S[args])
        print('DUx[%d] = '%(k) + str_coef(coef[args]) + ';')

    e = get_3d_strains(Ux=Ux)
    C = get_const(dim=6)
    s = C*e

    index = 0
    for k in range(6):
        for i in range(3):
            for j in range(3):
                print('dsdUx[%d] = '%(index) + str_coef(sym.simplify(sym.diff(s[k], Ux[i,j]))) + ';')
                index += 1

    # Compute the derivative of the coefficients w.r.t. Ux
    index = 3
    for arg in all:
        for i in range(3):
            for j in range(3):
                line = str(coef[arg])
                for a, ar in enumerate(voigt):
                    line = line.replace(str(S[ar]), 'dsdUx[%d]'%(9*a + 3*i + j))
                extra = sym.diff(coef[arg], Ux[i,j])
                if extra != 0:
                    line += ' + ' + str(extra)

                print('Jac[%d] = '%(index) + str_coef(line)  + ';')
                index += 1
    return

def print_3d_adjoint_coef():
    """Print coefficients for the adjoint-residual product"""
    voigt = ((0,0), (1,1), (2,2), (1,2), (0,2), (0,1))

    S = get_3d_stress()
    Ux = get_3d_disp_grad()
    Psix = get_3d_adjoint_grad()
    coef = S + Ux*S
    mat = Psix.T*coef
    expr = mat[0,0] + mat[1,1] + mat[2,2]

    for index, arg in enumerate(voigt):
        coef = sym.simplify(sym.diff(expr, S[arg]))
        print('phi[%d] = '%(index) + str_coef(coef) + ';')

def print_3d_xpt_adjoint_coef():
    voigt = ((0,0), (1,1), (2,2), (1,2), (0,2), (0,1))

    # Compute the derivatives of e^{T}*C*phi with respect to Psix and Ux.
    # This code is required for the derivative of the adjoint-residual
    # product w.r.t. the nodal locations.

    S = get_3d_stress()
    A = get_3d_stress(s='a')
    Ux = get_3d_disp_grad()
    e = get_3d_strains(Ux=Ux)
    Psix = get_3d_adjoint_grad()
    coef = S + Ux*S
    mat = Psix.T*coef
    expr = mat[0,0] + mat[1,1] + mat[2,2]

    phi = []
    for index, arg in enumerate(voigt):
        phi.append(sym.simplify(sym.diff(expr, S[arg])))
    phi = sym.Matrix(phi)

    # s = C*e
    result = (phi[0]*S[0,0] + phi[1]*S[1,1] + phi[2]*S[2,2] +
              phi[3]*S[1,2] + phi[4]*S[0,2] + phi[5]*S[0,1])
    index = 0
    for i in range(3):
        for j in range(3):
            coef = sym.simplify(sym.diff(result, Psix[i,j]))
            print('dfdPhix[%d] = '%(index) + str_coef(coef) + ';')
            index += 1

    # s = C*e
    result = (phi[0]*S[0,0] + phi[1]*S[1,1] + phi[2]*S[2,2] +
              phi[3]*S[1,2] + phi[4]*S[0,2] + phi[5]*S[0,1] +
              e[0]*A[0,0] + e[1]*A[1,1] + e[2]*A[2,2] +
              e[3]*A[1,2] + e[4]*A[0,2] + e[5]*A[0,1])
    index = 0
    for i in range(3):
        for j in range(3):
            coef = sym.simplify(sym.diff(result, Ux[i,j]))
            print('dfdUx[%d] = '%(index) + str_coef(coef) + ';')
            index += 1

def print_3d_geo_stiffness_coef():
    """Print coefficients for the adjoint-residual product"""
    voigt = ((0,0), (1,1), (2,2), (1,2), (0,2), (0,1))

    S = get_3d_stress()
    Ux = get_3d_disp_grad()
    mat = 0.5*S*(Ux + Ux.T + Ux.T*Ux)
    trace = mat[0,0] + mat[1,1] + mat[2,2]

    Psi = get_3d_adjoint_grad()
    Phi = get_3d_adjoint_grad(phi=True)

    # Compute the derivative of the coefficients w.r.t. Ux
    index = 0
    expr = None
    for arg in ((0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2)):
        for i in range(3):
            for j in range(3):
                deriv = sym.diff(sym.diff(trace, Ux[i,j]), Ux[arg])
                print('Jac[%d] = '%(index) + str_coef(deriv) + ';')

                index += 1
                if expr is None:
                    expr = deriv*Psi[arg]*Phi[i,j]
                else:
                    expr += deriv*Psi[arg]*Phi[i,j]

    for index, arg in enumerate(voigt):
        deriv = sym.diff(expr, S[arg])
        print('psi[%d] = '%(index) + str_coef(deriv) + ';')

def print_3d_strain_deriv():
    """Print coefficients for the adjoint-residual product"""
    voigt = ((0,0), (1,1), (2,2), (1,2), (0,2), (0,1))

    S = get_3d_stress(s='t1')
    Ux = get_3d_disp_grad()
    Psix = get_3d_adjoint_grad()
    mat = 0.5*S*(Ux + Ux.T + Ux.T*Ux)
    expr = mat[0,0] + mat[1,1] + mat[2,2]

    index = 0
    for i in range(3):
        for j in range(3):
            coef = sym.simplify(sym.diff(expr, Ux[i,j]))
            print('dfdUx[%d] = '%(index) + str_coef(coef) + ';')
            index += 1

print_3d_coef()
print_3d_linear_coef()
print_3d_adjoint_coef()
print_3d_xpt_adjoint_coef()
print_3d_geo_stiffness_coef()
print_3d_strain_deriv()
