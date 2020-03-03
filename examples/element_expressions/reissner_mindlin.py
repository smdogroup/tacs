"""
Generate the coefficients for the Reissner Mindlin plate
"""

import sympy as sym

def str_coef(coef):
    c = str(coef).replace('1.0*', '')
    c = c.replace('+ 1', '+ 1.0')
    c = c.replace('+ 2', '+ 2.0')
    return c

def print_coef(coef):
    print(str_coef(coef) + ';')

def get_disps():
    Ut = []
    for i in range(3*5):
        Ut.append(sym.Symbol('Ut[%d]'%(i)))

    Ux = []
    for i in range(6*2):
        Ux.append(sym.Symbol('Ux[%d]'%(i)))

    return Ut, Ux

def get_strains(Ut, Ux):
    e = []
    e.append(Ux[0]) # exx = u,x
    e.append(Ux[3]) # eyy = v,y
    e.append(Ux[1] + Ux[2]) # gxy = u,y + v,x

    e.append(Ux[8]) # kxx = roty,x
    e.append(-Ux[7]) # kyy = -rotx,y
    e.append(Ux[9] - Ux[6]) # kxy = roty,y - rotx,x

    e.append(Ux[5] - Ut[9]) # eyz = w,y - rotx
    e.append(Ux[4] + Ut[12]) # exz = w,x + roty

    return e

def get_stress_components():
    s = []
    for i in range(8):
        s.append(sym.Symbol('s[%d]'%(i)))
    return s

def get_stress(e):
    A = []
    B = []
    D = []
    As = []

    for i in range(6):
        A.append(sym.Symbol('A[%d]'%(i)))
        B.append(sym.Symbol('B[%d]'%(i)))
        D.append(sym.Symbol('D[%d]'%(i)))

    for i in range(3):
        As.append(sym.Symbol('As[%d]'%(i)))

    s = []
    s.append(A[0]*e[0]+A[1]*e[1]+A[2]*e[2] + B[0]*e[3]+B[1]*e[4]+B[2]*e[5])
    s.append(A[1]*e[0]+A[3]*e[1]+A[4]*e[2] + B[1]*e[3]+B[3]*e[4]+B[4]*e[5])
    s.append(A[2]*e[0]+A[4]*e[1]+A[5]*e[2] + B[2]*e[3]+B[4]*e[4]+B[5]*e[5])

    s.append(B[0]*e[0]+B[1]*e[1]+B[2]*e[2] + D[0]*e[3]+D[1]*e[4]+D[2]*e[5])
    s.append(B[1]*e[0]+B[3]*e[1]+B[4]*e[2] + D[1]*e[3]+D[3]*e[4]+D[4]*e[5])
    s.append(B[2]*e[0]+B[4]*e[1]+B[5]*e[2] + D[2]*e[3]+D[4]*e[4]+D[5]*e[5])

    s.append(As[0]*e[6]+As[1]*e[7])
    s.append(As[1]*e[6]+As[2]*e[7])

    return s

def print_plate_coef():
    Ut, Ux = get_disps()
    e = get_strains(Ut, Ux)
    s = get_stress(e)

    S = get_stress_components()

    U0 = e[0]*S[0]
    for i in range(1, 8):
        U0 = U0 + e[i]*S[i]

    for k in range(3*5):
        coef = sym.diff(U0, Ut[k])
        print('DUt[%d] = '%(k) + str_coef(coef) + ';')

    for k in range(10):
        coef = sym.diff(U0, Ux[k])
        print('DUx[%d] = '%(k) + str_coef(coef) + ';')

    U0 = e[0]*s[0]
    for i in range(1, 8):
        U0 = U0 + e[i]*s[i]

    U0 += sym.Symbol('rho')*(Ut[2]**2 + Ut[5]**2 + Ut[8]**2)

    for k in range(3*5):
        coef = 0.5*sym.diff(U0, Ut[k])
        print('DUt[%d] = '%(k) + str_coef(coef) + ';')

    for k in range(10):
        coef = 0.5*sym.diff(U0, Ux[k])
        print('DUx[%d] = '%(k) + str_coef(coef) + ';')

    # Arrange the variables in their order
    elem_vars = []
    for i in range(5):
        for j in range(3):
            elem_vars.append(Ut[3*i + j])
        for j in range(2):
            elem_vars.append(Ux[2*i + j])

    s = ''
    index = 0
    for i, ui in enumerate(elem_vars):
        for j, uj in enumerate(elem_vars):
            coef = 0.5*sym.diff(sym.diff(U0, ui), uj)
            if coef != 0:
                print('Jac[%d] = '%(index) + str_coef(coef) + ';')
                s += '%d, %d, '%(i, j)
                index += 1
        s += '\n'

    print('num_pairs = ', index)
    print('pairs = ', s)

    return


print_plate_coef()
