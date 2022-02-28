:orphan:



    The weak form consists of two groups of components, the coefficients
    of time-dependent terms (up to second-order), and the coefficients of
    the spatial derivative terms (only first-order).

    Note that we assume separability between the spatial derivatives and the
    temporal derivatives, so that DUt[] does not depend on Ux, and DUx does
    not depend on Udot or Uddot.

    The parameter *DUt* contains the time coefficients in the weak form in
    a matrix of size (vars_per_node x 3). Each column in the matrix represents
    the zero-th, first and second time derivatives with the rows representing
    each variable. Therefore, the weak form for a problem with the variable
    components U = (u, v) will have the following form:

    int_{Area} (DUt[0]*du + DUt[1]*d(dot{u}) + DUt[2]*d(ddot{u}) +
                DUt[3]*dv + DUt[4]*d(dot{v}) + DUt[5]*d(ddot{v}) +
                spatial terms) dA = 0

    The parameter *DUx* contains the spatial derivative components of the
    weak form in a matrix of size (vars_per_node x (spatial_dim + 1)).
    The first component represents the coefficient of the variable, while
    the second, third and possibly fourth component represent the remaining
    spatial derivative coefficients. A problem with the variable
    components U = (u, v) with a spatial dimension of two will have the
    following weak form:

    int_{Area} (time dependent terms +
                DUx[0]*du + DUx[1]*d(du/dx) + DUx[2]*d(du/dy)) +
                DUx[3]*dv + DUx[4]*d(dv/dx) + DUx[5]*d(dv/dy)) dA = 0

    Note that the coefficients DUt[0] and DUx[0], and DUt[3] and DUx[3],
    are both for the displacement u and v, respectively. This means that
    the implementation is not unique.



    The following code computes the weak form coefficients and their
    derivatives with respect to each of the input components. The
    descriptions of the terms DUt and DUx are the same as the
    evalWeakIntegrand() function described above.

    The parameter Jac contains a sparse matrix representation of the
    the derivatives of the coefficients in DUt and DUx. The matrix
    contains (3 + spatial_dim)*vars_per_node rows and columns.

    For instance, for the 2D problem (spatial_dim = 2) with the variables
    U = (u, v), the Jac matrix would contain 10 x 10 entries. The rows of the
    matrix (corresponding to DUt and DUx) are ordered first by variable, then
    by derivative. The columns of the matrix are ordered in a similar manner
    so that for this case:

    Index:     0;       1;      2;      3;      4;
    rows:  DUt[0]; DUt[1]; DUt[2]; DUx[0]; DUx[1];
    cols:      u;     u,t;   u,tt;    u,x;    u,y;

    Index:      5;      6;      7;      8;      9;
    rows:  DUt[3]; DUt[4]; DUt[5]; DUx[2]; DUx[3];
    cols:       v;    v,t;   v,tt;    v,x;    v,y;

    However, the Jacobian matrix of the terms DUt/DUx w.r.t. Ut and Ux is
    often sparse. For this reason, the sparsity pattern is returned in a
    pair-wise format with in Jac_pairs which stores the (row, column) entries
    that are non-zero.