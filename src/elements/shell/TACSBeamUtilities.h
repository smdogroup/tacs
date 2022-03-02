#ifndef TACS_BEAM_UTILITIES_H
#define TACS_BEAM_UTILITIES_H

#include "TACSElementAlgebra.h"

/*
  Compute the frame normals at each node location

  @param Xpts The node locations for the elements
  @param axis The coordinates of the reference axis
  @param fn1 The first normal direction
  @param fn2 The second normal direction
*/
template <class basis>
void TacsBeamComputeNodeNormals( const TacsScalar Xpts[],
                                 const TacsScalar axis[],
                                 TacsScalar fn1[],
                                 TacsScalar fn2[] ){
  for ( int i = 0; i < basis::NUM_NODES; i++ ){
    double pt[2];
    basis::getNodePoint(i, pt);

    // Compute the derivative X,xi at each node
    TacsScalar Xxi[3];
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);

    // Normalize the tangent direction
    TacsScalar t[3];
    TacsScalar tnorm = sqrt(vec3Dot(Xxi, Xxi));
    TacsScalar tinv = 1.0/tnorm;
    t[0] = tinv*Xxi[0];
    t[1] = tinv*Xxi[1];
    t[2] = tinv*Xxi[2];

    // Compute the first direction in the plane
    TacsScalar n1[3];
    TacsScalar tdot = vec3Dot(t, axis);
    n1[0] = axis[0] - tdot*t[0];
    n1[1] = axis[1] - tdot*t[1];
    n1[2] = axis[2] - tdot*t[2];

    // Compute the norm
    TacsScalar n1norm = sqrt(vec3Dot(n1, n1));
    TacsScalar n1inv = 1.0/n1norm;
    fn1[0] = n1inv*n1[0];
    fn1[1] = n1inv*n1[1];
    fn1[2] = n1inv*n1[2];

    // Compute the cross product
    TacsScalar n2[3];
    crossProduct(1.0, t, fn1, fn2);

    fn1 += 3;
    fn2 += 3;
  }
}




#endif // TACS_BEAM_UTILITIES_H
