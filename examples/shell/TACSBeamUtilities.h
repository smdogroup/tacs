#ifndef TACS_BEAM_UTILITIES_H
#define TACS_BEAM_UTILITIES_H

#include "TACSElementAlgebra.h"


template <class basis>
void getBeamNodeNormals( const TacsScalar Xpts[],
                         TacsScalar fn1[],
                         TacsScalar fn2[] ){
  for ( int i = 0; i < basis::NUM_NODES; i++ ){
    double pt[2];
    basis::getNodePoint(i, pt);

    // Compute the derivative X,xi at each node
    TacsScalar Xxi[3];
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);

    for
  }
}


#endif // TACS_BEAM_UTILITIES_H
