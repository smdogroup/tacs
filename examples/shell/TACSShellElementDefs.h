#ifndef TACS_SHELL_ELEMENT_DEFS_H
#define TACS_SHELL_ELEMENT_DEFS_H

#include "TACSShellElementBasis.h"
#include "TACSShellElementBasis.h"
#include "TACSDirector.h"
#include "TACSShellElement.h"
#include "TACSThermalShellElement.h"

/*
  Linear shell elements with appropriate quadrature schemes
*/
typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                         TACSLinearizedRotation, TACSShellLinearModel> TACSQuad2Shell;

typedef TACSShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadBasis<3>,
                         TACSLinearizedRotation, TACSShellLinearModel> TACSQuad3Shell;

typedef TACSShellElement<TACSQuadCubicQuadrature, TACSShellQuadBasis<4>,
                         TACSLinearizedRotation, TACSShellLinearModel> TACSQuad4Shell;

/*
  Thermal shell elements with appropriate quadrature schemes
*/
typedef TACSThermalShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                                TACSLinearizedRotation, TACSShellLinearModel> TACSQuad2ThermalShell;

typedef TACSThermalShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadBasis<3>,
                                TACSLinearizedRotation, TACSShellLinearModel> TACSQuad3ThermalShell;

typedef TACSThermalShellElement<TACSQuadCubicQuadrature, TACSShellQuadBasis<4>,
                                TACSLinearizedRotation, TACSShellLinearModel> TACSQuad4ThermalShell;

// typedef TACSShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadBasis<3>,
//                          TACSQuadraticRotation, TACSShellLinearModel> TACSQuadLinearShell;

// typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
//                          TACSQuaternionRotation, TACSShellLinearModel> TACSQuadLinearShell;

// typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
//                          TACSLinearizedRotation, TACSShellNonlinearModel> TACSQuadLinearShell;

// typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
//                          TACSQuadraticRotation, TACSShellNonlinearModel> TACSQuadLinearShell;

// typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
//                          TACSQuaternionRotation, TACSShellNoninearModel> TACSQuadLinearShell;


// typedef TACSShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadBasis<3>,
//                          TACSLinearizedRotation, TACSShellLinearModel> TACSQuadQuadraticShell;

// typedef TACSShellElement<TACSTriQuadraticQuadrature, TACSShellTriQuadraticBasis,
//     TACSLinearizedRotation, TACSShellLinearModel> TACSTriQuadraticShell;


#endif // TACS_SHELL_ELEMENT_DEFS_H
