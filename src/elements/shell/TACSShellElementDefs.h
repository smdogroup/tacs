#ifndef TACS_SHELL_ELEMENT_DEFS_H
#define TACS_SHELL_ELEMENT_DEFS_H

#include "TACSShellElementQuadBasis.h"
#include "TACSShellElementTriBasis.h"
#include "TACSDirector.h"
#include "TACSShellElementModel.h"
#include "TACSPlateElementModel.h"
#include "TACSShellElement.h"
#include "TACSThermalShellElement.h"

/*
  Linear shell elements with appropriate quadrature schemes
*/
typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                         TACSLinearizedRotation, TACSShellLinearModel> TACSQuad4Shell;

typedef TACSShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadBasis<3>,
                         TACSLinearizedRotation, TACSShellLinearModel> TACSQuad9Shell;

typedef TACSShellElement<TACSQuadCubicQuadrature, TACSShellQuadBasis<4>,
                         TACSLinearizedRotation, TACSShellLinearModel> TACSQuad16Shell;

typedef TACSShellElement<TACSTriLinearQuadrature, TACSShellTriLinearBasis,
                         TACSLinearizedRotation, TACSShellInplaneLinearModel> TACSTri3Shell;

// typedef TACSShellElement<TACSTriQuadraticQuadrature, TACSShellTriQuadraticBasis,
//                          TACSLinearizedRotation, TACSShellLinearModel> TACSTri6Shell;

/*
  Thermal shell elements with appropriate quadrature schemes
*/
typedef TACSThermalShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                                TACSLinearizedRotation, TACSShellLinearModel> TACSQuad4ThermalShell;

typedef TACSThermalShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadBasis<3>,
                                TACSLinearizedRotation, TACSShellLinearModel> TACSQuad9ThermalShell;

typedef TACSThermalShellElement<TACSQuadCubicQuadrature, TACSShellQuadBasis<4>,
                                TACSLinearizedRotation, TACSShellLinearModel> TACSQuad16ThermalShell;

/*
  Shell elements with a linearized rotation and nonlinear strain expressions
*/
typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                         TACSLinearizedRotation, TACSShellNonlinearModel> TACSQuad4NonlinearShell;

typedef TACSShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadBasis<3>,
                         TACSLinearizedRotation, TACSShellNonlinearModel> TACSQuad9NonlinearShell;

typedef TACSShellElement<TACSQuadCubicQuadrature, TACSShellQuadBasis<4>,
                         TACSLinearizedRotation, TACSShellNonlinearModel> TACSQuad16NonlinearShell;

#endif // TACS_SHELL_ELEMENT_DEFS_H
