#ifndef TACS_SHELL_ELEMENT_DEFS_H
#define TACS_SHELL_ELEMENT_DEFS_H

#include "TACSBeamElement.h"
#include "TACSBeamElementBasis.h"
#include "TACSBeamElementModel.h"
#include "TACSDirector.h"
#include "TACSShellElement.h"
#include "TACSShellElementModel.h"
#include "TACSShellElementQuadBasis.h"
#include "TACSShellElementTriBasis.h"
#include "TACSShellInplaneElementModel.h"
#include "TACSThermalShellElement.h"

/*
  Linear shell elements with appropriate quadrature schemes
*/
typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                         TACSLinearizedRotation, TACSShellLinearModel>
    TACSQuad4Shell;

typedef TACSShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadBasis<3>,
                         TACSLinearizedRotation, TACSShellLinearModel>
    TACSQuad9Shell;

typedef TACSShellElement<TACSQuadCubicQuadrature, TACSShellQuadBasis<4>,
                         TACSLinearizedRotation, TACSShellLinearModel>
    TACSQuad16Shell;

typedef TACSShellElement<TACSTriLinearQuadrature, TACSShellTriLinearBasis,
                         TACSLinearizedRotation, TACSShellInplaneLinearModel>
    TACSTri3Shell;

/*
  Thermal shell elements with appropriate quadrature schemes
*/
typedef TACSThermalShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                                TACSLinearizedRotation, TACSShellLinearModel>
    TACSQuad4ThermalShell;

typedef TACSThermalShellElement<TACSQuadQuadraticQuadrature,
                                TACSShellQuadBasis<3>, TACSLinearizedRotation,
                                TACSShellLinearModel>
    TACSQuad9ThermalShell;

typedef TACSThermalShellElement<TACSQuadCubicQuadrature, TACSShellQuadBasis<4>,
                                TACSLinearizedRotation, TACSShellLinearModel>
    TACSQuad16ThermalShell;

typedef TACSThermalShellElement<TACSTriLinearQuadrature,
                                TACSShellTriLinearBasis, TACSLinearizedRotation,
                                TACSShellInplaneLinearModel>
    TACSTri3ThermalShell;

/*
  Shell elements with a linearized rotation and nonlinear strain expressions
*/
typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                         TACSLinearizedRotation, TACSShellNonlinearModel>
    TACSQuad4NonlinearShell;

typedef TACSShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadBasis<3>,
                         TACSLinearizedRotation, TACSShellNonlinearModel>
    TACSQuad9NonlinearShell;

typedef TACSShellElement<TACSQuadCubicQuadrature, TACSShellQuadBasis<4>,
                         TACSLinearizedRotation, TACSShellNonlinearModel>
    TACSQuad16NonlinearShell;

typedef TACSShellElement<TACSQuadCubicQuadrature, TACSShellTriLinearBasis,
                         TACSLinearizedRotation, TACSShellInplaneNonlinearModel>
    TACSTri3NonlinearShell;

/*
  Thermal shell elements with appropriate quadrature schemes
*/
typedef TACSThermalShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                                TACSLinearizedRotation, TACSShellNonlinearModel>
    TACSQuad4NonlinearThermalShell;

typedef TACSThermalShellElement<TACSQuadQuadraticQuadrature,
                                TACSShellQuadBasis<3>, TACSLinearizedRotation,
                                TACSShellNonlinearModel>
    TACSQuad9NonlinearThermalShell;

typedef TACSThermalShellElement<TACSQuadCubicQuadrature, TACSShellQuadBasis<4>,
                                TACSLinearizedRotation, TACSShellNonlinearModel>
    TACSQuad16NonlinearThermalShell;

typedef TACSThermalShellElement<TACSTriLinearQuadrature,
                                TACSShellTriLinearBasis, TACSLinearizedRotation,
                                TACSShellInplaneNonlinearModel>
    TACSTri3NonlinearThermalShell;

/*
  Moderate rotation shell elements with appropriate quadrature schemes
*/
typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                         TACSQuadraticRotation, TACSShellLinearModel>
    TACSQuad4ShellModRot;

typedef TACSShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadBasis<3>,
                         TACSQuadraticRotation, TACSShellLinearModel>
    TACSQuad9ShellModRot;

typedef TACSShellElement<TACSQuadCubicQuadrature, TACSShellQuadBasis<4>,
                         TACSQuadraticRotation, TACSShellLinearModel>
    TACSQuad16ShellModRot;

typedef TACSShellElement<TACSTriLinearQuadrature, TACSShellTriLinearBasis,
                         TACSQuadraticRotation, TACSShellInplaneLinearModel>
    TACSTri3ShellModRot;

/*
  Quaternion shell elements
*/
typedef TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadBasis<2>,
                         TACSQuaternionRotation, TACSShellLinearModel>
    TACSQuad4ShellQuaternion;

typedef TACSShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadBasis<3>,
                         TACSQuaternionRotation, TACSShellLinearModel>
    TACSQuad9ShellQuaternion;

typedef TACSShellElement<TACSQuadCubicQuadrature, TACSShellQuadBasis<4>,
                         TACSQuaternionRotation, TACSShellLinearModel>
    TACSQuad16ShellQuaternion;

typedef TACSShellElement<TACSTriLinearQuadrature, TACSShellTriLinearBasis,
                         TACSQuaternionRotation, TACSShellInplaneLinearModel>
    TACSTri3ShellQuaternion;

typedef TACSBeamElement<TACSBeamLinearQuadrature, TACSBeamBasis<2>,
                        TACSLinearizedRotation, TACSBeamLinearModel>
    TACSBeam2;

typedef TACSBeamElement<TACSBeamQuadraticQuadrature, TACSBeamBasis<3>,
                        TACSLinearizedRotation, TACSBeamLinearModel>
    TACSBeam3;

typedef TACSBeamElement<TACSBeamLinearQuadrature, TACSBeamBasis<2>,
                        TACSQuadraticRotation, TACSBeamLinearModel>
    TACSBeam2ModRot;

typedef TACSBeamElement<TACSBeamQuadraticQuadrature, TACSBeamBasis<3>,
                        TACSQuadraticRotation, TACSBeamLinearModel>
    TACSBeam3ModRot;

typedef TACSBeamElement<TACSBeamLinearQuadrature, TACSBeamBasis<2>,
                        TACSQuaternionRotation, TACSBeamLinearModel>
    TACSBeam2Quaternion;

typedef TACSBeamElement<TACSBeamQuadraticQuadrature, TACSBeamBasis<3>,
                        TACSQuaternionRotation, TACSBeamLinearModel>
    TACSBeam3Quaternion;

/**
  Create a TACS shell element based on the name of the shell.

  @param name The name of the shell element
  @param transform The transformation used for the shell
  @param con The shell constitutive object
*/
inline TACSElement *TacsCreateShellByName(const char *name,
                                          TACSShellTransform *transform,
                                          TACSShellConstitutive *con) {
  TACSElement *shell = NULL;
  if (strcmp(name, "TACSQuad4ShellModRot") == 0) {
    shell = new TACSQuad4ShellModRot(transform, con);
  } else if (strcmp(name, "TACSQuad9ShellModRot") == 0) {
    shell = new TACSQuad9ShellModRot(transform, con);
  } else if (strcmp(name, "TACSQuad16ShellModRot") == 0) {
    shell = new TACSQuad16ShellModRot(transform, con);
  } else if (strcmp(name, "TACSTri3ShellModRot") == 0) {
    shell = new TACSTri3ShellModRot(transform, con);
  } else if (strcmp(name, "TACSQuad4ShellQuaternion") == 0) {
    shell = new TACSQuad4ShellQuaternion(transform, con);
  } else if (strcmp(name, "TACSQuad9ShellQuaternion") == 0) {
    shell = new TACSQuad9ShellQuaternion(transform, con);
  } else if (strcmp(name, "TACSQuad16ShellQuaternion") == 0) {
    shell = new TACSQuad16ShellQuaternion(transform, con);
  } else if (strcmp(name, "TACSTri3ShellQuaternion") == 0) {
    shell = new TACSTri3ShellQuaternion(transform, con);
  } else if (strcmp(name, "TACSQuad4Shell") == 0 ||
             strcmp(name, "CQUAD") == 0 || strcmp(name, "CQUAD4") == 0 ||
             strcmp(name, "CQUADR") == 0) {
    shell = new TACSQuad4Shell(transform, con);
  } else if (strcmp(name, "TACSQuad9Shell") == 0 ||
             strcmp(name, "CQUAD9") == 0) {
    shell = new TACSQuad9Shell(transform, con);
  } else if (strcmp(name, "TACSQuad16Shell") == 0) {
    shell = new TACSQuad16Shell(transform, con);
  } else if (strcmp(name, "TACSTri3Shell") == 0) {
    shell = new TACSTri3Shell(transform, con);
  } else if (strcmp(name, "TACSQuad4ThermalShell") == 0) {
    shell = new TACSQuad4ThermalShell(transform, con);
  } else if (strcmp(name, "TACSQuad9ThermalShell") == 0) {
    shell = new TACSQuad9ThermalShell(transform, con);
  } else if (strcmp(name, "TACSQuad16ThermalShell") == 0) {
    shell = new TACSQuad16ThermalShell(transform, con);
  } else if (strcmp(name, "TACSTri3ThermalShell") == 0) {
    shell = new TACSTri3ThermalShell(transform, con);
  } else if (strcmp(name, "TACSQuad4NonlinearThermalShell") == 0) {
    shell = new TACSQuad4NonlinearThermalShell(transform, con);
  } else if (strcmp(name, "TACSQuad9NonlinearThermalShell") == 0) {
    shell = new TACSQuad9NonlinearThermalShell(transform, con);
  } else if (strcmp(name, "TACSQuad16NonlinearThermalShell") == 0) {
    shell = new TACSQuad16NonlinearThermalShell(transform, con);
  } else if (strcmp(name, "TACSTri3NonlinearThermalShell") == 0) {
    shell = new TACSTri3NonlinearThermalShell(transform, con);
  } else if (strcmp(name, "TACSQuad4NonlinearShell") == 0) {
    shell = new TACSQuad4NonlinearShell(transform, con);
  } else if (strcmp(name, "TACSQuad9NonlinearShell") == 0) {
    shell = new TACSQuad9NonlinearShell(transform, con);
  } else if (strcmp(name, "TACSQuad16NonlinearShell") == 0) {
    shell = new TACSQuad16NonlinearShell(transform, con);
  } else if (strcmp(name, "TACSTri3NonlinearShell") == 0) {
    shell = new TACSTri3NonlinearShell(transform, con);
  }

  return shell;
}

#endif  // TACS_SHELL_ELEMENT_DEFS_H
