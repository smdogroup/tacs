/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_ELEMENT_TEMPLATES_H
#define TACS_ELEMENT_TEMPLATES_H

#include "PlaneStressQuad.h"
#include "PlaneStressTraction.h"
#include "TACSShellTraction.h"
#include "TACS3DTraction.h"
#include "MITCShell.h"
#include "Solid.h"
#include "PoissonElement.h"
#include "PlaneStressCoupledThermoQuad.h"
#include "PlaneStressCoupledThermoTraction.h"
#include "CoupledThermoSolid.h"
#include "TACS3DCoupledThermoTraction.h"

// Add define statements for the element types
#define PlaneStressQuad2 PlaneStressQuad<2>
#define PlaneStressQuad3 PlaneStressQuad<3>
#define PlaneStressQuad4 PlaneStressQuad<4>
#define PlaneStressQuad5 PlaneStressQuad<5>
#define PSQuadTraction2 PSQuadTraction<2>
#define PSQuadTraction3 PSQuadTraction<3>
#define PSQuadTraction4 PSQuadTraction<4>
#define PSQuadTraction5 PSQuadTraction<5>
#define MITCShell2 MITCShell<2>
#define MITCShell3 MITCShell<3>
#define MITCShell4 MITCShell<4>
#define MITCShell5 MITCShell<5>
#define MITCShell32 MITCShell<3,2>
#define MITCShell43 MITCShell<4,3>
#define MITCShell54 MITCShell<5,4>
#define Solid2 Solid<2>
#define Solid3 Solid<3>
#define Solid4 Solid<4>
#define Solid5 Solid<5>
#define TACS3DTraction2 TACS3DTraction<2>
#define TACS3DTraction3 TACS3DTraction<3>
#define TACS3DTraction4 TACS3DTraction<4>
#define TACS3DTraction5 TACS3DTraction<5>
#define TACS3DPressureTraction2 TACS3DPressureTraction<2>
#define TACS3DPressureTraction3 TACS3DPressureTraction<3>
#define TACS3DPressureTraction4 TACS3DPressureTraction<4>
#define TACS3DPressureTraction5 TACS3DPressureTraction<5>
#define TACS3DBoundingTraction2 TACS3DBoundingTraction<2>
#define TACSShellTraction2 TACSShellTraction<2>
#define TACSShellTraction3 TACSShellTraction<3>
#define TACSShellTraction4 TACSShellTraction<4>
#define TACSShellTraction5 TACSShellTraction<5>
#define PoissonQuad2 PoissonQuad<2>
#define PoissonQuad3 PoissonQuad<3>
#define PoissonQuad4 PoissonQuad<4>
#define PoissonQuad5 PoissonQuad<5>
#define PSThermoQuad2 PlaneStressCoupledThermoQuad<2>
#define PSThermoQuad3 PlaneStressCoupledThermoQuad<3>
#define PSThermoQuad4 PlaneStressCoupledThermoQuad<4>
#define PSThermoQuad5 PlaneStressCoupledThermoQuad<5>
#define PSThermoQuad6 PlaneStressCoupledThermoQuad<6>
#define PSThermoQuadTraction2 PSQuadThermoTraction<2>
#define PSThermoQuadTraction3 PSQuadThermoTraction<3>
#define PSThermoQuadTraction4 PSQuadThermoTraction<4>
#define PSThermoQuadTraction5 PSQuadThermoTraction<5>
#define PSThermoQuadTraction6 PSQuadThermoTraction<6>
#define PSThermoQuadHF2 PSQuadHeatFluxTraction<2>
#define PSThermoQuadHF3 PSQuadHeatFluxTraction<3>
#define PSThermoQuadHF4 PSQuadHeatFluxTraction<4>
#define PSThermoQuadHF5 PSQuadHeatFluxTraction<5>
#define SolidThermo2 CoupledThermoSolid<2>
#define SolidThermo3 CoupledThermoSolid<3>
#define SolidThermo4 CoupledThermoSolid<4>
#define SolidThermo5 CoupledThermoSolid<5>
#define SolidThermo6 CoupledThermoSolid<6>
#define TACS3DThermoTraction2 TACS3DThermoTraction<2>
#define TACS3DThermoTraction3 TACS3DThermoTraction<3>
#define TACS3DThermoTraction4 TACS3DThermoTraction<4>
#define TACS3DThermoTraction5 TACS3DThermoTraction<5>
#define TACS3DThermoTraction6 TACS3DThermoTraction<6>
#define TACS3DThermoPressureTraction2 TACS3DThermoPressureTraction<2>
#define TACS3DThermoPressureTraction3 TACS3DThermoPressureTraction<3>
#define TACS3DThermoPressureTraction4 TACS3DThermoPressureTraction<4>
#define TACS3DThermoPressureTraction5 TACS3DThermoPressureTraction<5>
#define TACS3DThermoPressureTraction6 TACS3DThermoPressureTraction<6>
#define TACS3DThermoHF2 TACS3DHeatFluxTraction<2>
#define TACS3DThermoHF3 TACS3DHeatFluxTraction<3>
#define TACS3DThermoHF4 TACS3DHeatFluxTraction<4>
#define TACS3DThermoHF5 TACS3DHeatFluxTraction<5>
#define TACS3DThermoNormalHF2 TACS3DNormalHeatFluxTraction<2>
#define TACS3DThermoNormalHF3 TACS3DNormalHeatFluxTraction<3>
#define TACS3DThermoNormalHF4 TACS3DNormalHeatFluxTraction<4>
#define TACS3DThermoNormalHF5 TACS3DNormalHeatFluxTraction<5>

#endif // TACS_ELEMENT_TEMPLATES_H
