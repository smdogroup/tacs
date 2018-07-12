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
#define Solid2 Solid<2>
#define Solid3 Solid<3>
#define Solid4 Solid<4>
#define Solid5 Solid<5>
#define TACS3DTraction2 TACS3DTraction<2>
#define TACS3DTraction3 TACS3DTraction<3>
#define TACS3DTraction4 TACS3DTraction<4>
#define TACS3DTraction5 TACS3DTraction<5>
#define TACSShellTraction2 TACSShellTraction<2>
#define TACSShellTraction3 TACSShellTraction<3>
#define TACSShellTraction4 TACSShellTraction<4>
#define TACSShellTraction5 TACSShellTraction<5>
#define PoissonQuad2 PoissonQuad<2>
#define PoissonQuad3 PoissonQuad<3>
#define PoissonQuad4 PoissonQuad<4>
#define PoissonQuad5 PoissonQuad<5>

#endif // TACS_ELEMENT_TEMPLATES_H
