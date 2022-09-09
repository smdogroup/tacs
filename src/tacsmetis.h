/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_INCLUDE_METIS_H
#define TACS_INCLUDE_METIS_H

#ifdef TACS_CPLUSPLUS_METIS
// Use this if METIS is compiled with C++
#include "metis.h"
#else
// Otherwise, assume metis is compiled with C
extern "C" {
#include "metis.h"
}
#endif  // TACS_CPLUSPLUS_METIS

#endif
