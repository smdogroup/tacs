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

#ifndef TACS_AMD_INTERFACE_H
#define TACS_AMD_INTERFACE_H

/*
  An approximate minimum degree reordering scheme for the schur
  complement method. Complete documentation is provided in the .c
  file. Note that the contents of rowp/cols are overwritten during
  this call so be careful!
*/

void amd_order_interface(int nvars, int *rowp, int *cols, int *perm,
                         int *interface_nodes, int ninterface_nodes,
                         int ndep_vars, const int *dep_vars,
                         const int *indep_ptr, const int *indep_vars,
                         int use_exact_degree);

#endif
