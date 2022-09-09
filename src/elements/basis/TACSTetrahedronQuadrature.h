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

#ifndef TACS_TETRAHEDRON_QUADRATURE_H
#define TACS_TETRAHEDRON_QUADRATURE_H

/*
  Quadrature points for tetrahedra on the [0, 1];[0, 1];[0, 1] intervals.
*/

const double TacsTetrahedronWts4[] = {1.0 / 24.0};
const double TacsTetrahedronPts4[] = {
    0.138196601125011, 0.138196601125011, 0.138196601125011, 0.585410196624969,
    0.138196601125011, 0.138196601125011, 0.138196601125011, 0.585410196624969,
    0.138196601125011, 0.138196601125011, 0.138196601125011, 0.585410196624969};

const double TacsTetrahedronWts5[] = {-2.0 / 15.0, 3.0 / 40.0};
const double TacsTetrahedronPts5[] = {
    0.25,      0.25,      0.25, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 0.5, 1.0 / 6.0,
    1.0 / 6.0, 1.0 / 6.0, 0.5,  1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 0.5};

/// Check this quadrature scheme for correctness
const double TacsTetrahedronWts11[] = {-0.013155555555556, 0.007622222222222,
                                       0.02488888888889};
const double TacsTetrahedronPts11[] = {0.25,
                                       0.25,
                                       0.25,
                                       0.071428571428571,
                                       0.071428571428571,
                                       0.071428571428571,
                                       0.785714285714286,
                                       0.071428571428571,
                                       0.071428571428571,
                                       0.071428571428571,
                                       0.785714285714286,
                                       0.071428571428571,
                                       0.071428571428571,
                                       0.071428571428571,
                                       0.785714285714286,
                                       0.399403576166799,
                                       0.100596423833201,
                                       0.100596423833201,
                                       0.100596423833201,
                                       0.399403576166799,
                                       0.100596423833201,
                                       0.100596423833201,
                                       0.100596423833201,
                                       0.399403576166799,
                                       0.100596423833201,
                                       0.399403576166799,
                                       0.399403576166799,
                                       0.399403576166799,
                                       0.100596423833201,
                                       0.399403576166799,
                                       0.399403576166799,
                                       0.399403576166799,
                                       0.100596423833201};

#endif  // TACS_TETRAHEDRON_QUADRATURE_H
