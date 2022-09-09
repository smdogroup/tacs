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

#ifndef TACS_TRIANGLE_QUADRATURE_H
#define TACS_TRIANGLE_QUADRATURE_H

/*
  Quadrature points for triangles on the [0, 1];[0, 1] intervals.
*/

const double TacsTriangleWts3[] = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};
const double TacsTrianglePts3[] = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0,
                                   1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};

const double TacsTriangleWts4[] = {-27.0 / 96.0, 25.0 / 96.0, 25.0 / 96.0,
                                   25.0 / 96.0};
const double TacsTrianglePts4[] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 5.0, 1.0 / 5.0,
                                   3.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 3.0 / 5.0};

const double TacsTriangleWts6[] = {0.054975871827661, 0.054975871827661,
                                   0.054975871827661, 0.111690794839006,
                                   0.111690794839006, 0.111690794839006};
const double TacsTrianglePts6[] = {
    0.091576213509771, 0.091576213509771, 0.816847572980459, 0.091576213509771,
    0.091576213509771, 0.816847572980459, 0.108103018168070, 0.108103018168070,
    0.445948490915965, 0.108103018168070, 0.108103018168070, 0.445948490915965};

#endif  // TACS_TRIANGLE_QUADRATURE_H
