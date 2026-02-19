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

#ifndef TACS_MARCHING_CUBES_H
#define TACS_MARCHING_CUBES_H

/*
  Marching cubes point
*/
class TACSMarchingCubesPoint {
 public:
  float x, y, z;
};

/*
  The data structures required for the marching cubes algorithm
*/
class TACSMarchingCubesCell {
 public:
  TACSMarchingCubesPoint p[8];
  float val[8];
};

/*
  Marching cube point
*/
class TACSMarchingCubesTriangle {
 public:
  TACSMarchingCubesPoint p[3];
};

/*
  Polygonalize
 */
int TacsPolygonizeCube(TACSMarchingCubesCell grid, float isolevel,
                       TACSMarchingCubesTriangle *triangles);

#endif  // TACS_MARCHING_CUBES_H
