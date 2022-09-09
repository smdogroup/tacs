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

#ifndef TACS_MATRIX_FREE_MAT_H
#define TACS_MATRIX_FREE_MAT_H

#include "TACSAssembler.h"

/*
  Wrapper for a matrix-free matrix class that utilizes
  matrix-vector products from the TACSAssembler interface.
*/
class TACSMatrixFreeMat : public TACSMat {
 public:
  TACSMatrixFreeMat(TACSAssembler *_assembler);
  ~TACSMatrixFreeMat();

  void assembleMatrixFreeData(ElementMatrixType _matType, double alpha,
                              double beta, double gamma);
  TACSVec *createVec();
  void mult(TACSVec *x, TACSVec *y);
  void multTranspose(TACSVec *x, TACSVec *y);
  const char *getObjectName();

 private:
  TACSAssembler *assembler;
  ElementMatrixType matType;
  TacsScalar *data, *temp;
  int data_size, temp_size;
};

#endif  // TACS_MATRIX_FREE_MAT_H
