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

#include "TACSMatrixFreeMat.h"

TACSMatrixFreeMat::TACSMatrixFreeMat(TACSAssembler *_assembler) {
  assembler = _assembler;
  assembler->incref();
  data_size = 0;
  temp_size = 0;
  data = NULL;
  temp = NULL;
  matType = TACS_JACOBIAN_MATRIX;
}

TACSMatrixFreeMat::~TACSMatrixFreeMat() {
  assembler->decref();
  if (data) {
    delete[] data;
  }
  if (temp) {
    delete[] temp;
  }
}

void TACSMatrixFreeMat::assembleMatrixFreeData(ElementMatrixType _matType,
                                               double alpha, double beta,
                                               double gamma) {
  if (data) {
    delete[] data;
  }
  if (temp) {
    delete[] temp;
  }
  matType = _matType;
  assembler->getMatrixFreeDataSize(matType, &data_size, &temp_size);
  data = new TacsScalar[data_size];
  temp = new TacsScalar[temp_size];
  assembler->assembleMatrixFreeData(matType, alpha, beta, gamma, data);
}

TACSVec *TACSMatrixFreeMat::createVec() { return assembler->createVec(); }

void TACSMatrixFreeMat::mult(TACSVec *tx, TACSVec *ty) {
  TACSBVec *x = dynamic_cast<TACSBVec *>(tx);
  TACSBVec *y = dynamic_cast<TACSBVec *>(ty);
  if (x && y) {
    y->zeroEntries();
    assembler->addMatrixFreeVecProduct(matType, data, temp, x, y,
                                       TACS_MAT_NORMAL);
  }
}

void TACSMatrixFreeMat::multTranspose(TACSVec *tx, TACSVec *ty) {
  TACSBVec *x = dynamic_cast<TACSBVec *>(tx);
  TACSBVec *y = dynamic_cast<TACSBVec *>(ty);
  if (x && y) {
    y->zeroEntries();
    assembler->addMatrixFreeVecProduct(matType, data, temp, x, y,
                                       TACS_MAT_TRANSPOSE);
  }
}

const char *TACSMatrixFreeMat::getObjectName() { return "TACSMatrixFreeMat"; }