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

#ifndef TACS_ELEMENT_VERIFICATION_H
#define TACS_ELEMENT_VERIFICATION_H

#include "TACSElement.h"

/**
  Assign variables randomly to an array. This is useful for
  testing various things.

 */
void TacsGenerateRandomArray( TacsScalar *array, int size,
                              TacsScalar lower=-1.0,
                              TacsScalar upper=1.0 );

/*
  Find the largest absolute value of the difference between the
  arrays a and b
*/
double TacsGetMaxError( TacsScalar *a, TacsScalar *b, int size,
                        int *max_index );

/*
  Find the maximum relative error between a and b and return the
*/
double TacsGetMaxRelError( TacsScalar *a, TacsScalar *b, int size,
                           int *max_index );

/*
  Print out the values and the relative errors
*/
void TacsPrintErrorComponents( FILE *fp, const char *descript,
                               TacsScalar *a, TacsScalar *b,
                               int size );

/*
  Perturb the input variables in the forward sense
*/
void TacsForwardDiffPerturb( TacsScalar *out, int size,
                             const TacsScalar *orig,
                             const TacsScalar *pert,
                             double dh );

/*
  Perturb the variables in the backward sense
*/
void TacsBackwardDiffPerturb( TacsScalar *out, int size,
                              const TacsScalar *orig,
                              const TacsScalar *pert,
                              double dh );

/*
  
*/
int TacsTestElementResidual( TACSElement *element,
                             double time,
                             const TacsScalar Xpts[],
                             const TacsScalar vars[],
                             const TacsScalar dvars[],
                             const TacsScalar ddvars[],
                             double dh,
                             int test_print_leve,
                             double test_fail_atol
                             double test_fail_rtol );

#endif TACS_ELEMENT_VERIFICATION_H