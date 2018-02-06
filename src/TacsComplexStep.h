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

#ifndef TACS_COMPLEX_STEP_H
#define TACS_COMPLEX_STEP_H

#include <complex>

// Define the real part function for the complex data type
inline double TacsRealPart( const std::complex<double>& c ){
  return real(c);
}

// Define the imaginary part function for the complex data type
inline double TacsImagPart( const std::complex<double>& c ){
  return imag(c);
}

// Dummy function for real part
inline double TacsRealPart( const double& r ){
  return r;
}

// Compute the absolute value
inline std::complex<double> fabs( const std::complex<double>& c ){
  if (real(c) < 0.0){
    return -c;
  }
  return c;
}

#endif // TACS_COMPLEX_STEP_H
