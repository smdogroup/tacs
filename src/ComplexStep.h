#ifndef COMPLEX_STEP_H
#define COMPLEX_STEP_H

#include <complex>

/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved
*/

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

#endif // COMPLEX_STEP_H
