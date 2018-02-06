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

#ifndef TACS_FE_LIBRARY_H
#define TACS_FE_LIBRARY_H

/*!  
  FElibrary.h contains many important functions and data that are
  repeatedly used in the formation of various finite element stiffness
  matricies. The intent of this code is to provide functions for
  calculating shape functions, Gauss points and functions for
  calculating the Jacobians for element stiffness matricies.
*/

#include <stdlib.h>
#include <math.h>
#include "TACSObject.h"

/*
  The following data defines the quadrature rules that can be used in
  the elements. The default is to use tensor-product Gauss quadrature
  rules, however, more sophisticated methods can be used.

  Currently, the options are Gauss quadrature or Lobatto (or
  Gauss-Lobatto) quadrature schemes that include the end points of the
  interval.
*/
enum QuadratureType { GAUSS_QUADRATURE, 
                      LOBATTO_QUADRATURE };

TACS_BEGIN_NAMESPACE(FElibrary)

int comparator( const void * a, const void * b );

// Design variable related operations

/*! 
  Sort an array and return the number of unique design variables within 
  that array start is the index of the first non-negative value and the 
  return value is the number of non-negative design variables
  &dvNums[start] is the start of the new array of unique design vars...
*/
int uniqueSort( int *dvNums, int numDVs );

/*!
  Merge two sorted arrays into a single sorted array, in place.

  Two part algorithm:
  1. Find the number of duplicates between a and b
  2. Run through the list backwards placing elements into a[]
  when appropriate.

  Memory requirements: note that len(a) >= na + nb
*/
int mergeArrays( int *a, int na, const int *b, int nb );

/*!
  Find the interval such that the given index satisfies:

  intv[k] <= index < intv[k+1]

  The intervals must be non-decreasing. Note that len is equal
  to the length of the intv array which is one more than the 
  total number of intervals.
*/
int findInterval( int index, const int intv[], int len );

/*!
  Match the intervals in a list of sorted variables.

  Given the intervals, the range of the intervals and 
  the list of sorted variables, local a point such that

  vars[ext_ptr[n]] <= ownerRange[n]
*/
void matchIntervals( int mpiSize, const int ownerRange[], 
                     int nvars, const int vars[], int ext_ptr[] );

/*!  
  Solve the quadratic equation and return the positive and negative
  roots.

  The code returns the roots of the equation:

  a*x^2 + b*x + c = 0
  
  This code avoids truncation error by testing the sign of the
  b coefficient and using the corresponding expression with the
  least susceptibility to truncation error.
*/
template <class ScalarType>
void solveQERoots( ScalarType * r1, ScalarType * r2,
                   ScalarType a, ScalarType b, ScalarType c ){
  ScalarType discrim = b*b - 4.0*a*c;
  if (TacsRealPart(discrim) < 0.0){
    *r1 = *r2 = 0.0;
    return;
  }

  if (TacsRealPart(a) == 0.0){
    if (TacsRealPart(b) == 0.0){
      *r1 = *r2 = 0.0;
      return;
    }

    // Solve b*x + c = 0
    *r1 = - c/b;
    *r2 = 0.0;
    return;
  }

  // Depending on the sign of b, use different expressions to 
  // avoid truncation error
  discrim = sqrt(discrim);

  if (TacsRealPart(b) > 0.0){
    *r1 = -(b + discrim)/(2.0*a);
    *r2 = c/((*r1)*a);
  }
  else { // b < 0.0
    *r1 = -(b - discrim)/(2.0*a);
    *r2 = c/((*r1)*a);
  }
}

/*
  Compute one-dimensional Lagrange shape functions on the interval
  (-1, 1) with evenly-spaced interpolation points. Note that 
  there is optimized code for order <= 4.

  input:
  a: the parametric location
  knots: the parametric locations of the interpolanting poitns

  output:
  sf: the shape functions
*/
inline void lagrangeSF( double sf[], 
                        const double a, int porder ){
  if (porder <= 1){
    sf[0] = 1.0;
  }
  else if (porder == 2){
    sf[0] = 0.5*(1.0 - a);
    sf[1] = 0.5*(a + 1.0);
  }
  else if (porder == 3){
    sf[0] = -0.5*a*(1.0 - a);
    sf[1] = (1.0 - a)*(1.0 + a);
    sf[2] = 0.5*(1.0 + a)*a;
  }
  else if (porder == 4){
    sf[0] = -(1.0/16.0)*(3.0*a + 1.0)*(3.0*a - 1.0)*(a - 1.0);
    sf[1] =  (9.0/16.0)*(a + 1.0)*(3.0*a - 1.0)*(a - 1.0);
    sf[2] = -(9.0/16.0)*(a + 1.0)*(3.0*a + 1.0)*(a - 1.0);
    sf[3] =  (1.0/16.0)*(a + 1.0)*(3.0*a + 1.0)*(3.0*a - 1.0);
  }
  else {
    // Loop over the shape function control points
    for ( int i = 0; i < porder; i++ ){
      double ki = -1.0 + 2.0*i/(porder - 1.0);
      sf[i] = 1.0;
        
      // Loop over each point again, except for the current control point, 
      // adding the contribution to the shape function
      for ( int j = 0; j < porder; j++ ){
        if (i != j){
          double kj = -1.0 + 2.0*j/(porder - 1.0);
            
          sf[i] = sf[i]*(a - kj)/(ki - kj);
        }
      }      
    }
  }
}

/*
  Compute one-dimensional Lagrange shape functions on the interval
  (-1, 1) with evenly-spaced interpolation points. Note that 
  there is optimized code for order <= 4.

  input:
  a: the parametric location
  knots: the parametric locations of the interpolanting poitns

  output:
  sf: the shape functions
  dsf: the derivative of the shape functions
*/
inline void lagrangeSF( double sf[], double dsf[], 
                        const double a, int porder ){
  if (porder <= 1){
    sf[0] = 1.0;
    dsf[0] = 0.0;
  }
  else if (porder == 2){
    sf[0] = 0.5*(1.0 - a);
    sf[1] = 0.5*(a + 1.0);

    dsf[0] = -0.5;
    dsf[1] = 0.5;
  }
  else if (porder == 3){
    sf[0] = -0.5*a*(1.0 - a);
    sf[1] = (1.0 - a)*(1.0 + a);
    sf[2] = 0.5*(1.0 + a)*a;

    dsf[0] = -0.5 + a;
    dsf[1] = - 2.0*a;
    dsf[2] = 0.5 + a;
  }
  else if (porder == 4){
    sf[0] = -(1.0/16.0)*(3.0*a + 1.0)*(3.0*a - 1.0)*(a - 1.0);
    sf[1] =  (9.0/16.0)*(a + 1.0)*(3.0*a - 1.0)*(a - 1.0);
    sf[2] = -(9.0/16.0)*(a + 1.0)*(3.0*a + 1.0)*(a - 1.0);
    sf[3] =  (1.0/16.0)*(a + 1.0)*(3.0*a + 1.0)*(3.0*a - 1.0);
    
    dsf[0] = - (1.0/16.0)*(27.0*a*a - 18.0*a - 1.0);
    dsf[1] =   (9.0/16.0)*(9.0*a*a - 2.0*a - 3.0);
    dsf[2] = - (9.0/16.0)*(9.0*a*a + 2.0*a - 3.0);
    dsf[3] =   (1.0/16.0)*(27.0*a*a + 18.0*a - 1.0);
  }
  else {
    double dn;
    
    for ( int i = 0; i < porder; i++ ){
      sf[i] = 1.0;
      dsf[i] = 0.0;
    }
    
    // Loop over the shape function control points
    for ( int i = 0; i < porder; i++ ){
      double ki = -1.0 + 2.0*i/(porder - 1.0);
      
      // Loop over each point again, except for the current control point, 
      // adding the contribution to the shape function
      for ( int j = 0; j < porder; j++ ){
        if ( i != j ){
          double kj = -1.0 + 2.0*j/(porder - 1.0);
          sf[i] = sf[i]*(a - kj)/(ki - kj);
          
          // Loop over the whole thing again to determine the contribution to
          // the derivative of the shape function
          dn = 1.0/(ki - kj);
          
          for ( int k = 0; k < porder; k++ ){
            if (k != i && k != j){
              double kk = -1.0 + 2.0*k/(porder - 1.0);
              dn = dn*(a - kk)/(ki - kk);
            }
          }
          
          dsf[i] += dn;
        }
      }
    }
  }
}

/*
  Compute one-dimensional Lagrange shape functions on the interval
  (-1, 1) with evenly-spaced interpolation points. Note that 
  there is optimized code for order <= 4.

  input:
  a: the parametric location
  knots: the parametric locations of the interpolanting poitns

  output:
  sf: the shape functions
  dsf: the derivative of the shape functions
  ddsf: the second derivative of the shape functions
*/
inline void lagrangeSF( double sf[], double dsf[], double ddsf[], 
                        const double a, int porder ){
  if (porder <= 1){
    sf[0] = 1.0;
    dsf[0] = 0.0;
    ddsf[0] = 0.0;
  }
  else if (porder == 2){
    sf[0] = 0.5*(1.0 - a);
    sf[1] = 0.5*(a + 1.0);

    dsf[0] = -0.5;
    dsf[1] = 0.5;

    ddsf[0] = 0.0;
    ddsf[1] = 0.0;
  }
  else if (porder == 3){
    sf[0] = -0.5*a*(1.0 - a);
    sf[1] = (1.0 - a)*(1.0 + a);
    sf[2] = 0.5*(1.0 + a)*a;

    dsf[0] = -0.5 + a;
    dsf[1] = - 2.0*a;
    dsf[2] = 0.5 + a;

    ddsf[0] = 1.0;
    ddsf[1] = -2.0;
    ddsf[2] = 1.0;
  }
  else if (porder == 4){
    sf[0] = -(1.0/16.0)*(3.0*a + 1.0)*(3.0*a - 1.0)*(a - 1.0);
    sf[1] =  (9.0/16.0)*(a + 1.0)*(3.0*a - 1.0)*(a - 1.0);
    sf[2] = -(9.0/16.0)*(a + 1.0)*(3.0*a + 1.0)*(a - 1.0);
    sf[3] =  (1.0/16.0)*(a + 1.0)*(3.0*a + 1.0)*(3.0*a - 1.0);
    
    dsf[0] = - (1.0/16.0)*(27.0*a*a - 18.0*a - 1.0);
    dsf[1] =   (9.0/16.0)*(9.0*a*a - 2.0*a - 3.0);
    dsf[2] = - (9.0/16.0)*(9.0*a*a + 2.0*a - 3.0);
    dsf[3] =   (1.0/16.0)*(27.0*a*a + 18.0*a - 1.0);

    ddsf[0] = - (1.0/16.0)*(54.0*a - 18.0);
    ddsf[1] =   (9.0/16.0)*(18.0*a - 2.0);
    ddsf[2] = - (9.0/16.0)*(18.0*a + 2.0);
    ddsf[3] =   (1.0/16.0)*(54.0*a + 18.0);
  }
  else {
    double dn, ddn;
    
    for ( int i = 0; i < porder; i++ ){
      sf[i] = 1.0;
      dsf[i] = 0.0;
      ddsf[i] = 0.0;
    }
    
    // Loop over the shape function control points
    for ( int i = 0; i < porder; i++ ){
      double ki = -1.0 + 2.0*i/(porder - 1.0);
      
      // Loop over each point again, except for the current control point, 
      // adding the contribution to the shape function
      for ( int j = 0; j < porder; j++ ){
        if ( i != j ){
          double kj = -1.0 + 2.0*j/(porder - 1.0);
          sf[i] = sf[i]*(a - kj)/(ki - kj);
          
          // Loop over the whole thing again to determine the contribution to
          // the derivative of the shape function
          dn = 1.0/(ki - kj);
          
          for ( int k = 0; k < porder; k++ ){
            if ( k != i && k != j ){
              double kk = -1.0 + 2.0*k/(porder - 1.0);
              
              dn = dn*(a - kk)/(ki - kk);             
              ddn = 1.0/((ki - kk)*(ki - kj));
              
              for ( int m = 0; m < porder; m++ ){
                if (((m != i) && (m != j)) && (m != k)){
                  double km = -1.0 + 2.0*m/(porder - 1.0);
                  
                  ddn = ddn*(a - km)/(ki - km);
                }
              }
              ddsf[i] += ddn;
            }
          }     
          dsf[i] += dn;
        }
      }
    }
  }
}

/*
  Compute one-dimensional Lagrange shape functions on an arbitrary
  interval with interpolant points as input

  input:
  a: the parametric location
  knots: the parametric locations of the interpolanting poitns

  output:
  sf: the shape functions
*/
inline void lagrangeSFKnots( double sf[], 
                             const double a, const double knots[], int porder ){
  if (porder <= 1){
    sf[0] = 1.0;
  }
  else if (porder == 2){
    const double k0 = knots[0], k1 = knots[1];
    sf[0] = (a - k1)/(k0 - k1);
    sf[1] = (a - k0)/(k1 - k0);
  }
  else if (porder == 3){
    const double k0 = knots[0], k1 = knots[1], k2 = knots[2];
    sf[0] = (a - k1)*(a - k2)/((k0 - k1)*(k0 - k2));
    sf[1] = (a - k0)*(a - k2)/((k1 - k0)*(k1 - k2));
    sf[2] = (a - k0)*(a - k1)/((k2 - k0)*(k2 - k1));
  }
  else {
    // Loop over the shape function control points
    for ( int i = 0; i < porder; i++ ){
      double ki = knots[i];
      sf[i] = 1.0;
        
      // Loop over each point again, except for the current control point, 
      // adding the contribution to the shape function
      for ( int j = 0; j < porder; j++ ){
        if (i != j){
          sf[i] = sf[i]*(a - knots[j])/(ki - knots[j]);
        }
      }      
    }
  }
}

/*
  Compute one-dimensional Lagrange shape functions on an arbitrary
  interval with interpolant points as input

  input:
  a: the parametric location
  knots: the parametric locations of the interpolanting poitns

  output:
  sf: the shape functions
  dsf: the derivatives of the shape functions
*/
inline void lagrangeSFKnots( double sf[], double dsf[], 
                             const double a, const double knots[], int porder ){
  if (porder <= 1){
    sf[0] = 1.0;
    dsf[0] = 0.0;
  }
  else if (porder == 2){
    const double k0 = knots[0], k1 = knots[1];
    sf[0] = (a - k1)/(k0 - k1);
    sf[1] = (a - k0)/(k1 - k0);

    dsf[0] = 1.0/(k0 - k1);
    dsf[1] = 1.0/(k1 - k0);
  }
  else if (porder == 3){
    const double k0 = knots[0], k1 = knots[1], k2 = knots[2];
    sf[0] = (a - k1)*(a - k2)/((k0 - k1)*(k0 - k2));
    sf[1] = (a - k0)*(a - k2)/((k1 - k0)*(k1 - k2));
    sf[2] = (a - k0)*(a - k1)/((k2 - k0)*(k2 - k1));

    dsf[0] = ((a - k1) + (a - k2))/((k0 - k1)*(k0 - k2));
    dsf[1] = ((a - k0) + (a - k2))/((k1 - k0)*(k1 - k2));
    dsf[2] = ((a - k0) + (a - k1))/((k2 - k0)*(k2 - k1));
  }
  else {
    double dn;
    
    for ( int i = 0; i < porder; i++ ){
      sf[i] = 1.0;
      dsf[i] = 0.0;
    }
    
    // Loop over the shape function control points
    for ( int i = 0; i < porder; i++ ){
      double ki = knots[i];
      
      // Loop over each point again, except for the current control
      // point, adding the contribution to the shape function
      for ( int j = 0; j < porder; j++ ){
        if ( i != j ){
          double kj = knots[j];
          sf[i] = sf[i]*(a - kj)/(ki - kj);
          
          // Loop over the whole thing again to determine the
          // contribution to the derivative of the shape function
          dn = 1.0/(ki - kj);
          
          for ( int k = 0; k < porder; k++ ){
            if (k != i && k != j){
              dn = dn*(a - knots[k])/(ki - knots[k]);
            }
          }
          
          dsf[i] += dn;
        }
      }
    }
  }
}
  
/*
  Compute one-dimensional Lagrange shape functions on an arbitrary
  interval with interpolant points as input

  input:
  a: the parametric location
  knots: the parametric locations of the interpolanting poitns

  output:
  sf: the shape functions
  dsf: the derivatives of the shape functions
  ddsf: the second derivatives of the shape functions
*/
inline void lagrangeSFKnots( double sf[], double dsf[], double ddsf[], 
                             const double a, const double knots[], 
                             int porder ){
  if (porder <= 1){
    sf[0] = 1.0;
    dsf[0] = 0.0;
    ddsf[0] = 0.0;
  }
  else if (porder == 2){
    const double k0 = knots[0], k1 = knots[1];
    sf[0] = (a - k1)/(k0 - k1);
    sf[1] = (a - k0)/(k1 - k0);

    dsf[0] = 1.0/(k0 - k1);
    dsf[1] = 1.0/(k1 - k0);

    ddsf[0] = ddsf[1] = 0.0;
  }
  else if (porder == 3){
    const double k0 = knots[0], k1 = knots[1], k2 = knots[2];
    sf[0] = (a - k1)*(a - k2)/((k0 - k1)*(k0 - k2));
    sf[1] = (a - k0)*(a - k2)/((k1 - k0)*(k1 - k2));
    sf[2] = (a - k0)*(a - k1)/((k2 - k0)*(k2 - k1));

    dsf[0] = ((a - k1) + (a - k2))/((k0 - k1)*(k0 - k2));
    dsf[1] = ((a - k0) + (a - k2))/((k1 - k0)*(k1 - k2));
    dsf[2] = ((a - k0) + (a - k1))/((k2 - k0)*(k2 - k1));

    ddsf[0] = 2.0/((k0 - k1)*(k0 - k2));
    ddsf[1] = 2.0/((k1 - k0)*(k1 - k2));
    ddsf[2] = 2.0/((k2 - k0)*(k2 - k1));
  }
  else {
    double dn, ddn;
    
    for ( int i = 0; i < porder; i++ ){
      sf[i] = 1.0;
      dsf[i] = 0.0;
      ddsf[i] = 0.0;
    }
    
    // Loop over the shape function control points
    for ( int i = 0; i < porder; i++ ){
      double ki = knots[i];
      
      // Loop over each point again, except for the current control
      // point, adding the contribution to the shape function
      for ( int j = 0; j < porder; j++ ){
        if (i != j){
          double kj = knots[j];
          sf[i] = sf[i]*(a - kj)/(ki - kj);
          
          // Loop over the whole thing again to determine the
          // contribution to the derivative of the shape function
          dn = 1.0/(ki - kj);
          
          for ( int k = 0; k < porder; k++ ){
            if (k != i && k != j){
              double kk = knots[k];
              
              dn = dn*(a - kk)/(ki - kk);             
              ddn = 1.0/((ki - kk)*(ki - kj));
              
              for ( int m = 0; m < porder; m++ ){
                if (((m != i) && (m != j)) && (m != k)){
                  ddn = ddn*(a - knots[m])/(ki - knots[m]);
                }
              }
              ddsf[i] += ddn;
            }
          }     
          dsf[i] += dn;
        }
      }
    }
  }
}

/*
  The following code computes tensor-product Lagrange polynomail shape
  functions up to 8th order.  
*/

/*
  Compute a Lagrange tensor-product interpolant for the
  two-dimensional case. 

  input:
  order: the order of the Lagrange polynomial >= 1
  gpt: the quadrature point on the interval (-1,1)^2

  output:
  N: the shape functions
*/
inline void biLagrangeSF( double N[], 
                          const double gpt[], const int order ){
  double na[8], nb[8];

  if (order == 1){
    N[0] = 1.0;
  }
  else if (order == 2){
    lagrangeSF(na, gpt[0], 2);
    lagrangeSF(nb, gpt[1], 2);

    N[0] = na[0]*nb[0];
    N[1] = na[1]*nb[0];
    N[2] = na[0]*nb[1];
    N[3] = na[1]*nb[1];
  }
  else if (order == 3){
    lagrangeSF(na, gpt[0], 3);
    lagrangeSF(nb, gpt[1], 3);

    N[0] = na[0]*nb[0];
    N[1] = na[1]*nb[0];
    N[2] = na[2]*nb[0];
    N[3] = na[0]*nb[1];
    N[4] = na[1]*nb[1];
    N[5] = na[2]*nb[1];
    N[6] = na[0]*nb[2];
    N[7] = na[1]*nb[2];
    N[8] = na[2]*nb[2];
  }
  else if (order == 4){
    FElibrary::lagrangeSF(na, gpt[0], 4);
    FElibrary::lagrangeSF(nb, gpt[1], 4);

    N[0] = na[0]*nb[0];
    N[1] = na[1]*nb[0];
    N[2] = na[2]*nb[0];
    N[3] = na[3]*nb[0];
    N[4] = na[0]*nb[1];
    N[5] = na[1]*nb[1];
    N[6] = na[2]*nb[1];
    N[7] = na[3]*nb[1];
    N[8] = na[0]*nb[2];
    N[9] = na[1]*nb[2];
    N[10] = na[2]*nb[2];
    N[11] = na[3]*nb[2];
    N[12] = na[0]*nb[3];
    N[13] = na[1]*nb[3];
    N[14] = na[2]*nb[3];
    N[15] = na[3]*nb[3];
  }
  else {
    FElibrary::lagrangeSF(na, gpt[0], order);
    FElibrary::lagrangeSF(nb, gpt[1], order);
    
    for ( int j = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++ ){
        N[0] = na[i]*nb[j];
        N++;
      }
    }
  }
}

/*
  Compute a Lagrange tensor-product interpolant for the
  two-dimensional case. 

  input:
  order: the order of the Lagrange polynomial >= 1
  gpt: the quadrature point on the interval (-1,1)^2

  output:
  N: the shape functions
  Na, Nb: the derivative of the shape functions along the
  parametric directions
*/
inline void biLagrangeSF( double N[], double Na[], double Nb[], 
                          const double gpt[], const int order ){
  double na[8], nb[8];
  double dna[8], dnb[8];

  if (order <= 1){
    N[0] = 1.0;
    Na[0] = Nb[0] = 0.0;
  }
  else if (order == 2){
    lagrangeSF(na, dna, gpt[0], 2);
    lagrangeSF(nb, dnb, gpt[1], 2);

    N[0] = na[0]*nb[0];
    N[1] = na[1]*nb[0];
    N[2] = na[0]*nb[1];
    N[3] = na[1]*nb[1];

    Na[0] = dna[0]*nb[0];
    Na[1] = dna[1]*nb[0];
    Na[2] = dna[0]*nb[1];
    Na[3] = dna[1]*nb[1];
    
    Nb[0] = na[0]*dnb[0];
    Nb[1] = na[1]*dnb[0];
    Nb[2] = na[0]*dnb[1];
    Nb[3] = na[1]*dnb[1];
  }
  else if (order == 3){
    lagrangeSF(na, dna, gpt[0], 3);
    lagrangeSF(nb, dnb, gpt[1], 3);

    N[0] = na[0]*nb[0];
    N[1] = na[1]*nb[0];
    N[2] = na[2]*nb[0];
    N[3] = na[0]*nb[1];
    N[4] = na[1]*nb[1];
    N[5] = na[2]*nb[1];
    N[6] = na[0]*nb[2];
    N[7] = na[1]*nb[2];
    N[8] = na[2]*nb[2];

    Na[0] = dna[0]*nb[0];
    Na[1] = dna[1]*nb[0];
    Na[2] = dna[2]*nb[0];
    Na[3] = dna[0]*nb[1];
    Na[4] = dna[1]*nb[1];
    Na[5] = dna[2]*nb[1];
    Na[6] = dna[0]*nb[2];
    Na[7] = dna[1]*nb[2];
    Na[8] = dna[2]*nb[2];
    
    Nb[0] = na[0]*dnb[0];
    Nb[1] = na[1]*dnb[0];
    Nb[2] = na[2]*dnb[0];
    Nb[3] = na[0]*dnb[1];
    Nb[4] = na[1]*dnb[1];
    Nb[5] = na[2]*dnb[1];
    Nb[6] = na[0]*dnb[2];
    Nb[7] = na[1]*dnb[2];
    Nb[8] = na[2]*dnb[2];
  }
  else if (order == 4){
    FElibrary::lagrangeSF(na, dna, gpt[0], 4);
    FElibrary::lagrangeSF(nb, dnb, gpt[1], 4);

    N[0] = na[0]*nb[0];
    N[1] = na[1]*nb[0];
    N[2] = na[2]*nb[0];
    N[3] = na[3]*nb[0];
    N[4] = na[0]*nb[1];
    N[5] = na[1]*nb[1];
    N[6] = na[2]*nb[1];
    N[7] = na[3]*nb[1];
    N[8] = na[0]*nb[2];
    N[9] = na[1]*nb[2];
    N[10] = na[2]*nb[2];
    N[11] = na[3]*nb[2];
    N[12] = na[0]*nb[3];
    N[13] = na[1]*nb[3];
    N[14] = na[2]*nb[3];
    N[15] = na[3]*nb[3]; 
  
    Na[0] = dna[0]*nb[0];
    Na[1] = dna[1]*nb[0];
    Na[2] = dna[2]*nb[0];
    Na[3] = dna[3]*nb[0];
    Na[4] = dna[0]*nb[1];
    Na[5] = dna[1]*nb[1];
    Na[6] = dna[2]*nb[1];
    Na[7] = dna[3]*nb[1];
    Na[8] = dna[0]*nb[2];
    Na[9] = dna[1]*nb[2];
    Na[10] = dna[2]*nb[2];
    Na[11] = dna[3]*nb[2];
    Na[12] = dna[0]*nb[3];
    Na[13] = dna[1]*nb[3];
    Na[14] = dna[2]*nb[3];
    Na[15] = dna[3]*nb[3];
    
    Nb[0] = na[0]*dnb[0];
    Nb[1] = na[1]*dnb[0];
    Nb[2] = na[2]*dnb[0];
    Nb[3] = na[3]*dnb[0];
    Nb[4] = na[0]*dnb[1];
    Nb[5] = na[1]*dnb[1];
    Nb[6] = na[2]*dnb[1];
    Nb[7] = na[3]*dnb[1];
    Nb[8] = na[0]*dnb[2];
    Nb[9] = na[1]*dnb[2];
    Nb[10] = na[2]*dnb[2];
    Nb[11] = na[3]*dnb[2];
    Nb[12] = na[0]*dnb[3];
    Nb[13] = na[1]*dnb[3];
    Nb[14] = na[2]*dnb[3];
    Nb[15] = na[3]*dnb[3];
  }
  else {
    FElibrary::lagrangeSF(na, dna, gpt[0], order);
    FElibrary::lagrangeSF(nb, dnb, gpt[1], order);
    
    for ( int j = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++ ){
        N[0] = na[i]*nb[j];
        Na[0] = dna[i]*nb[j];
        Nb[0] = na[i]*dnb[j];
        N++; Na++; Nb++;
      }
    }
  }
}

/*
  Compute a Lagrange tensor-product interpolant for the
  three-dimensional case. 

  input:
  order: the order of the Lagrange polynomial >= 1
  gpt: the quadrature point on the interval (-1,1)^3

  output:
  N: the shape functions
*/
inline void triLagrangeSF( double N[], 
                           const double gpt[], const int order ){
  double na[8], nb[8], nc[8];

  if (order <= 1){
    N[0] = 1.0;
  }
  else if (order == 2){
    lagrangeSF(na, gpt[0], 2);
    lagrangeSF(nb, gpt[1], 2);
    lagrangeSF(nc, gpt[2], 2);

    N[0] = na[0]*nb[0]*nc[0];
    N[1] = na[1]*nb[0]*nc[0];
    N[2] = na[0]*nb[1]*nc[0];
    N[3] = na[1]*nb[1]*nc[0];
    N[4] = na[0]*nb[0]*nc[1];
    N[5] = na[1]*nb[0]*nc[1];
    N[6] = na[0]*nb[1]*nc[1];
    N[7] = na[1]*nb[1]*nc[1];
  }
  else {
    lagrangeSF(na, gpt[0], order);
    lagrangeSF(nb, gpt[1], order);
    lagrangeSF(nc, gpt[2], order);

    for ( int k = 0; k < order; k++ ){
      for ( int j = 0; j < order; j++ ){
        for ( int i = 0; i < order; i++ ){
          N[0] = na[i]*nb[j]*nc[k];
          N++;
        }
      }
    }
  }
}

/*
  Compute a Lagrange tensor-product interpolant for the
  three-dimensional case. 

  input:
  order: the order of the Lagrange polynomial >= 1
  gpt: the quadrature point on the interval (-1,1)^3

  output:
  N: the shape functions
  Na, Nb, Nc: the derivative of the shape functions along the
  parametric directions
*/
inline void triLagrangeSF( double N[], double Na[], double Nb[], double Nc[], 
                           const double gpt[], const int order ){
  double na[8], nb[8], nc[8];
  double dna[8], dnb[8], dnc[8];

  if (order <= 1){
    N[0] = 1.0;
    Na[0] = Nb[0] = Nc[0] = 0.0;
  }
  else if (order == 2){
    lagrangeSF(na, dna, gpt[0], 2);
    lagrangeSF(nb, dnb, gpt[1], 2);
    lagrangeSF(nc, dnc, gpt[2], 2);

    N[0] = na[0]*nb[0]*nc[0];
    N[1] = na[1]*nb[0]*nc[0];
    N[2] = na[0]*nb[1]*nc[0];
    N[3] = na[1]*nb[1]*nc[0];
    N[4] = na[0]*nb[0]*nc[1];
    N[5] = na[1]*nb[0]*nc[1];
    N[6] = na[0]*nb[1]*nc[1];
    N[7] = na[1]*nb[1]*nc[1];

    Na[0] = dna[0]*nb[0]*nc[0];
    Na[1] = dna[1]*nb[0]*nc[0];
    Na[2] = dna[0]*nb[1]*nc[0];
    Na[3] = dna[1]*nb[1]*nc[0];
    Na[4] = dna[0]*nb[0]*nc[1];
    Na[5] = dna[1]*nb[0]*nc[1];
    Na[6] = dna[0]*nb[1]*nc[1];
    Na[7] = dna[1]*nb[1]*nc[1];

    Nb[0] = na[0]*dnb[0]*nc[0];
    Nb[1] = na[1]*dnb[0]*nc[0];
    Nb[2] = na[0]*dnb[1]*nc[0];
    Nb[3] = na[1]*dnb[1]*nc[0];
    Nb[4] = na[0]*dnb[0]*nc[1];
    Nb[5] = na[1]*dnb[0]*nc[1];
    Nb[6] = na[0]*dnb[1]*nc[1];
    Nb[7] = na[1]*dnb[1]*nc[1];

    Nc[0] = na[0]*nb[0]*dnc[0];
    Nc[1] = na[1]*nb[0]*dnc[0];
    Nc[2] = na[0]*nb[1]*dnc[0];
    Nc[3] = na[1]*nb[1]*dnc[0];
    Nc[4] = na[0]*nb[0]*dnc[1];
    Nc[5] = na[1]*nb[0]*dnc[1];
    Nc[6] = na[0]*nb[1]*dnc[1];
    Nc[7] = na[1]*nb[1]*dnc[1];
  }
  else {
    lagrangeSF(na, dna, gpt[0], order);
    lagrangeSF(nb, dnb, gpt[1], order);
    lagrangeSF(nc, dnc, gpt[2], order);

    for ( int k = 0; k < order; k++ ){
      for ( int j = 0; j < order; j++ ){
        for ( int i = 0; i < order; i++ ){
          N[0] = na[i]*nb[j]*nc[k];
          Na[0] = dna[i]*nb[j]*nc[k];
          Nb[0] = na[i]*dnb[j]*nc[k];
          Nc[0] = na[i]*nb[j]*dnc[k];
          N++;
          Na++; Nb++; Nc++;
        }
      }
    }
  }
}

/*
  Test the derivative of the b-spline basis
  
  input:
  k: the order of the basis to test
*/
void bspline_basis_test( int k );

/*
  Find the interval for the computing the basis function

  input:
  u: the parametric location
  T: the knot locations
  n: the number of control points
  k: the order of the b-spline
*/
int bspline_interval( double u, const double * T, int n, int k );

/*
  Evaluate the basis functions and optionally their derivatives

  output:
  N:  the basis functions - an array of size k

  input:
  idx:  the index of the interval for u
  u:    the parametric position
  Tu:   the knot vector
  ku:   the order of the basis functions
  work: a temporary array of size 2*k

  u is on the idx-th knot span such that u is in the interval
  u \in [Tu[idx], Tu[idx+1]) 
*/
void bspline_basis( double * N, const int idx, const double u, 
                    const double * Tu, 
                    const int ku, double * work );

/*
  Compute the derivative of the b-spline basis

  output:
  N:  the basis functions - an array of size k*(ideriv+1)

  input:
  idx:  the index of the interval for u
  u:    the parametric position
  Tu:   the knot vector
  ku:   the order of the basis functions
  work: a temporary array of size 2*k + k*k

  u is on the idx-th knot span such that u is in the interval
  u \in [Tu[idx], Tu[idx+1]) 
*/
void bspline_basis_derivative( double * N, const int idx, const double u,
                               int ideriv, const double * Tu,
                               const int ku, double * work );

/*
  Evaluate the one-dimensional b-spline

  input:
  u:     the parametric location for the spline evaluation
  idu:   the order of the derivative to use 
  T:     the knot vector of length n + k
  n:     the number of knots
  k:     the order of the spline to evaluate
  coef:  the spline coefficients
  work:  a working array for temporary storage

  the work array must be of size:
  if idu == 0: len = 3*ku   
  otherwise: len = (idu+3)*ku + ku*ku

  returns:
  the value of the interpolant (or its derivative) at u
*/
TacsScalar bspline1d( const double u, const int idu, const double * Tu, 
                      const int nu, const int ku, const TacsScalar * coef,
                      double * work );

/*
  Evaluate a two-dimensional tensor product b-spline

  input:
  u, v:       the parametric location for the spline evaluation
  idu, idv:   the order of the derivative to use 
  Tu, Tv:     the knot vector of length n + k
  nu, nv:     the number of knots
  ku, kv:     the order of the spline to evaluate
  coef:       the spline coefficients
  work:       a working array for temporary storage

  the work array must be of size:
  if idu == 0: len = ku + kv + 2*max(ku, kv)
  otherwise: len = (idu+1)*ku + (idv+1)*kv + max(2*ku + ku**2, 2*kv + kv**2)

  returns:
  the value of the interpolant (or its derivative) at u
*/
TacsScalar bspline2d( const double u, const double v, 
                      const int idu, const int idv,
                      const double * Tu, const double * Tv,
                      const int nu, const int nv, const int ku, const int kv, 
                      const TacsScalar * coef,
                      double * work );

/*
  Evaluate a three-dimensional tensor product b-spline

  input:
  u, v, w:        the parametric location for the spline evaluation
  idu, idv, idw:  the order of the derivative to use 
  Tu, Tv, Tw:     the knot vector of length n + k
  nu, nv, nw:     the number of knots
  ku, kv, kw:     the order of the spline to evaluate
  coef:           the spline coefficients
  work:           a working array for temporary storage

  the work array must be of size:
  if idu == 0: len = ku + kv + kw + 2*max(ku, kv, kw)
  otherwise: 
  len = (idu+1)*ku + (idv+1)*kv + (idw+1)*kw + 
  max(2*ku + ku**2, 2*kv + kv**2, 2*kw + kw**2)

  returns:
  the value of the interpolant (or its derivative) at u
*/
TacsScalar bspline3d( const double u, const double v, const double w, 
                      const int idu, const int idv, const int idw,
                      const double * Tu, const double * Tv, const double * Tw,
                      const int nu, const int nv, const int nw, 
                      const int ku, const int kv, const int kw, 
                      const TacsScalar * coef,
                      double * work );

/*
  C1 functions for one dimensional problems
*/
void cubicHP( double N[], double Na[], double Naa[], double a );
void quinticHP( double N[], double Na[], double Naa[], double a );

/*
  The following are the definitions of the Gauss (or Gauss-Legendre)
  quadrature points and weights for the interval [-1, 1]. Note that
  these schemes are exact for polynomials of degree 2n - 1.
*/
const double gaussPts1[] = { 0.0 };
const double gaussWts1[] = { 2.0 };
  
const double gaussPts2[] = { -0.577350269189626, 0.577350269189626 };
const double gaussWts2[] = {  1.0,               1.0 };
  
const double gaussPts3[] = { -0.774596669241483, 0.0,  0.774596669241483 };
const double gaussWts3[] = {  5.0/9.0,       8.0/9.0,            5.0/9.0 };

const double gaussPts4[] = { -0.861136311594053, -0.339981043584856, 
                             0.339981043584856, 0.861136311594053 };
const double gaussWts4[] = {  0.347854845137454,  0.652145154862546, 
                              0.652145154862546, 0.347854845137454 };

const double gaussPts5[] = { -0.906179845938664, -0.538469310105683, 
                             0.0,               
                             0.538469310105683, 0.906179845938664 };
const double gaussWts5[] = {  0.236926885056189,  0.478628670499366, 
                              0.568888888888889, 
                              0.478628670499366, 0.236926885056189 };
                                                          
const double gaussPts6[] = { -0.9324695142031520278123016, 
                             -0.6612093864662645136613996, 
                             -0.2386191860831969086305017, 
                             0.2386191860831969086305017, 
                             0.6612093864662645136613996, 
                             0.9324695142031520278123016  };
const double gaussWts6[] = {  0.1713244923791703450402961,  
                              0.3607615730481386075698335,  
                              0.4679139345726910473898703, 
                              0.4679139345726910473898703, 
                              0.3607615730481386075698335, 
                              0.1713244923791703450402961 };

const double gaussPts7[] = { -0.9491079123427585245261897, 
                             -0.7415311855993944398638648, 
                             -0.4058451513773971669066064, 
                             0.0,
                             0.4058451513773971669066064, 
                             0.7415311855993944398638648, 
                             0.9491079123427585245261897 };
const double gaussWts7[] = { 0.1294849661688696932706114, 
                             0.2797053914892766679014678, 
                             0.3818300505051189449503698, 
                             0.4179591836734693877551020,
                             0.3818300505051189449503698, 
                             0.2797053914892766679014678, 
                             0.1294849661688696932706114 };

const double gaussPts8[] = { -0.9602898564975362316835609, 
                             -0.7966664774136267395915539, 
                             -0.5255324099163289858177390, 
                             -0.1834346424956498049394761,
                             0.1834346424956498049394761, 
                             0.5255324099163289858177390, 
                             0.7966664774136267395915539, 
                             0.9602898564975362316835609 };
const double gaussWts8[] = { 0.1012285362903762591525314, 
                             0.2223810344533744705443560, 
                             0.3137066458778872873379622, 
                             0.3626837833783619829651504,
                             0.3626837833783619829651504, 
                             0.3137066458778872873379622, 
                             0.2223810344533744705443560, 
                             0.1012285362903762591525314 };

/*
  The following are the Gauss--Lobatto integration schemes. These
  quadrature schemes always include the starting/end points for the
  interval. These schemes are less accurate in that they integrate
  polynomials of degree 2n-3 exactly (compared to Gauss quadrature
  schemes which integrated 2n-1 exactly).
*/
const double lobattoPts2[] = { -1.0, 1.0 };
const double lobattoWts2[] = {  1.0, 1.0 };

const double lobattoPts3[] = { -1.0, 0.0, 1.0 };
const double lobattoWts3[] = { 1.0/3.0, 4.0/3.0, 1.0/3.0 };

const double lobattoPts4[] = { -1.0, -0.44721359549995793,
                               0.44721359549995793, 1.0 };
const double lobattoWts4[] = { 1.0/6.0, 5.0/6.0, 5.0/6.0, 1.0/6.0 };

const double lobattoPts5[] = { -1.0, -0.65465367070797709, 0.0,
                               0.65465367070797709, 1.0 };
const double lobattoWts5[] = { 1.0/10.0, 49.0/90.0, 32.0/45.0,
                               49.0/90.0, 1.0/10.0 };

const double lobattoPts6[] = { -1.0, -0.76505532392946474,
                               -0.2852315164806451, 
                               0.2852315164806451,
                               0.76505532392946474, 1.0 };
const double lobattoWts6[] = { 1.0/15.0, 
                               0.378474956297847,
                               0.55485837703548635,
                               0.55485837703548635,
                               0.378474956297847,
                               1.0/15.0 };

/*
  Retrieve the pointer to the specified Gauss quadrature weight and
  point scheme for up to the 8-point scheme written out above.

  input:
  numGauss: the desired number of points

  output:
  gaussPts: the Gauss quadrature points
  gaussWts: the Gauss quadrature weights
*/
int getGaussPtsWts( int _numGauss, 
                    const double ** gaussPts, 
                    const double ** gaussWts );

/*
  Retrieve a given Gauss-quadrature weight-point pair for
  the specified scheme.

  input:
  dim:    the dimension of the domain i.e. 1, 2 or 3
  npts:   the number of points in the scheme
  num:    the point within the domain

  output:
  pt:    the parametric point within the domain 
*/
double getGaussPtWt( const int dim, const int npts,
                     const int num, double pt[] );

/*
  Retrieve the pointer to the specified Gauss quadrature weight and
  point scheme for up to the 8-point scheme written out above.

  input:
  quad:     the type of quadrature scheme to return
  numGauss: the desired number of points

  output:
  gaussPts: the Gauss quadrature points
  gaussWts: the Gauss quadrature weights
*/
int getGaussPtsWts( enum QuadratureType quad, 
                    int _numGauss, 
                    const double ** gaussPts, 
                    const double ** gaussWts );

/*
  Retrieve a given Gauss-quadrature weight-point pair for
  the specified scheme.

  input:
  quad:   the type of quadrature scheme to return
  dim:    the dimension of the domain i.e. 1, 2 or 3
  npts:   the number of points in the scheme
  num:    the point within the domain

  output:
  pt:    the parametric point within the domain 
*/
double getGaussPtWt( enum QuadratureType quad, 
                     const int dim, const int npts,
                     const int num, double pt[] );


/*! 
  Calculate the inverse of the given 2x2 matrix
*/
template <class ScalarType> 
inline ScalarType jacobian2d( const ScalarType Xd[], 
                              ScalarType Jinv[] ){
  ScalarType h = Xd[0]*Xd[3] - Xd[1]*Xd[2];
  ScalarType hinv = 1.0/h;

  Jinv[0] =   Xd[3]*hinv;
  Jinv[1] = - Xd[1]*hinv;
    
  Jinv[2] = - Xd[2]*hinv;
  Jinv[3] =   Xd[0]*hinv;

  return h;
}

/*
  Calculate the determinant of a 2x2 matrix
*/
template <class ScalarType> 
inline ScalarType jacobian2d( const ScalarType Xd[] ){
  return Xd[0]*Xd[3] - Xd[1]*Xd[2];
}

/*
  Calculate the derivative of the inverse of the 2x2 matrix
*/
template <class ScalarType> 
inline ScalarType jacobian2dSens( const ScalarType Xd[], 
                                  const ScalarType sXd[], 
                                  ScalarType Jinv[], ScalarType sJinv[], 
                                  ScalarType * _sh ){
  ScalarType h = Xd[0]*Xd[3] - Xd[1]*Xd[2];
  ScalarType hinv = 1.0/h;

  ScalarType sh = (Xd[0]*sXd[3] + sXd[0]*Xd[3] -
                   Xd[1]*sXd[2] - sXd[1]*Xd[2]);

  Jinv[0] =   Xd[3]*hinv;
  Jinv[1] = - Xd[1]*hinv;
    
  Jinv[2] = - Xd[2]*hinv;
  Jinv[3] =   Xd[0]*hinv;

  sJinv[0] =   (sXd[3]*h - Xd[3]*sh)*hinv*hinv;
  sJinv[3] =   (sXd[0]*h - Xd[0]*sh)*hinv*hinv;
  sJinv[1] = - (sXd[1]*h - Xd[1]*sh)*hinv*hinv;
  sJinv[2] = - (sXd[2]*h - Xd[2]*sh)*hinv*hinv;

  *_sh = sh;
  return h;
}

/*
  Calculate the derivative of the determinant of a 2x2 matrix
*/
template <class ScalarType> 
inline ScalarType jacobian2dSens( const ScalarType Xd[], 
                                  const ScalarType sXd[], 
                                  ScalarType *sh ){
  ScalarType h = Xd[0]*Xd[3] - Xd[1]*Xd[2];

  *sh = (Xd[0]*sXd[3] + sXd[0]*Xd[3] -
         Xd[1]*sXd[2] - sXd[1]*Xd[2]);

  return h;
}

/*! 
  Invert the matrix Xd to obtain the inverse of the Jacobian: Jinv

  input: 
  Xd: the Jacobian matrix
  
  output:
  Jinv: the Jacobian inverse
  
  returns: the determinant of Xd
*/
template <class ScalarType> 
inline ScalarType jacobian3d( const ScalarType Xd[], 
                              ScalarType Jinv[] ){
  ScalarType h = (Xd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1]) 
                  - Xd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2]) 
                  + Xd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]));
  ScalarType hinv = 1.0/h;

  Jinv[0] =   (Xd[4]*Xd[8] - Xd[5]*Xd[7])*hinv;
  Jinv[1] = - (Xd[1]*Xd[8] - Xd[2]*Xd[7])*hinv;
  Jinv[2] =   (Xd[1]*Xd[5] - Xd[2]*Xd[4])*hinv;
    
  Jinv[3] = - (Xd[3]*Xd[8] - Xd[5]*Xd[6])*hinv;
  Jinv[4] =   (Xd[0]*Xd[8] - Xd[2]*Xd[6])*hinv;
  Jinv[5] = - (Xd[0]*Xd[5] - Xd[2]*Xd[3])*hinv;
    
  Jinv[6] =   (Xd[3]*Xd[7] - Xd[4]*Xd[6])*hinv;
  Jinv[7] = - (Xd[0]*Xd[7] - Xd[1]*Xd[6])*hinv;
  Jinv[8] =   (Xd[0]*Xd[4] - Xd[1]*Xd[3])*hinv;
    
  return h;
}

/*
  Compute the determinant of the Jacobian matrix Xd

  input: Xd the Jacobian
  returns the determinant
*/
template <class ScalarType> 
inline ScalarType jacobian3d( const ScalarType Xd[] ){
  ScalarType h;

  h = Xd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1]) 
    - Xd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2]) 
    + Xd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]);

  return h;
}

/*
  Compute the sensitivity of the determinant of the Jacobian matrix
  Xd.

  input: 
  Xd: the Jacobian
  sXd: the sensitivity of the Jacobian

  output:
  s: the sensitivity of the determinant of the Jacobian

  returns the determinant of the Jacobian
*/
template <class ScalarType> 
inline ScalarType jacobian3dSens( const ScalarType Xd[], 
                                  const ScalarType sXd[],
                                  ScalarType * s ){
  ScalarType h = (Xd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1]) 
                  - Xd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2]) 
                  + Xd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]));

  *s= sXd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1]) + 
    Xd[8]*(Xd[0]*sXd[4] + sXd[0]*Xd[4] - Xd[3]*sXd[1] - sXd[3]*Xd[1]) 
    - sXd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2]) - 
    Xd[7]*(Xd[0]*sXd[5] + sXd[0]*Xd[5] - Xd[3]*sXd[2] - sXd[3]*Xd[2]) 
    + sXd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]) + 
    Xd[6]*(Xd[1]*sXd[5] + sXd[1]*Xd[5] - Xd[2]*sXd[4] - sXd[2]*Xd[4]);

  return h;
}

/*
  Compute the sensitivity of the inverse of the Jacobian

  input:
  Xd: the input Jacobian
  sXd: the sensitivity of the input Jacobian

  output:
  Jinv: the inverse of the Jacobian
  sJinv: the sensitivity of the Jacobian
  sh: the sensitivity of the determinant of the Jacobian

  returns the determinant of Xd
*/
template <class ScalarType> 
inline ScalarType jacobian3dSens( const ScalarType Xd[], 
                                  const ScalarType sXd[], 
                                  ScalarType Jinv[], ScalarType sJinv[], 
                                  ScalarType * _sh ){

  ScalarType h = (Xd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1]) 
                  - Xd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2]) 
                  + Xd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]));
  ScalarType hinv = 1.0/h;

  ScalarType sh;
  sh= sXd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1]) + 
    Xd[8]*(Xd[0]*sXd[4] + sXd[0]*Xd[4] - Xd[3]*sXd[1] - sXd[3]*Xd[1]) 
    - sXd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2]) - 
    Xd[7]*(Xd[0]*sXd[5] + sXd[0]*Xd[5] - Xd[3]*sXd[2] - sXd[3]*Xd[2]) 
    + sXd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]) + 
    Xd[6]*(Xd[1]*sXd[5] + sXd[1]*Xd[5] - Xd[2]*sXd[4] - sXd[2]*Xd[4]);

  Jinv[0] =   (Xd[4]*Xd[8] - Xd[5]*Xd[7])*hinv;
  Jinv[1] = - (Xd[1]*Xd[8] - Xd[2]*Xd[7])*hinv;
  Jinv[2] =   (Xd[1]*Xd[5] - Xd[2]*Xd[4])*hinv;
    
  Jinv[3] = - (Xd[3]*Xd[8] - Xd[5]*Xd[6])*hinv;
  Jinv[4] =   (Xd[0]*Xd[8] - Xd[2]*Xd[6])*hinv;
  Jinv[5] = - (Xd[0]*Xd[5] - Xd[2]*Xd[3])*hinv;
    
  Jinv[6] =   (Xd[3]*Xd[7] - Xd[4]*Xd[6])*hinv;
  Jinv[7] = - (Xd[0]*Xd[7] - Xd[1]*Xd[6])*hinv;
  Jinv[8] =   (Xd[0]*Xd[4] - Xd[1]*Xd[3])*hinv;

  for ( int i = 0; i < 9; i++ ){
    sJinv[i] = - Jinv[i]*hinv*sh;
  }

  sJinv[0] +=  (Xd[4]*sXd[8] + sXd[4]*Xd[8] - Xd[5]*sXd[7] - sXd[5]*Xd[7])*hinv;
  sJinv[1] += -(Xd[1]*sXd[8] + sXd[1]*Xd[8] - Xd[2]*sXd[7] - sXd[2]*Xd[7])*hinv;
  sJinv[2] +=  (Xd[1]*sXd[5] + sXd[1]*Xd[5] - Xd[2]*sXd[4] - sXd[2]*Xd[4])*hinv;
    
  sJinv[3] += -(Xd[3]*sXd[8] + sXd[3]*Xd[8] - Xd[5]*sXd[6] - sXd[5]*Xd[6])*hinv;
  sJinv[4] +=  (Xd[0]*sXd[8] + sXd[0]*Xd[8] - Xd[2]*sXd[6] - sXd[2]*Xd[6])*hinv;
  sJinv[5] += -(Xd[0]*sXd[5] + sXd[0]*Xd[5] - Xd[2]*sXd[3] - sXd[2]*Xd[3])*hinv;
    
  sJinv[6] +=  (Xd[3]*sXd[7] + sXd[3]*Xd[7] - Xd[4]*sXd[6] - sXd[4]*Xd[6])*hinv;
  sJinv[7] += -(Xd[0]*sXd[7] + sXd[0]*Xd[7] - Xd[1]*sXd[6] - sXd[1]*Xd[6])*hinv;
  sJinv[8] +=  (Xd[0]*sXd[4] + sXd[0]*Xd[4] - Xd[1]*sXd[3] - sXd[1]*Xd[3])*hinv;

  *_sh = sh;    
  return h;
}

/*
  Compute the sensitivity of the inverse of the Jacobian
  given that the inverse has already been computed.

  input:
  Xd: the input Jacobian
  sXd: the sensitivity of the input Jacobian
  Jinv: the inverse of the Jacobian
  h: the determinant of the Jacobian

  output:
  sJinv: the sensitivity of the Jacobian
  sh: the sensitivity of the determinant of the Jacobian
*/
template <class ScalarType> 
inline void jacobian3dSens( const ScalarType Xd[], const ScalarType sXd[], 
                            ScalarType h,
                            const ScalarType Jinv[], ScalarType sJinv[], 
                            ScalarType * _sh ){
  ScalarType hinv = 1.0/h;

  ScalarType sh;
  sh= sXd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1]) + 
    Xd[8]*(Xd[0]*sXd[4] + sXd[0]*Xd[4] - Xd[3]*sXd[1] - sXd[3]*Xd[1]) 
    - sXd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2]) - 
    Xd[7]*(Xd[0]*sXd[5] + sXd[0]*Xd[5] - Xd[3]*sXd[2] - sXd[3]*Xd[2]) 
    + sXd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]) + 
    Xd[6]*(Xd[1]*sXd[5] + sXd[1]*Xd[5] - Xd[2]*sXd[4] - sXd[2]*Xd[4]);

  ScalarType scale = - sh*hinv;
  for ( int i = 0; i < 9; i++ ){
    sJinv[i] = scale*Jinv[i];
  }

  sJinv[0] +=  (Xd[4]*sXd[8] + sXd[4]*Xd[8] - Xd[5]*sXd[7] - sXd[5]*Xd[7])*hinv;
  sJinv[1] += -(Xd[1]*sXd[8] + sXd[1]*Xd[8] - Xd[2]*sXd[7] - sXd[2]*Xd[7])*hinv;
  sJinv[2] +=  (Xd[1]*sXd[5] + sXd[1]*Xd[5] - Xd[2]*sXd[4] - sXd[2]*Xd[4])*hinv;
    
  sJinv[3] += -(Xd[3]*sXd[8] + sXd[3]*Xd[8] - Xd[5]*sXd[6] - sXd[5]*Xd[6])*hinv;
  sJinv[4] +=  (Xd[0]*sXd[8] + sXd[0]*Xd[8] - Xd[2]*sXd[6] - sXd[2]*Xd[6])*hinv;
  sJinv[5] += -(Xd[0]*sXd[5] + sXd[0]*Xd[5] - Xd[2]*sXd[3] - sXd[2]*Xd[3])*hinv;
    
  sJinv[6] +=  (Xd[3]*sXd[7] + sXd[3]*Xd[7] - Xd[4]*sXd[6] - sXd[4]*Xd[6])*hinv;
  sJinv[7] += -(Xd[0]*sXd[7] + sXd[0]*Xd[7] - Xd[1]*sXd[6] - sXd[1]*Xd[6])*hinv;
  sJinv[8] +=  (Xd[0]*sXd[4] + sXd[0]*Xd[4] - Xd[1]*sXd[3] - sXd[1]*Xd[3])*hinv;

  *_sh = sh;
}

TACS_END_NAMESPACE

#endif
