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

#ifndef TACS_TENSOR_TOOLBOX_H
#define TACS_TENSOR_TOOLBOX_H

/*!
  This library of functions defines handy operations for 2 and 3
  dimensional tensors of various orders. Most importantly, tensor
  transformations are handled.
*/

#include <math.h>

#include "FElibrary.h"

TACS_BEGIN_NAMESPACE(Tensor)

// The following are some basic 2D and 3D vector operations that are useful

/*
  Find the Cartesian dot product of two 2D vectors
*/
template <class ScalarType>
inline ScalarType dot2D(const ScalarType A[], const ScalarType B[]) {
  return A[0] * B[0] + A[1] * B[1];
}

/*
  Find the Cartesian dot product of two 3D vectors
*/
template <class ScalarType>
inline ScalarType dot3D(const ScalarType A[], const ScalarType B[]) {
  return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}

/*
  Find the Cartesian cross product of two 3D vectors and return
  in the vector out
*/
template <class ScalarType>
inline void crossProduct3D(ScalarType out[], const ScalarType A[],
                           const ScalarType B[]) {
  //   (A2 * B3 - A3 * B2 ), ( A3 *B1 - A1 * B3 ), ( A1 * B2 - A2 * B1 )

  out[0] = A[1] * B[2] - A[2] * B[1];
  out[1] = A[2] * B[0] - A[0] * B[2];
  out[2] = A[0] * B[1] - A[1] * B[0];
}

/*
  Find the cross product and its sensitivity with respect to two
  3D input vectors and their sensitivity
*/
template <class ScalarType>
inline void crossProduct3DSens(ScalarType out[], ScalarType sout[],
                               const ScalarType A[], const ScalarType B[],
                               const ScalarType sA[], const ScalarType sB[]) {
  //   (A2 * B3 - A3 * B2 ), ( A3 *B1 - A1 * B3 ), ( A1 * B2 - A2 * B1 )

  out[0] = A[1] * B[2] - A[2] * B[1];
  out[1] = A[2] * B[0] - A[0] * B[2];
  out[2] = A[0] * B[1] - A[1] * B[0];

  sout[0] = (sA[1] * B[2] - sA[2] * B[1] + A[1] * sB[2] - A[2] * sB[1]);
  sout[1] = (sA[2] * B[0] - sA[0] * B[2] + A[2] * sB[0] - A[0] * sB[2]);
  sout[2] = (sA[0] * B[1] - sA[1] * B[0] + A[0] * sB[1] - A[1] * sB[0]);
}

/*
  Find the cross product and its sensitivity with respect to two
  3D input vectors and their sensitivity
*/
template <class ScalarType>
inline void crossProduct3DSens(ScalarType sout[], const ScalarType A[],
                               const ScalarType B[], const ScalarType sA[],
                               const ScalarType sB[]) {
  //   (A2 * B3 - A3 * B2 ), ( A3 *B1 - A1 * B3 ), ( A1 * B2 - A2 * B1 )
  sout[0] = (sA[1] * B[2] - sA[2] * B[1] + A[1] * sB[2] - A[2] * sB[1]);
  sout[1] = (sA[2] * B[0] - sA[0] * B[2] + A[2] * sB[0] - A[0] * sB[2]);
  sout[2] = (sA[0] * B[1] - sA[1] * B[0] + A[0] * sB[1] - A[1] * sB[0]);
}

/*
  Find the cross product and its sensitivity with respect to two
  3D input vectors and their sensitivity
*/
template <class ScalarType>
inline void crossProduct3DSens_d(ScalarType out[], ScalarType sout[],
                                 ScalarType outd[], ScalarType soutd[],
                                 const ScalarType A[], const ScalarType B[],
                                 const ScalarType sA[], const ScalarType sB[],
                                 const ScalarType Ad[], const ScalarType Bd[],
                                 const ScalarType sAd[],
                                 const ScalarType sBd[]) {
  //   (A2 * B3 - A3 * B2 ), ( A3 *B1 - A1 * B3 ), ( A1 * B2 - A2 * B1 )
  out[0] = A[1] * B[2] - A[2] * B[1];
  out[1] = A[2] * B[0] - A[0] * B[2];
  out[2] = A[0] * B[1] - A[1] * B[0];

  sout[0] = (sA[1] * B[2] - sA[2] * B[1] + A[1] * sB[2] - A[2] * sB[1]);
  sout[1] = (sA[2] * B[0] - sA[0] * B[2] + A[2] * sB[0] - A[0] * sB[2]);
  sout[2] = (sA[0] * B[1] - sA[1] * B[0] + A[0] * sB[1] - A[1] * sB[0]);

  outd[0] = (Ad[1] * B[2] - Ad[2] * B[1] + A[1] * Bd[2] - A[2] * Bd[1]);
  outd[1] = (Ad[2] * B[0] - Ad[0] * B[2] + A[2] * Bd[0] - A[0] * Bd[2]);
  outd[2] = (Ad[0] * B[1] - Ad[1] * B[0] + A[0] * Bd[1] - A[1] * Bd[0]);

  soutd[0] = (sAd[1] * B[2] - sAd[2] * B[1] + sA[1] * Bd[2] - sA[2] * Bd[1] +
              Ad[1] * sB[2] - Ad[2] * sB[1] + A[1] * sBd[2] - A[2] * sBd[1]);

  soutd[1] = (sAd[2] * B[0] - sAd[0] * B[2] + sA[2] * Bd[0] - sA[0] * Bd[2] +
              Ad[2] * sB[0] - Ad[0] * sB[2] + A[2] * sBd[0] - A[0] * sBd[2]);

  soutd[2] = (sAd[0] * B[1] - sAd[1] * B[0] + sA[0] * Bd[1] - sA[1] * Bd[0] +
              Ad[0] * sB[1] - Ad[1] * sB[0] + A[0] * sBd[1] - A[1] * sBd[0]);
}

/*
  Normalize a 2D vector in place
*/
template <class ScalarType>
inline ScalarType normalize2D(ScalarType A[]) {
  ScalarType Anrm = sqrt(A[0] * A[0] + A[1] * A[1]);
  ScalarType invAnrm = 1.0 / Anrm;

  if (Anrm != 0.0) {
    A[0] = A[0] * invAnrm;
    A[1] = A[1] * invAnrm;
  }

  return Anrm;
}

/*
  Find the sensitivity of the normalized 3D vector
*/
template <class ScalarType>
inline ScalarType normalize2DSens(ScalarType A[], ScalarType sA[],
                                  ScalarType* sAnrm) {
  ScalarType Anrm = sqrt(A[0] * A[0] + A[1] * A[1]);
  ScalarType invAnrm = 1.0 / Anrm;
  *sAnrm = invAnrm * (A[0] * sA[0] + A[1] * sA[1]);

  if (Anrm != 0.0) {
    A[0] = A[0] * invAnrm;
    A[1] = A[1] * invAnrm;

    sA[0] = (sA[0] - A[0] * (*sAnrm)) * invAnrm;
    sA[1] = (sA[1] - A[1] * (*sAnrm)) * invAnrm;
  }

  return Anrm;
}

/*
  Normalize a 3D vector and copy over the result
*/
template <class ScalarType>
inline ScalarType normalize3D(ScalarType out[], const ScalarType A[]) {
  ScalarType Anrm = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
  ScalarType invAnrm = 1.0 / Anrm;

  if (Anrm != 0.0) {
    out[0] = A[0] * invAnrm;
    out[1] = A[1] * invAnrm;
    out[2] = A[2] * invAnrm;
  } else {
    out[0] = out[1] = out[2] = ScalarType(0.0);
  }

  return Anrm;
}

/*
  Find the sensitivity of the normalized 3D vector and
  copy the result to sout
*/
template <class ScalarType>
inline ScalarType normalize3DSens(ScalarType sout[], const ScalarType A[],
                                  const ScalarType sA[],
                                  const ScalarType Anrm) {
  ScalarType sAnrm = (A[0] * sA[0] + A[1] * sA[1] + A[2] * sA[2]) / Anrm;
  ScalarType invAnrm = 1.0 / (Anrm * Anrm);

  if (Anrm != 0.0) {
    sout[0] = (sA[0] * Anrm - A[0] * sAnrm) * invAnrm;
    sout[1] = (sA[1] * Anrm - A[1] * sAnrm) * invAnrm;
    sout[2] = (sA[2] * Anrm - A[2] * sAnrm) * invAnrm;
  } else {
    sout[0] = sout[1] = sout[2] = ScalarType(0.0);
  }

  return sAnrm;
}

/*
  Normalize a 3D vector in place
*/
template <class ScalarType>
inline ScalarType normalize3D(ScalarType A[]) {
  ScalarType Anrm = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
  ScalarType invAnrm = 1.0 / Anrm;

  if (Anrm != 0.0) {
    A[0] = A[0] * invAnrm;
    A[1] = A[1] * invAnrm;
    A[2] = A[2] * invAnrm;
  }

  return Anrm;
}

/*
  Find the sensitivity of the normalized 3D vector
*/
template <class ScalarType>
inline ScalarType normalize3DSens(ScalarType* sAnrm, ScalarType A[],
                                  ScalarType sA[]) {
  ScalarType Anrm = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
  ScalarType invAnrm = 1.0 / Anrm;
  *sAnrm = (1.0 / Anrm) * (A[0] * sA[0] + A[1] * sA[1] + A[2] * sA[2]);

  if (Anrm != 0.0) {
    A[0] = A[0] * invAnrm;
    A[1] = A[1] * invAnrm;
    A[2] = A[2] * invAnrm;

    sA[0] = (sA[0] - (A[0]) * (*sAnrm)) * invAnrm;
    sA[1] = (sA[1] - (A[1]) * (*sAnrm)) * invAnrm;
    sA[2] = (sA[2] - (A[2]) * (*sAnrm)) * invAnrm;
  } else {
    sA[0] = sA[1] = sA[2] = ScalarType(0.0);
  }

  return Anrm;
}

/*
  Find the sensitivity of the normalized 3D vector - code produced by tapenade
*/
template <class ScalarType>
inline void normalize3DSens_d(ScalarType A[], ScalarType sA[], ScalarType Ad[],
                              ScalarType sAd[]) {
  ScalarType Anrm;
  ScalarType Anrmd;
  ScalarType arg1;
  ScalarType arg1d;
  arg1d = Ad[0] * A[0] + A[0] * Ad[0] + Ad[1] * A[1] + A[1] * Ad[1] +
          Ad[2] * A[2] + A[2] * Ad[2];
  arg1 = A[0] * A[0] + A[1] * A[1] + A[2] * A[2];

  if (arg1 == 0.0) {
    Anrmd = 0.0;
  } else {
    Anrmd = arg1d / (2.0 * sqrt(arg1));
  }

  Anrm = sqrt(arg1);
  ScalarType sAnrm = 1.0 / Anrm * (A[0] * sA[0] + A[1] * sA[1] + A[2] * sA[2]);
  ScalarType sAnrmd =
      (Ad[0] * sA[0] + A[0] * sAd[0] + Ad[1] * sA[1] + A[1] * sAd[1] +
       Ad[2] * sA[2] + A[2] * sAd[2]) /
          Anrm -
      Anrmd * (A[0] * sA[0] + A[1] * sA[1] + A[2] * sA[2]) / (Anrm * Anrm);

  Ad[0] = (Ad[0] * Anrm - A[0] * Anrmd) / (Anrm * Anrm);
  A[0] = A[0] / Anrm;
  Ad[1] = (Ad[1] * Anrm - A[1] * Anrmd) / (Anrm * Anrm);
  A[1] = A[1] / Anrm;
  Ad[2] = (Ad[2] * Anrm - A[2] * Anrmd) / (Anrm * Anrm);
  A[2] = A[2] / Anrm;
  sAd[0] = (sAd[0] * Anrm - sA[0] * Anrmd) / (Anrm * Anrm) -
           (Ad[0] * Anrm - A[0] * Anrmd) * sAnrm / (Anrm * Anrm) -
           A[0] * sAnrmd / Anrm;
  sA[0] = sA[0] / Anrm - A[0] / Anrm * sAnrm;
  sAd[1] = (sAd[1] * Anrm - sA[1] * Anrmd) / (Anrm * Anrm) -
           (Ad[1] * Anrm - A[1] * Anrmd) * sAnrm / (Anrm * Anrm) -
           A[1] * sAnrmd / Anrm;
  sA[1] = sA[1] / Anrm - A[1] / Anrm * sAnrm;
  sAd[2] = (sAd[2] * Anrm - sA[2] * Anrmd) / (Anrm * Anrm) -
           (Ad[2] * Anrm - A[2] * Anrmd) * sAnrm / (Anrm * Anrm) -
           A[2] * sAnrmd / Anrm;
  sA[2] = sA[2] / Anrm - A[2] / Anrm * sAnrm;
}

/*
  Some handy functions for computing the transformation of stress and strains
*/

template <class ScalarType>
inline void transform3DStress(ScalarType newStress[], const ScalarType stress[],
                              const ScalarType transfrom[]);

template <class ScalarType>
inline void invTransform3DStress(ScalarType output[], const ScalarType input[],
                                 const ScalarType trans[]);

template <class ScalarType>
inline void transform3DStrain(ScalarType newStrain[], const ScalarType strain[],
                              const ScalarType transform[]);

template <class ScalarType>
inline void C3transform(ScalarType transform[], const ScalarType init[],
                        const ScalarType cosTheta, const ScalarType sinTheta);

/*
  The following are some helpful functions for rotating tensors
*/

template <class ScalarType>
inline void transform3DStress(ScalarType output[], const ScalarType input[],
                              const ScalarType trans[]) {
  /*
    Transform the components of the symmetric input tensor to a different
    reference axis defined by the array transform.

    input: a second order 3x3 tensor
    transform: a 3x3 transformation matrix
    ouput: the 3x3 transformed tensor in the transformed coordinate system

    the tensor is symmetric:
    [ input[0], input[5], input[4] ]
    [ input[5], input[1], input[3] ]
    [ input[4], input[3], input[2] ]

    [ trans[0], trans[1], trans[2] ]
    [ trans[3], trans[4], trans[5] ]
    [ trans[6], trans[7], trans[8] ]

    as = trans * input;
    output = as * [ trans ]^{T}
  */

  ScalarType as[9];

  as[0] = trans[0] * input[0] + trans[1] * input[5] + trans[2] * input[4];
  as[3] = trans[3] * input[0] + trans[4] * input[5] + trans[5] * input[4];
  as[6] = trans[6] * input[0] + trans[7] * input[5] + trans[8] * input[4];

  as[1] = trans[0] * input[5] + trans[1] * input[1] + trans[2] * input[3];
  as[4] = trans[3] * input[5] + trans[4] * input[1] + trans[5] * input[3];
  as[7] = trans[6] * input[5] + trans[7] * input[1] + trans[8] * input[3];

  as[2] = trans[0] * input[4] + trans[1] * input[3] + trans[2] * input[2];
  as[5] = trans[3] * input[4] + trans[4] * input[3] + trans[5] * input[2];
  as[8] = trans[6] * input[4] + trans[7] * input[3] + trans[8] * input[2];

  output[0] = as[0] * trans[0] + as[1] * trans[1] + as[2] * trans[2];
  output[1] = as[3] * trans[3] + as[4] * trans[4] + as[5] * trans[5];
  output[2] = as[6] * trans[6] + as[7] * trans[7] + as[8] * trans[8];

  output[3] = as[3] * trans[6] + as[4] * trans[7] + as[5] * trans[8];
  output[4] = as[0] * trans[6] + as[1] * trans[7] + as[2] * trans[8];
  output[5] = as[0] * trans[3] + as[1] * trans[4] + as[2] * trans[5];
}

template <class ScalarType>
inline void invTransform3DStress(ScalarType output[], const ScalarType input[],
                                 const ScalarType trans[]) {
  /*
    Perform the inverse of the transformation provided by using the
    transpose of the transformation instead of the actual transformation
    matrix. This can be used to rotate a second order symmetric tensor
    (such as stress) back to an original reference frame.

    input: a second order 3x3 tensor
    transform: a 3x3 transformation matrix
    ouput: the 3x3 transformed tensor in the transformed coordinate system

    the tensor is symmetric:
    [ input[0], input[5], input[4] ]
    [ input[5], input[1], input[3] ]
    [ input[4], input[3], input[2] ]

    [ trans[0], trans[1], trans[2] ]
    [ trans[3], trans[4], trans[5] ]
    [ trans[6], trans[7], trans[8] ]

    as = trans^{T} * input;
    output = as * [ trans ]
  */

  ScalarType as[9];

  as[0] = trans[0] * input[0] + trans[3] * input[5] + trans[6] * input[4];
  as[3] = trans[1] * input[0] + trans[4] * input[5] + trans[7] * input[4];
  as[6] = trans[2] * input[0] + trans[5] * input[5] + trans[8] * input[4];

  as[1] = trans[0] * input[5] + trans[3] * input[1] + trans[6] * input[3];
  as[4] = trans[1] * input[5] + trans[4] * input[1] + trans[7] * input[3];
  as[7] = trans[2] * input[5] + trans[5] * input[1] + trans[8] * input[3];

  as[2] = trans[0] * input[4] + trans[3] * input[3] + trans[6] * input[2];
  as[5] = trans[1] * input[4] + trans[4] * input[3] + trans[7] * input[2];
  as[8] = trans[2] * input[4] + trans[5] * input[3] + trans[8] * input[2];

  output[0] = as[0] * trans[0] + as[1] * trans[3] + as[2] * trans[6];
  output[1] = as[3] * trans[1] + as[4] * trans[4] + as[5] * trans[7];
  output[2] = as[6] * trans[2] + as[7] * trans[5] + as[8] * trans[8];

  output[3] = as[3] * trans[2] + as[4] * trans[5] + as[5] * trans[8];
  output[4] = as[0] * trans[2] + as[1] * trans[5] + as[2] * trans[8];
  output[5] = as[0] * trans[1] + as[1] * trans[4] + as[2] * trans[7];
}

template <class ScalarType>
inline void transform3DStrain(ScalarType output[], const ScalarType input[],
                              const ScalarType trans[]) {
  /*
    input: a second order 3x3 tensor
    transform: a 3x3 transformation matrix
    ouput: the 3x3 transformed tensor in the transformed coordinate system

    the tensor is symmetric:
    [ input[0],     0.5*input[5], 0.5*input[4] ]
    [ 0.5*input[5],     input[1], 0.5*input[3] ]
    [ 0.5*input[4], 0.5*input[3], input[2] ]

    [ trans[0], trans[1], trans[2] ]
    [ trans[3], trans[4], trans[5] ]
    [ trans[6], trans[7], trans[8] ]

    as = trans * input;
    output = as * [ trans ]^{T}
  */

  ScalarType as[9];

  as[0] = trans[0] * input[0] + 0.5 * trans[1] * input[5] +
          0.5 * trans[2] * input[4];
  as[3] = trans[3] * input[0] + 0.5 * trans[4] * input[5] +
          0.5 * trans[5] * input[4];
  as[6] = trans[6] * input[0] + 0.5 * trans[7] * input[5] +
          0.5 * trans[8] * input[4];

  as[1] = 0.5 * trans[0] * input[5] + trans[1] * input[1] +
          0.5 * trans[2] * input[3];
  as[4] = 0.5 * trans[3] * input[5] + trans[4] * input[1] +
          0.5 * trans[5] * input[3];
  as[7] = 0.5 * trans[6] * input[5] + trans[7] * input[1] +
          0.5 * trans[8] * input[3];

  as[2] = 0.5 * trans[0] * input[4] + 0.5 * trans[1] * input[3] +
          trans[2] * input[2];
  as[5] = 0.5 * trans[3] * input[4] + 0.5 * trans[4] * input[3] +
          trans[5] * input[2];
  as[8] = 0.5 * trans[6] * input[4] + 0.5 * trans[7] * input[3] +
          trans[8] * input[2];

  output[0] = as[0] * trans[0] + as[1] * trans[1] + as[2] * trans[2];
  output[1] = as[3] * trans[3] + as[4] * trans[4] + as[5] * trans[5];
  output[2] = as[6] * trans[6] + as[7] * trans[7] + as[8] * trans[8];

  output[3] = 2.0 * (as[3] * trans[6] + as[4] * trans[7] + as[5] * trans[8]);
  output[4] = 2.0 * (as[0] * trans[6] + as[1] * trans[7] + as[2] * trans[8]);
  output[5] = 2.0 * (as[0] * trans[3] + as[1] * trans[4] + as[2] * trans[5]);
}

template <class ScalarType>
inline void transform3DStressSens(ScalarType output[], ScalarType doutput[],
                                  const ScalarType input[],
                                  const ScalarType dinput[],
                                  const ScalarType trans[],
                                  const ScalarType dtrans[]) {
  /*
    Transform the components of the symmetric input tensor to a different
    reference axis defined by the array transform.

    input: a second order 3x3 tensor
    transform: a 3x3 transformation matrix
    ouput: the 3x3 transformed tensor in the transformed coordinate system

    the tensor is symmetric:
    [ input[0], input[5], input[4] ]
    [ input[5], input[1], input[3] ]
    [ input[4], input[3], input[2] ]

    [ trans[0], trans[1], trans[2] ]
    [ trans[3], trans[4], trans[5] ]
    [ trans[6], trans[7], trans[8] ]

    as = trans * input;
    output = as * [ trans ]^{T}
  */

  ScalarType as[9], das[9];

  as[0] = trans[0] * input[0] + trans[1] * input[5] + trans[2] * input[4];
  as[3] = trans[3] * input[0] + trans[4] * input[5] + trans[5] * input[4];
  as[6] = trans[6] * input[0] + trans[7] * input[5] + trans[8] * input[4];

  as[1] = trans[0] * input[5] + trans[1] * input[1] + trans[2] * input[3];
  as[4] = trans[3] * input[5] + trans[4] * input[1] + trans[5] * input[3];
  as[7] = trans[6] * input[5] + trans[7] * input[1] + trans[8] * input[3];

  as[2] = trans[0] * input[4] + trans[1] * input[3] + trans[2] * input[2];
  as[5] = trans[3] * input[4] + trans[4] * input[3] + trans[5] * input[2];
  as[8] = trans[6] * input[4] + trans[7] * input[3] + trans[8] * input[2];

  das[0] = (trans[0] * dinput[0] + trans[1] * dinput[5] + trans[2] * dinput[4] +
            dtrans[0] * input[0] + dtrans[1] * input[5] + dtrans[2] * input[4]);
  das[3] = (trans[3] * dinput[0] + trans[4] * dinput[5] + trans[5] * dinput[4] +
            dtrans[3] * input[0] + dtrans[4] * input[5] + dtrans[5] * input[4]);
  das[6] = (trans[6] * dinput[0] + trans[7] * dinput[5] + trans[8] * dinput[4] +
            dtrans[6] * input[0] + dtrans[7] * input[5] + dtrans[8] * input[4]);

  das[1] = (trans[0] * dinput[5] + trans[1] * dinput[1] + trans[2] * dinput[3] +
            dtrans[0] * input[5] + dtrans[1] * input[1] + dtrans[2] * input[3]);
  das[4] = (trans[3] * dinput[5] + trans[4] * dinput[1] + trans[5] * dinput[3] +
            dtrans[3] * input[5] + dtrans[4] * input[1] + dtrans[5] * input[3]);
  das[7] = (trans[6] * dinput[5] + trans[7] * dinput[1] + trans[8] * dinput[3] +
            dtrans[6] * input[5] + dtrans[7] * input[1] + dtrans[8] * input[3]);

  das[2] = (trans[0] * dinput[4] + trans[1] * dinput[3] + trans[2] * dinput[2] +
            dtrans[0] * input[4] + dtrans[1] * input[3] + dtrans[2] * input[2]);
  das[5] = (trans[3] * dinput[4] + trans[4] * dinput[3] + trans[5] * dinput[2] +
            dtrans[3] * input[4] + dtrans[4] * input[3] + dtrans[5] * input[2]);
  das[8] = (trans[6] * dinput[4] + trans[7] * dinput[3] + trans[8] * dinput[2] +
            dtrans[6] * input[4] + dtrans[7] * input[3] + dtrans[8] * input[2]);

  output[0] = as[0] * trans[0] + as[1] * trans[1] + as[2] * trans[2];
  output[1] = as[3] * trans[3] + as[4] * trans[4] + as[5] * trans[5];
  output[2] = as[6] * trans[6] + as[7] * trans[7] + as[8] * trans[8];

  output[3] = as[3] * trans[6] + as[4] * trans[7] + as[5] * trans[8];
  output[4] = as[0] * trans[6] + as[1] * trans[7] + as[2] * trans[8];
  output[5] = as[0] * trans[3] + as[1] * trans[4] + as[2] * trans[5];

  doutput[0] = (as[0] * dtrans[0] + as[1] * dtrans[1] + as[2] * dtrans[2] +
                das[0] * trans[0] + das[1] * trans[1] + das[2] * trans[2]);
  doutput[1] = (as[3] * dtrans[3] + as[4] * dtrans[4] + as[5] * dtrans[5] +
                das[3] * trans[3] + das[4] * trans[4] + das[5] * trans[5]);
  doutput[2] = (as[6] * dtrans[6] + as[7] * dtrans[7] + as[8] * dtrans[8] +
                das[6] * trans[6] + das[7] * trans[7] + das[8] * trans[8]);

  doutput[3] = (as[3] * dtrans[6] + as[4] * dtrans[7] + as[5] * dtrans[8] +
                das[3] * trans[6] + das[4] * trans[7] + das[5] * trans[8]);
  doutput[4] = (as[0] * dtrans[6] + as[1] * dtrans[7] + as[2] * dtrans[8] +
                das[0] * trans[6] + das[1] * trans[7] + das[2] * trans[8]);
  doutput[5] = (as[0] * dtrans[3] + as[1] * dtrans[4] + as[2] * dtrans[5] +
                das[0] * trans[3] + das[1] * trans[4] + das[2] * trans[5]);
}

template <class ScalarType>
inline void C3transform(ScalarType transform[], const ScalarType init[],
                        const ScalarType cosTheta, const ScalarType sinTheta) {
  transform[0] = init[0] * cosTheta + init[3] * sinTheta;
  transform[1] = init[1] * cosTheta + init[4] * sinTheta;
  transform[2] = init[2] * cosTheta + init[5] * sinTheta;

  transform[3] = -init[0] * sinTheta + init[3] * cosTheta;
  transform[4] = -init[1] * sinTheta + init[4] * cosTheta;
  transform[5] = -init[2] * sinTheta + init[5] * cosTheta;

  transform[6] = init[6];
  transform[7] = init[7];
  transform[8] = init[8];
}

template <class ScalarType>
inline void d_transform3D(ScalarType input[], ScalarType output[],
                          ScalarType trans[], ScalarType dtrans[]) {
  ScalarType as[9];
  ScalarType das[9];

  as[0] = trans[0] * input[0] + trans[1] * input[5] + trans[2] * input[4];
  as[3] = trans[3] * input[0] + trans[4] * input[5] + trans[5] * input[4];
  as[6] = trans[6] * input[0] + trans[7] * input[5] + trans[8] * input[4];

  as[1] = trans[0] * input[5] + trans[1] * input[1] + trans[2] * input[3];
  as[4] = trans[3] * input[5] + trans[4] * input[1] + trans[5] * input[3];
  as[7] = trans[6] * input[5] + trans[7] * input[1] + trans[8] * input[3];

  as[2] = trans[0] * input[4] + trans[1] * input[3] + trans[2] * input[2];
  as[5] = trans[3] * input[4] + trans[4] * input[3] + trans[5] * input[2];
  as[8] = trans[6] * input[4] + trans[7] * input[3] + trans[8] * input[2];

  das[0] = dtrans[0] * input[0] + dtrans[1] * input[5] + dtrans[2] * input[4];
  das[3] = dtrans[3] * input[0] + dtrans[4] * input[5] + dtrans[5] * input[4];
  das[6] = dtrans[6] * input[0] + dtrans[7] * input[5] + dtrans[8] * input[4];

  das[1] = dtrans[0] * input[5] + dtrans[1] * input[1] + dtrans[2] * input[3];
  das[4] = dtrans[3] * input[5] + dtrans[4] * input[1] + dtrans[5] * input[3];
  das[7] = dtrans[6] * input[5] + dtrans[7] * input[1] + dtrans[8] * input[3];

  das[2] = dtrans[0] * input[4] + dtrans[1] * input[3] + dtrans[2] * input[2];
  das[5] = dtrans[3] * input[4] + dtrans[4] * input[3] + dtrans[5] * input[2];
  das[8] = dtrans[6] * input[4] + dtrans[7] * input[3] + dtrans[8] * input[2];

  output[0] = as[0] * dtrans[0] + as[1] * dtrans[1] + as[2] * dtrans[2] +
              das[0] * trans[0] + das[1] * trans[1] + das[2] * trans[2];

  output[1] = as[3] * dtrans[3] + as[4] * dtrans[4] + as[5] * dtrans[5] +
              das[3] * trans[3] + das[4] * trans[4] + das[5] * trans[5];

  output[2] = as[6] * dtrans[6] + as[7] * dtrans[7] + as[8] * dtrans[8] +
              das[6] * trans[6] + das[7] * trans[7] + das[8] * trans[8];

  output[3] = as[3] * dtrans[6] + as[4] * dtrans[7] + as[5] * dtrans[8] +
              das[3] * trans[6] + das[4] * trans[7] + das[5] * trans[8];

  output[4] = as[0] * dtrans[6] + as[1] * dtrans[7] + as[2] * dtrans[8] +
              das[0] * trans[6] + das[1] * trans[7] + das[2] * trans[8];

  output[5] = as[0] * dtrans[3] + as[1] * dtrans[4] + as[2] * dtrans[5] +
              das[0] * trans[3] + das[1] * trans[4] + das[2] * trans[5];
}

TACS_END_NAMESPACE

#endif
