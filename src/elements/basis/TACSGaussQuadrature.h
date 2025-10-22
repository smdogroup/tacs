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

#ifndef TACS_GAUSS_QUADRATURE_H
#define TACS_GAUSS_QUADRATURE_H

/*
  The following are the definitions of the Gauss (or Gauss-Legendre)
  quadrature points and weights for the interval [-1, 1]. Note that
  these schemes are exact for polynomials of degree 2n - 1.
*/
const double TacsGaussQuadPts1[] = {0.0};
const double TacsGaussQuadWts1[] = {2.0};

const double TacsGaussQuadPts2[] = {-0.577350269189626, 0.577350269189626};
const double TacsGaussQuadWts2[] = {1.0, 1.0};

const double TacsGaussQuadPts3[] = {-0.774596669241483, 0.0, 0.774596669241483};
const double TacsGaussQuadWts3[] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

const double TacsGaussQuadPts4[] = {-0.861136311594053, -0.339981043584856,
                                    0.339981043584856, 0.861136311594053};
const double TacsGaussQuadWts4[] = {0.347854845137454, 0.652145154862546,
                                    0.652145154862546, 0.347854845137454};

const double TacsGaussQuadPts5[] = {-0.906179845938664, -0.538469310105683, 0.0,
                                    0.538469310105683, 0.906179845938664};
const double TacsGaussQuadWts5[] = {0.236926885056189, 0.478628670499366,
                                    0.568888888888889, 0.478628670499366,
                                    0.236926885056189};

const double TacsGaussQuadPts6[] = {
    -0.9324695142031520278123016, -0.6612093864662645136613996,
    -0.2386191860831969086305017, 0.2386191860831969086305017,
    0.6612093864662645136613996,  0.9324695142031520278123016};
const double TacsGaussQuadWts6[] = {
    0.1713244923791703450402961, 0.3607615730481386075698335,
    0.4679139345726910473898703, 0.4679139345726910473898703,
    0.3607615730481386075698335, 0.1713244923791703450402961};

const double TacsGaussQuadPts7[] = {
    -0.9491079123427585245261897, -0.7415311855993944398638648,
    -0.4058451513773971669066064, 0.0,
    0.4058451513773971669066064,  0.7415311855993944398638648,
    0.9491079123427585245261897};
const double TacsGaussQuadWts7[] = {
    0.1294849661688696932706114, 0.2797053914892766679014678,
    0.3818300505051189449503698, 0.4179591836734693877551020,
    0.3818300505051189449503698, 0.2797053914892766679014678,
    0.1294849661688696932706114};

const double TacsGaussQuadPts8[] = {
    -0.9602898564975362316835609, -0.7966664774136267395915539,
    -0.5255324099163289858177390, -0.1834346424956498049394761,
    0.1834346424956498049394761,  0.5255324099163289858177390,
    0.7966664774136267395915539,  0.9602898564975362316835609};
const double TacsGaussQuadWts8[] = {
    0.1012285362903762591525314, 0.2223810344533744705443560,
    0.3137066458778872873379622, 0.3626837833783619829651504,
    0.3626837833783619829651504, 0.3137066458778872873379622,
    0.2223810344533744705443560, 0.1012285362903762591525314};

/*
  The following are the Gauss--Lobatto integration schemes. These
  quadrature schemes always include the starting/end points for the
  interval. These schemes are less accurate in that they integrate
  polynomials of degree 2n-3 exactly (compared to Gauss quadrature
  schemes which integrated 2n-1 exactly).
*/
const double TacsGaussLobattoPts2[] = {-1.0, 1.0};
const double TacsGaussLobattoWts2[] = {1.0, 1.0};

const double TacsGaussLobattoPts3[] = {-1.0, 0.0, 1.0};
const double TacsGaussLobattoWts3[] = {1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0};

const double TacsGaussLobattoPts4[] = {-1.0, -0.44721359549995793,
                                       0.44721359549995793, 1.0};
const double TacsGaussLobattoWts4[] = {1.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0,
                                       1.0 / 6.0};

const double TacsGaussLobattoPts5[] = {-1.0, -0.65465367070797709, 0.0,
                                       0.65465367070797709, 1.0};
const double TacsGaussLobattoWts5[] = {1.0 / 10.0, 49.0 / 90.0, 32.0 / 45.0,
                                       49.0 / 90.0, 1.0 / 10.0};

const double TacsGaussLobattoPts6[] = {-1.0,
                                       -0.76505532392946474,
                                       -0.2852315164806451,
                                       0.2852315164806451,
                                       0.76505532392946474,
                                       1.0};
const double TacsGaussLobattoWts6[] = {1.0 / 15.0,          0.378474956297847,
                                       0.55485837703548635, 0.55485837703548635,
                                       0.378474956297847,   1.0 / 15.0};

#endif  // TACS_GAUSS_QUADRATURE_H
