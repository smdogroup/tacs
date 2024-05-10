/*
========================================================================
Gaussian Process Model for TACS Buckling Constraints
========================================================================
@File   :   GaussianProcessModel.h
@Date   :   2024/05/10
@Author :   Sean Phillip Engelstad
@Description : Use Gaussian Processes for machine learning techniques to interpolate
learned buckling constraints on a training dataset to new test data points. This
approach is implemented into the TACSGPBladeStiffenedShellConstitutive class in TACS
for more physically accurate buckling constraints of stiffened panels.
*/

#pragma once

// =============================================================================
// Extension Includes
// =============================================================================
#include "TacsUtilities.h"

// =============================================================================
// Class Declaration
// =============================================================================

class GaussianProcessModel {
public:
 GaussianProcessModel();
 ~GaussianProcessModel();

private:
 void loadAxialData(filename);
 void loadShearData(filename);
}