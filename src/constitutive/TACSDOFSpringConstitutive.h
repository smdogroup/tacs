#ifndef TACS_DOF_SPRING_CONSTITUTIVE_H
#define TACS_DOF_SPRING_CONSTITUTIVE_H

/*
  The stiffness object for the variety of springs

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Copyright (C) 2020 Aerion Technologies Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.
*/

#include "TACSGeneralSpringConstitutive.h"

/**
  This is the base class for the traditional spring constitutive objects with no
  dof coupling. Assumes 6 dofs.
*/
class TACSDOFSpringConstitutive : public TACSGeneralSpringConstitutive {
 public:
  TACSDOFSpringConstitutive(TacsScalar _k[]);

  const char* constitutiveName();

 private:
  static const char* constName;
};

#endif  // TACS_DOF_SPRING_CONSTITUTIVE_H
