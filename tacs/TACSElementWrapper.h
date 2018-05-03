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

#ifndef TACS_ELEMENT_WRAPPER_H
#define TACS_ELEMENT_WRAPPER_H

/*
  The following file contains wrappers for element objects.

  These C++ classes are used to create python-level implementations of
  the underlying TACS constitutive objects. They are designed to use
  callbacks through the python layer. 
  
  Not much error checking is performed here, so beware.
*/

#include "TACSElement.h"
#include "Python.h"

/*
  The following class implements a basic wrapper for the
  PlaneStressStiffness type of constitutive object 
*/
class TACSElementWrapper : public TACSElement {
 public:
  TACSElementWrapper( PyObject *_self_ptr, 
                      int _num_nodes,
                      int _num_displacements ){ 
    self_ptr = _self_ptr;
    num_nodes = _num_nodes;
    num_displacements = _num_displacements;

    // This causes a circular reference so the object is never
    // deleted. This should be fixed properly using weak references,
    // but I'm not 100% sure how to do this yet...
    Py_INCREF(self_ptr);
    
    getinitconditions = NULL;
    addresidual = NULL;
    addjacobian = NULL;
    addadjresproduct = NULL;

  }
  ~TACSElementWrapper() {
    Py_DECREF(self_ptr);
  }

  // TACS Element member variables
  // -----------------------------
  int num_nodes;
  int num_displacements;

  // TACS Element member functions
  // -----------------------------
  int numDisplacements() {
    return num_displacements;
  }

  int numNodes() {
    return num_nodes;
  }

  // Retrieve the initial conditions and add the derivative
  // ------------------------------------------------------
  void getInitConditions( TacsScalar vars[],
                          TacsScalar dvars[],
                          TacsScalar ddvars[],
                          const TacsScalar Xpts[] ){
    memset(vars, 0, numVariables()*sizeof(TacsScalar));
    memset(dvars, 0, numVariables()*sizeof(TacsScalar));
    memset(ddvars, 0, numVariables()*sizeof(TacsScalar));
    if (self_ptr && getinitconditions) {
      int nvars = num_nodes*num_displacements;
      getinitconditions(self_ptr, nvars, num_nodes, 
                        vars, dvars, ddvars, Xpts);
    }
  }
  
  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ) {
    if (self_ptr && addresidual) {
      int nvars = num_nodes*num_displacements;
      addresidual(self_ptr, nvars, num_nodes, time, 
                  res, Xpts, vars, dvars, ddvars);
    }
  }

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian( double time, TacsScalar J[],
                    double alpha, double beta, double gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ) {
    if (self_ptr && addjacobian) {
      int nvars = num_nodes*num_displacements;
      addjacobian(self_ptr, nvars, num_nodes, time, 
                  J, alpha, beta, gamma, 
                  Xpts, vars, dvars, ddvars);
    }
  }

  // Add the product of the adjoint variables with the derivative of the residual
  // ----------------------------------------------------------------------------
  void addAdjResProduct( double time, double scale,
                         TacsScalar dvSens[], int dvLen,
                         const TacsScalar psi[],
                         const TacsScalar Xpts[],
                         const TacsScalar vars[],
                         const TacsScalar dvars[],
                         const TacsScalar ddvars[] ) {
    if (self_ptr && addadjresproduct) {
      int nvars = num_nodes*num_displacements;
      addadjresproduct(self_ptr, nvars, num_nodes, time, scale, dvSens, dvLen, 
                       psi, Xpts, vars, dvars, ddvars);
    }
  }

  // Define the object name 
  // ----------------------
  const char * elementName(){ 
    return "TACSElementWrapper";
  }

  // Function pointers
  // -----------------
  PyObject *self_ptr; // Pointer to the python object
  void (*getinitconditions)( void *, int, int, 
                             TacsScalar *, TacsScalar *, 
                             TacsScalar *, const TacsScalar * );
  
  void (*addresidual)( void *, int, int, double time, TacsScalar res[],
                       const TacsScalar Xpts[],
                       const TacsScalar vars[],
                       const TacsScalar dvars[],
                       const TacsScalar ddvars[] );
  
  void (*addjacobian)( void *, int, int, double time, TacsScalar J[],
                       double alpha, double beta, double gamma,
                       const TacsScalar Xpts[],
                       const TacsScalar vars[],
                       const TacsScalar dvars[],
                       const TacsScalar ddvars[] );

  void (*addadjresproduct)( void*, int, int, double, TacsScalar, TacsScalar *,
                            int, const TacsScalar *, const TacsScalar *, 
                            const TacsScalar *, const TacsScalar *, 
                            const TacsScalar * );
};

#endif
