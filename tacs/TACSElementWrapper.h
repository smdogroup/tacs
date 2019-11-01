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
                      int _vars_per_node,
                      int _num_nodes ){
    self_ptr = _self_ptr;
    vars_per_node = _vars_per_node;
    num_nodes = _num_nodes;

    // This causes a circular reference so the object is never
    // deleted. This should be fixed properly using weak references,
    // but I'm not 100% sure how to do this yet...
    Py_INCREF(self_ptr);
    
    getinitconditions = NULL;
    addresidual = NULL;
    addjacobian = NULL;

  }
  ~TACSElementWrapper() {
    Py_DECREF(self_ptr);
  }

  // TACS Element member variables
  // -----------------------------
  int vars_per_node;
  int num_nodes;

  // TACS Element member functions
  // -----------------------------
  int getVarsPerNode() {
    return vars_per_node;
  }

  int getNumNodes() {
    return num_nodes;
  }

  // Retrieve the initial conditions and add the derivative
  // ------------------------------------------------------
  void getInitConditions( int elem_index, const TacsScalar Xpts[],
                          TacsScalar vars[], TacsScalar dvars[],
                          TacsScalar ddvars[] ){
    memset(vars, 0, getNumVariables()*sizeof(TacsScalar));
    memset(dvars, 0, getNumVariables()*sizeof(TacsScalar));
    memset(ddvars, 0, getNumVariables()*sizeof(TacsScalar));
    if (self_ptr && getinitconditions) {
      int num_vars = num_nodes*vars_per_node;
      getinitconditions(self_ptr, elem_index, num_nodes, Xpts, num_vars,
                        vars, dvars, ddvars);
    }
  }
  
  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual( int elem_index, double time,
                    const TacsScalar Xpts[], const TacsScalar vars[],
                    const TacsScalar dvars[], const TacsScalar ddvars[],
                    TacsScalar res[] ) {
    if (self_ptr && addresidual){
      int num_vars = num_nodes*vars_per_node;
      addresidual(self_ptr, elem_index, time, num_nodes, Xpts,
                  num_vars, vars, dvars, ddvars, res);
    }
  }

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian( int elem_index, double time,
                    TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                    const TacsScalar Xpts[], const TacsScalar vars[],
                    const TacsScalar dvars[], const TacsScalar ddvars[],
                    TacsScalar res[], TacsScalar mat[] ) {
    if (self_ptr && addjacobian){
      int num_vars = num_nodes*vars_per_node;
      addjacobian(self_ptr, elem_index, time,
                  alpha, beta, gamma, num_nodes, Xpts,
                  num_vars, vars, dvars, ddvars, res, mat);
    }
  }

  // Define the object name 
  // ----------------------
  const char *getObjectName(){
    return "TACSElementWrapper";
  }

  // Function pointers
  // -----------------
  PyObject *self_ptr; // Pointer to the python object
  void (*getinitconditions)(void*, int, int, const TacsScalar*, int,
                            TacsScalar*, TacsScalar*, TacsScalar*);
  void (*addresidual)(void*, int, double, int, const TacsScalar*,
                      int, const TacsScalar*, const TacsScalar*,
                      const TacsScalar*, TacsScalar*);
  void (*addjacobian)(void*, int, double, TacsScalar, TacsScalar, TacsScalar,
                      int, const TacsScalar*, int, const TacsScalar*,
                      const TacsScalar*, const TacsScalar*,
                      TacsScalar*, TacsScalar*);
};

#endif
