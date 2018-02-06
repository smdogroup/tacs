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

#ifndef TACS_CONSTITUTIVE_WRAPPER_H
#define TACS_CONSTITUTIVE_WRAPPER_H

/*
  The following file contains wrappers for constitutive objects.

  These C++ classes are used to create python-level implementations of
  the underlying TACS constitutive objects. They are designed to use
  callbacks through the python layer. 
  
  Not much error checking is performed here, so beware.
*/

#include "PlaneStressStiffness.h"
#include "FSDTStiffness.h"
#include "Python.h"

/*
  The following class implements a basic wrapper for the
  PlaneStressStiffness type of constitutive object 
*/
class PSStiffnessWrapper : public PlaneStressStiffness {
 public:
  PSStiffnessWrapper( PyObject *_self_ptr ){
    self_ptr = _self_ptr;

    // This causes a circular reference so the object is never
    // deleted. This should be fixed properly using weak references,
    // but I'm not 100% sure how to do this yet...
    Py_INCREF(self_ptr);

    setdesignvars = NULL;
    getdesignvars = NULL;
    getdesignvarrange = NULL;
    calculatestress = NULL;
    addstressdvsens = NULL;
    getpointwisemass = NULL;
    addpointwisemassdvsens = NULL;
    fail = NULL;
    failstrainsens = NULL;
    addfaildvsens = NULL;
  }
  ~PSStiffnessWrapper(){
    Py_DECREF(self_ptr);
  }

  // Define the object name 
  // ----------------------
  const char * constitutiveName(){ 
    return "PSStiffnessWrapper";
  }

  // Function pointers
  // -----------------
  PyObject *self_ptr; // Pointer to the python object
  void (*setdesignvars)( void *, const TacsScalar *, int );
  void (*getdesignvars)( void *, TacsScalar *, int );
  void (*getdesignvarrange)( void *, TacsScalar *, TacsScalar *, int );
  void (*calculatestress)( void *, const double *, 
                           const TacsScalar *, TacsScalar* );
  void (*addstressdvsens)( void *, const double *, const TacsScalar *,
                           TacsScalar, const TacsScalar *,
                           TacsScalar *, int );
  void (*getpointwisemass)( void *, const double *, TacsScalar * );
  void (*addpointwisemassdvsens)( void *, const double *,
                                  const TacsScalar *, TacsScalar *, int );
  TacsScalar (*fail)( void *, const double *, const TacsScalar * );
  void (*failstrainsens)( void *, const double *,
                          const TacsScalar *, TacsScalar * );
  void (*addfaildvsens)( void *, const double *, const TacsScalar *, 
                         TacsScalar, TacsScalar *, int );

  // Design variable access functions
  // --------------------------------
  void setDesignVars( const TacsScalar *x, int dvLen ){
    if (self_ptr && setdesignvars){
      setdesignvars(self_ptr, x, dvLen);
    }
  }
  void getDesignVars( TacsScalar *x, int dvLen ){
    if (self_ptr && getdesignvars){
      getdesignvars(self_ptr, x, dvLen);
    }
  }
  void getDesignVarRange( TacsScalar *lb, TacsScalar *ub, int dvLen ){
    if (self_ptr && getdesignvarrange){
      getdesignvarrange(self_ptr, lb, ub, dvLen);
    }
  }

  // Stress member functions
  // -----------------------
  void calculateStress( const double pt[], 
                        const TacsScalar strain[], 
                        TacsScalar stress[] ){
    if (self_ptr && calculatestress){
      calculatestress(self_ptr, pt, strain, stress);
    }
  }
  void addStressDVSens( const double pt[], const TacsScalar strain[], 
                        TacsScalar alpha, const TacsScalar psi[], 
                        TacsScalar dvSens[], int dvLen ){
    if (self_ptr && addstressdvsens){
      addstressdvsens(self_ptr, pt, strain, alpha, psi, dvSens, dvLen);
    }
  }

  // Mass moment member functions
  // ----------------------------
  void getPointwiseMass( const double pt[], 
                         TacsScalar mass[] ){
    if (self_ptr && getpointwisemass){
      getpointwisemass(self_ptr, pt, mass);
    }
  }
  void addPointwiseMassDVSens( const double pt[], 
                               const TacsScalar alpha[],
                               TacsScalar dvSens[], int dvLen ){
    if (self_ptr && addpointwisemassdvsens){
      addpointwisemassdvsens(self_ptr, pt, alpha, dvSens, dvLen);
    }
  }

  // Evaluate the failure functions
  // ------------------------------
  void failure( const double pt[], const TacsScalar strain[],
                TacsScalar *fval ){
    if (self_ptr && fail){
      *fval = fail(self_ptr, pt, strain);
    }
  }
  void failureStrainSens( const double pt[], 
                          const TacsScalar strain[],
                          TacsScalar sens[] ){
    if (self_ptr && failstrainsens){
      failstrainsens(self_ptr, pt, strain, sens);
    }
  }
  void addFailureDVSens( const double pt[], const TacsScalar strain[],
                         TacsScalar alpha, TacsScalar dvSens[], int dvLen ){
    if (self_ptr && addfaildvsens){
      addfaildvsens(self_ptr, pt, strain, alpha, dvSens, dvLen);
    }
  }
};

/*
  The following class implements basic functionality for the
  FSDTStiffness type of class 
*/
class FSDTStiffnessWrapper : public FSDTStiffness {
 public:
  FSDTStiffnessWrapper( PyObject *_self_ptr ){
    self_ptr = _self_ptr;
    Py_INCREF(self_ptr);
    setdesignvars = NULL;
    getdesignvars = NULL;
    getdesignvarrange = NULL;
    getstiffness = NULL;
    addstiffnessdvsens = NULL;
    getpointwisemass = NULL;
    addpointwisemassdvsens = NULL;
    fail = NULL;
    failstrainsens = NULL;
    addfaildvsens = NULL;
  }
  ~FSDTStiffnessWrapper(){
    Py_DECREF(self_ptr);
  }

  // Define the object name 
  // ----------------------
  const char * constitutiveName(){ 
    return "FSDTStiffnessWrapper";
  }

  // Function pointers
  // -----------------
  PyObject *self_ptr; // Pointer to the python object
  void (*setdesignvars)( void *, const TacsScalar *, int );
  void (*getdesignvars)( void *, TacsScalar *, int );
  void (*getdesignvarrange)( void *, TacsScalar *, TacsScalar *, int );
  TacsScalar (*getstiffness)( void *, const double *, 
                              TacsScalar*, TacsScalar*, 
                              TacsScalar*, TacsScalar* );
  void (*addstiffnessdvsens)( void *, const double *, 
                              const TacsScalar *, const TacsScalar *,
                              TacsScalar, TacsScalar *, int );
  void (*getpointwisemass)( void *, const double *, TacsScalar * );
  void (*addpointwisemassdvsens)( void *, const double *,
                                  const TacsScalar *, TacsScalar *, int );
  TacsScalar (*fail)( void *, const double *, const TacsScalar * );
  void (*failstrainsens)( void *, const double *,
                          const TacsScalar *, TacsScalar * );
  void (*addfaildvsens)( void *, const double *, const TacsScalar *, 
                         TacsScalar, TacsScalar *, int );

  // Design variable access functions
  // --------------------------------
  void setDesignVars( const TacsScalar *x, int dvLen ){
    if (self_ptr && setdesignvars){
      setdesignvars(self_ptr, x, dvLen);
    }
  }
  void getDesignVars( TacsScalar *x, int dvLen ){
    if (self_ptr && getdesignvars){
      getdesignvars(self_ptr, x, dvLen);
    }
  }
  void getDesignVarRange( TacsScalar *lb, TacsScalar *ub, int dvLen ){
    if (self_ptr && getdesignvarrange){
      getdesignvarrange(self_ptr, lb, ub, dvLen);
    }
  }

  // Stress member functions
  // -----------------------
  TacsScalar getStiffness( const double pt[], 
                           TacsScalar A[], TacsScalar B[],
                           TacsScalar D[], TacsScalar As[] ){
    if (self_ptr && getstiffness){
      return getstiffness(self_ptr, pt, A, B, D, As);
    }
    return 0.0;
  }
  void addStiffnessDVSens( const double pt[],
                           const TacsScalar e[], const TacsScalar psi[],
                           TacsScalar rotPsi, 
                           TacsScalar fdvSens[], int dvLen ){
    if (self_ptr && addstiffnessdvsens){
      addstiffnessdvsens(self_ptr, pt, e, psi, rotPsi,
                         fdvSens, dvLen);
    }
  }

  // Mass moment member functions
  // ----------------------------
  void getPointwiseMass( const double pt[], 
                         TacsScalar mass[] ){
    if (self_ptr && getpointwisemass){
      getpointwisemass(self_ptr, pt, mass);
    }
  }
  void addPointwiseMassDVSens( const double pt[], 
                               const TacsScalar alpha[],
                               TacsScalar dvSens[], int dvLen ){
    if (self_ptr && addpointwisemassdvsens){
      addpointwisemassdvsens(self_ptr, pt, alpha, dvSens, dvLen);
    }
  }

  // Evaluate the failure functions
  // ------------------------------
  void failure( const double pt[], const TacsScalar strain[],
                TacsScalar *fval ){
    if (self_ptr && fail){
      *fval = fail(self_ptr, pt, strain);
    }
  }
  void failureStrainSens( const double pt[], 
                          const TacsScalar strain[],
                          TacsScalar sens[] ){
    if (self_ptr && failstrainsens){
      failstrainsens(self_ptr, pt, strain, sens);
    }
  }
  void addFailureDVSens( const double pt[], const TacsScalar strain[],
                         TacsScalar alpha, TacsScalar dvSens[], int dvLen ){
    if (self_ptr && addfaildvsens){
      addfaildvsens(self_ptr, pt, strain, alpha, dvSens, dvLen);
    }
  }
};


#endif
