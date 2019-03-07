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

#ifndef TACS_ELEMENT_H
#define TACS_ELEMENT_H

/*
  Basic TACSElement definition

  The purpose of this file is to provide an interface for creating and
  storing different instances of the finite elements that will be used
  by TACS.  This is what should be extended when including more
  elements and not the underlying TACS implementation itself.
*/

#include "TACSObject.h"
#include "TACSConstitutive.h"

/*
  TACSElement is the base class from which all other elements inherit.
  The functions in this class are broken into a few different groups:

  Functions inherited from TACSOptObject:
  ---------------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );
  void getDesignVarRange( TacsScalar lowerBound[],
                          TacsScalar upperBound[], int numDVs );

  These functions are used to set and retrieve variable information
  and set up internal data for sensitivity calculations. Further details
  can be found in TACSObject.h

  Required information about the sizes of stresses/nodes/variables
  ----------------------------------------------------------------
  int numDisplacements();
  int numStresses();
  int numNodes();
  int numVariables();

  These functions return the size of the displacements (number of
  degrees of freedom at each node), the number of nodes, the number of
  stresses and the number of variables required for this element.

  Information required for output visualization
  ---------------------------------------------

  const char *elementName();
  const char *displacementName( int i );
  const char *stressName( int i );
  const char *strainName( int i );
  const char *extraName( int i );
  int numExtras();
  ElementType getElementType();

  The functions returning char* provide a name for the
  displacement/strain etc.  components that are useful for writing
  output files with named fields. numDisplacements() etc. provides the
  number of components in each field. getElementType() provides an
  enumerated type that idicates the underlying type of element (BEAM,
  SHELL etc.). Note that different element classes can return the same
  ElementType.

  Visualization is a critical component of any finite-element analysis.
  The following functions are used to retrieve the element-level data
  for post-analysis visualization.

  addOutputCount(): Add the number of nodes and number of csr entries
  that are requried for this element

  getOutputData(): Get the output data of a given type from the element.

  getOutputConnectivity(): Get the connectivity of the data

  Functions for analysis:
  -----------------------

  getInitCondition(): Return the initial conditions associated with
  the element.

  computeEnergies(): Compute the kinetic and potential energies
  of this element at the provided time

  addResidual(): Add the residual associated with this element.  This
  includes all time-dependent components of the residual.

  addJacobian(): Add the Jacobian of the governing equations
  provided by the addResidual() call.

  getMatType(): Return a time-independent element matrix of a given
  type as specified by the ElementMatrixType enumeration. This can be
  used to evaluate mass and geometric stiffness matrices. Note that
  not all element classes implement all the matrix types.

  Functions for sensitivity analysis:
  -----------------------------------

  There are two main types of derivatives that are handled within the
  TACSElement class: material design variables and geometric design
  variables. The underlying assumption in TACS is that the geometry
  variables only change the nodal coordinates. Thus, all geometric
  sensivity information can be expanded through a derivative of the
  element residuals w.r.t. the nodal coordinates.

  For efficiency reasons, the TACSElement code computes the derivative
  of the product of the adjoint variables and the element residuals.
  The derivative of this product can be computed more efficiently than
  computing the derivative with respect to each design variable
  or node in sequence. The sensitivity functions are:

  addInitConditionAdjResProduct(): Add the derivative of the product
  of the adjoint variables and the initial conditions with respect
  to material (or other) design variables

  addInitConditionAdjResXptProduct(): Compute the derivative of the
  product of the adjoint variables and the initial conditions with
  respect to node locations.

  addAdjResProduct(): Add the derivative of the product of the adjoint
  variables and the residual to the material design variable vector

  addAdjResXptProduct(): Add the derivative of the product of the
  adjoint variables and the adjoint with respect to the node locations

  addMatDVSensInnerProduct(): Compute the derivative of the inner
  product of the given matrix type with respect to the material design
  variables.

  getMatSVSensInnerProduct(): Compute the derivative of the inner
  product of the given matrix type with respect to the state
  variables.  (This is used when the matrix type depends on the state
  variables e.g. the geometric stiffness matrix)

  Post-processing functions:
  --------------------------

  There are several post-processing calculations that are useful for
  evaluating functions of interest. These functions are primarily
  used in the TACSFunction classes for evaluating output functions
  of interest.

  TACSConstitutive * getConstitutive(): return the constitutive
  relationship object (this may return NULL if no constitutive class
  is provided).

  getNumGaussPts(): Get the number of quadrature points in the element

  getGaussWtsPts(): Get the quadrature locations and weights

  getShapeFunctions(): Return the shape functions for this element

  getDetJacobian(): Get the determinant of the Jacobian evaluated at the
  specified parametric location

  getDetJacobianXptSens(): Return the derivative of the determinant of
  the Jacobian w.r.t. all nodal coordinates.

  getStrain(): Retrieve the strain evaluated at a parametric location

  addStrainXptSens(): Add the derivative of the strain w.r.t. all
  nodal locations to the provided input

  addStrainSVSens(): Add the derivative of the strain w.r.t. all nodal
  displacements
*/

// The element types used for visualization
enum ElementType { TACS_ELEMENT_NONE,
                   TACS_POINT_ELEMENT,
                   TACS_EULER_BEAM,
                   TACS_TIMOSHENKO_BEAM,
                   TACS_PLANE_STRESS,
                   TACS_SHELL,
                   TACS_SOLID,
                   TACS_Q3D_ELEMENT,
                   TACS_RIGID,
                   TACS_POISSON_2D_ELEMENT,
                   TACS_POISSON_3D_ELEMENT };

// The different element matrix types
enum ElementMatrixType { STIFFNESS_MATRIX,
                         MASS_MATRIX,
                         GEOMETRIC_STIFFNESS_MATRIX,
                         STIFFNESS_PRODUCT_DERIVATIVE };

// Element behavior types
enum ElementBehaviorType{ LINEAR,
                          NONLINEAR,
                          LARGE_ROTATION };

// The TACSElement base class
class TACSElement : public TACSOptObject {
 public:
  TACSElement( int _componentNum=0 ){
    componentNum = _componentNum;
  }
  virtual ~TACSElement(){}

  // Get the number of displacements, stresses, nodes, etc.
  // ------------------------------------------------------
  virtual int numDisplacements() = 0; // Degrees of freedom per node
  virtual int numNodes() = 0; // Number of nodes for this element

  // Number of stresses (possibly zero)
  virtual int numStresses(){
    return 0;
  }

  // Number of variables for this element (nodes times dof/node)
  virtual int numVariables(){
    return numNodes()*numDisplacements();
  }

  // Identifies whether the nodes are associated with the multipliers
  virtual void getMultiplierIndex( int *multiplier ){
    *multiplier = -1;
  }

  // Retrieve the initial conditions and add the derivative
  // ------------------------------------------------------
  virtual void getInitConditions( TacsScalar vars[],
                                  TacsScalar dvars[],
                                  TacsScalar ddvars[],
                                  const TacsScalar Xpts[] ){
    memset(vars, 0, numVariables()*sizeof(TacsScalar));
    memset(dvars, 0, numVariables()*sizeof(TacsScalar));
    memset(ddvars, 0, numVariables()*sizeof(TacsScalar));
  }

  // Add the product of the initial condition with the adjoint variables
  // -------------------------------------------------------------------
  virtual void addInitConditionAdjResProduct( TacsScalar fdvSens[], int dvLen,
                                              const TacsScalar adjVars[],
                                              const TacsScalar adjDVars[],
                                              const TacsScalar Xpts[] ){}
  virtual void getInitConditionAdjResXptProduct( TacsScalar fXptSens[],
                                                 const TacsScalar adjVars[],
                                                 const TacsScalar adjDVars[],
                                                 const TacsScalar Xpts[] ){
    memset(fXptSens, 0, 3*numNodes()*sizeof(TacsScalar));
  }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  virtual void computeEnergies( double time,
                                TacsScalar *_Te,
                                TacsScalar *_Pe,
                                const TacsScalar Xpts[],
                                const TacsScalar vars[],
                                const TacsScalar dvars[] ){
    *_Te = 0.0;
    *_Pe = 0.0;
  }

  // Compute the residual of the governing equations
  // -----------------------------------------------
  virtual void addResidual( double time, TacsScalar res[],
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[] ) = 0;

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  virtual void addJacobian( double time, TacsScalar J[],
                            double alpha, double beta, double gamma,
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[] );

  // Add the product of the adjoint variables with the derivative of the residual
  // ----------------------------------------------------------------------------
  virtual void addAdjResProduct( double time, double scale,
                                 TacsScalar dvSens[], int dvLen,
                                 const TacsScalar psi[],
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[],
                                 const TacsScalar dvars[],
                                 const TacsScalar ddvars[] ){}
  virtual void addAdjResXptProduct( double time, double scale,
                                    TacsScalar fXptSens[],
                                    const TacsScalar psi[],
                                    const TacsScalar Xpts[],
                                    const TacsScalar vars[],
                                    const TacsScalar dvars[],
                                    const TacsScalar ddvars[] ){}

  // Retrieve a specific time-independent matrix from the element
  // ------------------------------------------------------------
  virtual void getMatType( ElementMatrixType matType,
                           TacsScalar mat[],
                           const TacsScalar Xpts[],
                           const TacsScalar vars[] ){
    int size = numVariables()*numVariables();
    memset(mat, 0, size*sizeof(TacsScalar));
  }

  // Compute the derivative of the inner product w.r.t. design variables
  // -------------------------------------------------------------------
  virtual void addMatDVSensInnerProduct( ElementMatrixType matType,
                                         double scale,
                                         TacsScalar dvSens[], int dvLen,
                                         const TacsScalar psi[],
                                         const TacsScalar phi[],
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[] ){}

  // Compute the derivative of the inner product w.r.t. vars[]
  // ---------------------------------------------------------
  virtual void getMatSVSensInnerProduct( ElementMatrixType matType,
                                         TacsScalar res[],
                                         const TacsScalar psi[],
                                         const TacsScalar phi[],
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[] ){
    memset(res, 0, numVariables()*sizeof(TacsScalar));
  }

  // Member functions for evaluating global functions of interest
  // ------------------------------------------------------------
  virtual TACSConstitutive *getConstitutive(){ return NULL; }

  // Get the number of Gauss quadrature points
  // -----------------------------------------
  virtual int getNumGaussPts(){ return 0; }

  // Get the quadrature points and weights
  // -------------------------------------
  virtual double getGaussWtsPts( const int num, double *pt ){
    return 0.0;
  }

  // Get the shape functions from the element
  // ----------------------------------------
  virtual void getShapeFunctions( const double pt[], double N[] ){}

  // Return the determinant of the Jacobian of the transformation
  // ------------------------------------------------------------
  virtual TacsScalar getDetJacobian( const double pt[],
                                     const TacsScalar Xpts[] ){
    return 0.0;
  }

  // Return the determinant of the Jacobian and its sensitivity at this point
  // ------------------------------------------------------------------------
  virtual TacsScalar getDetJacobianXptSens( TacsScalar *hXptSens,
                                            const double pt[],
                                            const TacsScalar Xpts[] ){
    memset(hXptSens, 0, 3*numNodes()*sizeof(TacsScalar));
    return getDetJacobian(pt, Xpts);
  }

  // This function returns the strain evaluated at pt
  // ------------------------------------------------
  virtual void getStrain( TacsScalar strain[],
                          const double pt[],
                          const TacsScalar Xpts[],
                          const TacsScalar vars[] ){}

  // This function adds the sensitivity of the strain w.r.t. Xpts
  // ------------------------------------------------------------
  virtual void addStrainXptSens( TacsScalar strainXptSens[],
                                 const double pt[],
                                 const TacsScalar scale,
                                 const TacsScalar strainSens[],
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[] ){}

  // This function adds the sensitivity of the strain to the state variables
  // -----------------------------------------------------------------------
  virtual void addStrainSVSens( TacsScalar strainSVSens[],
                                const double pt[],
                                const TacsScalar scale,
                                const TacsScalar strainSens[],
                                const TacsScalar Xpts[],
                                const TacsScalar vars[] ){}

  // Function used for localizing the error to nodes with PU-weights
  // ---------------------------------------------------------------
  virtual void addLocalizedError( double time, TacsScalar err[],
                                  const TacsScalar adjoint[],
                                  const TacsScalar Xpts[],
                                  const TacsScalar vars[] ){}

  // These constants are used to denote which output to obtain
  // ---------------------------------------------------------
  static const unsigned int OUTPUT_NODES = 1;
  static const unsigned int OUTPUT_DISPLACEMENTS = 2;
  static const unsigned int OUTPUT_STRAINS = 4;
  static const unsigned int OUTPUT_STRESSES = 8;
  static const unsigned int OUTPUT_EXTRAS = 16;
  static const unsigned int OUTPUT_COORDINATES = 32;

  // Retrieve information about the name and quantity of variables
  // -------------------------------------------------------------
  virtual const char *elementName(){ return NULL; }
  virtual const char *displacementName( int i ){ return NULL; }
  virtual const char *stressName( int i ){ return NULL; }
  virtual const char *strainName( int i ){ return NULL; }
  virtual const char *extraName( int i ){ return NULL; }

  // Return the name of the element
  // ------------------------------
  virtual const char *TACSObjectName(){ return this->elementName(); }

  // Get the number of extras and element type information
  // -----------------------------------------------------
  virtual int numExtras(){ return 0; }
  virtual enum ElementType getElementType(){ return TACS_ELEMENT_NONE; }

  // Functions for retrieving data from the element for visualization
  // ----------------------------------------------------------------
  void setComponentNum( int comp_num ){ componentNum = comp_num; }
  int getComponentNum(){ return componentNum; }
  virtual void addOutputCount( int *nelems, int *nnodes, int *ncsr ){}
  virtual void getOutputData( unsigned int out_type,
                              double * data, int ld_data,
                              const TacsScalar Xpts[],
                              const TacsScalar vars[] ){}
  virtual void getOutputConnectivity( int * con, int start_node ){}

  // Test functions used to test the derivative evaluation code
  // ----------------------------------------------------------
  static void setFailTolerances( double fail_rtol, double fail_atol );
  static void setPrintLevel( int flag );
  static void setStepSize( double dh );

  int testResidual( double time, const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );
  int testResidual( double time, const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[],
                    const int multipliers[],
                    int nmultipliers );
  int testJacobian( double time, const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[], int col=-1 );
  int testStrainSVSens( const TacsScalar Xpts[],
                        const TacsScalar vars[] );
  int testStrainXptSens( const TacsScalar Xpts[],
                         const TacsScalar vars[] );
  int testAdjResProduct( const TacsScalar *x, int dvLen,
                         double time, const TacsScalar Xpts[],
                         const TacsScalar vars[],
                         const TacsScalar dvars[],
                         const TacsScalar ddvars[] );
  int testAdjResXptProduct( double time, const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[] );
  int testJacobianXptSens( const TacsScalar Xpts[] );
  int testMatDVSensInnerProduct( ElementMatrixType matType,
                                 const TacsScalar *x, int dvLen,
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[] );

  int testMatSVSensInnerProduct( ElementMatrixType matType,
                                 const TacsScalar *x, int dvLen,
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[] );
 private:
  int componentNum;

  // Static information used in the test functions
  static int test_print_level;
  static double test_step_size;
  static double test_fail_rtol;
  static double test_fail_atol;
};

#endif // TACS_ELEMENT_H
