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
#include "TACSElementTypes.h"

// The TACSElement base class
class TACSElement : public TACSObject {
 public:
  TACSElement( int _componentNum=0 ){
    componentNum = _componentNum;
  }
  virtual ~TACSElement(){}

  /**
    Set the component number for this element.

    The component number can be used to identify groups of elements
    for visualization purposes

    @param comp_num The component number assigned to the element
  */
  void setComponentNum( int comp_num ){
    componentNum = comp_num;
  }
  
  /**
    Get the component number for this element
  
    @return The component number for the element
  */
  int getComponentNum(){
    return componentNum;
  }

  /**
    Get the number of degrees of freedom per node for this element

    @return The number of degrees of freedom per node
  */  
  virtual int getVarsPerNode() = 0;

  /**
    Get the number of nodes associated with this element

    @return The number of nodes for this element
  */
  virtual int getNumNodes() = 0; // Number of nodes for this element

  /**
    Get the node index where a Lagrange multiplier is defined.

    A negative index indicates that no multiplier is defined. The
    index is relative to the ordering in the element.
  
    @return Index of a Lagrange multiplier node
  */
  virtual int getMultiplierIndex(){
    return -1;
  }

  /**
    Get the type of element layout for visualization

    @return The layout type for this element
  */
  virtual ElementLayout getElementLayout(){
    return TACS_LAYOUT_NONE;
  }

  /**
    Set the element design variables from the design vector

    @param dvLen The length of the design array
    @param dvs The design variable values 
  */
  void setDesignVars( int dvLen, const TacsScalar dvs[] );

  /**
    Get the element design variables values

    @param dvLen The length of the design array
    @param dvs The design variable values 
  */
  void getDesignVars( int dvLen, TacsScalar dvs[] );

  /**
    Get the lower and upper bounds for the design variable values

    @param dvLen The length of the design array
    @param lowerBound The design variable lower bounds 
    @param lowerBound The design variable upper bounds
  */
  void getDesignVarRange( int dvLen,
                          TacsScalar lowerBound[], 
                          TacsScalar upperBound[] );
  
  /**
    Retrieve the initial conditions for time-dependent analysis

    By default, the initial displacements, velocities and accelerations
    are zero.

    @param Xpts The element node locations
    @param vars The values of the element degrees of freedom
    @param dvars The first time derivative of the element DOF
    @param ddvars The second time derivative of the element DOF
  */
  virtual void getInitConditions( const TacsScalar Xpts[],
                                  TacsScalar vars[],
                                  TacsScalar dvars[],
                                  TacsScalar ddvars[] ){
    int num_vars = getNumNodes()*getVarsPerNode();
    memset(vars, 0, num_vars*sizeof(TacsScalar));
    memset(dvars, 0, num_vars*sizeof(TacsScalar));
    memset(ddvars, 0, num_vars*sizeof(TacsScalar));
  }

  /**
    Add the contributions to the derivative from the initial conditions

    @param Xpts The element node locations
    @param adjVars The values of the element adjoint
    @param adjDVars The adjoint of the first time derivatives
    @param adjDDVars The adjoint of the first time derivatives
    @param ddvars The second time derivative of the element DOF
    @param dvLen The length of the design vector
  */
  virtual void addInitConditionAdjResProduct( const TacsScalar Xpts[],
                                              const TacsScalar adjVars[],
                                              const TacsScalar adjDVars[],
                                              const TacsScalar adjDDVars[],
                                              int dvLen,
                                              TacsScalar fdvSens[] ){}

  /**
    Get the contribution to the derivatives of the initial conditions w.r.t.
    the node locations

    @param Xpts The element node locations
    @param adjVars The values of the element adjoint
    @param adjDVars The adjoint of the first time derivatives
    @param adjDDVars The adjoint of the first time derivatives
    @param fXptSens Derivative w.r.t. the node locations
  */
  virtual void getInitConditionAdjResXptProduct( const TacsScalar Xpts[],
                                                 const TacsScalar adjVars[],
                                                 const TacsScalar adjDVars[],
                                                 const TacsScalar adjDDVars[],
                                                 TacsScalar fXptSens[] ){
    memset(fXptSens, 0, 3*getNumNodes()*sizeof(TacsScalar));
  }

  /**
    Compute the kinetic and potential energy within the element.

    This can be used to evaluate the Hamiltonian and test whether the
    element satisfies the Lagrangian equations of motion.

    @param time The simulation time
    @param Xpts The element node locations
    @param vars The values of the element degrees of freedom
    @param dvars The first time derivative of the element DOF
    @param Te The kinetic energy contributed by this element
    @param Pe the potential energy contributed by this element
  */
  virtual void computeEnergies( double time,
                                const TacsScalar Xpts[],
                                const TacsScalar vars[],
                                const TacsScalar dvars[],
                                TacsScalar *Te,
                                TacsScalar *Pe ){
    *Te = 0.0;
    *Pe = 0.0;
  }

  /**
    Add the contribution from this element to the residual.

    Note that this simply adds, and does not over-write the residual so
    that multiple contributions can be computed.

    @param time The simulation time
    @param Xpts The element node locations
    @param vars The values of the element degrees of freedom
    @param dvars The first time derivative of the element DOF
    @param ddvars The second time derivative of the element DOF
    @param res The element residual input/output
  */
  virtual void addResidual( double time, 
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[],
                            TacsScalar res[] ) = 0;

  /**
    Add the contribution from this element to the residual and Jacobian.

    Note that this simply adds, and does not over-write the Jacobian so
    that multiple contributions can be computed.

    The Jacobian contribution consists of a linear combination of the
    Jacobians with respect to the variables, and their first and second
    time derivatives as follows:

    mat += alpha*d(res)/d(vars) + beta*d(res)/d(dvars) + gamma*d(res)/d(ddvars)

    @param time The simulation time
    @param alpha The coefficient for the DOF Jacobian 
    @param beta The coefficient for the first time derivative DOF Jacobian
    @param gamma The coefficient for the second time derivative DOF Jacobian
    @param Xpts The element node locations
    @param vars The values of the element degrees of freedom
    @param dvars The first time derivative of the element DOF
    @param ddvars The second time derivative of the element DOF
    @param res The element residual input/output
    @param mat The element Jacobian input/output
  */
  virtual void addJacobian( double time,
                            double alpha, double beta, double gamma,
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[],
                            TacsScalar res[],
                            TacsScalar mat[] );

  /**
    Add the derivative of the adjoint-residual product to the output vector

    This adds the contribution scaled by an input factor as follows:

    dvSens += scale*d(psi^{T}*(res))/dx

    By default the code is not implemented, but is not required so that
    analysis can be performed. Correct derivatives require a specific
    implementation.

    @param time The simulation time
    @param scale The coefficient for the derivative result
    @param psi The element adjoint variables
    @param Xpts The element node locations
    @param vars The values of the element degrees of freedom
    @param dvars The first time derivative of the element DOF
    @param ddvars The second time derivative of the element DOF
    @param dvLen The length of the design variable vector
    @param dvSens The derivative vector
  */
  virtual void addAdjResProduct( double time, double scale,
                                 const TacsScalar psi[],
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[],
                                 const TacsScalar dvars[],
                                 const TacsScalar ddvars[],
                                 int dvLen,
                                 TacsScalar dvSens[] ){}

  /**
    Add the derivative of the adjoint-residual product to the output vector

    This adds the contribution scaled by an input factor as follows:

    dvSens += scale*d(psi^{T}*(res))/d(Xpts)

    By default the code is not implemented, but is not required so that
    analysis can be performed. Correct derivatives require a specific
    implementation.

    @param time The simulation time
    @param scale The coefficient for the derivative result
    @param psi The element adjoint variables
    @param Xpts The element node locations
    @param vars The values of the element degrees of freedom
    @param dvars The first time derivative of the element DOF
    @param ddvars The second time derivative of the element DOF
    @param dvLen The length of the design variable vector
    @param dvSens The derivative vector
  */
  virtual void addAdjResXptProduct( double time, double scale,
                                    const TacsScalar psi[],
                                    const TacsScalar Xpts[],
                                    const TacsScalar vars[],
                                    const TacsScalar dvars[],
                                    const TacsScalar ddvars[],
                                    TacsScalar fXptSens[] ){}

  /**
    Compute a specific type of element matrix (mass, stiffness, geometric
    stiffness, etc.)
    
    @param matType The type of element matrix to compute
    @param Xpts The element node locations
    @param vars The values of element degrees of freedom
    @param mat The element matrix output
  */
  virtual void getMatType( ElementMatrixType matType,
                           const TacsScalar Xpts[],
                           const TacsScalar vars[],
                           TacsScalar mat[] ){
    int size = getNumNodes()*getVarsPerNode();
    memset(mat, 0, size*size*sizeof(TacsScalar));
  }

  /**
    Add the derivative of the product of a specific matrix w.r.t.
    the design variables

    dvSens += scale*d(psi^{T}*(mat)*phi)/d(x)
    
    where mat is computed via the getMatType().

    @param matType The type of element matrix to compute
    @param scale The scalar value that multiplies the derivative
    @param psi The left-hand vector
    @param phi The right-hand vector
    @param Xpts The element node locations
    @param vars The values of element degrees of freedom
    @param mat The element matrix output
  */
  virtual void addMatDVSensInnerProduct( ElementMatrixType matType,
                                         double scale,
                                         const TacsScalar psi[],
                                         const TacsScalar phi[],
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[],
                                         int dvLen,
                                         TacsScalar dvSens[] ){}

  /**
    Compute the derivative of the product of a specific matrix w.r.t.
    the input variables (vars).

    dvSens = d(psi^{T}*(mat)*phi)/d(vars)
    
    where mat is computed via the getMatType().

    @param matType The type of element matrix to compute
    @param psi The left-hand vector
    @param phi The right-hand vector
    @param Xpts The element node locations
    @param vars The values of element degrees of freedom
    @param res The residual output The element matrix output
  */
  virtual void getMatSVSensInnerProduct( ElementMatrixType matType,
                                         const TacsScalar psi[],
                                         const TacsScalar phi[],
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[],
                                         TacsScalar res[] ){
    int size = getNumNodes()*getVarsPerNode();
    memset(res, 0, size*sizeof(TacsScalar));
  }

  /**
    Compute the output data for visualization

    @param etype The type of element data to be output
    @param write_flag The type of data to be output
    @param Xpts The element node locations
    @param vars The values of the element degrees of freedom
    @param dvars The first time derivative of the element DOF
    @param ddvars The second time derivative of the element DOF
    @param ld_data The dimension of the data
    @param data The data to be created
  */
  virtual void getOutputData( ElementType etype, int write_flag,
                              const TacsScalar Xpts[],
                              const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[],
                              int ld_data, double *data ){}

 private:
  int componentNum;
};

#endif // TACS_ELEMENT_H
