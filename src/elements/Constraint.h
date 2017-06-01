#ifndef TACS_CONSTRAINT_H
#define TACS_CONSTRAINT_H

/*
  Constraints for Rigid-body dynamics in TACS.
  
  Copyright (c) 2015-2017 Graeme Kennedy. All rights reserved.
*/

#include "TACSElement.h"
#include "TACSGibbsVector.h"

/*
  Cylindrical constraint between rigid bodies
*/
class TACSCylindricalConstraint : public TACSElement {
 public:
  TACSCylindricalConstraint( TACSRigidBody *_bodyA,
                             TACSRigidBody *_bodyB,
                             TACSGibbsVector *_point,
                             TACSGibbsVector *_eAVec );
  TACSCylindricalConstraint( TACSRigidBody *_bodyA,
                             TACSGibbsVector *_point,
                             TACSGibbsVector *_eAVec );
  ~TACSCylindricalConstraint();

  // Set and retrieve design variable values
  // ---------------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );

  // Return the number of displacements and nodes
  // --------------------------------------------
  int numDisplacements(){ return 8; }
  int numNodes();
  const char* elementName(){ return elem_name; }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies( double time,
                        TacsScalar *_Te, 
                        TacsScalar *_Pe,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[] );

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );
 private:
  void updatePoints( int init_e=0 );  // Update the local data
  TACSRigidBody *bodyA, *bodyB;   // The rigid bodies involved in the joint
  TACSGibbsVector *point;           // The point where the joint is located in global frame
  TACSGibbsVector *eAVec;           // The revolute direction in global frame
  TACSGibbsVector *xAVec, *xBVec;   // The positions of joint from each body in global frame
  TACSGibbsVector *eB1Vec, *eB2Vec; // The positions of joint from each body in global fram
  TACSGibbsVector *eVec;            // The coordinate direction in global frame with minimal dot product with eAVec
  static const char *elem_name;       // The name of the element
};

#endif // TACS_CONSTRAINT_H
