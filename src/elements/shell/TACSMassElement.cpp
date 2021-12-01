#include "TACSMassElement.h"

/*
  A 6 DOF point mass element
*/

/*
  Create the MassElement.
*/
TACSMassElement::TACSMassElement( TACSGeneralMassConstitutive *_con ){

  con = _con;
  con->incref();

}

TACSMassElement::~TACSMassElement(){
    con->decref();
}

/*
  Retrieve information about the names of the element variables
*/
const char * TACSMassElement::getObjectName(){
  return elemName;
}

/*
  Retrieve the numbers of displacements, nodes, stress and variables
*/
int TACSMassElement::getVarsPerNode() { return NUM_DISPS; }

int TACSMassElement::getNumNodes() { return NUM_NODES; }

enum ElementType TACSMassElement::getElementType(){ return TACS_MASS_ELEMENT; }

/*
  The element name, variable, stress and strain names.
*/
const char * TACSMassElement::elemName = "TACSMassElement";

/*
  Assemble the element residual associated with the given design
  variables and elements.
*/
void TACSMassElement::addResidual( int elemIndex, double time,
                                   const TacsScalar Xpts[],
                                   const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar ddvars[],
                                   TacsScalar res[] ){

  TacsScalar f[6];
  double pt[3] = {0.0, 0.0, 0.0};
  con->evalInertia(elemIndex, pt, Xpts, ddvars, f);
  for (int i = 0; i < NUM_DISPS; i++){
    res[i] += f[i];
  }

}

/*
  Assemble the stiffness matrix for the mass element.
*/
void TACSMassElement::addJacobian( int elemIndex, double time,
                            TacsScalar alpha,
                            TacsScalar beta,
                            TacsScalar gamma,
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[],
                            TacsScalar res[],
                            TacsScalar J[] ){
  double pt[3] = {0.0, 0.0, 0.0};
  for (int j = 0; j < NUM_DISPS; j++){
    TacsScalar N[6], f[6];
    // Shape functions
    memset(N, 0, 6*sizeof(TacsScalar));
    N[j] = 1.0;
    con->evalInertia(elemIndex, pt, Xpts, N, f);
    for (int i = 0; i < NUM_DISPS; i++){
      J[j + i*NUM_VARIABLES] += gamma * f[i];
    }
  }
}

int TACSMassElement::evalPointQuantity( int elemIndex, int quantityType,
                                        double time,
                                        int n, double pt[],
                                        const TacsScalar Xpts[],
                                        const TacsScalar vars[],
                                        const TacsScalar dvars[],
                                        const TacsScalar ddvars[],
                                        TacsScalar *detXd,
                                        TacsScalar *quantity ){
  *detXd = 1.0;
  if (quantityType == TACS_ELEMENT_DENSITY){
    if (quantity){
      *quantity = con->evalDensity(elemIndex, pt, Xpts);
    }
    return 1;
  }
  else if (quantityType == TACS_ELEMENT_DISPLACEMENT){
    if (quantity){
      quantity[0] = vars[0];
      quantity[1] = vars[1];
      quantity[2] = vars[2];
    }
    return 3;
  }

  return 0;
}
