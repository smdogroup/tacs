#include "CyTACSMat.h"
/*
  Constructor for the TACSMat wrapper
*/
CyTACSMat::CyTACSMat():
TACSMat(){
  // Set the initial values for the callbacks
  self = NULL;
  zeroentries = NULL;
  applybcs = NULL;
  mult = NULL;
}

CyTACSMat::~CyTACSMat(){}

/*
  Set the member callback functions that are required
*/
void CyTACSMat::setSelfPointer( void *_self){
  self = _self;
}
