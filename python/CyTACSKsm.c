#include "CyTACSKsm.h"

CyTACSKsm::CyTACSKsm():
TACSKsm(){
  self = NULL;
}

CyTACSKsm::~CyTACSKsm(){}

/*
  Set the member callback functions that are required
*/
void CyTACSKsm::setSelfPointer(void *_self){
  self = _self;
}
