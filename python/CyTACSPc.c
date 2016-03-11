#include "CyTACSPc.h"

CyTACSPc::CyTACSPc():
TACSPc(){
  self = NULL;
}

CyTACSPc::~CyTACSPc(){}

/*
  Set the member callback functions that are required
*/
void CyTACSPc::setSelfPointer(void *_self){
  self = _self;
}
