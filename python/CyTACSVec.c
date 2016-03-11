#include "CyTACSVec.h"

CyTACSVec::CyTACSVec():
TACSVec(){

  self = NULL;

}

CyTACSVec::~CyTACSVec(){}

/*
  Set the member callback functions that are required
*/
void CyTACSVec::setSelfPointer(void *_self){
  self = _self;
}
