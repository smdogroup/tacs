#include "CyKSMPrint.h"

CyKSMPrint::CyKSMPrint():
KSMPrint(){
  self = NULL;
}

CyKSMPrint::~CyKSMPrint(){}

/*
  Set the member callback functions that are required
*/
void CyKSMPrint::setSelfPointer(void *_self){
  self = _self;
}
