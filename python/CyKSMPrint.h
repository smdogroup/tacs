#ifndef CY_TACS_KSM_PRINT_H
#define CY_TACS_KSM_PRINT_H

#include "KSM.h"

class CyKSMPrint : public KSMPrint {
 public:
  KSMPrint();

  // Set the member callback functions that are required
  // ---------------------------------------------------
  void setSelfPointer( void *_self );

 private:
  // Public member function pointers to the callbacks that are
  // required before class can be used
  // ---------------------------------------------------------
  void *self;
}

#endif // CY_TACS_KSM_PRINT
