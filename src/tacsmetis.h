#ifndef TACS_INCLUDE_METIS_H
#define TACS_INCLUDE_METIS_H

#ifdef TACS_CPLUSPLUS_METIS
// Use this if METIS is compiled with C++
#include "metis.h"
#else
// Otherwise, assume metis is compiled with C
extern "C" {
  #include "metis.h"
}
#endif // TACS_CPLUSPLUS_METIS

#endif
