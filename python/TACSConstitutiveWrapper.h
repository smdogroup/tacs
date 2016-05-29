#ifndef TACS_CONSTITUTIVE_WRAPPER_H
#define TACS_CONSTITUTIVE_WRAPPER_H

#include "PlaneStressStiffness.h"

class PSStiffnessWrapper : public PlaneStressStiffness {
 public:
  PSStiffnessWrapper(){
    self_ptr = NULL;
    calculatestress = NULL;
    addstressdvsens = NULL;
    getpointwisemass = NULL;
    addpointwisemassdvsens = NULL;
    fail = NULL;
    failstrainsens = NULL;
    addstressdvsens = NULL;
  }
  ~PSStiffnessWrapper(){}

  // Define the object name 
  // ----------------------
  const char * constitutiveName(){ 
    return "PSStiffnessWrapper";
  }

  // Function pointers
  // -----------------
  void *self_ptr; // Pointer to the python object
  void (*calculatestress)( void *, const double *, 
                           const TacsScalar *, TacsScalar* );
  void (*addstressdvsens)( void *, const double *, const TacsScalar *,
                           TacsScalar, const TacsScalar *,
                           TacsScalar *, int );
  void (*getpointwisemass)( void *, const double *, TacsScalar * );
  void (*addpointwisemassdvsens)( void *, const double *,
                                  const TacsScalar *, TacsScalar *, int );
  TacsScalar (*fail)( void *, const double *, const TacsScalar * );
  void (*failstrainsens)( void *, const double *,
                          const TacsScalar *, TacsScalar * );
  void (*addfaildvsens)( void *, const double *, const TacsScalar *, 
                         TacsScalar, TacsScalar *, int );

  // Stress member functions
  // -----------------------
  void calculateStress( const double pt[], 
                        const TacsScalar strain[], 
                        TacsScalar stress[] ){
    calculatestress(self_ptr, pt, strain, stress);
  }
  void addStressDVSens( const double pt[], const TacsScalar strain[], 
                        TacsScalar alpha, const TacsScalar psi[], 
                        TacsScalar dvSens[], int dvLen ){
    addstressdvsens(self_ptr, pt, strain, alpha, psi, dvSens, dvLen);
  }

  // Mass moment member functions
  // ----------------------------
  void getPointwiseMass( const double pt[], 
                         TacsScalar mass[] ){
    getpointwisemass(self_ptr, pt, mass);
  }
  void addPointwiseMassDVSens( const double pt[], 
                               const TacsScalar alpha[],
                               TacsScalar dvSens[], int dvLen ){
    addpointwisemassdvsens(self_ptr, pt, alpha, dvSens, dvLen);
  }

  // Evaluate the failure functions
  // ------------------------------
  void failure( const double pt[], const TacsScalar strain[],
                TacsScalar *fval ){
    *fval = fail(self_ptr, pt, strain);
  }
  void failureStrainSens( const double pt[], 
                          const TacsScalar strain[],
                          TacsScalar sens[] ){
    failstrainsens(self_ptr, pt, strain, sens);
  }
  void addFailureDVSens( const double pt[], const TacsScalar strain[],
                         TacsScalar alpha, TacsScalar dvSens[], int dvLen ){
    addfaildvsens(self_ptr, pt, strain, alpha, dvSens, dvLen);
  }
};

#endif
