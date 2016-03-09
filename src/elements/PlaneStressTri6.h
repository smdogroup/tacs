#ifndef TACS_PLANE_STRESS_TRI6_H
#define TACS_PLANE_STRESS_TRI6_H

#include "TACS2DElement.h"

class PlaneStressTri6 : public TACS2DElement<6> {
 public:
  enum PlaneStressElementType { LINEAR, NONLINEAR };

  PlaneStressTri6( PlaneStressStiffness * _stiff, 
		   PlaneStressElementType type = LINEAR, 
		   int _componentNum = 0 );
  ~PlaneStressTri6();
  
    // Return the name of this element
  // -------------------------------
  const char * elementName() const { return elemName; }

  // Retrieve the shape functions
  // ----------------------------
  void getShapeFunctions( const double pt[], double N[]);
  void getShapeFunctions( const double pt[], double N[],
			  double Na[], double Nb[] );

  // Retrieve the Gauss points/weights
  // ---------------------------------
  int getNumGaussPts();
  double getGaussWtsPts( const int num, double pt[] ); 

  // Functions for post-processing
  // -----------------------------
  void addOutputCount( int * nelems, int * nnodes, int * ncsr );
  void getOutputData( unsigned int out_type, 
		      double * data, int ld_data, 
		      const TacsScalar Xpts[],
		      const TacsScalar vars[] );
  void getOutputConnectivity( int * con, int node );

 private:
  static const int NUM_NODES = 6;
  static const char *elemName; 
};

#endif // TACS_PLANE_STRESS_TRI6_H
