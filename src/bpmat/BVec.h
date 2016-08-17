#ifndef TACS_BVEC_H
#define TACS_BVEC_H

/*
  Block vector classes

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include "KSM.h"

// Declare the classes that are required within this header
class TACSBVecDistribute;
class TACSBVecDistCtx;

#include "BVecDist.h"

/*!
  Variable map for the parallel distribution of a vector

  This class defines the mapping between the variables and processors
  and should be instantiated once for each analysis model.
*/
class TACSVarMap : public TACSObject {
 public:
  TACSVarMap( MPI_Comm _comm, int _N );
  ~TACSVarMap();

  int getDim();
  MPI_Comm getMPIComm();
  void getOwnerRange( const int ** _ownerRange );

 private:
  MPI_Comm comm; // The MPI communicator
  int *ownerRange; // The ownership range of the variables
  int N; // Number of nodes on this processor
};

/*
  Define boundary conditions that are applied after all the
  matrix/vector values have been set.  

  BCMap should only be instantiated once for a given analysis.  All
  other classes in this file should point that instance.
*/
class TACSBcMap : public TACSObject {
 public:
  TACSBcMap( int num_bcs );
  ~TACSBcMap();

  // Add/get the boundary conditions stored in this object
  // -----------------------------------------------------
  void addBC( int local_var, int global_var, 
	      const int bc_nums[], const TacsScalar bc_vals[], 
	      int _nvals );
  int getBCs( const int ** _local, const int ** _global,
	      const int ** _var_ptr, 
	      const int ** _vars, const TacsScalar ** _values );

 private:
  // Information used to apply boundary conditions
  int * local;
  int * global;
  int * var_ptr;
  int * vars;
  TacsScalar * values;
  int nbcs;

  int max_size;
  int max_var_ptr_size;

  int bc_increment;
  int var_ptr_increment;
};

/*
  Parallel block-based vector class, distributed based on the variable
  map provided.

  This class inherits from TACSVec and is used for TACSAssembler
  operations. There are two flavours of this class: One with the
  capability to set/add and get values from anywhere within the array,
  and the plain vanilla variety. The conversion from one to the other
  takes place behind the scenes and is transparent to the user.
*/
class TACSBVec : public TACSVec {
 public: 
  TACSBVec( TACSVarMap *map, int bsize, 
            TACSBcMap *bcs=NULL, TACSBVecDistribute *ext_dist=NULL );
  TACSBVec( MPI_Comm _comm, int size, int bsize );
  ~TACSBVec();  

  // The basic vector operations
  // ---------------------------
  void getSize( int *_size );                // Number of local entries
  int getBlockSize(){ return bsize; }        // Get the block size
  TacsScalar norm();                         // Compute the Cartesian 2 norm
  void scale( TacsScalar alpha );            // Scale the vector by a value
  TacsScalar dot( TACSVec *x );              // Compute x^{T}*y
  void mdot( TACSVec **x, 
             TacsScalar *ans, int m );       // Multiple dot product
  void axpy( TacsScalar alpha, TACSVec *x ); // y <- y + alpha*x
  void copyValues( TACSVec *x );             // Copy values from x to this
  void axpby( TacsScalar alpha, 
	      TacsScalar beta, TACSVec *x ); // y <- alpha*x + beta*y 

  // Get/set the vector elements
  // ---------------------------
  void set( TacsScalar val );         // Set all values of the vector
  void zeroEntries();                 // Zero all the entries
  int getArray( TacsScalar **vals );  // Get the local values
  void applyBCs();                    // Zero rows corresponding to BCs
  void initRand();                    // Init random number generator
  void setRand( double lower, double upper ); // Set random values

  // Read/write the vector to a binary file -- the same on all procs
  // ---------------------------------------------------------------
  int writeToFile( const char *filename );
  int readFromFile( const char *filename );

  // These functions are sometimes required
  // --------------------------------------
  TACSBcMap *getBCMap(){ return bcs; }

  // Convert the global ordering into ids
  // ------------------------------------
  int convertIndex( int n, const int *index, int *ids );
  int convertIndex( int n, int *ids );

  // Add/set/get the values from the array
  // -------------------------------------
  void addValues( int n, const int *ids, const TacsScalar *vals );
  void setValues( int n, const int *ids, const TacsScalar *vals );
  void getValues( int n, const int *ids, TacsScalar *vals );


  const char *TACSObjectName();

 private:
  // The MPI communicator
  MPI_Comm comm;

  // The variable map that defines the global distribution of nodes
  TACSVarMap *var_map;

  // The boundary conditions which may or may not be set.
  TACSBcMap *bcs;

  // The vector block size
  int bsize;

  // On-processor owned part of the array
  int size; // The size of the array
  TacsScalar *x; // The entries allocated for the array

  // External off-processor part (not always allocated)
  int ext_size;
  TacsScalar *x_ext;

  // The data used to communicate with the external part of the
  // vector. This is not always allocated so be careful.
  TACSBVecDistribute *ext_dist;
  TACSBVecDistCtx *ext_ctx;

  // Name for the vector
  static const char *vecName;
};

#endif
