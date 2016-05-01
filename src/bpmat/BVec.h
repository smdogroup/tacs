#ifndef TACS_BVEC_H
#define TACS_BVEC_H

/*
  Block vector classes

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include <stdio.h>
#include "mpi.h"
#include "KSM.h"

/*
  BCMap should only be instantiated once for a given analysis. 
  All other classes in this file should then point that instance.
*/

/*!
  Variable map for the parallel distribution of a vector
*/
class VarMap : public TACSObject {
 public:
  VarMap( MPI_Comm _comm, int _N );
  ~VarMap();

  int getDim(){ return N; }
  MPI_Comm getMPIComm(){ return comm; }
  void getOwnerRange( const int ** _ownerRange, 
		      int * _mpiRank, int * _mpiSize ){
    *_ownerRange = ownerRange;
    *_mpiRank = mpiRank; *_mpiSize = mpiSize;
  }  

 private:
  MPI_Comm comm; // The MPI communicator
  int *ownerRange; // The ownership range of the variables
  int mpiSize; // Number of MPI procs
  int mpiRank; // MPI rank
  int N; // Number of nodes on this processor
};

/*
  Define boundary conditions that are applied after all 
  the matrix/vector values have been set.
*/
class BCMap : public TACSObject {
 public:
  BCMap( int num_bcs );
  ~BCMap();

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

/*!
  The basic parallel block-based vector class, distributed
  based on the variable map provided. The boundary conditions
  zero the rows.
*/
class BVec : public TACSVec {
 public: 
  BVec( VarMap * _rmap, int bsize, BCMap * bcs=NULL );
  BVec( MPI_Comm _comm, int size, int bsize );
  ~BVec();  

  // The basic vector operations
  // ---------------------------
  void getSize( int * _size ); // Get the number of local entries
  int getBlockSize(){ return bsize; } // Get the block size
  TacsScalar norm(); // Compute the Cartesian 2 norm
  void scale( TacsScalar alpha ); // Scale the vector by a value
  TacsScalar dot( TACSVec * x ); // Compute x^{T} * y
  void mdot( TACSVec ** x, TacsScalar * ans, int m ); // Multiple dot product
  void axpy( TacsScalar alpha, TACSVec * x ); // Compute y <- y + alpha * x
  void copyValues( TACSVec * x ); // Copy values from x to this
  void axpby( TacsScalar alpha, 
	      TacsScalar beta, TACSVec * x ); // Compute y <- alpha * x + beta * y 

  // Get/set the vector elements
  // ---------------------------
  void set( TacsScalar val );         // Set all values of the vector
  void zeroEntries();                 // Zero all the entries
  int getArray( TacsScalar ** vals ); // Get the local values
  void applyBCs();                    // Zero rows corresponding to BCs
  void placeArray( TacsScalar * _x ); // Place the array x[] into the vector
  void restoreArray();                // Restore the original array

  void initRand(); // Initialize the random number generator
  void setRand( double lower, double upper ); // Set values to a random number

  // Copy values to or from a sequential version of the vector
  // ---------------------------------------------------------
  void copySeqValues( BVec * x ); // Copy values from x to this.
  void setValuesSeq( BVec * x );  // Set values from this into x.

  // Read/write the vector to a binary file -- the same on all procs
  // ---------------------------------------------------------------
  int writeToFile( const char * filename );
  int readFromFile( const char * filename );

  // These functions are sometimes required
  // --------------------------------------
  BCMap *getBCMap(){ return bcs; }    // This may return NULL

  virtual const char * TACSObjectName();

  void getOwnership( int *_mpiRank, int *_mpiSize, const int **_ownerRange ){
    *_mpiRank = mpiRank;
    *_mpiSize = mpiSize;
    *_ownerRange = ownerRange;
  }  

 private:
  VarMap *rmap;

  MPI_Comm comm;
  const int *ownerRange;
  int mpiSize;
  int mpiRank;

  TacsScalar *x;
  TacsScalar *displaced; // The displaced array x 
  int bsize; // The block size
  int size; // The size of the array

  // These may be defined - if they're not, applyBCs has no effect
  // Local boundary conditions -- already in local numbering
  BCMap *bcs;

  static const char *vecName;
};

#endif
