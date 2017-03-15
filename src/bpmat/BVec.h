#ifndef TACS_BVEC_H
#define TACS_BVEC_H

/*
  Block vector classes

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include "KSM.h"
#include "BVecDist.h"

/*
  The following class defines the dependent node information.

  Note that this class steals the ownership of the data. 
*/
class TACSBVecDepNodes : public TACSObject {
 public:
  TACSBVecDepNodes( int _ndep_nodes, 
                    int **_dep_ptr, int **_dep_conn,
                    double **_dep_weights ){
    ndep_nodes = _ndep_nodes;
    dep_ptr = *_dep_ptr;  *_dep_ptr = NULL; 
    dep_conn = *_dep_conn;  *_dep_conn = NULL; 
    dep_weights = *_dep_weights;  *_dep_weights = NULL; 
  }
  ~TACSBVecDepNodes(){
    delete [] dep_ptr;
    delete [] dep_conn;
    delete [] dep_weights;
  }
  
  // Get the dependent connectivity and weights
  // ------------------------------------------
  int getDepNodes( const int **_dep_ptr, const int **_dep_conn,
                   const double **_dep_weights ){
    if (_dep_ptr){ *_dep_ptr = dep_ptr; }
    if (_dep_conn){ *_dep_conn = dep_conn; }
    if (_dep_weights){ *_dep_weights = dep_weights; }
    return ndep_nodes;
  }

  // Get the dependent node connectivity for reordering
  // --------------------------------------------------
  int getDepNodeReorder( const int **_dep_ptr, int **_dep_conn ){
    if (_dep_ptr){ *_dep_ptr = dep_ptr; }
    if (_dep_conn){ *_dep_conn = dep_conn; }
    return ndep_nodes;
  }

 private:
  int ndep_nodes;
  int *dep_ptr, *dep_conn;
  double *dep_weights;
};

/*
  Parallel block-based vector class, distributed based on the variable
  map provided.

  This class inherits from TACSVec and is used for TACSAssembler
  operations. There are two flavours of this class: One with the
  capability to set/add and get values from both local and pre-defined
  global locations within the vector, and the plain vanilla variety
  which only has local variables. 

  The communication of the local/global variable value locations takes
  place behind the scenes and is transparent to the user.
*/
class TACSBVec : public TACSVec {
 public: 
  TACSBVec( TACSVarMap *map, int bsize,
            TACSBVecDistribute *ext_dist=NULL,
            TACSBVecDepNodes *dep_nodes=NULL );
  TACSBVec( MPI_Comm _comm, int size, int bsize );
  ~TACSBVec();  

  // The basic vector operations
  // ---------------------------
  MPI_Comm getMPIComm(){ return comm; }
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
  void applyBCs( TACSBcMap *map, TACSVec *vec=NULL );

  // Get/set the vector elements
  // ---------------------------
  void set( TacsScalar val );           // Set all values of the vector
  void zeroEntries();                   // Zero all the entries
  int getArray( TacsScalar **vals );    // Get the local values
  int getExtArray( TacsScalar **vals ); // Get the external values
  void initRand();                      // Init random number generator
  void setRand( double lower, double upper ); // Set random values

  // Read/write the vector to a binary file -- the same on all procs
  // ---------------------------------------------------------------
  int writeToFile( const char *filename );
  int readFromFile( const char *filename );

  // Retrieve objects stored within the vector class
  // -----------------------------------------------
  TACSBcMap *getBcMap();
  TACSBVecIndices *getBVecIndices();
  TACSBVecDistribute *getBVecDistribute();

  // Add/set the values from the array
  // ---------------------------------
  void setValues( int n, const int *index, const TacsScalar *vals,
                  TACSBVecOperation op );
  void beginSetValues( TACSBVecOperation op );
  void endSetValues( TACSBVecOperation op );

  // Retrieve the values that have been set
  // --------------------------------------
  void beginDistributeValues();
  void endDistributeValues();
  void getValues( int n, const int *index, TacsScalar *vals );

  // Get the name of this object
  // ---------------------------
  const char *TACSObjectName();

 private:
  // The MPI communicator
  MPI_Comm comm;

  // The variable map that defines the global distribution of nodes
  TACSVarMap *var_map;

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
  TACSBVecIndices *ext_indices;
  TACSBVecDistCtx *ext_ctx;

  // Dependent node information
  int dep_size;
  TacsScalar *x_dep;

  // Pointer to the dependent node data
  TACSBVecDepNodes *dep_nodes;

  // Name for the vector
  static const char *vecName;
};

#endif
