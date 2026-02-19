/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_BVEC_H
#define TACS_BVEC_H

/*
  Block vector classes
*/

#include "KSM.h"
#include "TACSBVecDistribute.h"

/*
  The following class defines the dependent node information.

  Note that this class steals the ownership of the data.
*/
class TACSBVecDepNodes : public TACSObject {
 public:
  TACSBVecDepNodes(int _ndep_nodes, int **_dep_ptr, int **_dep_conn,
                   double **_dep_weights);
  ~TACSBVecDepNodes();

  // Get the dependent connectivity and weights
  int getDepNodes(const int **_dep_ptr, const int **_dep_conn,
                  const double **_dep_weights);

  // Get the dependent node connectivity for reordering
  // --------------------------------------------------
  int getDepNodeReorder(const int **_dep_ptr, int **_dep_conn);

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
  TACSBVec(TACSNodeMap *map, int bsize, TACSBVecDistribute *ext_dist = NULL,
           TACSBVecDepNodes *dep_nodes = NULL);
  TACSBVec(MPI_Comm _comm, int size, int bsize);
  ~TACSBVec();

  // The basic vector operations
  // ---------------------------
  MPI_Comm getMPIComm() { return comm; }
  void getSize(int *_size);             // Number of local entries
  int getBlockSize() { return bsize; }  // Get the block size
  TacsScalar norm();                    // Compute the Cartesian 2 norm
  void scale(TacsScalar alpha);         // Scale the vector by a value
  TacsScalar dot(TACSVec *x);           // Compute x^{T}*y
  void mdot(TACSVec **x, TacsScalar *ans, int m);  // Multiple dot product
  void axpy(TacsScalar alpha, TACSVec *x);         // y <- y + alpha*x
  void copyValues(TACSVec *x);                     // Copy values from x to this
  void axpby(TacsScalar alpha, TacsScalar beta,
             TACSVec *x);  // y <- alpha*x + beta*y
  void applyBCs(TACSBcMap *map, TACSVec *vec = NULL);
  void setBCs(TACSBcMap *map);

  // Get/set the vector elements
  // ---------------------------
  void set(TacsScalar val);                  // Set all values of the vector
  void zeroEntries();                        // Zero all the entries
  int getArray(TacsScalar **vals);           // Get the local values
  int getDepArray(TacsScalar **vals);        // Get the dependent values
  int getExtArray(TacsScalar **vals);        // Get the external values
  void initRand();                           // Init random number generator
  void setRand(double lower, double upper);  // Set random values

  // Read/write the vector to a binary file -- the same on all procs
  // ---------------------------------------------------------------
  int writeToFile(const char *filename);
  int readFromFile(const char *filename);

  // Retrieve objects stored within the vector class
  // -----------------------------------------------
  TACSNodeMap *getNodeMap();
  TACSBVecIndices *getBVecIndices();
  TACSBVecDistribute *getBVecDistribute();
  TACSBVecDepNodes *getBVecDepNodes();

  // Add/set the values from the array
  // ---------------------------------
  void setValues(int n, const int *index, const TacsScalar *vals,
                 TACSBVecOperation op);
  void beginSetValues(TACSBVecOperation op);
  void endSetValues(TACSBVecOperation op);

  // Retrieve the values that have been set
  // --------------------------------------
  void beginDistributeValues();
  void endDistributeValues();
  int getValues(int n, const int *index, TacsScalar *vals);

  // Get the name of this object
  // ---------------------------
  const char *getObjectName();

 private:
  // The MPI communicator
  MPI_Comm comm;

  // The variable map that defines the global distribution of nodes
  TACSNodeMap *node_map;

  // The vector block size
  int bsize;

  // On-processor owned part of the array
  int size;       // The size of the array
  TacsScalar *x;  // The entries allocated for the array

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

#endif  // TACS_BVEC_H
