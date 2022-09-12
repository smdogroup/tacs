/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_NODE_MAP_H
#define TACS_NODE_MAP_H

#include "TACSObject.h"

/**
  Variable map for the parallel distribution of a vector

  This class defines the mapping between the variables and processors
  and should be instantiated once for each analysis model.
*/
class TACSNodeMap : public TACSObject {
 public:
  // Contro the type of variable map
  enum NodeMapType { DISTRIBUTED, GLOBAL };

  TACSNodeMap(MPI_Comm _comm, int _N, NodeMapType vtype = DISTRIBUTED);
  ~TACSNodeMap();

  int getNumNodes();
  MPI_Comm getMPIComm();
  void getOwnerRange(const int **_ownerRange);
  int getNodeOwner(int node);

 private:
  NodeMapType ntype;     // The type of variable map
  MPI_Comm comm;         // The MPI communicator
  int mpiSize, mpiRank;  // The size/rank of the processor
  int *ownerRange;       // The ownership range of the variables
  int N;                 // Number of nodes on this processor
};

#endif  // TACS_NODE_MAP_H
