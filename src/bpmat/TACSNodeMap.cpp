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

#include "TACSNodeMap.h"

/**
  Defines the mapping from nodes to processor owners.

  The mapping may be either distributed, in which case the node
  variable values are distributed across all processors, or global, in
  which the node values are duplicated across all procs.

  @param comm Communicator for this object
  @param N The number of variables defined (either locally or globally)
  @param vtype The type of node map
*/
TACSNodeMap::TACSNodeMap(MPI_Comm _comm, int _N, NodeMapType _ntype) {
  comm = _comm;
  ntype = _ntype;

  // Get the communicator size
  MPI_Comm_size(comm, &mpiSize);
  MPI_Comm_rank(comm, &mpiRank);

  // The ownership ranges for all processes
  N = _N;

  if (ntype == DISTRIBUTED) {
    ownerRange = new int[mpiSize + 1];
    memset(ownerRange, 0, (mpiSize + 1) * sizeof(int));

    // Get the number of variables
    ownerRange[0] = 0;
    MPI_Allgather(&N, 1, MPI_INT, &ownerRange[1], 1, MPI_INT, comm);

    // Set the ownership values so that they range over the owned
    // unknown node numbers
    for (int i = 0; i < mpiSize; i++) {
      ownerRange[i + 1] += ownerRange[i];
    }
  } else {
    ownerRange = NULL;
  }
}

TACSNodeMap::~TACSNodeMap() {
  if (ownerRange) {
    delete[] ownerRange;
  }
}

/*
  Get the number of nodes on this processor
*/
int TACSNodeMap::getNumNodes() { return N; }

/*
  Get the MPI communicator
*/
MPI_Comm TACSNodeMap::getMPIComm() { return comm; }

/*
  Get the ownership range for this processor
*/
void TACSNodeMap::getOwnerRange(const int **_ownerRange) {
  *_ownerRange = ownerRange;
}

/*
  Get the owner of this processor. If the node number is out of range,
  then return -1;
*/
int TACSNodeMap::getNodeOwner(int node) {
  // If the node is out of range, return immediately
  if (node < 0 || node >= ownerRange[mpiSize]) {
    return -1;
  }

  // Check whether the node is on the current processor
  if (node >= ownerRange[mpiRank] && node < ownerRange[mpiRank + 1]) {
    return mpiRank;
  }

  // Perform a binary search starting from the mid-processor
  int low = 0;
  int high = mpiSize - 1;

  while (high - low >= 0) {
    // Compute the mid-point
    int mid = (high + low) / 2;
    if (node >= ownerRange[mid] && node < ownerRange[mid + 1]) {
      return mid;
    } else if (node < ownerRange[mid]) {
      high = mid - 1;
    } else if (node >= ownerRange[mid + 1]) {
      low = mid + 1;
    }
  }

  return (high + low) / 2;
}
