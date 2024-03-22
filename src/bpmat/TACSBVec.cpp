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

#include "TACSBVec.h"

#include "tacslapack.h"

/*
  The following class defines the dependent node information.

  Note that this class steals the ownership of the data.
*/
TACSBVecDepNodes::TACSBVecDepNodes(int _ndep_nodes, int **_dep_ptr,
                                   int **_dep_conn, double **_dep_weights) {
  ndep_nodes = _ndep_nodes;
  dep_ptr = *_dep_ptr;
  *_dep_ptr = NULL;
  dep_conn = *_dep_conn;
  *_dep_conn = NULL;
  dep_weights = *_dep_weights;
  *_dep_weights = NULL;
}

TACSBVecDepNodes::~TACSBVecDepNodes() {
  delete[] dep_ptr;
  delete[] dep_conn;
  delete[] dep_weights;
}

/*
  Get the dependent connectivity and weights
*/
int TACSBVecDepNodes::getDepNodes(const int **_dep_ptr, const int **_dep_conn,
                                  const double **_dep_weights) {
  if (_dep_ptr) {
    *_dep_ptr = dep_ptr;
  }
  if (_dep_conn) {
    *_dep_conn = dep_conn;
  }
  if (_dep_weights) {
    *_dep_weights = dep_weights;
  }
  return ndep_nodes;
}

/*
  Get the dependent node connectivity for reordering
*/
int TACSBVecDepNodes::getDepNodeReorder(const int **_dep_ptr, int **_dep_conn) {
  if (_dep_ptr) {
    *_dep_ptr = dep_ptr;
  }
  if (_dep_conn) {
    *_dep_conn = dep_conn;
  }
  return ndep_nodes;
}

/*!
  Create a block-based parallel vector

  input:
  rmap:   the variable->processor map for the unknowns
  bcs:    the boundary conditions associated with this vector
*/
TACSBVec::TACSBVec(TACSNodeMap *map, int _bsize, TACSBVecDistribute *_ext_dist,
                   TACSBVecDepNodes *_dep_nodes) {
  node_map = map;
  node_map->incref();

  // Get the MPI communicator
  comm = node_map->getMPIComm();

  // Set the block size
  bsize = _bsize;
  size = bsize * node_map->getNumNodes();

  // Allocate the array of owned unknowns
  x = new TacsScalar[size];
  memset(x, 0, size * sizeof(TacsScalar));

  // Set the external data
  ext_dist = _ext_dist;
  if (ext_dist) {
    ext_dist->incref();
    ext_indices = ext_dist->getIndices();
    ext_indices->incref();
    ext_size = bsize * ext_dist->getNumNodes();
    x_ext = new TacsScalar[ext_size];
    memset(x_ext, 0, ext_size * sizeof(TacsScalar));

    // Create the communicator context
    ext_ctx = ext_dist->createCtx(bsize);
    ext_ctx->incref();
  } else {
    ext_size = 0;
    ext_indices = NULL;
    x_ext = NULL;
    ext_ctx = NULL;
  }

  // Set the dependent node data (if defined)
  dep_nodes = _dep_nodes;
  if (dep_nodes) {
    dep_nodes->incref();
    dep_size = bsize * dep_nodes->getDepNodes(NULL, NULL, NULL);
    x_dep = new TacsScalar[dep_size];
    memset(x_dep, 0, dep_size * sizeof(TacsScalar));
  } else {
    dep_size = 0;
    x_dep = NULL;
  }
}

/*
  Create the block-based parallel vector without a TACSNodeMap object
  or boundary conditions - this is required for some parallel matrix
  objects, or other situtations in which neither object exists.

  input:
  comm:  the communicator that this object is distributed over
  size:  the number of local entries stored by this vector
*/
TACSBVec::TACSBVec(MPI_Comm _comm, int _size, int _bsize) {
  bsize = _bsize;
  size = _size;
  comm = _comm;
  node_map = NULL;

  x = new TacsScalar[size];
  memset(x, 0, size * sizeof(TacsScalar));

  // Zero/NULL the external data
  ext_size = 0;
  x_ext = NULL;
  ext_dist = NULL;
  ext_ctx = NULL;

  // Zero/NULL the dependent node data
  dep_size = 0;
  x_dep = NULL;
  dep_nodes = NULL;
}

TACSBVec::~TACSBVec() {
  if (node_map) {
    node_map->decref();
  }
  if (x) {
    delete[] x;
  }
  if (x_ext) {
    delete[] x_ext;
  }
  if (ext_dist) {
    ext_dist->decref();
  }
  if (ext_indices) {
    ext_indices->decref();
  }
  if (ext_ctx) {
    ext_ctx->decref();
  }
  if (x_dep) {
    delete[] x_dep;
  }
  if (dep_nodes) {
    dep_nodes->decref();
  }
}

/*
  Get the local size of the vector on this processor
*/
void TACSBVec::getSize(int *_size) { *_size = size; }

/*
  Compute the norm of the vector
*/
TacsScalar TACSBVec::norm() {
  // Compute the norm for each processor
  TacsScalar res, sum;
#if defined(TACS_USE_COMPLEX)
  res = 0.0;
  int i = 0;
  int rem = size % 4;
  TacsScalar *y = x;
  for (; i < rem; i++) {
    res += y[0] * y[0];
    y++;
  }

  for (; i < size; i += 4) {
    res += y[0] * y[0] + y[1] * y[1] + y[2] * y[2] + y[3] * y[3];
    y += 4;
  }
#else
  int one = 1;
  res = BLASnrm2(&size, x, &one);
  res *= res;
#endif
  TacsAddFlops(2 * size);

  MPI_Allreduce(&res, &sum, 1, TACS_MPI_TYPE, MPI_SUM, comm);

  return sqrt(sum);
}

/*
  Scale the vector by a scalar
*/
void TACSBVec::scale(TacsScalar alpha) {
  int one = 1;
  BLASscal(&size, &alpha, x, &one);
  TacsAddFlops(size);
}

/*
  Compute the dot product of two vectors
*/
TacsScalar TACSBVec::dot(TACSVec *tvec) {
  TacsScalar sum = 0.0;
  TACSBVec *vec = dynamic_cast<TACSBVec *>(tvec);
  if (vec) {
    if (vec->size != size) {
      fprintf(stderr, "TACSBVec::dot Error, the sizes must be the same\n");
      return 0.0;
    }

    TacsScalar res;
#if defined(TACS_USE_COMPLEX)
    res = 0.0;
    int i = 0;
    int rem = size % 4;
    TacsScalar *y = x;
    TacsScalar *z = vec->x;
    for (; i < rem; i++) {
      res += y[0] * z[0];
      y++;
      z++;
    }

    for (; i < size; i += 4) {
      res += y[0] * z[0] + y[1] * z[1] + y[2] * z[2] + y[3] * z[3];
      y += 4;
      z += 4;
    }
#else
    int one = 1;
    res = BLASdot(&size, x, &one, vec->x, &one);
#endif
    MPI_Allreduce(&res, &sum, 1, TACS_MPI_TYPE, MPI_SUM, comm);
  } else {
    fprintf(stderr, "TACSBVec type error: Input must be TACSBVec\n");
  }

  TacsAddFlops(2 * size);

  return sum;
}

/*
  Compute multiple dot products. This is more efficient for parallel
  computations since there are fewer gather operations for the same
  number of dot products.
*/
void TACSBVec::mdot(TACSVec **tvec, TacsScalar *ans, int nvecs) {
  for (int k = 0; k < nvecs; k++) {
    ans[k] = 0.0;

    TACSBVec *vec = dynamic_cast<TACSBVec *>(tvec[k]);
    if (vec) {
      if (vec->size != size) {
        fprintf(stderr, "TACSBVec::dot Error, the sizes must be the same\n");
        continue;
      }

#if defined(TACS_USE_COMPLEX)
      TacsScalar res = 0.0;
      int i = 0;
      int rem = size % 4;
      TacsScalar *y = x;
      TacsScalar *z = vec->x;
      for (; i < rem; i++) {
        res += y[0] * z[0];
        y++;
        z++;
      }

      for (; i < size; i += 4) {
        res += y[0] * z[0] + y[1] * z[1] + y[2] * z[2] + y[3] * z[3];
        y += 4;
        z += 4;
      }

      ans[k] = res;
#else
      int one = 1;
      ans[k] = BLASdot(&size, x, &one, vec->x, &one);
#endif
    } else {
      fprintf(stderr, "TACSBVec type error: Input must be TACSBVec\n");
    }
  }

  TacsAddFlops(2 * nvecs * size);

  MPI_Allreduce(MPI_IN_PLACE, ans, nvecs, TACS_MPI_TYPE, MPI_SUM, comm);
}

/*
  Compute y = alpha*x + y
*/
void TACSBVec::axpy(TacsScalar alpha, TACSVec *tvec) {
  TACSBVec *vec = dynamic_cast<TACSBVec *>(tvec);

  if (vec) {
    if (vec->size != size) {
      fprintf(stderr, "TACSBVec::axpy Error, the sizes must be the same\n");
      return;
    }

    int one = 1;
    BLASaxpy(&size, &alpha, vec->x, &one, x, &one);
  } else {
    fprintf(stderr, "TACSBVec type error: Input must be TACSBVec\n");
  }

  TacsAddFlops(2 * size);
}

/*
  Compute x <- alpha *vec + beta *x
*/
void TACSBVec::axpby(TacsScalar alpha, TacsScalar beta, TACSVec *tvec) {
  TACSBVec *vec = dynamic_cast<TACSBVec *>(tvec);

  if (vec) {
    if (vec->size != size) {
      fprintf(stderr, "TACSBVec::axpby Error sizes must be the same\n");
      return;
    }

    int i = 0;
    int rem = size % 4;
    TacsScalar *y = x;
    TacsScalar *z = vec->x;

    for (; i < rem; i++) {
      y[0] = beta * y[0] + alpha * z[0];
      y++;
      z++;
    }

    for (; i < size; i += 4) {
      y[0] = beta * y[0] + alpha * z[0];
      y[1] = beta * y[1] + alpha * z[1];
      y[2] = beta * y[2] + alpha * z[2];
      y[3] = beta * y[3] + alpha * z[3];
      y += 4;
      z += 4;
    }
  } else {
    fprintf(stderr, "TACSBVec type error: Input must be TACSBVec\n");
  }

  TacsAddFlops(3 * size);
}

/*
  Copy the values x <- vec->x
*/
void TACSBVec::copyValues(TACSVec *tvec) {
  TACSBVec *vec = dynamic_cast<TACSBVec *>(tvec);
  if (vec) {
    if (vec->size != size) {
      fprintf(stderr,
              "TACSBVec::copyValues error, sizes %d and %d must be the same\n",
              size, vec->size);
      return;
    }

    int one = 1;
    BLAScopy(&size, vec->x, &one, x, &one);

    // Copy the external nodes and dependent nodes only if this
    // vector and the source share the same objects
    if (x_ext && vec->x_ext && ext_dist == vec->ext_dist) {
      BLAScopy(&ext_size, vec->x_ext, &one, x_ext, &one);
    }
    if (x_dep && vec->x_dep && dep_nodes == vec->dep_nodes) {
      BLAScopy(&dep_size, vec->x_dep, &one, x_dep, &one);
    }
  } else {
    fprintf(stderr, "TACSBVec type error: Input must be TACSBVec\n");
  }
}

/*
  Zero all the entries in the vector
*/
void TACSBVec::zeroEntries() {
  memset(x, 0, size * sizeof(TacsScalar));

  // Also zero the external and dependent nodes
  if (x_ext) {
    memset(x_ext, 0, ext_size * sizeof(TacsScalar));
  }
  if (x_dep) {
    memset(x_dep, 0, dep_size * sizeof(TacsScalar));
  }
}

/*
  Set all the entries in the vector to val
*/
void TACSBVec::set(TacsScalar val) {
  int i = 0;
  int rem = size % 4;
  TacsScalar *y = x;

  for (; i < rem; i++) {
    y[0] = val;
    y++;
  }

  for (; i < size; i += 4) {
    y[0] = y[1] = y[2] = y[3] = val;
    y += 4;
  }
}

/*
  Initialize the random value generator
*/
void TACSBVec::initRand() {
  unsigned int t = time(NULL);
  MPI_Bcast(&t, 1, MPI_INT, 0, comm);
  srand(t);
}

/*
  Set all the values in the vector using a uniform pseudo-random
  distribution over an interval between lower/upper.
*/
void TACSBVec::setRand(double lower, double upper) {
  if (!node_map) {
    for (int i = 0; i < size; i++) {
      x[i] = lower + ((upper - lower) * rand()) / (1.0 * RAND_MAX);
    }
  } else {
    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    const int *owner_range;
    node_map->getOwnerRange(&owner_range);

    // Generate random values for each processor sequentially.
    // This should result in the same number of calls on all
    // processors.
    for (int k = 0; k < mpi_size; k++) {
      if (k != mpi_rank) {
        int end = bsize * owner_range[k + 1];
        for (int i = bsize * owner_range[k]; i < end; i++) {
          rand();
        }
      } else {
        for (int i = 0; i < size; i++) {
          x[i] = lower + ((upper - lower) * rand()) / (1.0 * RAND_MAX);
        }
      }
    }
  }
}

/*
  Retrieve the locally stored values from the array
*/
int TACSBVec::getArray(TacsScalar **array) {
  if (array) {
    *array = x;
  }
  return size;
}

/*
  Retrieve the array of dependent node values
*/
int TACSBVec::getDepArray(TacsScalar **array) {
  if (array) {
    *array = x_dep;
  }
  return dep_size;
}

/*
  Retrieve the external values stored on this processor
*/
int TACSBVec::getExtArray(TacsScalar **array) {
  if (array) {
    *array = x_ext;
  }
  return ext_size;
}

/*
  Access the var map
*/
TACSNodeMap *TACSBVec::getNodeMap() { return node_map; }

/*
  Retrieve the external index values
*/
TACSBVecIndices *TACSBVec::getBVecIndices() { return ext_indices; }

/*
  Retrieve the external value distribution object
*/
TACSBVecDistribute *TACSBVec::getBVecDistribute() { return ext_dist; }

/*
  Retrieve the dependent node data
*/
TACSBVecDepNodes *TACSBVec::getBVecDepNodes() { return dep_nodes; }

/*
  Apply the Dirichlet boundary conditions to the vector
*/
void TACSBVec::applyBCs(TACSBcMap *bcmap, TACSVec *tvec) {
  TacsScalar *uvals = NULL;
  if (tvec) {
    TACSBVec *vec = dynamic_cast<TACSBVec *>(tvec);
    vec->getArray(&uvals);
  }

  // apply the boundary conditions
  if (x) {
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);

    // Get ownership range
    const int *owner_range;
    node_map->getOwnerRange(&owner_range);

    // Get the values from the boundary condition arrays
    const int *nodes, *vars;
    TacsScalar *values;
    int nbcs = bcmap->getBCs(&nodes, &vars, &values);

    if (uvals) {
      for (int i = 0; i < nbcs; i++) {
        if (nodes[i] >= owner_range[mpi_rank] &&
            nodes[i] < owner_range[mpi_rank + 1]) {
          int var = bsize * (nodes[i] - owner_range[mpi_rank]);
          for (int k = 0; k < bsize; k++) {
            if (vars[i] & (1 << k)) {
              // Scan through the rows to be zeroed
              x[var + k] = uvals[var + k] - values[bsize * i + k];
            }
          }
        }
      }
    } else {
      for (int i = 0; i < nbcs; i++) {
        if (nodes[i] >= owner_range[mpi_rank] &&
            nodes[i] < owner_range[mpi_rank + 1]) {
          int var = bsize * (nodes[i] - owner_range[mpi_rank]);
          for (int k = 0; k < bsize; k++) {
            if (vars[i] & (1 << k)) {
              // Scan through the rows to be zeroed
              x[var + k] = 0.0;
            }
          }
        }
      }
    }
  }
}

/*
  Set the boundary conditions values (both zero and non-zero values)
*/
void TACSBVec::setBCs(TACSBcMap *bcmap) {
  // apply the boundary conditions
  if (x) {
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);

    // Get ownership range
    const int *owner_range;
    node_map->getOwnerRange(&owner_range);

    // Get the values from the boundary condition arrays
    const int *nodes, *vars;
    TacsScalar *values;
    int nbcs = bcmap->getBCs(&nodes, &vars, &values);

    for (int i = 0; i < nbcs; i++) {
      if (nodes[i] >= owner_range[mpi_rank] &&
          nodes[i] < owner_range[mpi_rank + 1]) {
        int var = bsize * (nodes[i] - owner_range[mpi_rank]);
        // Scan through the rows to be set
        for (int k = 0; k < bsize; k++) {
          if (vars[i] & (1 << k)) {
            x[var + k] = values[bsize * i + k];
          }
        }
      }
    }
  }
}

/*
  Return the object name
*/
const char *TACSBVec::getObjectName() { return vecName; }

/*!
  Write the values to a file.

  This uses MPI file I/O. The filenames must be the same on all
  processors. The format is independent of the number of processors.

  The file format is as follows:
  int                      The length of the vector
  len *sizeof(TacsScalar)  The vector entries
*/
int TACSBVec::writeToFile(const char *filename) {
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Copy the filename
  char *fname = new char[strlen(filename) + 1];
  strcpy(fname, filename);

  // Get the range of variable numbers
  int *range = new int[mpi_size + 1];
  range[0] = 0;
  MPI_Allgather(&size, 1, MPI_INT, &range[1], 1, MPI_INT, comm);
  for (int i = 0; i < mpi_size; i++) {
    range[i + 1] += range[i];
  }

  // Open the MPI file
  int fail = 0;
  MPI_File fp = NULL;
  MPI_File_open(comm, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,
                &fp);

  if (fp) {
    if (mpi_rank == 0) {
      MPI_File_write(fp, &range[mpi_size], 1, MPI_INT, MPI_STATUS_IGNORE);
    }

    char datarep[] = "native";
    MPI_File_set_view(fp, sizeof(int), TACS_MPI_TYPE, TACS_MPI_TYPE, datarep,
                      MPI_INFO_NULL);
    MPI_File_write_at_all(fp, range[mpi_rank], x, size, TACS_MPI_TYPE,
                          MPI_STATUS_IGNORE);
    MPI_File_close(&fp);
  } else {
    fail = 1;
  }

  delete[] range;
  delete[] fname;

  return fail;
}

/*!
  Read values from a binary data file.

  The size of this vector must be the size of the vector originally
  stored in the file otherwise nothing is read in.

  The file format is as follows:
  int                      The length of the vector
  len *sizeof(TacsScalar)  The vector entries
*/
int TACSBVec::readFromFile(const char *filename) {
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Copy the filename
  char *fname = new char[strlen(filename) + 1];
  strcpy(fname, filename);

  // Get the range of variable numbers
  int *range = new int[mpi_size + 1];
  range[0] = 0;
  MPI_Allgather(&size, 1, MPI_INT, &range[1], 1, MPI_INT, comm);
  for (int i = 0; i < mpi_size; i++) {
    range[i + 1] += range[i];
  }

  int fail = 0;
  MPI_File fp = NULL;
  MPI_File_open(comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);

  if (fp) {
    int len = 0;
    if (mpi_rank == 0) {
      MPI_File_read(fp, &len, 1, MPI_INT, MPI_STATUS_IGNORE);
    }
    MPI_Bcast(&len, 1, MPI_INT, 0, comm);
    if (len != range[mpi_size]) {
      fprintf(stderr,
              "[%d] Cannot read TACSBVec from file, incorrect "
              "size %d != %d\n",
              mpi_rank, range[mpi_size], len);
      memset(x, 0, size * sizeof(TacsScalar));

      // Mark this as a failure
      fail = 1;
    }

    char datarep[] = "native";
    MPI_File_set_view(fp, sizeof(int), TACS_MPI_TYPE, TACS_MPI_TYPE, datarep,
                      MPI_INFO_NULL);
    MPI_File_read_at_all(fp, range[mpi_rank], x, size, TACS_MPI_TYPE,
                         MPI_STATUS_IGNORE);
    MPI_File_close(&fp);
  } else {
    fail = 1;
  }

  delete[] range;
  delete[] fname;

  return fail;
}

/*
  Add values to the local entries of the array
*/
void TACSBVec::setValues(int n, const int *index, const TacsScalar *vals,
                         TACSBVecOperation op) {
  // Get the MPI rank
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Get the ownership range
  const int *owner_range;
  node_map->getOwnerRange(&owner_range);

  // Loop over the values
  for (int i = 0; i < n; i++) {
    // Three different options: 1. The variable is locally
    // owned, 2. the variable is dependent, 3. the variable
    // is external
    if (index[i] >= owner_range[rank] && index[i] < owner_range[rank + 1]) {
      // Find the offset into the local array
      int x_index = bsize * (index[i] - owner_range[rank]);
      TacsScalar *y = &x[x_index];

      if (op == TACS_INSERT_VALUES) {
        // Set the values into the array
        for (int k = 0; k < bsize; k++, vals++, y++) {
          y[0] = vals[0];
        }
      } else if (op == TACS_ADD_VALUES) {
        // Add the values to the array
        for (int k = 0; k < bsize; k++, vals++, y++) {
          y[0] += vals[0];
        }
      } else {
        // Insert only the non-zero values
        for (int k = 0; k < bsize; k++, vals++, y++) {
          if (TacsRealPart(vals[0]) != 0.0) {
            y[0] = vals[0];
          }
        }
      }
    } else if (index[i] < 0) {
      // Compute the dependent node
      int dep_index = -bsize * (index[i] + 1);
      TacsScalar *y = &x_dep[dep_index];

      if (op == TACS_INSERT_VALUES) {
        // Insert the values
        for (int k = 0; k < bsize; k++, vals++, y++) {
          y[0] = vals[0];
        }
      } else if (op == TACS_ADD_VALUES) {
        // Add the values to the array
        for (int k = 0; k < bsize; k++, vals++, y++) {
          y[0] += vals[0];
        }
      } else {
        // Insert only the non-zero values
        for (int k = 0; k < bsize; k++, vals++, y++) {
          if (TacsRealPart(vals[0]) != 0.0) {
            y[0] = vals[0];
          }
        }
      }
    } else {
      int ext_index = bsize * ext_indices->findIndex(index[i]);
      TacsScalar *y = &x_ext[ext_index];

      if (op == TACS_INSERT_VALUES) {
        // Insert the values into the array
        for (int k = 0; k < bsize; k++, vals++, y++) {
          y[0] = vals[0];
        }
      } else if (op == TACS_ADD_VALUES) {
        // Add the values into the array
        for (int k = 0; k < bsize; k++, vals++, y++) {
          y[0] += vals[0];
        }
      } else {
        // Insert only the non-zero values
        for (int k = 0; k < bsize; k++, vals++, y++) {
          if (TacsRealPart(vals[0]) != 0.0) {
            y[0] = vals[0];
          }
        }
      }
    }
  }
}

/*
  Begin collecting the vector values to their owners
*/
void TACSBVec::beginSetValues(TACSBVecOperation op) {
  // Get the MPI rank
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Get the ownership range
  const int *owner_range;
  node_map->getOwnerRange(&owner_range);

  if (dep_nodes && (op == TACS_ADD_VALUES)) {
    const int *dep_ptr, *dep_conn;
    const double *dep_weights;
    int ndep = dep_nodes->getDepNodes(&dep_ptr, &dep_conn, &dep_weights);

    const TacsScalar *z = x_dep;
    for (int i = 0; i < ndep; i++, z += bsize) {
      for (int jp = dep_ptr[i]; jp < dep_ptr[i + 1]; jp++) {
        // Check if the dependent node is locally owned
        if (dep_conn[jp] >= owner_range[rank] &&
            dep_conn[jp] < owner_range[rank + 1]) {
          // Find the offset into the local array
          int x_index = bsize * (dep_conn[jp] - owner_range[rank]);
          TacsScalar *y = &x[x_index];

          // Add the values to the array
          for (int k = 0; k < bsize; k++, y++) {
            y[0] += dep_weights[jp] * z[k];
          }
        } else {
          // Add the dependent values to external array
          int ext_index = bsize * ext_indices->findIndex(dep_conn[jp]);
          TacsScalar *y = &x_ext[ext_index];

          // Add the values to the array
          for (int k = 0; k < bsize; k++, y++) {
            y[0] += dep_weights[jp] * z[k];
          }
        }
      }
    }
  }

  // Now initiate the assembly of the residual
  if (ext_dist) {
    ext_dist->beginReverse(ext_ctx, x_ext, x, op);
  }
}

/*
  Finish adding values from the external contributions
*/
void TACSBVec::endSetValues(TACSBVecOperation op) {
  if (ext_dist) {
    ext_dist->endReverse(ext_ctx, x_ext, x, op);
  }

  // Zero the external part and dependent parts
  if (x_ext) {
    memset(x_ext, 0, ext_size * sizeof(TacsScalar));
  }
  if (x_dep) {
    memset(x_dep, 0, dep_size * sizeof(TacsScalar));
  }
}

/*
  Initiate sending values to their destinations
*/
void TACSBVec::beginDistributeValues() {
  if (ext_dist) {
    ext_dist->beginForward(ext_ctx, x, x_ext);
  }
}

/*
  End the distribution of values to their destination processes.

  This must be called before getValues
*/
void TACSBVec::endDistributeValues() {
  if (ext_dist) {
    ext_dist->endForward(ext_ctx, x, x_ext);
  }

  if (dep_nodes) {
    // Get the MPI rank
    int rank;
    MPI_Comm_rank(comm, &rank);

    // Get the ownership range
    const int *owner_range;
    node_map->getOwnerRange(&owner_range);

    // Get the dependent node information
    const int *dep_ptr, *dep_conn;
    const double *dep_weights;
    int ndep = dep_nodes->getDepNodes(&dep_ptr, &dep_conn, &dep_weights);

    // Set a pointer into the dependent variables
    TacsScalar *z = x_dep;
    for (int i = 0; i < ndep; i++, z += bsize) {
      // Zero the variables
      for (int k = 0; k < bsize; k++) {
        z[k] = 0.0;
      }

      // Compute the weighted value of the dependent node
      for (int jp = dep_ptr[i]; jp < dep_ptr[i + 1]; jp++) {
        // Check if the dependent node is locally owned
        if (dep_conn[jp] >= owner_range[rank] &&
            dep_conn[jp] < owner_range[rank + 1]) {
          // Find the offset into the local array
          int x_index = bsize * (dep_conn[jp] - owner_range[rank]);
          TacsScalar *y = &x[x_index];

          // Add the values to the array
          for (int k = 0; k < bsize; k++, y++) {
            z[k] += dep_weights[jp] * y[0];
          }
        } else {
          int ext_index = bsize * ext_indices->findIndex(dep_conn[jp]);
          TacsScalar *y = &x_ext[ext_index];

          // Add the values to the array
          for (int k = 0; k < bsize; k++, y++) {
            z[k] += dep_weights[jp] * y[0];
          }
        }
      }
    }
  }
}

/*
  Get the values from the vector
*/
int TACSBVec::getValues(int n, const int *index, TacsScalar *vals) {
  // Get the MPI rank
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Get the ownership range
  const int *owner_range;
  node_map->getOwnerRange(&owner_range);

  // Fail flag
  int fail = 0;

  // Loop over the values
  for (int i = 0; i < n; i++) {
    // Three different options: 1. The variable is locally
    // owned, 2. the variable is dependent, 3. the variable
    // is external
    if (index[i] >= owner_range[rank] && index[i] < owner_range[rank + 1]) {
      // Find the offset into the local array
      int x_index = bsize * (index[i] - owner_range[rank]);
      TacsScalar *y = &x[x_index];

      // Add the values to the array
      for (int k = 0; k < bsize; k++, vals++, y++) {
        vals[0] = y[0];
      }
    } else if (index[i] < 0) {
      // Compute the dependent node
      int dep_index = -bsize * (index[i] + 1);
      TacsScalar *y = &x_dep[dep_index];

      // Add the values to the array
      if (dep_index < dep_size) {
        for (int k = 0; k < bsize; k++, vals++, y++) {
          vals[0] = y[0];
        }
      } else {
        fail = 1;
      }
    } else {
      int ext_index = bsize * ext_indices->findIndex(index[i]);
      TacsScalar *y = &x_ext[ext_index];

      if (ext_index >= 0) {
        // Add the values to the array
        for (int k = 0; k < bsize; k++, vals++, y++) {
          vals[0] = y[0];
        }
      } else {
        fail = 1;
      }
    }
  }

  return fail;
}

const char *TACSBVec::vecName = "TACSBVec";
