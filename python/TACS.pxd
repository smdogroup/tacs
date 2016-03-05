# For MPI capabilities
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
import numpy as np
cimport numpy as np
from libc.string cimport const_char

cdef extern from "TACSObject.h":
    cdef cppclass TACSObject:
       TACSObject()

cdef extern from "TACSAssembler.h":
    cdef cppclass TACSAssembler:
       TACSAssembler(MPI_Comm _tacs_comm, int numOwnedNodes,
                     int _varsPerNode, int _numElements, 
                     int _numNodes, int numDependentNodes,
                     int _nodeMaxCSRsize)
       # Return information about the TACSObject
       int getNumNodes(void *_self)
       int getNumDependentNodes(void *_self)
       int getNumElements(void *_self)
       #VarMap getVarMap(void *_self)
       
       # Add nodes to TACS
       void addNode(int localNodeNum, int tacsNodeNum)
