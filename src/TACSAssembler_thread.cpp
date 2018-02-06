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

#include "TACSAssembler.h"
#include "tacslapack.h"

/*!
  Schedule the parts of the matrix/residual to assemble
*/
static pthread_mutex_t sched_mutex = PTHREAD_MUTEX_INITIALIZER;

void TACSAssembler::schedPthreadJob( TACSAssembler * tacs,
                                     int * index, int total_size ){
  pthread_mutex_lock(&sched_mutex);

  if (tacs->numCompletedElements < total_size){
    *index = tacs->numCompletedElements;
    tacs->numCompletedElements += 1;
  }
  else {
    *index = -1;
  }
  
  pthread_mutex_unlock(&sched_mutex);
}

/*!
  The threaded-implementation of the residual assembly 

  Note that the input must be the TACSAssemblerPthreadInfo data type.
  This function only uses the following data members:
  
  tacs:     the pointer to the TACSAssembler object
*/
void *TACSAssembler::assembleRes_thread( void *t ){
  TACSAssemblerPthreadInfo *pinfo = 
    static_cast<TACSAssemblerPthreadInfo*>(t);

  // Un-pack information for this computation
  TACSAssembler *tacs = pinfo->tacs;
  TACSBVec *res = pinfo->res;

  // Allocate a temporary array large enough to store everything required  
  int s = tacs->maxElementSize;
  int sx = 3*tacs->maxElementNodes;
  int dataSize = 4*s + sx;
  TacsScalar *data = new TacsScalar[ dataSize ];
  
  // Set pointers to the allocate memory
  TacsScalar *vars = &data[0];
  TacsScalar *dvars = &data[s];
  TacsScalar *ddvars = &data[2*s];
  TacsScalar *elemRes = &data[3*s];
  TacsScalar *elemXpts = &data[4*s];

  // Set the data for the auxiliary elements - if there are any
  int naux = 0, aux_count = 0;
  TACSAuxElem *aux = NULL;
  if (tacs->auxElements){
    naux = tacs->auxElements->getAuxElements(&aux);
  }

  while (tacs->numCompletedElements < tacs->numElements){
    int elemIndex = -1;
    TACSAssembler::schedPthreadJob(tacs, &elemIndex, tacs->numElements);
    
    if (elemIndex >= 0){
      // Get the element object
      TACSElement *element = tacs->elements[elemIndex];

      // Retrieve the variable values
      int ptr = tacs->elementNodeIndex[elemIndex];
      int len = tacs->elementNodeIndex[elemIndex+1] - ptr;
      const int *nodes = &tacs->elementTacsNodes[ptr];
      tacs->xptVec->getValues(len, nodes, elemXpts);
      tacs->varsVec->getValues(len, nodes, vars);
      tacs->dvarsVec->getValues(len, nodes, dvars);
      tacs->ddvarsVec->getValues(len, nodes, ddvars);
  
      // Generate the Jacobian of the element
      element->addResidual(tacs->time, elemRes, elemXpts, 
                           vars, dvars, ddvars);

      // Increment the aux counter until we possibly have
      // aux[aux_count].num == elemIndex
      while (aux_count < naux && aux[aux_count].num < elemIndex){
        aux_count++;
      }

      // Add the residual from the auxiliary elements 
      while (aux_count < naux && aux[aux_count].num == elemIndex){
        aux[aux_count].elem->addResidual(tacs->time, elemRes, elemXpts, 
                                         vars, dvars, ddvars);
        aux_count++;
      }

      // Add the values to the residual when the memory unlocks
      pthread_mutex_lock(&tacs->tacs_mutex);
      res->setValues(len, nodes, elemRes, TACS_ADD_VALUES);
      pthread_mutex_unlock(&tacs->tacs_mutex);
    }
  }

  delete [] data;

  pthread_exit(NULL);
}

/*!
  The threaded-implementation of the matrix assembly 

  Note that the input must be the TACSAssemblerPthreadInfo data type.
  This function only uses the following data members:
  
  tacs:     the pointer to the TACSAssembler object
  A:        the generic TACSMat base class
*/
void *TACSAssembler::assembleJacobian_thread( void *t ){
  TACSAssemblerPthreadInfo *pinfo = 
    static_cast<TACSAssemblerPthreadInfo*>(t);

  // Un-pack information for this computation
  TACSAssembler *tacs = pinfo->tacs;
  TACSBVec *res = pinfo->res;
  TACSMat *A = pinfo->mat;
  double alpha = pinfo->alpha;
  double beta = pinfo->beta;
  double gamma = pinfo->gamma;
  MatrixOrientation matOr = pinfo->matOr;

  // Allocate a temporary array large enough to store everything
  // required
  int s = tacs->maxElementSize;
  int sx = 3*tacs->maxElementNodes;
  int sw = tacs->maxElementIndepNodes;
  int dataSize = 4*s + sx + s*s + sw;
  TacsScalar *data = new TacsScalar[ dataSize ];
  int *idata = new int[ sw + tacs->maxElementNodes + 1 ];

  // Set pointers to the allocate memory
  TacsScalar *vars = &data[0];
  TacsScalar *dvars = &data[s];
  TacsScalar *ddvars = &data[2*s];
  TacsScalar *elemRes = &data[3*s];
  TacsScalar *elemXpts = &data[4*s];
  TacsScalar *elemWeights = &data[4*s + sx];
  TacsScalar *elemMat = &data[4*s + sx + sw];

  // Set the data for the auxiliary elements - if there are any
  int naux = 0, aux_count = 0;
  TACSAuxElem *aux = NULL;
  if (tacs->auxElements){
    naux = tacs->auxElements->getAuxElements(&aux);
  }

  while (tacs->numCompletedElements < tacs->numElements){
    int elemIndex = -1;
    TACSAssembler::schedPthreadJob(tacs, &elemIndex, tacs->numElements);
    
    if (elemIndex >= 0){
      // Get the element object
      TACSElement *element = tacs->elements[elemIndex];

      // Retrieve the variable values
      int ptr = tacs->elementNodeIndex[elemIndex];
      int len = tacs->elementNodeIndex[elemIndex+1] - ptr;
      const int *nodes = &tacs->elementTacsNodes[ptr];
      tacs->xptVec->getValues(len, nodes, elemXpts);
      tacs->varsVec->getValues(len, nodes, vars);
      tacs->dvarsVec->getValues(len, nodes, dvars);
      tacs->ddvarsVec->getValues(len, nodes, ddvars);
  
      // Retrieve the number of element variables
      int nvars = element->numVariables();
      memset(elemRes, 0, nvars*sizeof(TacsScalar));
      memset(elemMat, 0, nvars*nvars*sizeof(TacsScalar));
  
      // Generate the Jacobian of the element
      if (res){
        element->addResidual(tacs->time, elemRes, elemXpts, 
                             vars, dvars, ddvars);
      }
      element->addJacobian(tacs->time, elemMat, 
			   alpha, beta, gamma,
			   elemXpts, vars, dvars, ddvars);

      // Increment the aux counter until we possibly have
      // aux[aux_count].num == elemIndex
      while (aux_count < naux && aux[aux_count].num < elemIndex){
        aux_count++;
      }

      // Add the residual from the auxiliary elements 
      while (aux_count < naux && aux[aux_count].num == elemIndex){
        if (res){
          aux[aux_count].elem->addResidual(tacs->time, elemRes, elemXpts, 
                                           vars, dvars, ddvars);
        }
        aux[aux_count].elem->addJacobian(tacs->time, elemMat, 
                                         alpha, beta, gamma,
                                         elemXpts, vars, dvars, ddvars);
        aux_count++;
      }
      
      pthread_mutex_lock(&tacs->tacs_mutex);
      // Add values to the residual
      if (res){ res->setValues(len, nodes, elemRes, TACS_ADD_VALUES); }

      // Add values to the matrix
      tacs->addMatValues(A, elemIndex, elemMat, idata, elemWeights, matOr);
      pthread_mutex_unlock(&tacs->tacs_mutex);
    }
  }

  delete [] data;
  delete [] idata;

  pthread_exit(NULL);
}

/*!
  The threaded-implementation of the matrix-type assembly 

  This function uses the following information from the
  TACSAssemblerPthreadInfo class:

  A:            the matrix to assemble (output)
  scaleFactor:  scaling factor applied to the matrix
  matType:      the matrix type defined in Element.h
  matOr:        the matrix orientation: NORMAL or TRANSPOSE
*/
void *TACSAssembler::assembleMatType_thread( void *t ){
  TACSAssemblerPthreadInfo *pinfo = 
    static_cast<TACSAssemblerPthreadInfo*>(t);

  // Un-pack information for this computation
  TACSAssembler *tacs = pinfo->tacs;
  TACSMat *A = pinfo->mat;
  ElementMatrixType matType = pinfo->matType;
  MatrixOrientation matOr = pinfo->matOr;

  // Allocate a temporary array large enough to store everything required  
  int s = tacs->maxElementSize;
  int sx = 3*tacs->maxElementNodes;
  int sw = tacs->maxElementIndepNodes;
  int dataSize = s + sx + s*s + sw;
  TacsScalar *data = new TacsScalar[ dataSize ];
  int *idata = new int[ sw + tacs->maxElementNodes + 1 ];
  
  TacsScalar *vars = &data[0];
  TacsScalar *elemXpts = &data[s];
  TacsScalar *elemWeights = &data[s + sx];
  TacsScalar *elemMat = &data[s + sx + sw];
  
  while (tacs->numCompletedElements < tacs->numElements){
    int elemIndex = -1;
    TACSAssembler::schedPthreadJob(tacs, &elemIndex, tacs->numElements);
    
    if (elemIndex >= 0){
      // Get the element
      TACSElement *element = tacs->elements[elemIndex];
      
      // Retrieve the variable values
      // Retrieve the variable values
      int ptr = tacs->elementNodeIndex[elemIndex];
      int len = tacs->elementNodeIndex[elemIndex+1] - ptr;
      const int *nodes = &tacs->elementTacsNodes[ptr];
      tacs->xptVec->getValues(len, nodes, elemXpts);
      tacs->varsVec->getValues(len, nodes, vars);

      // Retrieve the type of the matrix
      element->getMatType(matType, elemMat, elemXpts, vars);
      
      pthread_mutex_lock(&tacs->tacs_mutex);
      // Add values to the matrix
      tacs->addMatValues(A, elemIndex, elemMat, idata, elemWeights, matOr);
      pthread_mutex_unlock(&tacs->tacs_mutex);
    }
  }

  delete [] data;

  pthread_exit(NULL);
}

