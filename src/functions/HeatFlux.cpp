/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2019 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "HeatFlux.h"
#include "TACSAssembler.h"
#include  "ThermoElements.h"
#include "CoupledThermoSolidStiffness.h"
#include "CoupledThermoPlaneStressStiffness.h"
#include "TensorToolbox.h"

/*
  The context for the TMRHeatFlux function
*/
class HeatFluxIntCtx : public TACSFunctionCtx {
 public:
  HeatFluxIntCtx( TACSFunction *func,
                  int maxNodes ){
    // Allocate the working array
    value = 0.0;
    int maxStrains = 3;
    // Allocate the working array
    work = new TacsScalar[2*maxStrains + 3*maxNodes];

    // Set the pointers into the work array
    strain = &work[0];
    failSens = &work[maxStrains];
    hXptSens = &work[2*maxStrains];
    
  }
  ~HeatFluxIntCtx(){
    delete [] work;
    delete [] hXptSens;
    delete [] failSens;
    delete [] strain;
  }

  // Data to be used for the function computation
  TacsScalar value;
  TacsScalar maxFail;
  TacsScalar ksFailSum;
  TacsScalar *strain;
  TacsScalar *failSens;
  TacsScalar *hXptSens;
  TacsScalar *work;
};

HeatFluxIntegral::HeatFluxIntegral( TACSAssembler *_tacs,
                                    int *_elem_index,
                                    int *_surface_index, 
                                    int _num_elems ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN,
             TACSFunction::SINGLE_STAGE, 0){
  maxNumNodes = _tacs->getMaxElementNodes();
  value = 0.0;
  // Create map
  num_elems = _num_elems;
  for (int i = 0; i < num_elems; i++){
    elem_to_surf.insert(std::make_pair(_elem_index[i], _surface_index[i]));
  }  
}
HeatFluxIntegral::~HeatFluxIntegral(){}

/*
  HeatFluxIntegral function name
*/
const char * HeatFluxIntegral::funcName = "HeatFluxIntegral";

/*
  Return the function name
*/
const char *HeatFluxIntegral::functionName(){
  return funcName;
}

/*
  Retrieve the function value
*/
TacsScalar HeatFluxIntegral::getFunctionValue(){
  return value;
}

/*
  Allocate and return the function-specific context
*/
TACSFunctionCtx *HeatFluxIntegral::createFunctionCtx(){
  return new HeatFluxIntCtx(this, maxNumNodes);
}

/*
  Initialize the internal values stored within the KS function
*/
void HeatFluxIntegral::initEvaluation( EvaluationType ftype ){
  value = 0.0;
}

/*
  Reduce the function values across all MPI processes
*/
void HeatFluxIntegral::finalEvaluation( EvaluationType ftype ){
  TacsScalar temp = value;
  MPI_Allreduce(&temp, &value, 1, TACS_MPI_TYPE,
                MPI_SUM, tacs->getMPIComm());
}

/*
  Initialize the context for either integration or initialization
*/
void HeatFluxIntegral::initThread( const double tcoef,
                                   EvaluationType ftype,
                                   TACSFunctionCtx *fctx ){
  HeatFluxIntCtx *ctx = dynamic_cast<HeatFluxIntCtx*>(fctx);
  if (ctx){
    ctx->value = 0.0;
  }  
}

/*
  Perform the element-wise evaluation of the TACSDisplacementIntegral function.
*/
void HeatFluxIntegral::elementWiseEval( EvaluationType ftype,
                                        TACSElement *element,
                                        int elemNum,
                                        const TacsScalar Xpts[],
                                        const TacsScalar vars[],
                                        const TacsScalar dvars[],
                                        const TacsScalar ddvars[],
                                        TACSFunctionCtx *fctx ){
  HeatFluxIntCtx *ctx = dynamic_cast<HeatFluxIntCtx*>(fctx);
  elem_to_surf_it = elem_to_surf.find(elemNum);
  int surface = elem_to_surf_it->second;
  if (ctx){
    // Get the number of quadrature points for this element
    const int numGauss = element->getNumGaussPts();
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();
        
    double N[numNodes], Na[numNodes], Nb[numNodes];
    TacsScalar q[numDisps-1], strain[numDisps-1];
    int order = 0;
    // Get the constitutive object for this element
    TACSConstitutive *constitutive = element->getConstitutive();
    // Direction vector of the surface/edge
    TacsScalar dir[3];
    if (numDisps == 4){      
      order = cbrt(numNodes);
    }
    else {      
      order = sqrt(numNodes);
    }
    if (constitutive){
      // With the first iteration, find the maximum over the domain
      for ( int i = 0; i < numGauss; i++ ){
        // Get the Gauss points one at a time
        double pt[3];
        double weight = element->getGaussWtsPts(i, pt);
                
        dir[0] = dir[1] = dir[2] = 0.0;
        
        // Get the strain B*u and temperature dT
        // If 3D structure
        if (numDisps == 4){
          ThermoSolid* elem = dynamic_cast<ThermoSolid*>(element);
          if (elem){
            elem->getShapeFunctions(pt, N, Na, Nb);
            computeDirections(dir, surface, order, numDisps,
                              Xpts, Na, Nb);
            elem->getBT(strain, pt, Xpts, vars);
          }         
          CoupledThermoSolidStiffness *con =
            dynamic_cast<CoupledThermoSolidStiffness*>(element->getConstitutive());
          // Compute the heat flux to the surface
          if (con){
            con->calculateConduction(pt, strain, q);
          }
          value += dir[0]*q[0] + dir[1]*q[1] + dir[2]*q[2];
        }
        else {
          ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
          if (elem){
            elem->getShapeFunctions(pt, N, Na, Nb);
            computeDirections(dir, surface, order, numDisps,
                              Xpts, Na, Nb);
            elem->getBT(strain, pt, Xpts, vars);
          }
          CoupledThermoPlaneStressStiffness *con =
            dynamic_cast<CoupledThermoPlaneStressStiffness*>(element->getConstitutive());
          if (con){
            con->calculateConduction(pt, strain, q);
          }
          value += dir[0]*q[0] + dir[1]*q[1];
        }          
        // Add up the contribution from the quadrature
        TacsScalar h = element->getDetJacobian(pt, Xpts);
        h *= weight;        
        ctx->value += h*value;
      }      
    } // end if constitutive    
  }  
}
/*
  For each thread used to evaluate the function, call the
  post-evaluation code once.
*/
void HeatFluxIntegral::finalThread( const double tcoef,
                                    EvaluationType ftype,
                                    TACSFunctionCtx *fctx ){
  HeatFluxIntCtx *ctx = dynamic_cast<HeatFluxIntCtx*>(fctx);
  if (ctx){
    value += ctx->value;
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void HeatFluxIntegral::getElementSVSens( double alpha,
                                         double beta,
                                         double gamma,
                                         TacsScalar *elemSVSens,
                                         TACSElement *element,
                                         int elemNum,
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[],
                                         const TacsScalar dvars[],
                                         const TacsScalar ddvars[],
                                         TACSFunctionCtx *fctx ){
  HeatFluxIntCtx *ctx = dynamic_cast<HeatFluxIntCtx*>(fctx);
  elem_to_surf_it = elem_to_surf.find(elemNum);
  int surface = elem_to_surf_it->second;
  if (ctx){
    int numGauss = element->getNumGaussPts();
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();
    double N[numNodes], Na[numNodes], Nb[numNodes];

    TacsScalar dir[3];
    dir[0] = dir[1] = dir[2] = 0.0;
    // Set the stress/strain arrays
    TacsScalar stress[numDisps-1];
    for ( int i = 0; i < numGauss; i++ ){
      // Get the quadrature point
      double pt[3];
      TacsScalar weight = element->getGaussWtsPts(i, pt);
      TacsScalar h = weight*element->getDetJacobian(pt, Xpts);
            
      dir[0] = dir[1] = dir[2] = 0.0;
      // Add the sensitivity of the heat flux wrt to dT
      // q = {1,1}* H(x) * B * dT
      // dq/d(dT) = {1,1}* H(x) * B
      // Get the derivative of dT at the current point
      if (numDisps == 4){
        int order = cbrt(numNodes);
        ThermoSolid* elem = dynamic_cast<ThermoSolid*>(element);
        if (elem){
          elem->getShapeFunctions(pt, N, Na, Nb);
          computeDirections(dir, surface, order, numDisps, 
                            Xpts, Na, Nb);
        }        
        CoupledThermoSolidStiffness *con =
          dynamic_cast<CoupledThermoSolidStiffness*>(element->getConstitutive());
        // Compute the heat flux to the surface
        if (con){
          con->calculateConduction(pt, dir, stress);
        }        
        if (elem){
          elem->addBTSVSens(elemSVSens, pt, h*alpha, stress, 
                            Xpts, vars);
        }
      }
      else {
        int order = sqrt(numNodes);
        ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
        if (elem){
          elem->getShapeFunctions(pt, N, Na, Nb);
          computeDirections(dir, surface, order, numDisps, Xpts, 
                            Na, Nb);
        }
        CoupledThermoPlaneStressStiffness *con =
          dynamic_cast<CoupledThermoPlaneStressStiffness*>(element->getConstitutive());
        if (con){
          con->calculateConduction(pt, dir, stress); 
        }
        if (elem){
          elem->addBTSVSens(elemSVSens, pt, h*alpha, stress, 
                            Xpts, vars);
        }
      }
    }
  }
}
/*
  Determine the derivative of the function with respect to
  the element nodal locations
*/
void HeatFluxIntegral::getElementXptSens( const double tcoef,
                                          TacsScalar fXptSens[],
                                          TACSElement *element,
                                          int elemNum,
                                          const TacsScalar Xpts[],
                                          const TacsScalar vars[],
                                          const TacsScalar dvars[],
                                          const TacsScalar ddvars[],
                                          TACSFunctionCtx *fctx ){}

/*
  Determine the derivative of the function with respect to
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void HeatFluxIntegral::addElementDVSens( const double tcoef,
                                         TacsScalar *fdvSens,
                                         int numDVs,
                                         TACSElement *element,
                                         int elemNum,
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[],
                                         const TacsScalar dvars[],
                                         const TacsScalar ddvars[],
                                         TACSFunctionCtx *fctx ){
  HeatFluxIntCtx *ctx = dynamic_cast<HeatFluxIntCtx*>(fctx);
  elem_to_surf_it = elem_to_surf.find(elemNum);
  int surface = elem_to_surf_it->second;
  if (ctx){
    int numGauss = element->getNumGaussPts();
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();

    double N[numNodes], Na[numNodes], Nb[numNodes];
    TacsScalar dir[3];
    dir[0] = dir[1] = dir[2] = 0.0;
    
    // Set the stress/strain arrays
    TacsScalar strain[numDisps-1];
    for ( int i = 0; i < numGauss; i++ ){
      // Get the quadrature point
      double pt[3];
      TacsScalar weight = element->getGaussWtsPts(i, pt);
      TacsScalar h = weight*element->getDetJacobian(pt, Xpts);
      
      dir[0] = dir[1] = dir[2] = 0.0;
      // Get the derivative of dT at the current point
      if (numDisps == 4){
        int order = cbrt(numNodes);
        ThermoSolid* elem = dynamic_cast<ThermoSolid*>(element);
        if (elem){
          elem->getShapeFunctions(pt, N, Na, Nb);
          computeDirections(dir, surface, order, numDisps,
                            Xpts, Na, Nb);
          elem->getBT(strain, pt, Xpts, vars);
        }
        CoupledThermoSolidStiffness *con =
          dynamic_cast<CoupledThermoSolidStiffness*>(element->getConstitutive());
        // Compute the heat flux to the surface
        if (con){
          con->addConductionDVSens(pt, strain, tcoef*h, dir,
                                   fdvSens, numDVs);
        }
      }
      else {
        int order = sqrt(numNodes);        
        ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
        if (elem){
          elem->getShapeFunctions(pt, N, Na, Nb);
          computeDirections(dir, surface, order, numDisps, 
                            Xpts, Na, Nb);
          elem->getBT(strain, pt, Xpts, vars);
        }
        CoupledThermoPlaneStressStiffness *con =
          dynamic_cast<CoupledThermoPlaneStressStiffness*>(element->getConstitutive());
        if (con){
          con->addConductionDVSens(pt, strain, tcoef*h, dir,
                                    fdvSens, numDVs);
        }
      }      
    }    
  }
}
/*
  Compute the surface or edge normal
*/
void HeatFluxIntegral::computeDirections( double dir[],
                                          const int surface,
                                          const int order,
                                          const int numDisps,
                                          const TacsScalar Xpts[],
                                          const double Na[],
                                          const double Nb[] ){
  TacsScalar Xa[3], Xb[3];
  Xa[0] = Xa[1] = Xa[2] = 0.0;
  Xb[0] = Xb[1] = Xb[2] = 0.0;

  if (numDisps-1 == 4){
    if (surface < 2){
      const int ii = (order-1)*(surface % 2);
      for ( int kk = 0; kk < order; kk++ ){
        for ( int jj = 0; jj < order; jj++ ){
          const int node = ii + jj*order + kk*order*order;
              
          Xa[0] += Na[jj + kk*order]*Xpts[3*node];
          Xa[1] += Na[jj + kk*order]*Xpts[3*node+1];
          Xa[2] += Na[jj + kk*order]*Xpts[3*node+2];

          Xb[0] += Nb[jj + kk*order]*Xpts[3*node];
          Xb[1] += Nb[jj + kk*order]*Xpts[3*node+1];
          Xb[2] += Nb[jj + kk*order]*Xpts[3*node+2];
        }
      }
    }
    else if (surface < 4){
      const int jj = (order-1)*(surface % 2);
      for ( int kk = 0; kk < order; kk++ ){
        for ( int ii = 0; ii < order; ii++ ){
          const int node = ii + jj*order + kk*order*order;
              
          Xa[0] += Na[ii + kk*order]*Xpts[3*node];
          Xa[1] += Na[ii + kk*order]*Xpts[3*node+1];
          Xa[2] += Na[ii + kk*order]*Xpts[3*node+2];

          Xb[0] += Nb[ii + kk*order]*Xpts[3*node];
          Xb[1] += Nb[ii + kk*order]*Xpts[3*node+1];
          Xb[2] += Nb[ii + kk*order]*Xpts[3*node+2];
        }
      }
    }
    else {
      const int kk = (order-1)*(surface % 2);
      for ( int jj = 0; jj < order; jj++ ){
        for ( int ii = 0; ii < order; ii++ ){
          const int node = ii + jj*order + kk*order*order;
              
          Xa[0] += Na[ii + jj*order]*Xpts[3*node];
          Xa[1] += Na[ii + jj*order]*Xpts[3*node+1];
          Xa[2] += Na[ii + jj*order]*Xpts[3*node+2];

          Xb[0] += Nb[ii + jj*order]*Xpts[3*node];
          Xb[1] += Nb[ii + jj*order]*Xpts[3*node+1];
          Xb[2] += Nb[ii + jj*order]*Xpts[3*node+2];
        }
      }
    }
    // Compute the normal to the element
    Tensor::crossProduct3D(dir, Xa, Xb);
  }
  else {
    if (surface == 0 || surface == 1){
      dir[0] = 0.0;
      dir[1] = 1.0;
    }
    else {
      dir[0] = 1.0;
      dir[1] = 0.0;
    }
  }
}
