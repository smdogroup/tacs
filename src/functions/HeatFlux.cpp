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
#include "ThermoElements.h"
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
  }
  ~HeatFluxIntCtx(){}

  // Data to be used for the function computation
  TacsScalar value;
};

HeatFluxIntegral::HeatFluxIntegral( TACSAssembler *_tacs,
                                    int *_elem_index,
                                    int *_surface_index,
                                    int _num_elems ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN,
             TACSFunction::SINGLE_STAGE, 0){
  maxNumNodes = _tacs->getMaxElementNodes();
  MPI_Comm_rank(_tacs->getMPIComm(), &mpi_rank);
  
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
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();

    TacsScalar q[numDisps-1], strain[numDisps-1];
    int order = 0;
    // Get the constitutive object for this element
    TACSConstitutive *constitutive = element->getConstitutive();

    // Direction vector of the surface/edge
    double dir1[3], dir2[3], base[3];
    if (numDisps == 4){
      order = cbrt(numNodes);
    }
    else {
      order = sqrt(numNodes);
    }

    // The knot locations for the basis functions
    double knots[order];
    // Set the knot locations
    if (order == 2){
      knots[0] = -1.0;
      knots[1] = 1.0;
    }
    else if (order == 3){
      knots[0] = -1.0;
      knots[1] = 0.0;
      knots[2] = 1.0;
    }
    else {
      // Set a co-sine spacing for the knot locations
      for ( int k = 0; k < order; k++ ){
        knots[k] = -cos(M_PI*k/(order-1));
      }
    }

    if (constitutive){
      // Get the quadrature points and weights
      const double *gaussPts, *gaussWts;
      FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);
      // With the first iteration, find the maximum over the domain
      if (numDisps == 4){
        for ( int m = 0; m < order; m++ ){
          for ( int n = 0; n < order; n++ ){
            // Set the quadrature point
            double pt[3];

            // Compute the base point and direction of the surface tangent in
            // parameter space
            dir1[0] = dir1[1] = dir1[2] = 0.0;
            dir2[0] = dir2[1] = dir2[2] = 0.0;
            base[0] = base[1] = base[2] = 0.0;
            
            if (surface == 0 || surface == 1){
              base[0] = -1.0 + 2.0*(surface % 2);
              dir1[1] = -1.0 + 2.0*(surface % 2);
              dir2[2] = 1.0;
            }
            else if (surface == 2 || surface == 3) {
              base[1] = -1.0 + 2.0*(surface % 2);
              dir1[0] = -1.0 + 2.0*(surface % 2);
              dir2[2] = -1.0;
            }
            else {
              base[2] = -1.0 + 2.0*(surface % 2);
              dir1[0] = -1.0 + 2.0*(surface % 2);
              dir2[1] = 1.0;
            }

            pt[0] = gaussPts[n]*dir1[0] + gaussPts[m]*dir2[0] + base[0];
            pt[1] = gaussPts[n]*dir1[1] + gaussPts[m]*dir2[1] + base[1];
            pt[2] = gaussPts[n]*dir1[2] + gaussPts[m]*dir2[2] + base[2];
            
            // Get magnitude and direction of the normal to surface
            TacsScalar normal[3];
            TacsScalar tn = computeDirections3D(pt, knots, dir1, dir2, surface,
                                                order, Xpts, normal); 
            // Multiply by the quadrature weight
            tn *= gaussWts[n]*gaussWts[m];
            
            ThermoSolid *elem = dynamic_cast<ThermoSolid*>(element);
            if (elem){
              elem->getBT(strain, pt, Xpts, vars);
            }
            CoupledThermoSolidStiffness *con =
              dynamic_cast<CoupledThermoSolidStiffness*>(element->getConstitutive());
            // Compute the heat flux to the surface
            if (con){
              //con->calculateConduction(pt, strain, q);
              con->heatflux(pt, normal, strain, &value);
            }
            // value = (normal[0]*q[0] + normal[1]*q[1] + normal[2]*q[2]);
            // Add up the contribution from the quadrature
            ctx->value += tn*value;
          } // end int n
        } // end int m
      }
      else {
        for ( int n = 0; n < order; n++ ){
          double pt[2];
          
          // Compute the base point and direction of the surface tangent in
          // parameter space
          dir1[0] = dir1[1] = 0.0;
          base[0] = base[1] = 0.0;
          if (surface == 0 || surface == 1){
            dir1[1] = -1.0 + 2.0*(surface % 2);
            base[0] = -1.0 + 2.0*(surface % 2);
          }
          else {
            dir1[0] = -1.0 + 2.0*(surface % 2);
            base[1] = -1.0 + 2.0*(surface % 2);
          }
         
          // Get the integration gauss points on the element surface
          pt[0] = gaussPts[n]*dir1[0] + base[0];
          pt[1] = gaussPts[n]*dir1[1] + base[1];

          // Compute the normal direction and element tangent magnitude
          TacsScalar normal[3];
          TacsScalar tn = computeDirections2D(pt, knots, dir1, surface, order,
                                              Xpts, normal);
          
          // Multiply by the quadrature weight
          tn *= gaussWts[n];

          ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
          if (elem){            
            elem->getBT(strain, pt, Xpts, vars);
          }
          CoupledThermoPlaneStressStiffness *con =
            dynamic_cast<CoupledThermoPlaneStressStiffness*>(element->getConstitutive());
          // Compute the heat flux to the edge
          if (con){
            //con->calculateConduction(pt, strain, q);
            con->heatflux(pt, normal, strain, &value);
          }
          // value = (normal[0]*q[0] + normal[1]*q[1]);
          // Add up the contribution from the quadrature
          ctx->value += tn*value;
        } // end int n
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

  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->numVariables();
  memset(elemSVSens, 0, numVars*sizeof(TacsScalar));
  
  if (ctx){
    int numGauss = element->getNumGaussPts();
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();
    
    int order = 0;
    // Get the constitutive object for this element
    TACSConstitutive *constitutive = element->getConstitutive();
    // Direction vector of the surface/edge
    double dir1[3], dir2[3], base[3];
    if (numDisps == 4){
      order = cbrt(numNodes);
    }
    else {
      order = sqrt(numNodes);
    }
    // The knot locations for the basis functions
    double knots[order];
    // Set the knot locations
    if (order == 2){
      knots[0] = -1.0;
      knots[1] = 1.0;
    }
    else if (order == 3){
      knots[0] = -1.0;
      knots[1] = 0.0;
      knots[2] = 1.0;
    }
    else {
      // Set a co-sine spacing for the knot locations
      for ( int k = 0; k < order; k++ ){
        knots[k] = -cos(M_PI*k/(order-1));
      }
    }

    TacsScalar q[numDisps-1];    
    if (constitutive){
      // Get the quadrature points and weights
      const double *gaussPts, *gaussWts;
      FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);
      // With the first iteration, find the maximum over the domain
      if (numDisps == 4){
        for ( int m = 0; m < order; m++ ){
          for ( int n = 0; n < order; n++ ){
            // Set the quadrature point
            double pt[3];
            // Compute the base point and direction of the surface tangent in
            // parameter space
            dir1[0] = dir1[1] = dir1[2] = 0.0;
            dir2[0] = dir2[1] = dir2[2] = 0.0;
            base[0] = base[1] = base[2] = 0.0;
            
            if (surface == 0 || surface == 1){
              base[0] = -1.0 + 2.0*(surface % 2);
              dir1[1] = -1.0 + 2.0*(surface % 2);
              dir2[2] = 1.0;
            }
            else if (surface == 2 || surface == 3) {
              base[1] = -1.0 + 2.0*(surface % 2);
              dir1[0] = -1.0 + 2.0*(surface % 2);
              dir2[2] = -1.0;
            }
            else {
              base[2] = -1.0 + 2.0*(surface % 2);
              dir1[0] = -1.0 + 2.0*(surface % 2);
              dir2[1] = 1.0;
            }

            pt[0] = gaussPts[n]*dir1[0] + gaussPts[m]*dir2[0] + base[0];
            pt[1] = gaussPts[n]*dir1[1] + gaussPts[m]*dir2[1] + base[1];
            pt[2] = gaussPts[n]*dir1[2] + gaussPts[m]*dir2[2] + base[2];
            // Get magnitude and direction of the normal to surface
            TacsScalar normal[3];
            TacsScalar tn = computeDirections3D(pt, knots, dir1, dir2, surface,
                                                order, Xpts, normal); 
            // Multiply by the quadrature weight
            tn *= gaussWts[n]*gaussWts[m];

            // Add the sensitivity of the heat flux wrt to dT
            // q = {nx,ny,nz}* H(x) * B * dT
            // dq/d(dT) = {nx,ny,nz}* H(x) * B

            CoupledThermoSolidStiffness *con =
              dynamic_cast<CoupledThermoSolidStiffness*>(element->getConstitutive());
            // Compute the heat flux to the surface
            if (con){
              //con->calculateConduction(pt, normal, q);
              con->heatfluxStrainSens(pt, normal, q);
            }
            ThermoSolid *elem = dynamic_cast<ThermoSolid*>(element);
            if (elem){
              elem->addBTSVSens(elemSVSens, pt, tn*alpha, q,
                                Xpts, vars);
            }
          }
        }
      }
      else {
        for ( int n = 0; n < order; n++ ){
          double pt[2];
          
          // Compute the base point and direction of the surface tangent in
          // parameter space
          dir1[0] = dir1[1] = 0.0;
          base[0] = base[1] = 0.0;
          if (surface == 0 || surface == 1){
            dir1[1] = -1.0 + 2.0*(surface % 2);
            base[0] = -1.0 + 2.0*(surface % 2);
          }
          else {
            dir1[0] = -1.0 + 2.0*(surface % 2);
            base[1] = -1.0 + 2.0*(surface % 2);
          }
         
          // Get the integration gauss points on the element surface
          pt[0] = gaussPts[n]*dir1[0] + base[0];
          pt[1] = gaussPts[n]*dir1[1] + base[1];

          // Compute the normal direction and element tangent magnitude
          TacsScalar normal[3];
          TacsScalar tn = computeDirections2D(pt, knots, dir1, surface, order,
                                              Xpts, normal);

          // Multiply by the quadrature weight
          tn *= gaussWts[n];
          
          CoupledThermoPlaneStressStiffness *con =
            dynamic_cast<CoupledThermoPlaneStressStiffness*>(element->getConstitutive());
          // Compute the heat flux to the edge
          if (con){
            //con->calculateConduction(pt, normal, q);
            con->heatfluxStrainSens(pt, normal, q);
          }
          ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
          if (elem){
            elem->addBTSVSens(elemSVSens, pt, tn*alpha, q,
                              Xpts, vars);
          }
        }
      }
    } // end if constitutive
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
                                          TACSFunctionCtx *fctx ){
  int numNodes = element->numNodes();
  memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));
}

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
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();

    TacsScalar q[numDisps-1], strain[numDisps-1];
    int order = 0;
    // Get the constitutive object for this element
    TACSConstitutive *constitutive = element->getConstitutive();

    // Direction vector of the surface/edge
    double dir1[3], dir2[3], base[3];
    if (numDisps == 4){
      order = cbrt(numNodes);
    }
    else {
      order = sqrt(numNodes);
    }

    // The knot locations for the basis functions
    double knots[order];
    // Set the knot locations
    if (order == 2){
      knots[0] = -1.0;
      knots[1] = 1.0;
    }
    else if (order == 3){
      knots[0] = -1.0;
      knots[1] = 0.0;
      knots[2] = 1.0;
    }
    else {
      // Set a co-sine spacing for the knot locations
      for ( int k = 0; k < order; k++ ){
        knots[k] = -cos(M_PI*k/(order-1));
      }
    }
    if (constitutive){
      // Get the quadrature points and weights
      const double *gaussPts, *gaussWts;
      FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);
      // With the first iteration, find the maximum over the domain
      if (numDisps == 4){
        for ( int m = 0; m < order; m++ ){
          for ( int n = 0; n < order; n++ ){
            // Set the quadrature point
            double pt[3];
            // Compute the base point and direction of the surface tangent in
            // parameter space
            dir1[0] = dir1[1] = dir1[2] = 0.0;
            dir2[0] = dir2[1] = dir2[2] = 0.0;
            base[0] = base[1] = base[2] = 0.0;
            
            if (surface == 0 || surface == 1){
              base[0] = -1.0 + 2.0*(surface % 2);
              dir1[1] = -1.0 + 2.0*(surface % 2);
              dir2[2] = 1.0;
            }
            else if (surface == 2 || surface == 3) {
              base[1] = -1.0 + 2.0*(surface % 2);
              dir1[0] = -1.0 + 2.0*(surface % 2);
              dir2[2] = -1.0;
            }
            else {
              base[2] = -1.0 + 2.0*(surface % 2);
              dir1[0] = -1.0 + 2.0*(surface % 2);
              dir2[1] = 1.0;
            }

            pt[0] = gaussPts[n]*dir1[0] + gaussPts[m]*dir2[0] + base[0];
            pt[1] = gaussPts[n]*dir1[1] + gaussPts[m]*dir2[1] + base[1];
            pt[2] = gaussPts[n]*dir1[2] + gaussPts[m]*dir2[2] + base[2];

            // Get magnitude and direction of the normal to surface
            TacsScalar normal[3];
            TacsScalar tn = computeDirections3D(pt, knots, dir1, dir2, surface,
                                                order, Xpts, normal); 
            // Multiply by the quadrature weight
            tn *= gaussWts[n]*gaussWts[m];
            
            // Get the derivative of dT at the current point
            ThermoSolid* elem = dynamic_cast<ThermoSolid*>(element);
            if (elem){
              elem->getBT(strain, pt, Xpts, vars);
            }
            CoupledThermoSolidStiffness *con =
              dynamic_cast<CoupledThermoSolidStiffness*>(element->getConstitutive());
            // Compute the heat flux to the surface
            if (con){
              // con->addConductionDVSens(pt, strain, tcoef*tn, normal,
              //                          fdvSens, numDVs);
              con->addHeatFluxDVSens(pt, normal, strain, tcoef*tn,
                                     fdvSens, numDVs);
            }
          }
        }
      }
      else{
        for ( int n = 0; n < order; n++ ){
          double pt[2];
          // Compute the base point and direction of the surface tangent in
          // parameter space
          dir1[0] = dir1[1] = 0.0;
          base[0] = base[1] = 0.0;
          if (surface == 0 || surface == 1){
            dir1[1] = -1.0 + 2.0*(surface % 2);
            base[0] = -1.0 + 2.0*(surface % 2);
          }
          else {
            dir1[0] = -1.0 + 2.0*(surface % 2);
            base[1] = -1.0 + 2.0*(surface % 2);
          }
         
          // Get the integration gauss points on the element surface
          pt[0] = gaussPts[n]*dir1[0] + base[0];
          pt[1] = gaussPts[n]*dir1[1] + base[1];

          // Compute the normal direction and element tangent magnitude
          TacsScalar normal[3];
          TacsScalar tn = computeDirections2D(pt, knots, dir1, surface, order,
                                              Xpts, normal);
          // Multiply by the quadrature weight
          tn *= gaussWts[n];
         
          ThermoQuad* elem = dynamic_cast<ThermoQuad*>(element);
          if (elem){            
            elem->getBT(strain, pt, Xpts, vars);
          }
          CoupledThermoPlaneStressStiffness *con =
            dynamic_cast<CoupledThermoPlaneStressStiffness*>(element->getConstitutive());
          if (con){
            // con->addConductionDVSens(pt, strain, tcoef*tn, normal,
            //                          fdvSens, numDVs);
            con->addHeatFluxDVSens(pt, normal, strain, tcoef*tn,
                                     fdvSens, numDVs);
          }
        }
      }
    }    
  }
}

/*
  Compute the surface or edge normal
*/
TacsScalar HeatFluxIntegral::computeDirections2D( const double pt[],
                                                  const double knots[],
                                                  const double dir[],
                                                  const int surface,
                                                  const int order,
                                                  const TacsScalar Xpts[],
                                                  TacsScalar n[] ){
  // Evaluate the Lagrange basis in each direction
  double na[order], nb[order];
  double dna[order], dnb[order];
  FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
  FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);
  
  // Calcualte the Jacobian at the current point
  const TacsScalar *x = Xpts;
  TacsScalar Xd[4] = {0.0, 0.0, 0.0, 0.0};
  for ( int j = 0; j < order; j++ ){
    for ( int i = 0; i < order; i++ ){
      Xd[0] += x[0]*dna[i]*nb[j];
      Xd[1] += x[0]*na[i]*dnb[j];
      
      Xd[2] += x[1]*dna[i]*nb[j];
      Xd[3] += x[1]*na[i]*dnb[j];
      x += 3;
    }
  }
  
  // Compute the derivative along each direction
  TacsScalar tx = Xd[0]*dir[0] + Xd[2]*dir[1];
  TacsScalar ty = Xd[1]*dir[0] + Xd[3]*dir[1];
  
  // Compute the magnitude of the tangent vector
  TacsScalar tn = sqrt(tx*tx + ty*ty);

  // Compute the normal vector (outward facing from the element edge)
  n[0] = ty/tn;
  n[1] = -tx/tn;
  
  return tn;
}

TacsScalar HeatFluxIntegral::computeDirections3D( const double pt[],
                                                  const double knots[],
                                                  const double dir1[],
                                                  const double dir2[],
                                                  const int surface,
                                                  const int order,
                                                  const TacsScalar Xpts[],
                                                  TacsScalar n[] ){
    
  // Evaluate the Lagrange basis in each direction
  double na[order], nb[order], nc[order];
  double dna[order], dnb[order], dnc[order];
  FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
  FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);
  FElibrary::lagrangeSFKnots(nc, dnc, pt[2], knots, order);
  
  // Calcualte the Jacobian at the current point
  const TacsScalar *x = Xpts;
  TacsScalar Xd[9] = {0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0};
  
  for ( int k  = 0; k < order; k++ ){
    for ( int j = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++ ){
        Xd[0] += x[0]*dna[i]*nb[j]*nc[k];
        Xd[1] += x[0]*na[i]*dnb[j]*nc[k];
        Xd[2] += x[0]*na[i]*nb[j]*dnc[k];
        
        Xd[3] += x[1]*dna[i]*nb[j]*nc[k];
        Xd[4] += x[1]*na[i]*dnb[j]*nc[k];
        Xd[5] += x[1]*na[i]*nb[j]*dnc[k];

        Xd[6] += x[2]*dna[i]*nb[j]*nc[k];
        Xd[7] += x[2]*na[i]*dnb[j]*nc[k];
        Xd[8] += x[2]*na[i]*nb[j]*dnc[k];

        x += 3;
      }
    }
  }
  // Compute the first tangent direction
  TacsScalar t1[3];
  t1[0] = Xd[0]*dir1[0] + Xd[1]*dir1[1] + Xd[2]*dir1[2];
  t1[1] = Xd[3]*dir1[0] + Xd[4]*dir1[1] + Xd[5]*dir1[2];
  t1[2] = Xd[6]*dir1[0] + Xd[7]*dir1[1] + Xd[8]*dir1[2];

  // Compute the second tangent direction
  TacsScalar t2[3];
  t2[0] = Xd[0]*dir2[0] + Xd[1]*dir2[1] + Xd[2]*dir2[2];
  t2[1] = Xd[3]*dir2[0] + Xd[4]*dir2[1] + Xd[5]*dir2[2];
  t2[2] = Xd[6]*dir2[0] + Xd[7]*dir2[1] + Xd[8]*dir2[2];

  // Compute the normal to the element
  Tensor::crossProduct3D(n, t1, t2);

  TacsScalar tn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

  n[0] /= tn;
  n[1] /= tn;
  n[2] /= tn;

  return tn;
}

/*
  Compute the shape functions and their derivatives w.r.t. the
  parametric element location
*/
void HeatFluxIntegral::getShapeFunctions( const double pt[],
                                          const double knots[],
                                          const int order,
                                          double N[],
                                          double Na[], double Nb[] ){
    double na[order], nb[order];
    double dna[order], dnb[order];
    FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
    FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);
    for ( int j = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++ ){
        N[0] = na[i]*nb[j];
        Na[0] = dna[i]*nb[j];
        Nb[0] = na[i]*dnb[j];
        N++;
        Na++;  Nb++;
      }
    }
  }
