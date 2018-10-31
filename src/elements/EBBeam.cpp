#include "FElibrary.h"
#include "TensorToolbox.h"
#include "EBBeam.h"

/*
  An Euler-Bernoulli beam element

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

/*
  Create the Euler-Bernoulli beam. Set the reference axis.
*/
EBBeam::EBBeam( EBStiffness * _beamStiff){
  beamStiff = _beamStiff;
  beamStiff->incref();

  ref_dir_type = beamStiff->ref_type;

  // Initialize the transformation
  ref_dir[0] = beamStiff->ref_dir[0];
  ref_dir[1] = beamStiff->ref_dir[1];
  ref_dir[2] = beamStiff->ref_dir[2];
  Tensor::normalize3D(ref_dir);

  numGauss = 2;
  gaussWts = FElibrary::gaussWts2;
  gaussPts = FElibrary::gaussPts2;
}

EBBeam::~EBBeam(){
  beamStiff->decref();
}

/*
  Retrieve information about the names of the element variables
*/
const char * EBBeam:: elementName() const {
  return elemName;
}

const char * EBBeam::displacementName( int i ) const {
  if (i >= 0 && i < NUM_DISPS ){
    return dispNames[i];
  }
  return NULL;
}

const char * EBBeam::stressName( int i ) const{
  if (i >= 0 && i < NUM_STRESSES ){
    return stressNames[i];
  }
  return NULL;
}

const char * EBBeam::strainName( int i ) const{
  if (i >= 0 && i < NUM_STRESSES ){
    return strainNames[i];
  }
  return NULL;
}

const char * EBBeam::extraName( int i ) const{
  if (i >= 0 && i < NUM_EXTRAS ){
    return extraNames[i];
  }
  return NULL;
}

/*
  Retrieve the numbers of displacements, nodes, stress and variables
*/
int EBBeam::numDisplacements() { return NUM_DISPS; }
int EBBeam::numStresses() const { return NUM_STRESSES; }
int EBBeam::numNodes() { return NUM_NODES; }
int EBBeam::numVariables() const { return NUM_VARIABLES; }
int EBBeam::numExtras() const { return NUM_EXTRAS; }

enum ElementType EBBeam::getElementType(){ return TACS_EULER_BEAM; }

/*
  The element name, variable, stress and strain names.
*/
const char * EBBeam::elemName = "EBBeam";

const char * EBBeam::dispNames[] = { "u0", "v0", "w0",
				     "rotx", "roty", "rotz" };

const char * EBBeam::stressNames[] = { "sx0", "sx1", "sy1", "st" };

const char * EBBeam::strainNames[] = { "ex0", "ex1", "ey1", "et" };

const char * EBBeam::extraNames[] = { "lambda", "buckling" };

// Return true if the element has this design variable
int EBBeam::ownsDesignVar( const int dvNum ) const {
  return beamStiff->ownsDesignVar(dvNum);
}

// Return the number of design variables owned by this element
int EBBeam::getNumDesignVars() const {
  return beamStiff->getNumDesignVars();
}

// Return a constant array of the design variable numbers
int EBBeam::getDesignVarNums( int * dvNums,
			      int * dvIndex, int dvLen ) const {
  return beamStiff->getDesignVarNums(dvNums, dvIndex, dvLen);
}

// Set the values of the design variables
void EBBeam::setDesignVars( const TacsScalar dvs[], int numDVs ){
  beamStiff->setDesignVars(dvs, numDVs);
}

// Get the design variable values -- populate the array dvs[]
void EBBeam::getDesignVars( TacsScalar dvs[], int numDVs ) const {
  beamStiff->getDesignVars(dvs, numDVs);
}

// Populate the arrays lowerBound[] and upperBound[]
void EBBeam::getDesignVarRange( TacsScalar lowerBound[],
				TacsScalar upperBound[],
				int numDVs ) const {
  beamStiff->getDesignVarRange(lowerBound, upperBound, numDVs);
}

/*
  Compute the local transformation from the global axis to the local
  reference frame for the beam element.

  If the reference direction is the strong axis:
  ----------------------------------------------

  The first direction is along the axis of the beam:

  t[0] = (Xpts[3:6] - Xpts[0:3])/||Xpts[3:6] - Xpts[0:3]||

  The second direction is perpendicular to the reference axis and
  the axial direction:

  t[3] = (ref_dir x t[0])/||ref_dir x t[0]||

  The third direction is perpendicular to both t[0] and t[3]

  t[6] = t[0] x t[3]

  If the reference direction is the weak axis:
  --------------------------------------------

  The first direction is along the axis of the beam:

  t[0] = (Xpts[3:6] - Xpts[0:3])/||Xpts[3:6] - Xpts[0:3]||

  The third direction is perpendicular to the reference axis and the
  direction along the length of the beam

  t[6] = (t[0] x ref_dir)/||t[0] x ref_dir||

  The second direction is perpendicular to both the t[6] and t[0]

  t[3] = t[6] x t[0]
*/
TacsScalar EBBeam::computeTransform( TacsScalar t[],
                                     const TacsScalar Xpts[] ){
  t[0] = Xpts[3] - Xpts[0];
  t[1] = Xpts[4] - Xpts[1];
  t[2] = Xpts[5] - Xpts[2];
  TacsScalar L = Tensor::normalize3D(t);

  if (Tensor::dot3D(ref_dir, &t[0]) == 1.0){
    fprintf(stderr,
            "EBBeam: Error, reference axis and beam axis are parallel\n");
  }
  else {
    if (ref_dir_type == STRONG_AXIS){
      Tensor::crossProduct3D(&t[3], ref_dir, &t[0]);
      Tensor::normalize3D(&t[3]);
      Tensor::crossProduct3D(&t[6], &t[0], &t[3]);
    }
    else { // ref_dir_type == WEAK_AXIS
      Tensor::crossProduct3D(&t[6], &t[0], ref_dir);
      Tensor::normalize3D(&t[6]);
      Tensor::crossProduct3D(&t[3], &t[6], &t[0]);
    }
  }
  
 // fprintf(stderr, "Transformation Matrix\n");
 // fprintf(stderr, "\t%f %f %f\n", t[0], t[1], t[2]);
 // fprintf(stderr, "\t%f %f %f\n", t[3], t[4], t[5]);
 // fprintf(stderr, "\t%f %f %f\n", t[6], t[7], t[8]);
  return L;
}

/*
  Compute the derivative of the transform with respect to the given
  nodal coordinate.
*/
TacsScalar EBBeam::computeTransformXptSens( TacsScalar t[],
                                            TacsScalar tSens[],
                                            TacsScalar * LSens,
                                            int component,
                                            const TacsScalar Xpts[] ){
  t[0] = Xpts[3] - Xpts[0];
  t[1] = Xpts[4] - Xpts[1];
  t[2] = Xpts[5] - Xpts[2];

  tSens[0] = tSens[1] = tSens[2] = 0.0;
  if (component < 3){
    tSens[component] = -1.0;
  }
  else {
    tSens[component-3] = 1.0;
  }

  TacsScalar L = Tensor::normalize3DSens(LSens, t, tSens);

  if (ref_dir_type == STRONG_AXIS){
    Tensor::crossProduct3D(&t[3], ref_dir, &t[0]);
    Tensor::crossProduct3D(&tSens[3], ref_dir, &tSens[0]);
    TacsScalar temp;
    Tensor::normalize3DSens(&temp, &t[3], &tSens[3]);
    Tensor::crossProduct3DSens(&t[6], &tSens[6],
                               &t[0], &t[3], &tSens[0], &tSens[3]);
  }
  else { // ref_dir_type == WEAK_AXIS
    Tensor::crossProduct3D(&t[6], &t[0], ref_dir);
    Tensor::crossProduct3D(&tSens[6], &tSens[0], ref_dir);
    TacsScalar temp;
    Tensor::normalize3DSens(&temp, &t[6], &tSens[6]);
    Tensor::crossProduct3DSens(&t[3], &tSens[3],
                               &t[6], &t[0], &tSens[6], &tSens[0]);
  }

  return L;
}

/*
  Transform the global element variables gvars to the set of local
  element variables vars.
*/
void EBBeam::transformVarsGlobalToLocal( TacsScalar vars[],
                                         const TacsScalar gvars[],
                                         const TacsScalar t[] ){
  for ( int i = 0; i < 2*NUM_NODES; i++ ){
    vars[0] = t[0]*gvars[0] + t[1]*gvars[1] + t[2]*gvars[2];
    vars[1] = t[3]*gvars[0] + t[4]*gvars[1] + t[5]*gvars[2];
    vars[2] = t[6]*gvars[0] + t[7]*gvars[1] + t[8]*gvars[2];

    vars += 3;
    gvars += 3;
  }
}

/*
  Compute the local variables in-place such that:

  local vars = transform * global vars
*/
void EBBeam::transformVarsGlobalToLocal( TacsScalar vars[],
                                         const TacsScalar t[] ){
  TacsScalar a[3];
  for ( int i = 0; i < 2*NUM_NODES; i++ ){
    a[0] = t[0]*vars[0] + t[1]*vars[1] + t[2]*vars[2];
    a[1] = t[3]*vars[0] + t[4]*vars[1] + t[5]*vars[2];
    a[2] = t[6]*vars[0] + t[7]*vars[1] + t[8]*vars[2];

    vars[0] = a[0];
    vars[1] = a[1];
    vars[2] = a[2];

    vars += 3;
  }
}

/*
  Compute the global residual in-place such that:

  global residual = transform^{T} * local residual
*/
void EBBeam::transformResLocalToGlobal( TacsScalar res[],
                                        const TacsScalar t[] ){
  TacsScalar a[3];
  for ( int i = 0; i < 2*NUM_NODES; i++ ){
    a[0] = t[0]*res[0] + t[3]*res[1] + t[6]*res[2];
    a[1] = t[1]*res[0] + t[4]*res[1] + t[7]*res[2];
    a[2] = t[2]*res[0] + t[5]*res[1] + t[8]*res[2];

    res[0] = a[0];
    res[1] = a[1];
    res[2] = a[2];

    res += 3;
  }
}

/*
  Compute the derivative of the residual in place

  global residual sensitivity = transformSens^{T} local residual +
  trasnform^{T} * local residual sensitivity
*/
void EBBeam::transformResLocalToGlobalSens( TacsScalar res[],
                                            const TacsScalar resSens[],
                                            const TacsScalar t[],
                                            const TacsScalar tSens[] ){
  TacsScalar a[3];
  for ( int i = 0; i < 2*NUM_NODES; i++ ){
    a[0] = tSens[0]*res[0] + tSens[3]*res[1] + tSens[6]*res[2];
    a[1] = tSens[1]*res[0] + tSens[4]*res[1] + tSens[7]*res[2];
    a[2] = tSens[2]*res[0] + tSens[5]*res[1] + tSens[8]*res[2];

    res[0] = a[0] + t[0]*resSens[0] + t[3]*resSens[1] + t[6]*resSens[2];
    res[1] = a[1] + t[1]*resSens[0] + t[4]*resSens[1] + t[7]*resSens[2];
    res[2] = a[2] + t[2]*resSens[0] + t[5]*resSens[1] + t[8]*resSens[2];

    res += 3;
    resSens += 3;
  }
}

/*
  Compute the global stiffness matrix in-place such that:

  global mat = transform^{T} * mat * transform
*/
void EBBeam::transformStiffnessMat( TacsScalar mat[],
                                    const TacsScalar t[] ){

  for ( int i = 0; i < 2*NUM_NODES; i++ ){
    for ( int j = 0; j < 2*NUM_NODES; j++ ){
      TacsScalar a[9], b[9];

      for ( int ii = 0; ii < 3; ii++ ){
        for ( int jj = 0; jj < 3; jj++ ){
          a[3*ii + jj] = mat[(3*i+ii)*NUM_VARIABLES + 3*j+jj];
        }
      }

      // Compute b = tramsform^{T} * a
      // b = [t[0]  t[3]  t[6]][a[0]  a[1]  a[2]]
      // .   [t[1]  t[4]  t[7]][a[3]  a[4]  a[5]]
      // .   [t[2]  t[5]  t[8]][a[6]  a[7]  a[8]]

      b[0] = t[0]*a[0] + t[3]*a[3] + t[6]*a[6];
      b[1] = t[0]*a[1] + t[3]*a[4] + t[6]*a[7];
      b[2] = t[0]*a[2] + t[3]*a[5] + t[6]*a[8];

      b[3] = t[1]*a[0] + t[4]*a[3] + t[7]*a[6];
      b[4] = t[1]*a[1] + t[4]*a[4] + t[7]*a[7];
      b[5] = t[1]*a[2] + t[4]*a[5] + t[7]*a[8];

      b[6] = t[2]*a[0] + t[5]*a[3] + t[8]*a[6];
      b[7] = t[2]*a[1] + t[5]*a[4] + t[8]*a[7];
      b[8] = t[2]*a[2] + t[5]*a[5] + t[8]*a[8];

      // Compute a = b * transform
      // a = [b[0]  b[1]  b[2]][t[0]  t[1]  t[2]]
      // .   [b[3]  b[4]  b[5]][t[3]  t[4]  t[5]]
      // .   [b[6]  b[7]  b[8]][t[6]  t[7]  t[8]]

      a[0] = b[0]*t[0] + b[1]*t[3] + b[2]*t[6];
      a[3] = b[3]*t[0] + b[4]*t[3] + b[5]*t[6];
      a[6] = b[6]*t[0] + b[7]*t[3] + b[8]*t[6];

      a[1] = b[0]*t[1] + b[1]*t[4] + b[2]*t[7];
      a[4] = b[3]*t[1] + b[4]*t[4] + b[5]*t[7];
      a[7] = b[6]*t[1] + b[7]*t[4] + b[8]*t[7];

      a[2] = b[0]*t[2] + b[1]*t[5] + b[2]*t[8];
      a[5] = b[3]*t[2] + b[4]*t[5] + b[5]*t[8];
      a[8] = b[6]*t[2] + b[7]*t[5] + b[8]*t[8];

      for ( int ii = 0; ii < 3; ii++ ){
        for ( int jj = 0; jj < 3; jj++ ){
          mat[(3*i+ii)*NUM_VARIABLES + 3*j+jj] = a[3*ii + jj];
        }
      }
    }
  }
}

/*
  Assemble the element residual associated with the given design
  variables and elements.
*/
void EBBeam::getRes( TacsScalar * res,
		     const TacsScalar vars[],
		     const TacsScalar Xpts[] ){
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2*NUM_NODES], Nahp[2*NUM_NODES], Naahp[2*NUM_NODES];
  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  TacsScalar transform[9];
  TacsScalar elemVars[NUM_VARIABLES];

  TacsScalar L = computeTransform(transform, Xpts);
  TacsScalar h = 2.0/L;
  TacsScalar hinv = 0.5*L;

  memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));
  transformVarsGlobalToLocal(elemVars, vars, transform);

  for ( int n = 0; n < numGauss; n++ ){
    FElibrary::lagrangeSF(N, Na, gaussPts[n], 2);
    FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[n]);

    evalStrain(h, Na, Nahp, Naahp, elemVars, strain);
    evalBvec(h, Na, Nahp, Naahp, B);

    beamStiff->calculateStress(&gaussPts[n], strain, stress);

    for ( int i = 0; i < NUM_NODES; i++ ){
      for ( int ii = 0; ii < NUM_DISPS; ii++ ){
      	int row = ii + i*NUM_DISPS;

      	res[row] +=
                hinv*gaussWts[n]*strain_product(&B[NUM_STRESSES*row], stress);
      }
    }
  }

  transformResLocalToGlobal(res, transform);
}

/*
  Assemble the stiffness matrix for the Euler-Bernoulli beam element.
*/
void EBBeam::getMat( TacsScalar * mat, TacsScalar * res,
		     const TacsScalar vars[],
		     const TacsScalar Xpts[],
		     MatrixOrientation matOr ){

  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2*NUM_NODES], Nahp[2*NUM_NODES], Naahp[2*NUM_NODES];

  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];
  TacsScalar BiStress[NUM_STRESSES];

  TacsScalar transform[9];
  TacsScalar elemVars[NUM_VARIABLES];

  TacsScalar L = computeTransform(transform, Xpts);
  TacsScalar h = 2.0/L;
  TacsScalar hinv = 0.5*L;

  memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));
  memset(mat, 0, NUM_VARIABLES*NUM_VARIABLES*sizeof(TacsScalar));
  transformVarsGlobalToLocal(elemVars, vars, transform);

  for ( int n = 0; n < numGauss; n++ ){
    FElibrary::lagrangeSF(N, Na, gaussPts[n], 2);
    FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[n]);

    evalStrain(h, Na, Nahp, Naahp, elemVars, strain);
    evalBvec(h, Na, Nahp, Naahp, B);

    beamStiff->calculateStress(&gaussPts[n], strain, stress);

    for ( int i = 0; i < NUM_NODES; i++ ){
      for ( int ii = 0; ii < NUM_DISPS; ii++ ){
	int row = ii + i*NUM_DISPS;

	res[row] +=
          hinv*gaussWts[n]*strain_product(&B[NUM_STRESSES*row], stress);

        beamStiff->calculateStress(&gaussPts[n],
                                   &B[NUM_STRESSES*row], BiStress);

        for ( int j = 0; j < NUM_NODES; j++ ){
          for ( int jj = 0; jj < NUM_DISPS; jj++ ){
	    // generate another B vector
	    int col = jj + j*NUM_DISPS;

	    // The regular element matrix
	    mat[col + row*NUM_VARIABLES] +=
              hinv*gaussWts[n]*strain_product(BiStress, &B[NUM_STRESSES*col]);
	  }
	}
      }
    }
  }

  transformStiffnessMat(mat, transform);
  transformResLocalToGlobal(res, transform);
}

void EBBeam::getMatType( ElementMatrixType matType, TacsScalar timeFactor,
			 TacsScalar * mat, const TacsScalar vars[],
			 const TacsScalar Xpts[],
                         MatrixOrientation matOr ){
  memset(mat, 0, NUM_VARIABLES*NUM_VARIABLES*sizeof(TacsScalar));
}

void EBBeam::addResidual( double time, TacsScalar res[],
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[] ){
  // fprintf(stderr, "addResidual is being used\n");
  // The shape functions associated with the element  
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2*NUM_NODES], Nahp[2*NUM_NODES], Naahp[2*NUM_NODES];
  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];

  // The derivative of the stress with respect to the strain
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  TacsScalar transform[9];
  TacsScalar elemVars[NUM_VARIABLES];

  TacsScalar L = computeTransform(transform, Xpts);
  TacsScalar h = 2.0/L;
  TacsScalar hinv = 0.5*L;

  memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));
  transformVarsGlobalToLocal(elemVars, vars, transform);

  for ( int n = 0; n < numGauss; n++ ){
    FElibrary::lagrangeSF(N, Na, gaussPts[n], 2);
    FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[n]);

    evalStrain(h, Na, Nahp, Naahp, elemVars, strain);
    evalBvec(h, Na, Nahp, Naahp, B);

    beamStiff->calculateStress(&gaussPts[n], strain, stress);

    for ( int i = 0; i < NUM_NODES; i++ ){
      for ( int ii = 0; ii < NUM_DISPS; ii++ ){
        int row = ii + i*NUM_DISPS;

        res[row] +=
                hinv*gaussWts[n]*strain_product(&B[NUM_STRESSES*row], stress);
      }
    }

    // TacsScalar mass;
    // beamStiff->getPointwiseMass(gaussPts[n], &mass);
  }

  transformResLocalToGlobal(res, transform);
  // res = NULL;
}

/*
  Get the derivative of the residual w.r.t. the nodal coordinates.
*/
void EBBeam::getResXptSens( TacsScalar * res,
			    const TacsScalar vars[],
			    const TacsScalar Xpts[] ){

  memset(res, 0, 3*NUM_NODES*NUM_VARIABLES*sizeof(TacsScalar));

  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2*NUM_NODES], Nahp[2*NUM_NODES], Naahp[2*NUM_NODES];

  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
  TacsScalar stressSens[NUM_STRESSES], strainSens[NUM_STRESSES];

  TacsScalar B[NUM_STRESSES*NUM_VARIABLES], BSens[NUM_STRESSES*NUM_VARIABLES];

  TacsScalar elemVars[NUM_VARIABLES], elemVarsSens[NUM_VARIABLES];
  TacsScalar transform[9], transformSens[9];

  computeTransform(transform, Xpts);
  transformVarsGlobalToLocal(elemVars, vars, transform);

  for ( int k = 0; k < 3*NUM_NODES; k++ ){
    TacsScalar resSens[NUM_VARIABLES];
    memset(resSens, 0, NUM_VARIABLES*sizeof(TacsScalar));

    TacsScalar LSens;
    TacsScalar L = computeTransformXptSens(transform, transformSens,
                                           &LSens, k, Xpts);

    TacsScalar h = 2.0/L;
    TacsScalar hSens = - 2.0*LSens/(L*L);
    TacsScalar hinv = 0.5*L;
    TacsScalar hinvSens = 0.5*LSens;

    transformVarsGlobalToLocal(elemVarsSens, vars, transformSens);

    for ( int n = 0; n < numGauss; n++ ){
      FElibrary::lagrangeSF(N, Na, gaussPts[n], 2);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[n]);

      evalStrain(h, Na, Nahp, Naahp, elemVars, strain);
      evalBvec(h, Na, Nahp, Naahp, B);

      evalStrainSens(h, hSens, Na, Nahp, Naahp,
                     elemVars, elemVarsSens, strainSens);
      evalBvecSens(h, hSens, Na, Nahp, Naahp, BSens);

      beamStiff->calculateStress(&gaussPts[n], strain, stress);
      beamStiff->calculateStress(&gaussPts[n], strainSens, stressSens);

      for ( int i = 0; i < NUM_NODES; i++ ){
        for ( int ii = 0; ii < NUM_DISPS; ii++ ){
          int row = ii + i*NUM_DISPS;

          res[row] +=
            gaussWts[n]*hinv*strain_product(&B[NUM_STRESSES*row], stress);

          resSens[row] +=
            gaussWts[n]*(hinvSens*strain_product(&B[NUM_STRESSES*row], stress) +
                         hinv*(strain_product(&B[NUM_STRESSES*row],
                                              stressSens) +
                               strain_product(&BSens[NUM_STRESSES*row],
                                              stress)));

        }
      }
    }

    transformResLocalToGlobalSens(res, resSens,
                                  transform, transformSens);
    res += NUM_VARIABLES;
  }
}

void EBBeam::getResDVSens( int dvNum, TacsScalar * res,
			   const TacsScalar vars[],
			   const TacsScalar Xpts[] ){
  if (dvNum >= 0 && ownsDesignVar(dvNum)){
    double N[NUM_NODES], Na[NUM_NODES];
    double Nhp[2*NUM_NODES], Nahp[2*NUM_NODES], Naahp[2*NUM_NODES];

    TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
    TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

    TacsScalar transform[9];
    TacsScalar elemVars[NUM_VARIABLES];

    TacsScalar L = computeTransform(transform, Xpts);
    TacsScalar h = 2.0/L;
    TacsScalar hinv = 0.5*L;

    memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));
    transformVarsGlobalToLocal(elemVars, vars, transform);

    for ( int n = 0; n < numGauss; n++ ){
      FElibrary::lagrangeSF(N, Na, gaussPts[n], 2);
      FElibrary::cubicHP(Nhp, Nahp, Naahp, gaussPts[n]);

      evalStrain(h, Na, Nahp, Naahp, elemVars, strain);
      evalBvec(h, Na, Nahp, Naahp, B);

      beamStiff->calculateStressDVSens(dvNum, &gaussPts[n], strain, stress);

      for ( int ii = 0; ii < NUM_DISPS; ii++ ){
        for ( int i = 0; i < NUM_NODES; i++ ){
          int row = ii + i*NUM_DISPS;

          res[row] +=
            hinv*gaussWts[n]*strain_product(&B[NUM_STRESSES*row], stress);
        }
      }
    }

    transformResLocalToGlobal(res, transform);
  }
}

/*
  Get the rigid displacement at point Xa, eminating from the
  parametric location pt, with element variables vars
*/
void EBBeam::getRD( TacsScalar Uaero[],
                    const TacsScalar Xa[], const double pt[],
                    const TacsScalar vars[], const TacsScalar Xpts[] ){
  // Calculate the linear shape functions
  double N1, N2;
  N1 = 0.5*(1.0 - pt[0]);
  N2 = 0.5*(1.0 + pt[0]);

  // Find the rigid link vector in three-space
  TacsScalar Rd[3];
  Rd[0] = Xa[0] - N1*Xpts[0] - N2*Xpts[3];
  Rd[1] = Xa[1] - N1*Xpts[1] - N2*Xpts[4];
  Rd[2] = Xa[2] - N1*Xpts[2] - N2*Xpts[5];

  TacsScalar U[6];
  for ( int k = 0; k < 6; k++ ){
    U[k] = N1*vars[k] + N2*vars[6+k];
  }

  Uaero[0] = U[0] + U[4]*Rd[2] - U[5]*Rd[1];
  Uaero[1] = U[1] + U[5]*Rd[0] - U[3]*Rd[2];
  Uaero[2] = U[2] + U[3]*Rd[1] - U[4]*Rd[0];
}

/*
  Calculate the derivative of the aerodynamic point 'Uaero' w.r.t. the
  input perturbation XptSens
*/
void EBBeam::getRDXptSens( TacsScalar Uaero[],
                           const TacsScalar Xa[], const TacsScalar XaSens[],
                           const double pt[], const TacsScalar vars[],
                           const TacsScalar Xpts[],
                           const TacsScalar XptSens[] ){
  // Calculate the linear shape functions
  double N1, N2;
  N1 = 0.5*(1.0 - pt[0]);
  N2 = 0.5*(1.0 + pt[0]);

  // Find the rigid link vector in three-space
  TacsScalar Rd[3];
  Rd[0] = XaSens[0] - N1*XptSens[0] - N2*XptSens[3];
  Rd[1] = XaSens[1] - N1*XptSens[1] - N2*XptSens[4];
  Rd[2] = XaSens[2] - N1*XptSens[2] - N2*XptSens[5];

  TacsScalar U[6];
  for ( int k = 0; k < 6; k++ ){
    U[k] = N1*vars[k] + N2*vars[6+k];
  }

  Uaero[0] = U[4]*Rd[2] - U[5]*Rd[1];
  Uaero[1] = U[5]*Rd[0] - U[3]*Rd[2];
  Uaero[2] = U[3]*Rd[1] - U[4]*Rd[0];
}

void EBBeam::getRDTranspose( TacsScalar elemAdj[], const TacsScalar Uaero[],
                             const TacsScalar Xa[], const double pt[],
                             const TacsScalar Xpts[] ){
  // Calculate the linear shape functions
  double N1, N2;
  N1 = 0.5*(1.0 - pt[0]);
  N2 = 0.5*(1.0 + pt[0]);

  // Find the rigid link vector in three-space
  TacsScalar Rd[3];
  Rd[0] = Xa[0] - N1*Xpts[0] - N2*Xpts[3];
  Rd[1] = Xa[1] - N1*Xpts[1] - N2*Xpts[4];
  Rd[2] = Xa[2] - N1*Xpts[2] - N2*Xpts[5];

  elemAdj[0] = N1*Uaero[0];
  elemAdj[1] = N1*Uaero[1];
  elemAdj[2] = N1*Uaero[2];
  elemAdj[3] = N1*(Rd[1]*Uaero[2] - Rd[2]*Uaero[1]);
  elemAdj[4] = N1*(Rd[2]*Uaero[0] - Rd[0]*Uaero[2]);
  elemAdj[5] = N1*(Rd[0]*Uaero[1] - Rd[1]*Uaero[0]);

  elemAdj[6] = N2*Uaero[0];
  elemAdj[7] = N2*Uaero[1];
  elemAdj[8] = N2*Uaero[2];
  elemAdj[9]  = N2*(Rd[1]*Uaero[2] - Rd[2]*Uaero[1]);
  elemAdj[10] = N2*(Rd[2]*Uaero[0] - Rd[0]*Uaero[2]);
  elemAdj[11] = N2*(Rd[0]*Uaero[1] - Rd[1]*Uaero[0]);
}

/*
  Get the additional contribution to the residual from the aerodynamic
  forces.
*/
void EBBeam::getRF( TacsScalar res[], const TacsScalar Fnode[],
                    const TacsScalar Xa[], const double pt[],
                    const TacsScalar Xpts[] ){
  // Calculate the linear shape functions
  double N1, N2;
  N1 = 0.5*(1.0 - pt[0]);
  N2 = 0.5*(1.0 + pt[0]);

  // Find the rigid link vector in three-space
  TacsScalar Rd[3];
  Rd[0] = Xa[0] - N1*Xpts[0] - N2*Xpts[3];
  Rd[1] = Xa[1] - N1*Xpts[1] - N2*Xpts[4];
  Rd[2] = Xa[2] - N1*Xpts[2] - N2*Xpts[5];

  res[0] = -Fnode[0]*N1;
  res[1] = -Fnode[1]*N1;
  res[2] = -Fnode[2]*N1;
  res[3] = - (Fnode[2]*Rd[1] - Fnode[1]*Rd[2])*N1;
  res[4] = - (Fnode[0]*Rd[2] - Fnode[2]*Rd[0])*N1;
  res[5] = - (Fnode[1]*Rd[0] - Fnode[0]*Rd[1])*N1;

  res[6] = -Fnode[0]*N2;
  res[7] = -Fnode[1]*N2;
  res[8] = -Fnode[2]*N2;
  res[9]  = - (Fnode[2]*Rd[1] - Fnode[1]*Rd[2])*N2;
  res[10] = - (Fnode[0]*Rd[2] - Fnode[2]*Rd[0])*N2;
  res[11] = - (Fnode[1]*Rd[0] - Fnode[0]*Rd[1])*N2;
}

void EBBeam::getRFXptSens( TacsScalar res[], const TacsScalar Fnode[],
                           const TacsScalar Xa[], const TacsScalar XaSens[],
                           const double pt[], const TacsScalar Xpts[],
                           const TacsScalar XptSens[] ){
  // Calculate the linear shape functions
  double N1, N2;
  N1 = 0.5*(1.0 - pt[0]);
  N2 = 0.5*(1.0 + pt[0]);

  // Find the rigid link vector in three-space
  TacsScalar Rd[3];
  Rd[0] = XaSens[0] - N1*XptSens[0] - N2*XptSens[3];
  Rd[1] = XaSens[1] - N1*XptSens[1] - N2*XptSens[4];
  Rd[2] = XaSens[2] - N1*XptSens[2] - N2*XptSens[5];

  res[0] = 0.0;
  res[1] = 0.0;
  res[2] = 0.0;
  res[3] = - (Fnode[2]*Rd[1] - Fnode[1]*Rd[2])*N1;
  res[4] = - (Fnode[0]*Rd[2] - Fnode[2]*Rd[0])*N1;
  res[5] = - (Fnode[1]*Rd[0] - Fnode[0]*Rd[1])*N1;

  res[6] = 0.0;
  res[7] = 0.0;
  res[8] = 0.0;
  res[9]  = - (Fnode[2]*Rd[1] - Fnode[1]*Rd[2])*N2;
  res[10] = - (Fnode[0]*Rd[2] - Fnode[2]*Rd[0])*N2;
  res[11] = - (Fnode[1]*Rd[0] - Fnode[0]*Rd[1])*N2;
}

void EBBeam::getRFTranspose( TacsScalar Fnode[], const TacsScalar Xa[],
                             const double pt[], const TacsScalar res[],
                             const TacsScalar Xpts[] ){
  // Calculate the linear shape functions
  double N1, N2;
  N1 = 0.5*(1.0 - pt[0]);
  N2 = 0.5*(1.0 + pt[0]);

  // Find the rigid link vector in three-space
  TacsScalar Rd[3];
  Rd[0] = Xa[0] - N1*Xpts[0] - N2*Xpts[3];
  Rd[1] = Xa[1] - N1*Xpts[1] - N2*Xpts[4];
  Rd[2] = Xa[2] - N1*Xpts[2] - N2*Xpts[5];

  Fnode[0] = -(res[0]*N1 + res[6]*N2);
  Fnode[1] = -(res[1]*N1 + res[7]*N2);
  Fnode[2] = -(res[2]*N1 + res[8]*N2);

  Fnode[0] -= (res[4]*Rd[2] - res[5]*Rd[1])*N1;
  Fnode[1] -= (res[3]*Rd[2] + res[5]*Rd[0])*N1;
  Fnode[2] -= (res[3]*Rd[1] - res[4]*Rd[0])*N1;

  Fnode[0] -= (res[10]*Rd[2] - res[11]*Rd[1])*N2;
  Fnode[1] -= ( res[9]*Rd[2] + res[11]*Rd[0])*N2;
  Fnode[2] -= ( res[9]*Rd[1] - res[10]*Rd[0])*N2;
}

/*
  Get the number of Gauss quadrature points
*/
int EBBeam::getGaussPtScheme(){
  return numGauss;
}

/*
  Get the number of Gauss points associated with this element.
*/
int EBBeam::getNumGaussPts( const int scheme ){
  if (scheme <= 0){
    return numGauss;
  }
  return scheme;
}

/*
  Get the 'num'th Gauss point and weight
*/
TacsScalar EBBeam::getGaussWtsPts( const int scheme,
				   const int num, double * pt ){
  if (scheme <= 0){
    pt[0] = gaussPts[num];
    return gaussWts[num];
  }

  return FElibrary::getGaussPtWt(1, scheme, num, pt);
}

/*
  Locate the closest point in physical space to the provided point.

  Here, dist is an output value and is equal to the distance to the
  point. Only overwrite arrays the minimum distance is less than
  the value that exists in dist
*/
TacsScalar EBBeam::getClosestPt( double pt[], int * fflag,
                                 const TacsScalar Xp[],
                                 const TacsScalar Xpts[] ){

  // The square of the distance is:
  // dist^2 = dot(X0 - Xu*u, X0 - Xu*u)
  //        = dot(X0, X0) - 2.0*u*dot(X0, Xu) + u^2*dot(Xu, Xu)

  // The parametric location with the minimum distance is:
  // u_min = dot(X0, Xu)/dot(Xu, Xu)

  TacsScalar X0[3], Xu[3];
  X0[0] = Xp[0] - 0.5*(Xpts[0] + Xpts[3]);
  X0[1] = Xp[1] - 0.5*(Xpts[1] + Xpts[4]);
  X0[2] = Xp[2] - 0.5*(Xpts[2] + Xpts[5]);

  Xu[0] = 0.5*(Xpts[3] - Xpts[0]);
  Xu[1] = 0.5*(Xpts[4] - Xpts[1]);
  Xu[2] = 0.5*(Xpts[5] - Xpts[2]);
  TacsScalar u_min = Tensor::dot3D(X0, Xu)/Tensor::dot3D(Xu, Xu);

  if (u_min < -1.0){
    u_min = -1.0;
  }
  else if (u_min > 1.0){
    u_min = 1.0;
  }

  X0[0] = X0[0] + u_min*Xu[0];
  X0[1] = X0[1] + u_min*Xu[1];
  X0[2] = X0[2] + u_min*Xu[2];

  TacsScalar dist = sqrt(Tensor::dot3D(X0, X0));

  pt[0] = u_min;
  pt[1] = pt[2] = 0.0;

  return dist;
}

void EBBeam::getPoint( TacsScalar Xp[], const double pt[],
                       const TacsScalar Xpts[] ){
  Xp[0] = 0.5*((1.0 - pt[0])*Xpts[0] + (1.0 + pt[0])*Xpts[3]);
  Xp[1] = 0.5*((1.0 - pt[0])*Xpts[1] + (1.0 + pt[0])*Xpts[4]);
  Xp[2] = 0.5*((1.0 - pt[0])*Xpts[2] + (1.0 + pt[0])*Xpts[5]);
}

/*
  Just provide the shape functions for the position
*/
void EBBeam::getShapeFunctions( const double pt[], double N[] ){
  N[0] = 0.5*(1.0 - pt[0]);
  N[1] = 0.5*(1.0 + pt[0]);
}
/*
  Calculate the value of determinant of the Jacobian transform.  In
  this case, it is half the length of the beam.
*/
TacsScalar EBBeam::getJacobian( const double * pt, const TacsScalar Xpts[] ){
  TacsScalar L = sqrt((Xpts[3]-Xpts[0])*(Xpts[3]-Xpts[0]) +
                      (Xpts[4]-Xpts[1])*(Xpts[4]-Xpts[1]) +
                      (Xpts[5]-Xpts[2])*(Xpts[5]-Xpts[2]));
  return 0.5*L;
}

/*
  Calculate the derivative of half the length of the beam with respect
  to the end points.
*/
TacsScalar EBBeam::getJacobianXptSens( TacsScalar hXptSens[],
                                       const double * pt,
				       const TacsScalar Xpts[] ){

  TacsScalar L = sqrt((Xpts[3]-Xpts[0])*(Xpts[3]-Xpts[0]) +
                      (Xpts[4]-Xpts[1])*(Xpts[4]-Xpts[1]) +
                      (Xpts[5]-Xpts[2])*(Xpts[5]-Xpts[2]));

  for ( int j = 0; j < 3; j++ ){
    hXptSens[j] = -0.5*(Xpts[3+j]-Xpts[j])/L;
  }

  for ( int j = 0; j < 3; j++ ){
    hXptSens[j+3] = 0.5*(Xpts[3+j]-Xpts[j])/L;
  }

  return 0.5*L;
}

/*
  Calculate the strain at the given location within the beam
*/
void EBBeam::getPtwiseStrain( TacsScalar strain[], const double * pt,
                              const TacsScalar vars[],
                              const TacsScalar Xpts[] ){
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2*NUM_NODES], Nahp[2*NUM_NODES], Naahp[2*NUM_NODES];

  TacsScalar transform[9];
  TacsScalar elemVars[NUM_VARIABLES];

  TacsScalar L = computeTransform(transform, Xpts);
  transformVarsGlobalToLocal(elemVars, vars, transform);
  TacsScalar h = 2.0/L;

  FElibrary::lagrangeSF(N, Na, pt[0], 2);
  FElibrary::cubicHP(Nhp, Nahp, Naahp, pt[0]);

  evalStrain(h, Na, Nahp, Naahp, elemVars, strain);
}

/*
  Calculate the derivative of the strain w.r.t. the nodal coordinates.
*/
void EBBeam::getPtwiseStrainXptSens( TacsScalar strain[],
				     TacsScalar strainSens[],
				     const double * pt, const TacsScalar vars[],
				     const TacsScalar Xpts[] ){
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2*NUM_NODES], Nahp[2*NUM_NODES], Naahp[2*NUM_NODES];

  TacsScalar transform[9], transformSens[9];
  TacsScalar elemVars[NUM_VARIABLES], elemVarsSens[NUM_VARIABLES];

  TacsScalar L = computeTransform(transform, Xpts);
  transformVarsGlobalToLocal(elemVars, vars, transform);
  TacsScalar h = 2.0/L;

  FElibrary::lagrangeSF(N, Na, pt[0], 2);
  FElibrary::cubicHP(Nhp, Nahp, Naahp, pt[0]);

  evalStrain(h, Na, Nahp, Naahp, elemVars, strain);

  for ( int k = 0; k < 3*NUM_NODES; k++ ){
    TacsScalar LSens;
    computeTransformXptSens(transform, transformSens,
                            &LSens, k, Xpts);
    transformVarsGlobalToLocal(elemVarsSens, vars, transformSens);
    TacsScalar hSens = - 2.0*LSens/(L*L);

    evalStrainSens(h, hSens, Na, Nahp, Naahp,
                   elemVars, elemVarsSens, &strainSens[k*NUM_STRESSES]);
  }
}

/*
  Calculate the derivative of the product of the strain and a vector.
  Add the value of this derivative, times the scalar 'scaleFactor' to
  the array 'elementSens'.
*/
void EBBeam::addPtwiseStrainSVSens( TacsScalar elementSens[],
				    const double * pt,
				    const TacsScalar scaleFactor,
				    const TacsScalar strainSens[],
				    const TacsScalar vars[],
				    const TacsScalar Xpts[] ){

  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2*NUM_NODES], Nahp[2*NUM_NODES], Naahp[2*NUM_NODES];
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  TacsScalar transform[9];
  TacsScalar elemVars[NUM_VARIABLES], elemSens[NUM_VARIABLES];
  memset(elemSens, 0, NUM_VARIABLES*sizeof(TacsScalar));

  TacsScalar L = computeTransform(transform, Xpts);
  TacsScalar h = 2.0/L;
  transformVarsGlobalToLocal(elemVars, vars, transform);

  FElibrary::lagrangeSF(N, Na, pt[0], 2);
  FElibrary::cubicHP(Nhp, Nahp, Naahp, pt[0]);
  evalBvec(h, Na, Nahp, Naahp, B);

  for ( int ii = 0; ii < NUM_DISPS; ii++ ){
    for ( int i = 0; i < NUM_NODES; i++ ){
      int row = i + NUM_NODES*ii;
      elemSens[row] +=
        scaleFactor*strain_product(strainSens, &B[row*NUM_STRESSES]);
    }
  }

  transformResLocalToGlobal(elemSens, transform);

  for ( int k = 0; k < NUM_VARIABLES; k++ ){
    elementSens[k] += elemSens[k];
  }
}

/*
  The linear strain expressions, and their derivatives for the beam
  element

  The displacement function is as follows:

  v = (v0[0]*Nhp[0] + rotz[0]*Nhp[1]
  .    v0[1]*Nhp[2] + rotz[1]*Nhp[3])

  w = (w0[0]*Nhp[0] - roty[0]*Nhp[1]
  .    w0[1]*Nhp[2] - roty[1]*Nhp[3])
*/

void EBBeam::evalStrain( TacsScalar h,
                         const double Na[],
                         const double Nahp[],
                         const double Naahp[],
                         const TacsScalar vars[],
                         TacsScalar strain[] ){
  TacsScalar u0_x = h*(Na[0]*vars[0] + Na[1]*vars[6]);

  TacsScalar w0_xx = h*(h*(Naahp[0]*vars[2] + Naahp[2]*vars[8]) -
                        (Naahp[1]*vars[4] + Naahp[3]*vars[10]));

  TacsScalar v0_xx = h*(h*(Naahp[0]*vars[1] + Naahp[2]*vars[7]) +
                        (Naahp[1]*vars[5] + Naahp[3]*vars[11]));

  TacsScalar rotx_x = h*(Na[0]*vars[3] + Na[1]*vars[9]);

  strain[0] = u0_x;
  strain[1] = -w0_xx;
  strain[2] = -v0_xx;
  strain[3] = rotx_x;
}

/*
  Evaluate the derivative of the strain w.r.t. the nodal degrees of
  freedom.
*/
void EBBeam::evalBvec( TacsScalar h,
                       const double Na[],
                       const double Nahp[],
                       const double Naahp[],
                       TacsScalar B[] ){
  memset(B, 0, NUM_STRESSES*NUM_VARIABLES*sizeof(TacsScalar));

  // u0
  B[0] = h*Na[0];          B += NUM_STRESSES;

  // v0
  B[2] = -h*h*Naahp[0];    B += NUM_STRESSES;

  // w0
  B[1] = -h*h*Naahp[0];    B += NUM_STRESSES;

  // rotx
  B[3] = h*Na[0];          B += NUM_STRESSES;

  // roty
  B[1] = h*Naahp[1];       B += NUM_STRESSES;

  // rotz
  B[2] = -h*Naahp[1];      B += NUM_STRESSES;

  // ----------------------------------------

  // u0
  B[0] = h*Na[1];          B += NUM_STRESSES;

  // v0
  B[2] = -h*h*Naahp[2];    B += NUM_STRESSES;

  // w0
  B[1] = -h*h*Naahp[2];    B += NUM_STRESSES;

  // rotx
  B[3] = h*Na[1];          B += NUM_STRESSES;

  // roty
  B[1] = h*Naahp[3];       B += NUM_STRESSES;

  // rotz
  B[2] = -h*Naahp[3];
}

void EBBeam::evalStrainSens( TacsScalar h, TacsScalar hSens,
                             const double Na[],
                             const double Nahp[],
                             const double Naahp[],
                             const TacsScalar vars[],
                             const TacsScalar varsSens[],
                             TacsScalar strainSens[] ){

  TacsScalar u0_x = (hSens*(Na[0]*vars[0] + Na[1]*vars[6]) +
                     h*(Na[0]*varsSens[0] + Na[1]*varsSens[6]));

  TacsScalar w0_xx = (hSens*(2.0*h*(Naahp[0]*vars[2] + Naahp[2]*vars[8]) -
                             (Naahp[1]*vars[4] + Naahp[3]*vars[10])) +
                      h*(h*(Naahp[0]*varsSens[2] + Naahp[2]*varsSens[8]) -
                         (Naahp[1]*varsSens[4] + Naahp[3]*varsSens[10])));

  TacsScalar v0_xx = (hSens*(2.0*h*(Naahp[0]*vars[1] + Naahp[2]*vars[7]) +
                             (Naahp[1]*vars[5] + Naahp[3]*vars[11])) +
                      h*(h*(Naahp[0]*varsSens[1] + Naahp[2]*varsSens[7]) +
                         (Naahp[1]*varsSens[5] + Naahp[3]*varsSens[11])));

  TacsScalar rotx_x = (hSens*(Na[0]*vars[3] + Na[1]*vars[9]) +
                       h*(Na[0]*varsSens[3] + Na[1]*varsSens[9]));

  strainSens[0] = u0_x;
  strainSens[1] = -w0_xx;
  strainSens[2] = -v0_xx;
  strainSens[3] = rotx_x;
}

void EBBeam::evalBvecSens( TacsScalar h, TacsScalar hSens,
                           const double Na[],
                           const double Nahp[],
                           const double Naahp[],
                           TacsScalar B[] ){
  memset(B, 0, NUM_STRESSES*NUM_VARIABLES*sizeof(TacsScalar));

  // u0
  B[0] = hSens*Na[0];              B += NUM_STRESSES;

  // v0
  B[2] = -2.0*h*hSens*Naahp[0];    B += NUM_STRESSES;

  // w0
  B[1] = -2.0*h*hSens*Naahp[0];    B += NUM_STRESSES;

  // rotx
  B[3] = hSens*Na[0];              B += NUM_STRESSES;

  // roty
  B[1] = hSens*Naahp[1];           B += NUM_STRESSES;

  // rotz
  B[2] = -hSens*Naahp[1];          B += NUM_STRESSES;

  // ----------------------------------------

  // u0
  B[0] = hSens*Na[1];              B += NUM_STRESSES;

  // v0
  B[2] = -2.0*h*hSens*Naahp[2];    B += NUM_STRESSES;

  // w0
  B[1] = -2.0*h*hSens*Naahp[2];    B += NUM_STRESSES;

  // rotx
  B[3] = hSens*Na[1];              B += NUM_STRESSES;

  // roty
  B[1] = hSens*Naahp[3];           B += NUM_STRESSES;

  // rotz
  B[2] = -hSens*Naahp[3];
}

void EBBeam::addOutputCount( int * nelems, int * nnodes, int * ncsr ){
  *nelems += NUM_NODES-1;
  *nnodes += NUM_NODES;
  *ncsr += NUM_NODES;
}

void EBBeam::getOutputConnectivity( int * con, int node ){
  for ( int n = 0; n < NUM_NODES-1; n++ ){
    con[2*n]   = node;
    con[2*n+1] = node+1;
  }
}

/*
  Retrieve the output data for this element
*/
void EBBeam::getOutputData( unsigned int out_type,
			    double * data, int ld_data,
			    const TacsScalar vars[],
                            const TacsScalar Xpts[] ){
  double N[NUM_NODES], Na[NUM_NODES];
  double Nhp[2*NUM_NODES], Nahp[2*NUM_NODES], Naahp[2*NUM_NODES];
  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];

  TacsScalar transform[9];
  TacsScalar elemVars[NUM_VARIABLES];

  TacsScalar L = computeTransform(transform, Xpts);
  TacsScalar h = 2.0/L;

  transformVarsGlobalToLocal(elemVars, vars, transform);

  for ( int n = 0; n < NUM_NODES; n++ ){
    double pt = -1.0 + (2.0*n)/NUM_NODES;
    int index = 0;
    if (out_type & TACSElement::OUTPUT_NODES){
      for ( int k = 0; k < 3; k++ ){
	data[index+k] = TacsRealPart(Xpts[3*n+k]);
      }
      index += 3;
    }
    if (out_type & TACSElement::OUTPUT_DISPLACEMENTS){
      for ( int k = 0; k < NUM_DISPS; k++ ){
	data[index+k] = TacsRealPart(vars[NUM_DISPS*n+k]);
      }
      index += NUM_DISPS;
    }

    FElibrary::lagrangeSF(N, Na, pt, 2);
    FElibrary::cubicHP(Nhp, Nahp, Naahp, pt);

    evalStrain(h, Na, Nahp, Naahp, elemVars, strain);
    if (out_type & TACSElement::OUTPUT_STRAINS){
      for ( int k = 0; k < NUM_STRESSES; k++ ){
        data[index+k] = TacsRealPart(strain[k]);
      }
      index += NUM_STRESSES;
    }
    if (out_type & TACSElement::OUTPUT_STRESSES){
      beamStiff->calculateStress(&pt, strain, stress);
      for ( int k = 0; k < NUM_STRESSES; k++ ){
        data[index+k] = TacsRealPart(stress[k]);
      }
      index += NUM_STRESSES;
    }

    data += ld_data;
  }
}
