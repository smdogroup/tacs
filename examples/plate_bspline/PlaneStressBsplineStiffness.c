#include "PlaneStressBsplineStiffness.h"
#include "YSlibrary.h"
/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/
const char * PlaneStressBsplineStiffness::constName = "PlaneStressBsplineStiffness";

const char * PlaneStressBsplineStiffness::constitutiveName(){ 
  return constName; 
}

PlaneStressBsplineStiffness::PlaneStressBsplineStiffness( TacsScalar _rho, 
                                                          TacsScalar _E,
                                                          TacsScalar _nu,
                                                          TacsScalar _ys,
                                                          TacsScalar _epsilon,
                                                          TacsScalar *_x,
                                                          double _qpenalty,
                                                          double _lower,
                                                          double *_Tu, double *_Tv,
                                                          int _Lu, int _Lv, 
                                                          int _pNum,
                                                          int _order ){
  // Number of knot intervals
  Lu = _Lu;
  Lv = _Lv;
  // 
  // Copy over the patch number and order
  pNum = _pNum;
  order = _order;

  Tu = new double[Lu+1+2*(order-1)];
  Tv = new double[Lv+1+2*(order-1)];
  memcpy(Tu, _Tu, (Lu+1+2*(order-1))*sizeof(double));
  memcpy(Tv, _Tv, (Lv+1+2*(order-1))*sizeof(double));

  // Compute the interval at which the patch is in
  t2 = int(floor(pNum/(Lu))); // y direction
  t1 = int(floor(pNum-t2*(Lu))); // x direction  
  
  // Initialize the constitutive properties
  E = _E;
  nu = _nu;
  ys = _ys;
  rho = _rho;
  epsilon = _epsilon;
  // Initialize the topology variables
  dvNum = t1+t2*Lu;
  q = _qpenalty;
  lowerBound = _lower;
  x = new TacsScalar[Lu*Lv];
  xw = NULL;
}
int PlaneStressBsplineStiffness::getNumStresses(){
  return NUM_STRESSES;
}
// Set the design variables
void PlaneStressBsplineStiffness::setDesignVars( const TacsScalar dvs[], 
                                                 int numDVs ){
  if (dvNum < numDVs){
    memcpy(x, dvs, numDVs*sizeof(TacsScalar));
  }    
  computeIndex(&index, numDVs);
}
// Get the design variables
void PlaneStressBsplineStiffness::getDesignVars( TacsScalar dvs[],
                                                 int numDVs ){
  if (dvNum < numDVs){
    dvs[dvNum] = rho[dvNum];
  }
}
// Get the lower and upper bound of design variables
void PlaneStressBsplineStiffness::getDesignVarRange( TacsScalar lb[],
                                                     TacsScalar ub[],
                                                     int numDVs ){
  if (dvNum < numDVs){
    lb[dvNum] = lowerBound;
    ub[dvNum] = 1.0;
  }
}
// Evaluate the mass of the domain
void PlaneStressBsplineStiffness::getPointwiseMass( const double pt[],
                                                    TacsScalar mass[] ){
  double N[16];
  getShapeFunctions(pt, N);
  // Compute the topology variable
  xw = 0.0;
  for (int i = 0; i < order*order; i++){
    xw += N[i]*x[index[i]];
  }
  mass[0] = xw*rho;
}
// Evaluate the derivative of the mass w.r.t. the design variable and 
void PlaneStressBsplineStiffness::addPointwiseMassDVSens( const double pt[], 
                                                          const TacsScalar alpha[],
                                                          TacsScalar dvSens[], 
                                                          int dvLen ){
  if (dvNum < dvLen){
    double N[16];
    getShapeFunctions(pt, N);
    int ind = findPatch();
    // Multiply by corresponding shape function for pNum
    dvSens[dvNum] += alpha[0]*rho*N[ind];
  }
}
// Calculate the stress at the integration point
void PlaneStressBsplineStiffness::calculateStress( const double pt[], 
                                                   const TacsScalar strain[],
                                                   TacsScalar stress[]){
  // Compute the material parameters
  TacsScalar D = E/(1.0-nu*nu);
  double N[16];
  getShapeFunctions(pt, N);
  // Compute the topology variable
  xw = 0.0;
  for (int i = 0; i < order*order; i++){
    xw += N[i]*x[index[i]];
  }
  TacsScalar w = xw/(1.0+q*(1.0-xw));
  
  D11 = D*(1.0-nu)*w;
  D12 = D*nu*w;
  G = 0.5*(1.0-2.0*nu)*D*w;  
  
  stress[0] = D11*strain[0]+D12*(strain[1]+strain[2]);
  stress[0] = D11*strain[1]+D12*(strain[0]+strain[2]);
  stress[0] = D11*strain[2]+D12*(strain[0]+strain[1]);

  stress[3] = G*strain[3];
  stress[4] = G*strain[4];
  stress[5] = G*strain[5];
  
}
// Compute the derivative of the stress w.r.t. the design variable
void PlaneStressBsplineStiffness::addStressDVSens( const double pt[], 
                                                   const TacsScalar strain[], 
                                                   TacsScalar alpha, 
                                                   const TacsScalar psi[], 
                                                   TacsScalar dvSens[], 
                                                   int dvLen ){
  if (dvNum < dvLen){
    double N[16];
    getShapeFunctions(pt, N);
    xw = 0.0;
    for (int i = 0; i < order*order; i++){
      xw += N[i]*x[index[i]];
    }
    TacsScalar s[6];
    calculateStress(pt, strain,s);
    TacsScalar dxw = 1.0+q*(1.0-xw);
    TacsScalar w = ((1.0+q)/(dxw*dxw));
    TacsScalar inner = alpha*w*(psi[0]*s[0]+psi[1]*s[1]+
                                psi[2]*s[2]+psi[3]*s[3]+ 
                                psi[4]*s[4]+psi[5]*s[5]);
    // Add the shape function corresponding to pNum
    int ind = findPatch();
    dvSens[dvNum] += N[ind]*inner;
    
  }
}

// Evaluate the failure criteria
void PlaneStressBsplineStiffness::failure( const double pt[], 
                                           const TacsScalar strain[],
                                           TacsScalar * fail ){
  // Use the von Mises failure criterion
  // Compute the relaxation factor
  TacsScalar r_factor = 1.0;
  // Compute the topology variable
  xw = 0.0;
  for (int i = 0; i < order*order; i++){
    xw += N[i]*x[index[i]];
  }
  
  if (epsilon > 0.0){
    r_factor = xw/(epsilon*(1.0-xw)+xw);
  }
  TacsScalar s[6];
  calculateStress(pt, strain, s);
  *fail = r_factor*VonMisesFailurePlaneStress(s,ys);
}
// Evaluate the failure criteria w.r.t. design variables
void PlaneStressBsplineStiffness::addFailureDVSens( const double pt[], 
                                                    const TacsScalar strain[],
                                                    TacsScalar alpha,
                                                    TacsScalar dvSens[], 
                                                    int dvLen ){
  // Compute the relaxation factor
  TacsScalar r_factor_sens = 0.0;
  // Compute the topology variable
  xw = 0.0;
  for (int i = 0; i < order*order; i++){
    xw += N[i]*x[index[i]];
  }
  if (epsilon > 0.0){
    TacsScalar d = 1.0/(epsilon*(1.0-xw)+xw);
    r_factor_sens = epsilon*d*d;
  }
  TacsScalar s[6];
  calculateStress(pt, strain, s);
  TacsScalar dxw = 1.0+q*(1.0-xw);
  TacsScalar w = (1.0+q)/(dxw*dxw);
  TacsScalar fail = VonMisesFailurePlaneStress(s,ys);
  TacsScalar inner = alpha*r_factor_sens*fail;
  // Add the shape function product corresponding to pNum
  int ind = findPatch();
  dvSens[dvNum] += N[ind]*inner;
  
}
// Find the index in the shape function tensor that belongs to the
// current patch
int PlaneStressBsplineStiffness::findPatch(){
  for (int i = 0; i < order*order; i++){
    if (pNum == index[i]){
      return i;
    }
  }
  return -1;
}

void PlaneStressBsplineStiffness::getShapeFunctions( double pt[], double N[]){
  // Do the filtering
  int intu = t1+order-1;
  int intv = t2+order-1;
  double u = Tu[intu]+Tu[intu+1]+pt[0]*(Tu[intu+1]-Tu[intu]);
  u *= 0.5;
  double v = Tv[intv]+Tv[intv+1]+pt[1]*(Tv[intv+1]-Tv[intv]);
  v *= 0.5;
  double Nu[4], Nv[4], work[8], N[16];
  // Compute the basis values in each direction
  FElibrary::bspline_basis_(Nu, intu, u, Tu, order, work);
  FElibrary::bspline_basis(Nv, intv, v, Tv, order, work);
  
  if (order == 3){
    N[0] = Nu[0]*Nv[0];
    N[1] = Nu[1]*Nv[0];
    N[2] = Nu[2]*Nv[0];
    
    N[3] = Nu[0]*Nv[1];
    N[4] = Nu[1]*Nv[1];
    N[5] = Nu[2]*Nv[1];
    
    N[6] = Nu[0]*Nv[2];
    N[7] = Nu[1]*Nv[2];
    N[8] = Nu[2]*Nv[2];
  }
  else if (order == 4){
    N[0] = Nu[0]*Nv[0];
    N[1] = Nu[1]*Nv[0];
    N[2] = Nu[2]*Nv[0];
    N[3] = Nu[3]*Nv[0];
    
    N[4] = Nu[0]*Nv[1];
    N[5] = Nu[1]*Nv[1];
    N[6] = Nu[2]*Nv[1];
    N[7] = Nu[3]*Nv[1];
    
    N[8] = Nu[0]*Nv[2];
    N[9] = Nu[1]*Nv[2];
    N[10] = Nu[2]*Nv[2];
    N[11] = Nu[3]*Nv[2];

    N[12] = Nu[0]*Nv[3];
    N[13] = Nu[1]*Nv[3];
    N[14] = Nu[2]*Nv[3];
    N[15] = Nu[3]*Nv[3];
  }
}

// Compute the patches affected by filter
void PlaneStressBsplineStiffness::computeIndexList( int **index, 
                                                    int numDVs ){
  int count = 0;
  int *_index = new int[order*order];
  // Check if patch is on the boundary (bottom, left, right, top)
  if (pNum < (order-2)*Lu || pNum % Lu == 0 || 
      (pNum-(order-3)) % Lu == 0 || 
      (pNum-Lv) % Lu == 0 ||
      (pNum-Lv+(order-3)) % Lu == 0 ||
      pNum >= numDVs-(order-2)*Lu){    
    
    int start = 0;
    // Check if it is at one of the four corner
    // Lower left corner
    if (pNum < 2 || pNum == Lu || pNum == Lu+order-3){
      start = 0;
    }
    // Lower right corner
    else if ((pNum >= Lu-2 && pNum < Lu) || pNum == Lu+Lv ||
             pNum == Lu+Lv-(order-3)){
      start = Lu-order;
    }
    // Top left corner
    else if ((pNum >= numDVs-Lu && pNum <= numDVs-Lu+1) ||
             pNum == numDVs-2*Lu ||
             pNum == numDVs-2*Lu+order-3){
      start = numDVs-order*Lu;
    }
    // Top right corner
    else if ((pNum >= numDVs-2 && pNum < numDVs) || 
             pNum == numDVs-Lu-1 ||
             pNum == numDVs-Lu-1-(order-3)){
      start = numDVs-order-(order-1)*Lu;
    }
    // Left with the interior boundary patches
    else {
      // For the bottom interior boundary patches
      if (pNum < (order-2)*Lu){
        if (pNum < Lu){
          start = pNum-1;
        }
        else{
          start = pNum-Lu-1;
        }
      }
      // For the left vertical interior boundary patches
      else if (pNum % Lu == 0 || 
               (pNum-(order-3)) % Lu == 0 ){
        if (pNum % Lu == 0){
          start = pNum-Lu;
        }
        else {
          start = pNum-Lu-1;
        }
      }
      // For the right vertical interior boundary patches
      else if ((pNum-Lv) % Lu == 0 ||
               (pNum-Lv+(order-3)) % Lu == 0 ){
        if ((pNum-Lv) % Lu == 0){
          start = pNum-Lu-(order-1);
        }
        else{
          start = pNum-Lu-(order-2);
        }
      }
      // For the top interior boundary patches
      else {
        if (pNum >= numDVs-Lu){
          start = pNum-(order-1)*Lu-1;
        }
        else {
          start = pNum-(order-2)*Lu-1;
        }
      }
    }
    // Assign the connectivity
    for (int j = 0; j < order; j++){
      for (int i = 0; i < order; i++){
        _index[count] = start+i+j*Lu;
        count++;
      }
    }
  }
  else{
    // Indices for interior patches
    for (int j = -1; j < -1+order; j++){
      for (int i = -1; i < -1+order; i++){
        _index[count] = pNum+j*Lu+i;
        count++;
      }      
    }    
  }
  *index = _index;
}

