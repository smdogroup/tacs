#ifndef TACS_PLANE_STRESS_BSPLINE_ALL_H
#define TACS_PLANE_STRESS_BSPLINE_ALL_H

#include "TACS2DElement.h"
#include "FElibrary.h"

template <int order>
class PlaneStressBsplineAll : public TACS2DElement<order*order> {
 public:
  PlaneStressBsplineAll( PlaneStressStiffness *_stiff,
                         double _Tu[], double _Tv[], 
                         int _Lu, int _Lv, 
                         ElementBehaviorType type=LINEAR,
                         int _componentNum=0,
                         int _pNum=0 );
  ~PlaneStressBsplineAll();
  
  // Return the name of this element
  // -------------------------------
  const char *elementName(){ return elemName; }

  // Retrieve the shape functions
  // ----------------------------
  void getShapeFunctions( const double pt[], double N[] );
  void getShapeFunctions( const double pt[], double N[],
			  double Na[], double Nb[] );

  // Retrieve the Gauss points/weights
  // ---------------------------------
  int getNumGaussPts();
  double getGaussWtsPts( const int num, double pt[] ); 

  // Functions for post-processing
  // -----------------------------
  void addOutputCount( int *nelems, int *nnodes, int *ncsr );
  void getOutputData( unsigned int out_type, 
		      double *data, int ld_data, 
		      const TacsScalar Xpts[],
		      const TacsScalar vars[] );
  void getOutputConnectivity( int *con, int node );
  void getXpts( TacsScalar U[], 
                const double N[],
                const TacsScalar Xpts[] );


 private:
  static const int NUM_NODES = order*order;

  // The Gauss quadrature scheme
  int numGauss;
  const double *gaussWts, *gaussPts;

  // Store the name of the element
  static const char *elemName;
  
  double Tu[11], Tv[11];
  // Patch number
  int pNum;
  // Interval number in x,y
  int t1, t2;
  int Lu, Lv;
};

template <int order>
PlaneStressBsplineAll<order>::PlaneStressBsplineAll( PlaneStressStiffness *_stiff,
                                                     double _Tu[], double _Tv[], 
                                                     int _Lu, int _Lv, 
                                                     ElementBehaviorType type,
                                                     int _componentNum,
                                                     int _pNum ):
TACS2DElement<order*order>(_stiff, type, _componentNum){
  // Number of knot intervals
  Lu = _Lu;
  Lv = _Lv;
  
  // Copy over the patch number
  pNum = _pNum;
  memcpy(Tu, _Tu, (Lu+1+2*(order-1))*sizeof(double));
  memcpy(Tv, _Tv, (Lv+1+2*(order-1))*sizeof(double));
  numGauss = FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);
  // Compute the interval at which the patch is in
  t2 = int(floor(pNum/(Lu))); // y direction
  t1 = int(floor(pNum-t2*(Lu))); // x direction
}

template <int order>
PlaneStressBsplineAll<order>::~PlaneStressBsplineAll(){}

template <int order>
const char *PlaneStressBsplineAll<order>::elemName = "PlaneStressBsplineAll";

/*
  Get the number of Gauss points in the Gauss quadrature scheme
*/
template <int order>
int PlaneStressBsplineAll<order>::getNumGaussPts(){
  return numGauss*numGauss;
}

/*
  Get the Gauss points
*/
template <int order>
double PlaneStressBsplineAll<order>::getGaussWtsPts( const int num,
                                                     double pt[] ){
  int m = (int)((num)/numGauss);
  int n = num-numGauss*m;
  
  pt[0] = gaussPts[n];
  pt[1] = gaussPts[m];
  return gaussWts[n]*gaussWts[m];
}

/*
  Evaluate the shape functions and their derivatives
*/
template <int order>
void PlaneStressBsplineAll<order>::getShapeFunctions( const double pt[], 
                                                      double N[] ){
  // Compute the knot intervals
  int intu = t1+order-1;
  int intv = t2+order-1;
  // Convert parametric point to knot vector point
  double Nu[4], Nv[4], work[8];
  double u = Tu[intu]+Tu[intu+1]+pt[0]*(Tu[intu+1]-Tu[intu]);
  u *= 0.5;
  double v = Tv[intv]+Tv[intv+1]+pt[1]*(Tv[intv+1]-Tv[intv]);
  v *= 0.5;
  
  // Compute the basis values in each direction
  FElibrary::bspline_basis(Nu, intu, u, Tu, order, work);
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

/*
  Compute the shape functions and their derivatives w.r.t. the
  parametric element location 
*/
template <int order>
void PlaneStressBsplineAll<order>::getShapeFunctions( const double pt[], 
                                                      double N[],
                                                      double Na[],
                                                      double Nb[] ){
  // Compute the knot intervals
  int intu = t1+order-1;
  int intv = t2+order-1;
  
  // Convert parametric point to knot vector point
  double dNu[8], dNv[8], work[30];
  double u = Tu[intu]+Tu[intu+1]+pt[0]*(Tu[intu+1]-Tu[intu]);
  u *= 0.5;
  double v = Tv[intv]+Tv[intv+1]+pt[1]*(Tv[intv+1]-Tv[intv]);
  v *= 0.5;
  
  // Compute the basis values and its derivative in each direction
  FElibrary::bspline_basis_derivative(dNu, intu, u, 1, Tu, order, work);
  FElibrary::bspline_basis_derivative(dNv, intv, v, 1, Tv, order, work);

  // Compute the 2d shape functions and the partial derivative along
  // each direction
  if (order == 3){
    /* printf("u: %e Nu: %e %e %e \n", u, dNu[0], dNu[1], */
    /*        dNu[2]); */
    /* printf("v: %e Nv: %e %e %e \n", v, dNv[0], dNv[1], */
    /*        dNv[2]); */

    N[0] = dNu[0]*dNv[0];
    N[1] = dNu[1]*dNv[0];
    N[2] = dNu[2]*dNv[0];

    N[3] = dNu[0]*dNv[1];
    N[4] = dNu[1]*dNv[1];
    N[5] = dNu[2]*dNv[1];
    
    N[6] = dNu[0]*dNv[2];
    N[7] = dNu[1]*dNv[2];
    N[8] = dNu[2]*dNv[2];
    /* for (int i = 0; i < 9; i++){ */
    /*   printf("N[%d]: %e \n", i, N[i]); */
    /* } */
    
    // Compute the partial derivative wrt u
    Na[0] = dNu[3]*dNv[0]*(Tu[intu+1]-Tu[intu])*.5;
    Na[1] = dNu[4]*dNv[0]*(Tu[intu+1]-Tu[intu])*.5;
    Na[2] = dNu[5]*dNv[0]*(Tu[intu+1]-Tu[intu])*.5;

    Na[3] = dNu[3]*dNv[1]*(Tu[intu+1]-Tu[intu])*.5;
    Na[4] = dNu[4]*dNv[1]*(Tu[intu+1]-Tu[intu])*.5;
    Na[5] = dNu[5]*dNv[1]*(Tu[intu+1]-Tu[intu])*.5;
    
    Na[6] = dNu[3]*dNv[2]*(Tu[intu+1]-Tu[intu])*.5;
    Na[7] = dNu[4]*dNv[2]*(Tu[intu+1]-Tu[intu])*.5;
    Na[8] = dNu[5]*dNv[2]*(Tu[intu+1]-Tu[intu])*.5;
    // Compute the derivative wrt v
    Nb[0] = dNu[0]*dNv[3]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[1] = dNu[1]*dNv[3]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[2] = dNu[2]*dNv[3]*(Tv[intv+1]-Tv[intv])*.5;

    Nb[3] = dNu[0]*dNv[4]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[4] = dNu[1]*dNv[4]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[5] = dNu[2]*dNv[4]*(Tv[intv+1]-Tv[intv])*.5;
    
    Nb[6] = dNu[0]*dNv[5]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[7] = dNu[1]*dNv[5]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[8] = dNu[2]*dNv[5]*(Tv[intv+1]-Tv[intv])*.5;    
  }
  else if (order == 4){
    /* printf("u: %e Nu: %e %e %e %e \n", u, dNu[0], dNu[1], */
    /*        dNu[2],dNu[3]); */
    /* printf("v: %e Nv: %e %e %e %e \n", v, dNv[0], dNv[1], */
    /*        dNv[2],dNv[3]); */
    N[0] = dNu[0]*dNv[0];
    N[1] = dNu[1]*dNv[0];
    N[2] = dNu[2]*dNv[0];
    N[3] = dNu[3]*dNv[0];
    
    N[4] = dNu[0]*dNv[1];
    N[5] = dNu[1]*dNv[1];
    N[6] = dNu[2]*dNv[1];
    N[7] = dNu[3]*dNv[1];
    
    N[8] = dNu[0]*dNv[2];
    N[9] = dNu[1]*dNv[2];
    N[10] = dNu[2]*dNv[2];
    N[11] = dNu[3]*dNv[2];
    
    N[12] = dNu[0]*dNv[3];
    N[13] = dNu[1]*dNv[3];
    N[14] = dNu[2]*dNv[3];
    N[15] = dNu[3]*dNv[3];
    
    /* for (int i = 0; i < 16; i++){ */
    /*   printf("N[%d]: %e \n", i, N[i]); */
    /* } */
    // Compute the partial derivative wrt u
    Na[0] = dNu[4]*dNv[0]*(Tu[intu+1]-Tu[intu])*.5;
    Na[1] = dNu[5]*dNv[0]*(Tu[intu+1]-Tu[intu])*.5;
    Na[2] = dNu[6]*dNv[0]*(Tu[intu+1]-Tu[intu])*.5;
    Na[3] = dNu[7]*dNv[0]*(Tu[intu+1]-Tu[intu])*.5;
    
    Na[4] = dNu[4]*dNv[1]*(Tu[intu+1]-Tu[intu])*.5;
    Na[5] = dNu[5]*dNv[1]*(Tu[intu+1]-Tu[intu])*.5;
    Na[6] = dNu[6]*dNv[1]*(Tu[intu+1]-Tu[intu])*.5;
    Na[7] = dNu[7]*dNv[1]*(Tu[intu+1]-Tu[intu])*.5;
    
    Na[8] = dNu[4]*dNv[2]*(Tu[intu+1]-Tu[intu])*.5;
    Na[9] = dNu[5]*dNv[2]*(Tu[intu+1]-Tu[intu])*.5;
    Na[10] = dNu[6]*dNv[2]*(Tu[intu+1]-Tu[intu])*.5;
    Na[11] = dNu[7]*dNv[2]*(Tu[intu+1]-Tu[intu])*.5;
    
    Na[12] = dNu[4]*dNv[3]*(Tu[intu+1]-Tu[intu])*.5;
    Na[13] = dNu[5]*dNv[3]*(Tu[intu+1]-Tu[intu])*.5;
    Na[14] = dNu[6]*dNv[3]*(Tu[intu+1]-Tu[intu])*.5;
    Na[15] = dNu[7]*dNv[3]*(Tu[intu+1]-Tu[intu])*.5;
    
    // Compute the derivative wrt v
    Nb[0] = dNu[0]*dNv[4]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[1] = dNu[1]*dNv[4]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[2] = dNu[2]*dNv[4]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[3] = dNu[3]*dNv[4]*(Tv[intv+1]-Tv[intv])*.5;

    Nb[4] = dNu[0]*dNv[5]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[5] = dNu[1]*dNv[5]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[6] = dNu[2]*dNv[5]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[7] = dNu[3]*dNv[5]*(Tv[intv+1]-Tv[intv])*.5;
    
    Nb[8] = dNu[0]*dNv[6]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[9] = dNu[1]*dNv[6]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[10] = dNu[2]*dNv[6]*(Tv[intv+1]-Tv[intv])*.5;   
    Nb[11] = dNu[3]*dNv[6]*(Tv[intv+1]-Tv[intv])*.5;   
    
    Nb[12] = dNu[0]*dNv[7]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[13] = dNu[1]*dNv[7]*(Tv[intv+1]-Tv[intv])*.5;
    Nb[14] = dNu[2]*dNv[7]*(Tv[intv+1]-Tv[intv])*.5;   
    Nb[15] = dNu[3]*dNv[7]*(Tv[intv+1]-Tv[intv])*.5;   
  }  
}

/*
  Get the number of elemens/nodes and CSR size of the contributed by
  this element.  
*/
template <int order>
void PlaneStressBsplineAll<order>::addOutputCount( int *nelems, 
                                                   int *nnodes, 
                                                   int *ncsr ){  
  *nelems += (order-1)*(order-1);
  *nnodes += order*order;
  *ncsr += 4*(order-1)*(order-1);
  /*
   *nelems += (order-2)*(order-2);
  *nnodes += (order-1)*(order-1);
  *ncsr += 4*(order-2)*(order-2);
  */
}

/*
  Get the output data from this element and place it in a real
  array for visualization later. The values generated for visualization
  are determined by a bit-wise selection variable 'out_type' which is 
  can be used to simultaneously write out different data. Note that this
  is why the bitwise operation & is used below. 

  The output may consist of the following:
  - the nodal locations
  - the displacements and rotations
  - the strains or strains within the element
  - extra variables that are used for optimization
  
  output:
  data:     the data to write to the file (eventually)

  input:
  out_type: the bit-wise variable used to specify what data to generate
  vars:     the element variables
  Xpts:     the element nodal locations
*/
template <int order>
void PlaneStressBsplineAll<order>::getOutputData( unsigned int out_type, 
                                                  double *data, int ld_data,
                                                  const TacsScalar Xpts[],
                                                  const TacsScalar vars[] ){
  /*  
  for (int m = 0; m < order-1; m++){
    for (int n = 0; n < order-1; n++){
      int p = n+(order-1)*m;
      int index = 0;
      // Compute the integration points in each direction
      double pt[2];

      pt[0] = -1.0 + 2.0*n/(order-2.0);
      pt[1] = -1.0 + 2.0*m/(order-2.0);
      */

  for (int m = 0; m < order; m++){
    for (int n = 0; n < order; n++){
      int p = n+(order-1)*m;
      int index = 0;
      // Compute the integration points in each direction
      double pt[2];

      pt[0] = -1.0 + 2.0*n/(order-1.0);
      pt[1] = -1.0 + 2.0*m/(order-1.0);

      // Compute the shape functions
      double N[NUM_NODES];
      double Na[NUM_NODES], Nb[NUM_NODES];
      getShapeFunctions(pt, N, Na, Nb);

      TacsScalar Xpts_new[3];
      if (out_type & TACSElement::OUTPUT_NODES){
        getXpts(Xpts_new, N, Xpts);
        for (int k = 0; k < 3; k++){
          data[index+k] = RealPart(Xpts_new[k]);
        }
        index += 3;
      }
      TacsScalar vars_new[2];
      if (out_type & TACSElement::OUTPUT_DISPLACEMENTS){
        this->getDisplacement(vars_new, N, vars);
        for (int k = 0; k < 2; k++){
          data[index+k] = RealPart(vars_new[k]);
        }
        index += 2;
      }
      // Compute the directive of X wrt to coordinate directions
      TacsScalar X[3], Xa[4];
      this->planeJacobian(X, Xa, N, Na, Nb, Xpts);
      // Compute the determinant of Xa and its transformation
      TacsScalar J[4];
      FElibrary::jacobian2d(Xa, J);

      // Compute the strains
      TacsScalar strain[3];
      this->evalStrain(strain, J, Na, Nb, vars);
        
      if (out_type & TACSElement::OUTPUT_STRAINS){
        for (int k = 0; k < 3; k++){
          data[index+k] = RealPart(strain[k]);
        }
        index += 3;
      }
      // Compute the stresses
      if (out_type & TACSElement::OUTPUT_STRESSES){
        TacsScalar stress[3];
        this->stiff->calculateStress(pt, strain, stress);
        for (int k = 0; k < 3; k++){
          data[index+k] = RealPart(stress[k]);
        }
        index += 3;
      }
      // Compute failure and buckling
      if (out_type & TACSElement::OUTPUT_EXTRAS){
        TacsScalar lambda;
        this->stiff->failure(pt, strain, &lambda);
        data[index] = RealPart(lambda);
          
        this->stiff->buckling(strain, &lambda);
        data[index+1] = RealPart(lambda);

        data[index+2] = RealPart(this->stiff->getDVOutputValue(0, pt));
        data[index+3] = RealPart(this->stiff->getDVOutputValue(1, pt));
          
        index += this->NUM_EXTRAS;
      }
      data += ld_data;
    }
  }
}
/*
  Get the element connectivity for visualization purposes. Since each
  element consists of a series of sub-elements used for visualization,
  we also need the connectivity of these visualization elements.

  output:
  con:  the connectivity of the local visualization elements contributed
  by this finite-element

  input:
  node:  the node offset number - so that this connectivity is more or 
  less global
*/
template <int order>
void PlaneStressBsplineAll<order>::getOutputConnectivity( int *con, 
                                                          int node ){
  int p = 0;
  for ( int m = 0; m < order-1; m++ ){
    for ( int n = 0; n < order-1; n++ ){
      con[4*p]   = node + n   + m*order;
      con[4*p+1] = node + n+1 + m*order; 
      con[4*p+2] = node + n+1 + (m+1)*order;
      con[4*p+3] = node + n   + (m+1)*order;
      p++;
    }
  }
  /*
  int p = 0;
  for ( int m = 0; m < order-2; m++ ){
    for ( int n = 0; n < order-2; n++ ){
      con[4*p] = (order-1)*(order-1)*pNum+n+(order-1)*m;
      con[4*p+1] = (order-1)*(order-1)*pNum+n+1+(order-1)*m;
      con[4*p+2] = (order-1)*(order-1)*pNum+n+1+(order-1)*(m+1);
      con[4*p+3] = (order-1)*(order-1)*pNum+n+(order-1)*(m+1);
      printf("pNum: %d %d, %d %d %d %d\n", pNum, node, 
             con[4*p], con[4*p+1], con[4*p+2],con[4*p+3]);
      p++;
    }
  }
  */
}

/*
  Get the interpolated positions
*/
template <int order>
void PlaneStressBsplineAll<order>::getXpts( TacsScalar U[], 
                                            const double N[],
                                            const TacsScalar Xpts[] ){
  // Zero the position values
  U[0] = U[1] = U[2] = 0.0;
  for (int i = 0; i < NUM_NODES; i++){
    U[0] += N[0]*Xpts[0];
    U[1] += N[0]*Xpts[1];
    U[2] += N[0]*Xpts[2];
    
    N++;
    Xpts+=3;
  }
}

#endif //TACS_PLANE_STRESS_BSPLINE_ALL_H 
