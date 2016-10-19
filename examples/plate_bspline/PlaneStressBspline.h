#ifndef TACS_PLANE_STRESS_BSPLINE_H
#define TACS_PLANE_STRESS_BSPLINE_H

#include "TACS2DElement.h"
#include "FElibrary.h"

class PlaneStressBspline : public TACS2DElement<9> {
public:
 PlaneStressBspline( PlaneStressStiffness *_stiff,
                     double _Tu[], double _Tv[], 
                     int _Lu, int _Lv, 
                     ElementBehaviorType type=LINEAR,
                     int _componentNum=0,
                     int _pNum=0):
  TACS2DElement<9>(_stiff, type, _componentNum){
    // Number of knot intervals
    Lu = _Lu;
    Lv = _Lv;
    
    // Copy over the patch number
    pNum = _pNum;
    memcpy(Tu, _Tu, (Lu+1+4)*sizeof(double));
    memcpy(Tv, _Tv, (Lv+1+4)*sizeof(double));
    numGauss = FElibrary::getGaussPtsWts(3, &gaussPts, &gaussWts);
    // Compute the interval at which the patch is in
    t2 = int(floor(pNum/(Lu))); // y direction
    t1 = int(floor(pNum-t2*(Lu))); // x direction
  }  
  ~PlaneStressBspline(){}
 
  // Return the name of this element
  // -------------------------------
  const char *elementName(){ return elemName; }  // Retrieve the shape functions
  // ----------------------------
  void getShapeFunctions( const double pt[], double N[]){
    // Compute the knot intervals
    int intu = t1+2;
    int intv = t2+2;
    // Convert parametric point to knot vector point    
    double Nu[3], Nv[3], work[15];
    double u = Tu[intu]+Tu[intu+1]+pt[0]*(Tu[intu+1]-Tu[intu]);
    u *= 0.5;
    double v = Tv[intv]+Tv[intv+1]+pt[1]*(Tv[intv+1]-Tv[intv]);
    v *= 0.5;
    
    // Compute the basis values in each direction
    FElibrary::bspline_basis(Nu, intu, u, Tu, 3, work);
    FElibrary::bspline_basis(Nv, intv, v, Tv, 3, work);

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
  void getShapeFunctions( const double pt[], double N[],
                          double Na[], double Nb[] ){
    // Compute the knot intervals
    int intu = t1+2;
    int intv = t2+2;

    // Convert parametric point to knot vector point
    double dNu[6], dNv[6], work[15];
    double u = Tu[intu]+Tu[intu+1]+pt[0]*(Tu[intu+1]-Tu[intu]);
    u *= 0.5;
    double v = Tv[intv]+Tv[intv+1]+pt[1]*(Tv[intv+1]-Tv[intv]);
    v *= 0.5;
        
    // Compute the basis values and its derivative in each direction
    FElibrary::bspline_basis_derivative(dNu, intu, u, 1, Tu, 3, work);
    FElibrary::bspline_basis_derivative(dNv, intv, v, 1, Tv, 3, work);
    
    // Compute the 2d shape functions
    N[0] = dNu[0]*dNv[0];
    N[1] = dNu[1]*dNv[0];
    N[2] = dNu[2]*dNv[0];

    N[3] = dNu[0]*dNv[1];
    N[4] = dNu[1]*dNv[1];
    N[5] = dNu[2]*dNv[1];
    
    N[6] = dNu[0]*dNv[2];
    N[7] = dNu[1]*dNv[2];
    N[8] = dNu[2]*dNv[2];
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
  };
  // Retrieve the Gauss points/weights
  // ---------------------------------
  int getNumGaussPts(){
    return 9;
  }
  double getGaussWtsPts( const int num, double pt[] ){
    int m = (int)((num)/numGauss);
    int n = num-numGauss*m;

    pt[0] = gaussPts[n];
    pt[1] = gaussPts[m];

    return gaussWts[n]*gaussWts[m];
  };  
  // Functions for post-processing
  // -----------------------------
  void addOutputCount( int *nelems, int *nnodes, int *ncsr ){
    *nelems += 1;
    *nnodes += 4;
    *ncsr += 4;
    
  }
  void getOutputData( unsigned int out_type,
                      double *data, int ld_data,
                      const TacsScalar Xpts[],
                      const TacsScalar vars[] ){
    
    for (int m = 0; m < 2; m++){
      for (int n = 0; n < 2; n++){
        int p = n+2*m;
        int index = 0;
        // Compute the integration points in each direction
        double pt[2];
        pt[0] = -1.0+2.0*n;
        pt[1] = -1.0+2.0*m;
        
        // Compute the shape functions
        double N[9];
        double Na[9], Nb[9];
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
  void getOutputConnectivity( int *con, int node ){
    int p = 0;
    for (int m = 0; m < 1; m++){
      for (int n = 0; n < 1; n++){
        con[4*p] = 4*pNum+n+2*m;
        con[4*p+1] = 4*pNum+n+1+2*m;
        con[4*p+2] = 4*pNum+n+1+2*(m+1);
        con[4*p+3] = 4*pNum+n+2*(m+1);
        p++;
      }
    }
  } 
  void getXpts( TacsScalar U[], 
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
 private:
  double Tu[8], Tv[8];
  static const int NUM_NODES = 9;
  static const char *elemName;
  // The Gauss quadrature scheme
  int numGauss;
  const double *gaussWts, *gaussPts;
  // Patch number
  int pNum;
  // Interval number in x,y
  int t1, t2;
  int Lu, Lv;
};

const char *PlaneStressBspline::elemName = "PlaneStressBspline";

#endif // TACS_PLANE_STRESS_BSPLINE_H
