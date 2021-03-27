
#ifndef TACS_SHELL_TRACTION_H
#define TACS_SHELL_TRACTION_H

#include "TACSShellElementBasis.h"
#include "TACSElementAlgebra.h"

template <int vars_per_node, class quadrature, class basis>
class TACSShellTraction : public TACSElement {
 public:
  TACSShellTraction( const TacsScalar _t[] ){
    for ( int i = 0; i < 3*basis::NUM_NODES; i++ ){
      t[i] = _t[i];
    }
  }

  int getVarsPerNode(){
    return vars_per_node;
  }
  int getNumNodes(){
    return basis::NUM_NODES;
  }

  ElementLayout getLayoutType(){
    return basis::getLayoutType();
  }

  int getNumQuadraturePoints(){
    return quadrature::getNumQuadraturePoints();
  }

  double getQuadratureWeight( int n ){
    return quadrature::getQuadratureWeight(n);
  }

  double getQuadraturePoint( int n, double pt[] ){
    return quadrature::getQuadraturePoint(n, pt);
  }

  int getNumElementFaces(){
    return quadrature::getNumElementFaces();
  }

  int getNumFaceQuadraturePoints( int face ){
    return quadrature::getNumFaceQuadraturePoints(face);
  }

  double getFaceQuadraturePoint( int face, int n, double pt[],
                                 double tangent[] ){
    return quadrature::getFaceQuadraturePoint(face, n, pt, tangent);
  }

  void addResidual( int elemIndex,
                    double time,
                    const TacsScalar *Xpts,
                    const TacsScalar *vars,
                    const TacsScalar *dvars,
                    const TacsScalar *ddvars,
                    TacsScalar *res ){
    // Compute the number of quadrature points
    const int nquad = quadrature::getNumQuadraturePoints();

    // Compute the node normal directions
    TacsScalar fn[3*basis::NUM_NODES];
    getNodeNormals(Xpts, fn);

    // Loop over each quadrature point and add the residual contribution
    for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
      // Get the quadrature weight
      double pt[3];
      double weight = quadrature::getQuadraturePoint(quad_index, pt);

      TacsScalar Xxi[6], n[3];
      basis::interpFields(pt, 3, fn, 3, n);
      basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);

      // Assemble the terms Xd = [Xxi; n] and Xdz
      TacsScalar Xd[9];
      assembleFrame(Xxi, n, Xd);

      // Compute the inverse of the 3x3 Jacobian transformation
      TacsScalar detXd = det3x3(Xd);
      detXd *= weight;

      // Interpolate the traction
      TacsScalar tr[3];
      basis::interpFields(pt, 3, t, 3, tr);

      // Scale the traction
      tr[0] *= -detXd;
      tr[1] *= -detXd;
      tr[2] *= -detXd;

      basis::addInterpFieldsTranspose(pt, 3, tr, vars_per_node, res);
    }
  }

  void getNodeNormals( const TacsScalar Xpts[],
                       TacsScalar fn[],
                       TacsScalar fnorm[]=NULL ){
    for ( int i = 0; i < basis::NUM_NODES; i++ ){
      double pt[2];
      basis::getNodePoint(i, pt);

      // Compute the derivative X,xi at each node
      TacsScalar Xxi[6];
      basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);

      TacsScalar a[3], b[3];
      a[0] = Xxi[0];
      a[1] = Xxi[2];
      a[2] = Xxi[4];

      b[0] = Xxi[1];
      b[1] = Xxi[3];
      b[2] = Xxi[5];

      // Compute the normal direction at the point
      crossProduct(a, b, &fn[3*i]);

      // Compute the 2-norm of the vector in the normal direction
      TacsScalar norm = sqrt(vec3Dot(&fn[3*i], &fn[3*i]));

      // Save the 2-norm value if the fnorm argument is not NULL
      if (fnorm){
        fnorm[i] = norm;
      }

      // Scale the normal direction
      if (norm != 0.0){
        vec3Scale(1.0/norm, &fn[3*i]);
      }
    }
  }

  static void assembleFrame( const TacsScalar Xxi[],
                             const TacsScalar n[],
                             TacsScalar Xd[] ){
    Xd[0] = Xxi[0];
    Xd[1] = Xxi[1];
    Xd[2] = n[0];

    Xd[3] = Xxi[2];
    Xd[4] = Xxi[3];
    Xd[5] = n[1];

    Xd[6] = Xxi[4];
    Xd[7] = Xxi[5];
    Xd[8] = n[2];
  }

 private:
  TacsScalar t[3*basis::NUM_NODES];
};

#endif // TACS_SHELL_TRACTION_H