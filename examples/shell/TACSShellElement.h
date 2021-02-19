#ifndef TACS_SHELL_ELEMENT_H
#define TACS_SHELL_ELEMENT_H

#include "TACSShellElementModel.h"
#include "TACSGaussQuadrature.h"
#include "TACSElementAlgebra.h"
#include "TACSShellConstitutive.h"
#include "TACSElement.h"
#include "TACSElementTypes.h"

/*
  Compute the transformation from the local coordinates
  to
*/
class TACSShellTransform : public TACSObject {
 public:
  /*
    Given the local shell element reference frame Xf, compute the
    transformation from the global coordinates to the shell-aligned local axis.
  */
  virtual void computeTransform( const TacsScalar Xxi[], TacsScalar T[] ) = 0;
};

class TACSShellNaturalTransform : public TACSShellTransform {
 public:
  TACSShellNaturalTransform(){}

  void computeTransform( const TacsScalar Xxi[], TacsScalar T[] ){
    // Compute the transformation
    TacsScalar a[3], b[3];
    a[0] = Xxi[0];
    a[1] = Xxi[2];
    a[2] = Xxi[4];

    b[0] = Xxi[1];
    b[1] = Xxi[3];
    b[2] = Xxi[5];

    // Compute the normal direction
    TacsScalar n[3];
    crossProduct(a, b, n);

    // Normalize the normal direction
    TacsScalar invNorm = 1.0/sqrt(vec3Dot(n, n));
    vec3Scale(invNorm, n);

    // Normalize the 1-direction of the element
    TacsScalar inv = 1.0/sqrt(vec3Dot(a, a));
    vec3Scale(inv, a);

    // Take the cross product to determine the 2-direction
    crossProduct(n, a, b);

    // Set the components of the transformation
    T[0] = a[0];
    T[3] = a[1];
    T[6] = a[2];

    T[1] = b[0];
    T[4] = b[1];
    T[7] = b[2];

    T[2] = n[0];
    T[5] = n[1];
    T[8] = n[2];
  }
};

/*
  The director class.

  Given a reference vector, t, from the element geometry, the director computes
  the exact or approximate rate of change of the displacement t.
*/
class TACSLinearizedRotation {
 public:

  static const int NUM_PARAMETERS = 3;

  /**
    Compute the director at a point.

    d = Q(q)*t = (C(q)^{T} - I)*t

    @param q The input rotation parametrization
    @param t The reference direction
    @param d The director values
  */
  static void computeDirector( const TacsScalar q[],
                               const TacsScalar t[],
                               TacsScalar d[] ){
    crossProduct(q, t, d);
  }

  /*
    Compute the director and rates at a point.

    d = Q(q)*t = (C(q)^{T} - I)*t
    ddot = d/dt(Q(q))*t
    dddot = d^2/dt^2(Q(q))*t

    @param q The input rotation parametrization
    @param t The reference direction
    @param d The director values
  */
  static void computeDirectorRates( const TacsScalar q[],
                                    const TacsScalar qdot[],
                                    const TacsScalar t[],
                                    TacsScalar d[],
                                    TacsScalar ddot[] ){
    crossProduct(q, t, d);
    crossProduct(qdot, t, ddot);
  }

  /*
    Compute the director and rates at a point.

    d = Q(q)*t = (C(q)^{T} - I)*t
    ddot = d/dt(Q(q))*t
    dddot = d^2/dt^2(Q(q))*t

    @param q The input rotation parametrization
    @param t The reference direction
    @param d The director values
  */
  static void computeDirectorRates( const TacsScalar q[],
                                    const TacsScalar qdot[],
                                    const TacsScalar qddot[],
                                    const TacsScalar t[],
                                    TacsScalar d[],
                                    TacsScalar ddot[],
                                    TacsScalar dddot[] ){
    crossProduct(q, t, d);
    crossProduct(qdot, t, ddot);
    crossProduct(qddot, t, dddot);
  }

  /*
    Given the derivatives of the kinetic energy expression with respect to time,
    add the contributions to the derivative of the

    Given the partial derivatives of the Lagrangian with respect to the
    director and the time derivative of the vector, compute

    ddtdTddot = d/dt(dT/d(ddot))
    dTddot = dT/d(dot)

    Compute:

    res += scale*(d/dt(dT/d(ddot))*d(ddot)/d(qdot) + dT/d(ddot)*d/dt(d(ddot)/d(qdot)))

  */
  static void addDirectorResidual( const TacsScalar q[],
                                    const TacsScalar qdot[],
                                    const TacsScalar qddot[],
                                    const TacsScalar t[],
                                    const TacsScalar dLdd[],
                                    TacsScalar res[] ){
    crossProductAdd(1.0, t, dLdd, res);
  }
};

/**
  Defines the quadrature over both the face and quadrature
*/
class TACSQuadQuadrature {
 public:
  static int getNumParameters(){
    return 2;
  }
  static int getNumQuadraturePoints(){
    return 9;
  }
  static double getQuadratureWeight( int n ){
    return TacsGaussQuadWts3[n % 3]*TacsGaussQuadWts3[n / 3];
  }
  static double getQuadraturePoint( int n, double pt[] ){
    pt[0] = TacsGaussQuadPts3[n % 3];
    pt[1] = TacsGaussQuadPts3[n / 3];

    return TacsGaussQuadWts3[n % 3]*TacsGaussQuadWts3[n/3];
  }
  static int getNumFaces(){
    return 4;
  }
  static int getNumFaceQuadraturePoints( int face ){
    return 3;
  }
  static double getFaceQuadraturePoint( int face, int n,
                                        double pt[],
                                        double t[] ){
    if (face/2 == 0){
      pt[0] = -1.0 + 2.0*(face % 2);
      pt[1] = TacsGaussQuadPts3[n];
    }
    else {
      pt[0] = TacsGaussQuadPts3[n];
      pt[1] = -1.0 + 2.0*(face % 2);
    }

    if (face == 0){
      // -X edge
      t[0] = 0.0;  t[1] = -1.0;
    }
    else if (face == 1){
      // +X edge
      t[0] = 0.0;  t[1] = 1.0;
    }
    else if (face == 2){
      // -Y edge
      t[0] = 1.0;  t[1] = 0.0;
    }
    else if (face == 3){
      // +Y edge
      t[0] = -1.0;  t[1] = 0.0;
    }

    return TacsGaussQuadWts3[n];
  }
};

class TACSShellQuadQuadraticBasis {
 public:
  static const int NUM_NODES = 9;

  // Set the number of tying points for each of the 5 components
  // of the tying strain
  static const int NUM_G11_TYING_POINTS = 6;
  static const int NUM_G22_TYING_POINTS = 6;
  static const int NUM_G12_TYING_POINTS = 4;
  static const int NUM_G13_TYING_POINTS = 6;
  static const int NUM_G23_TYING_POINTS = 6;

  static const int NUM_TYING_POINTS =
    NUM_G11_TYING_POINTS +
    NUM_G22_TYING_POINTS +
    NUM_G12_TYING_POINTS +
    NUM_G13_TYING_POINTS +
    NUM_G23_TYING_POINTS;

  static void getNodePoint( const int n, double pt[] ){
    pt[0] = -1.0 + 1.0*(n % 3);
    pt[1] = -1.0 + 1.0*(n / 3);
  }
  static ElementLayout getLayoutType(){
    return TACS_QUAD_QUADRATIC_ELEMENT;
  }

  static void interpFields( const double pt[],
                            const int vars_per_node,
                            const TacsScalar values[],
                            const int m,
                            TacsScalar field[] ){
    double na[3];
    na[0] = -0.5*pt[0]*(1.0 - pt[0]);
    na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
    na[2] = 0.5*(1.0 + pt[0])*pt[0];

    double nb[3];
    nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
    nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
    nb[2] = 0.5*(1.0 + pt[1])*pt[1];

    for ( int k = 0; k < m; k++ ){
      field[k] = 0.0;
    }

    for ( int j = 0; j < 3; j++ ){
      for ( int i = 0; i < 3; i++ ){
        for ( int k = 0; k < m; k++ ){
          field[k] += na[i]*nb[j]*values[k];
        }
        values += vars_per_node;
      }
    }
  }

  static void addInterpFieldsTranspose( const double pt[],
                                        const int m,
                                        const TacsScalar field[],
                                        const int vars_per_node,
                                        TacsScalar values[] ){
    double na[3];
    na[0] = -0.5*pt[0]*(1.0 - pt[0]);
    na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
    na[2] = 0.5*(1.0 + pt[0])*pt[0];

    double nb[3];
    nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
    nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
    nb[2] = 0.5*(1.0 + pt[1])*pt[1];

    for ( int j = 0; j < 3; j++ ){
      for ( int i = 0; i < 3; i++ ){
        for ( int k = 0; k < m; k++ ){
          values[k] += na[i]*nb[j]*field[k];
        }
        values += vars_per_node;
      }
    }
  }

  static void interpFieldsGrad( const double pt[],
                                const int vars_per_node,
                                const TacsScalar values[],
                                const int m,
                                TacsScalar grad[] ){
    double na[3];
    na[0] = -0.5*pt[0]*(1.0 - pt[0]);
    na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
    na[2] = 0.5*(1.0 + pt[0])*pt[0];

    double nb[3];
    nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
    nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
    nb[2] = 0.5*(1.0 + pt[1])*pt[1];

    double dna[3];
    dna[0] = -0.5 + pt[0];
    dna[1] = -2.0*pt[0];
    dna[2] = 0.5 + pt[0];

    double dnb[3];
    dnb[0] = -0.5 + pt[1];
    dnb[1] = -2.0*pt[1];
    dnb[2] = 0.5 + pt[1];

    for ( int k = 0; k < m; k++ ){
      grad[2*k] = 0.0;
      grad[2*k+1] = 0.0;
    }

    for ( int j = 0; j < 3; j++ ){
      for ( int i = 0; i < 3; i++ ){
        for ( int k = 0; k < m; k++ ){
          grad[2*k]   += dna[i]*nb[j]*values[k];
          grad[2*k+1] += na[i]*dnb[j]*values[k];
        }
        values += vars_per_node;
      }
    }
  }

  static void addInterpFieldsGradTranspose( const double pt[],
                                            const int m,
                                            TacsScalar grad[],
                                            const int vars_per_node,
                                            TacsScalar values[] ){
    double na[3];
    na[0] = -0.5*pt[0]*(1.0 - pt[0]);
    na[1] = (1.0 - pt[0])*(1.0 + pt[0]);
    na[2] = 0.5*(1.0 + pt[0])*pt[0];

    double nb[3];
    nb[0] = -0.5*pt[1]*(1.0 - pt[1]);
    nb[1] = (1.0 - pt[1])*(1.0 + pt[1]);
    nb[2] = 0.5*(1.0 + pt[1])*pt[1];

    double dna[3];
    dna[0] = -0.5 + pt[0];
    dna[1] = -2.0*pt[0];
    dna[2] = 0.5 + pt[0];

    double dnb[3];
    dnb[0] = -0.5 + pt[1];
    dnb[1] = -2.0*pt[1];
    dnb[2] = 0.5 + pt[1];

    for ( int j = 0; j < 3; j++ ){
      for ( int i = 0; i < 3; i++ ){
        for ( int k = 0; k < m; k++ ){
          values[k] += (dna[i]*nb[j]*grad[2*k] + na[i]*dnb[j]*grad[2*k+1]);
        }
        values += vars_per_node;
      }
    }
  }

  static int getNumTyingFields(){
    return 5;
  }
  static int getNumTyingPoints( const int field ){
    if (field == 0){ return 6; }
    else if (field == 1){ return 6; }
    else if (field == 2){ return 4; }
    else if (field == 3){ return 6; }
    else if (field == 4){ return 6; }
    return 0;
  }
  static void getTyingPoint( const int field,
                             const int ty,
                             double pt[] ){
    const double s = 0.774596669241483;
    const double t = 0.577350269189626;

    if (field == 0 || field == 4){ // g11 or g13
      if (ty % 2 == 0){
        pt[0] = -t;
      }
      else {
        pt[0] = t;
      }
      if (ty / 2 == 0){
        pt[1] = -s;
      }
      else if (ty / 2 == 1){
        pt[1] = 0.0;
      }
      else {
        pt[1] = s;
      }
    }
    else if (field == 1 || field == 3){ // g22 or g23
      if (ty / 2 == 0){
        pt[1] = -s;
      }
      else if (ty / 2 == 1){
        pt[1] = 0.0;
      }
      else {
        pt[1] = s;
      }
      if (ty % 2 == 0){
        pt[1] = -t;
      }
      else {
        pt[1] = t;
      }
    }
    else { // (field == 2) g12
      if (ty % 2 == 0){
        pt[0] = -t;
      }
      else {
        pt[0] = t;
      }
      if (ty / 2 == 0){
        pt[1] = -t;
      }
      else {
        pt[1] = t;
      }
    }
  }
  static TacsScalar interpTying( const int field,
                                 const double pt[],
                                 const TacsScalar ety[] ){
    const double s = 0.774596669241483;
    const double t = 0.577350269189626;
    const double tinv = 1.0/t;
    const double sinv = 1.0/(s*s);

    TacsScalar value = 0.0;
    if (field == 0 || field == 4){
      double ntu[2];
      ntu[0] = 0.5*tinv*(t - pt[0]);
      ntu[1] = 0.5*tinv*(t + pt[0]);

      double nv[3];
      nv[0] = 0.5*sinv*pt[1]*(pt[1] - s);
      nv[1] = sinv*(s - pt[1])*(s + pt[1]);
      nv[2] = 0.5*sinv*pt[1]*(s + pt[1]);

      for ( int j = 0; j < 3; j++ ){
        for ( int i = 0; i < 2; i++ ){
          value += ntu[i]*nv[j]*ety[i + 2*j];
        }
      }
    }
    else if (field == 1 || field == 3){
      double nu[3];
      nu[0] = 0.5*sinv*pt[0]*(pt[0] - s);
      nu[1] = sinv*(s - pt[0])*(s + pt[0]);
      nu[2] = 0.5*sinv*pt[0]*(s + pt[0]);

      double ntv[2];
      ntv[0] = 0.5*tinv*(t - pt[1]);
      ntv[1] = 0.5*tinv*(t + pt[1]);

      for ( int j = 0; j < 2; j++ ){
        for ( int i = 0; i < 3; i++ ){
          value += nu[i]*ntv[j]*ety[i + 3*j];
        }
      }
    }
    else { // field == 2
      double ntu[2];
      ntu[0] = 0.5*tinv*(t - pt[0]);
      ntu[1] = 0.5*tinv*(t + pt[0]);

      double ntv[2];
      ntv[0] = 0.5*tinv*(t - pt[1]);
      ntv[1] = 0.5*tinv*(t + pt[1]);

      for ( int j = 0; j < 2; j++ ){
        for ( int i = 0; i < 2; i++ ){
          value += ntu[i]*ntv[j]*ety[i + 2*j];
        }
      }
    }

    return value;
  }

  static void addInterpTyingTranspose( const int field,
                                       const double pt[],
                                       const TacsScalar value,
                                       TacsScalar ety[] ){
    const double s = 0.774596669241483;
    const double t = 0.577350269189626;
    const double tinv = 1.0/t;
    const double sinv = 1.0/(s*s);

    if (field == 0 || field == 4){
      double ntu[2];
      ntu[0] = 0.5*tinv*(t - pt[0]);
      ntu[1] = 0.5*tinv*(t + pt[0]);

      double nv[3];
      nv[0] = 0.5*sinv*pt[1]*(pt[1] - s);
      nv[1] = sinv*(s - pt[1])*(s + pt[1]);
      nv[2] = 0.5*sinv*pt[1]*(s + pt[1]);

      for ( int j = 0; j < 3; j++ ){
        for ( int i = 0; i < 2; i++ ){
          ety[i + 2*j] += ntu[i]*nv[j]*value;
        }
      }
    }
    else if (field == 1 || field == 3){
      double nu[3];
      nu[0] = 0.5*sinv*pt[0]*(pt[0] - s);
      nu[1] = sinv*(s - pt[0])*(s + pt[0]);
      nu[2] = 0.5*sinv*pt[0]*(s + pt[0]);

      double ntv[2];
      ntv[0] = 0.5*tinv*(t - pt[1]);
      ntv[1] = 0.5*tinv*(t + pt[1]);

      for ( int j = 0; j < 2; j++ ){
        for ( int i = 0; i < 3; i++ ){
          ety[i + 3*j] += nu[i]*ntv[j]*value;
        }
      }
    }
    else { // field == 2
      double ntu[2];
      ntu[0] = 0.5*tinv*(t - pt[0]);
      ntu[1] = 0.5*tinv*(t + pt[0]);

      double ntv[2];
      ntv[0] = 0.5*tinv*(t - pt[1]);
      ntv[1] = 0.5*tinv*(t + pt[1]);

      for ( int j = 0; j < 2; j++ ){
        for ( int i = 0; i < 2; i++ ){
          ety[i + 2*j] += ntu[i]*ntv[j]*value;
        }
      }
    }
  }
};


template <class quadrature, class basis, class director, class model>
class TACSShellElement : public TACSElement {
 public:
  TACSShellElement( TACSShellTransform *_transform,
                    TACSShellConstitutive *_con ){
    transform = _transform;
    transform->incref();

    con = _con;
    con->incref();
  }

  int getVarsPerNode(){
    return 3 + director::NUM_PARAMETERS;
  }
  int getNumNodes(){
    return basis::NUM_NODES;
  }

  ElementLayout getLayoutType(){
    return basis::getLayoutType();
  }


  void computeEnergies( int elemIndex,
                        double time,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[],
                        TacsScalar *Te,
                        TacsScalar *Pe );

  void addResidual( int elemIndex,
                    double time,
                    const TacsScalar *Xpts,
                    const TacsScalar *vars,
                    const TacsScalar *dvars,
                    const TacsScalar *ddvars,
                    TacsScalar *res );

  void getOutputData( int elemIndex,
                      ElementType etype,
                      int write_flag,
                      const TacsScalar Xpts[],
                      const TacsScalar vars[],
                      const TacsScalar dvars[],
                      const TacsScalar ddvars[],
                      int ld_data,
                      TacsScalar *data );

 private:
  TACSShellTransform *transform;
  TACSShellConstitutive *con;

  void getNodeNormals( const TacsScalar Xpts[],
                       TacsScalar fn[],
                       TacsScalar fnorm[]=NULL );

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

  static void assembleFrame( const TacsScalar nxi[],
                             TacsScalar Xdz[] ){
    Xdz[0] = nxi[0];
    Xdz[1] = nxi[1];
    Xdz[2] = 0.0;

    Xdz[3] = nxi[2];
    Xdz[4] = nxi[3];
    Xdz[5] = 0.0;

    Xdz[6] = nxi[4];
    Xdz[7] = nxi[5];
    Xdz[8] = 0.0;
  }

  static void extractFrame( const TacsScalar Xd[],
                            TacsScalar Xxi[],
                            TacsScalar n[] ){
    Xxi[0] = Xd[0];
    Xxi[1] = Xd[1];
    n[0] = Xd[2];

    Xxi[2] = Xd[3];
    Xxi[3] = Xd[4];
    n[1] = Xd[5];

    Xxi[4] = Xd[6];
    Xxi[5] = Xd[7];
    n[2] = Xd[8];
  }

  static void extractFrame( const TacsScalar Xd[],
                            TacsScalar Xxi[] ){
    Xxi[0] = Xd[0];
    Xxi[1] = Xd[1];

    Xxi[2] = Xd[3];
    Xxi[3] = Xd[4];

    Xxi[4] = Xd[6];
    Xxi[5] = Xd[7];
  }
};

template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  getNodeNormals( const TacsScalar Xpts[],
                  TacsScalar fn[],
                  TacsScalar fnorm[] ){
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

/*
  Add the residual to the provided vector
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  computeEnergies( int elemIndex,
                   double time,
                   const TacsScalar *Xpts,
                   const TacsScalar *vars,
                   const TacsScalar *dvars,
                   TacsScalar *_Te, TacsScalar *_Ue ){
  // Zero the kinetic and potential energies
  TacsScalar Te = 0.0;
  TacsScalar Ue = 0.0;

  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();
  const int vars_per_node = 3 + director::NUM_PARAMETERS;

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals(Xpts, fn);

  // Compute the frame normal and directors at each node
  TacsScalar d[3*basis::NUM_NODES];
  TacsScalar ddot[3*basis::NUM_NODES];
  for ( int i = 0, offset = 3; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::computeDirectorRates(&vars[offset], &dvars[offset],
                                   &fn[3*i], &d[3*i], &ddot[3*i]);
  }

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<basis>(Xpts, fn, vars_per_node, vars, d, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X,xi = [dX/dxi1 ; dX/dxi2] and n,xi = [dn/dxi1; dn/dxi2]
    TacsScalar X[3], Xxi[6], n[3], nxi[6];
    basis::interpFields(pt, 3, Xpts, 3, X);
    basis::interpFields(pt, 3, fn, 3, n);
    basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
    basis::interpFieldsGrad(pt, 3, fn, 3, nxi);

    // Compute the transformation at the quadrature point
    TacsScalar T[9];
    transform->computeTransform(Xxi, T);

    // Assemble the terms Xd = [Xxi; n] and Xdz
    TacsScalar Xd[9], Xdz[9];
    assembleFrame(Xxi, n, Xd);
    assembleFrame(nxi, Xdz);

    // Compute the inverse of the 3x3 Jacobian transformation
    TacsScalar Xdinv[9];
    TacsScalar detXd = inv3x3(Xd, Xdinv);
    detXd *= weight;

    // Compute negXdinvXdz = -Xdinv*Xdz
    TacsScalar negXdinvXdz[9];
    mat3x3MatMult(Xdinv, Xdz, negXdinvXdz);
    for ( int i = 0; i < 9; i++ ){
      negXdinvXdz[i] *= -1.0;
    }

    // Compute XdinvTT = Xdinv*T^{T}
    TacsScalar XdinvTT[9];
    mat3x3MatTransMult(Xdinv, T, XdinvTT);

    // Compute the director field and the gradient of the director
    // field at the specified point
    TacsScalar d0[3], d0xi[6];
    basis::interpFields(pt, 3, d, 3, d0);
    basis::interpFieldsGrad(pt, 3, d, 3, d0xi);

    // Compute the gradient of the displacement solution at the quadrature points
    TacsScalar u0xi[6];
    basis::interpFieldsGrad(pt, vars_per_node, vars, 3, u0xi);

    // Input: u0xi, d0, d0xi, T, negXdinvXdz, XdinvTT
    // Output: u0x, u1x

    // Compute the transformation u0x = T*ueta*Xdinv*T^{T}
    // u0x = T*u0d*Xdinv*T^{T}
    TacsScalar u0x[9];
    assembleFrame(u0xi, d0, u0x); // Use u0x to store [u0,xi; d0]

    // u1x = T*(u0d*(-Xdinv*Xdz) + u1d)*Xdinv*T^{T}
    TacsScalar u1x[9], tmp[9];
    assembleFrame(d0xi, u1x); // Use u1x to store [d0,xi; 0]
    mat3x3MatMultAdd(u0x, negXdinvXdz, u1x);
    mat3x3MatMult(T, u1x, tmp);
    mat3x3MatMult(tmp, XdinvTT, u1x);

    // Compute u0x = T*[u0,xi ; d]*Xdinv*T^{T}
    mat3x3MatMult(u0x, XdinvTT, tmp);
    mat3x3MatMult(T, tmp, u0x);

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    model::template interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvTT^{T}*gty*XdinvTT
    mat3x3SymmTransformTranspose(XdinvTT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

    Ue += 0.5*detXd*(s[0]*e[0] + s[1]*e[1] + s[2]*e[2] +
                     s[3]*e[3] + s[4]*e[4] + s[5]*e[5] +
                     s[6]*e[6] + s[7]*e[7] + s[8]*e[8]);
  }

  for ( int i = 0; i < basis::NUM_NODES; i++ ){
    Ue += 0.5*vars[5 + 6*i]*vars[5 + 6*i];
  }

  *_Te = Te;
  *_Ue = Ue;
}

/*
  Add the residual to the provided vector
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  addResidual( int elemIndex,
               double time,
               const TacsScalar *Xpts,
               const TacsScalar *vars,
               const TacsScalar *dvars,
               const TacsScalar *ddvars,
               TacsScalar *res ){
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();
  const int vars_per_node = 3 + director::NUM_PARAMETERS;

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals(Xpts, fn);

  // Compute the frame normal and directors at each node
  TacsScalar d[3*basis::NUM_NODES];
  TacsScalar ddot[3*basis::NUM_NODES];
  TacsScalar dddot[3*basis::NUM_NODES];
  for ( int i = 0, offset = 3; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::computeDirectorRates(&vars[offset], &dvars[offset], &ddvars[offset],
                                   &fn[3*i], &d[3*i], &ddot[3*i], &dddot[3*i]);
  }

  // Derivative of the director field
  TacsScalar dd[3*basis::NUM_NODES];
  memset(dd, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<basis>(Xpts, fn, vars_per_node, vars, d, ety);

  // Zero the contributions to the
  TacsScalar dety[basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X,xi = [dX/dxi1 ; dX/dxi2] and n,xi = [dn/dxi1; dn/dxi2]
    TacsScalar X[3], Xxi[6], n[3], nxi[6];
    basis::interpFields(pt, 3, Xpts, 3, X);
    basis::interpFields(pt, 3, fn, 3, n);
    basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
    basis::interpFieldsGrad(pt, 3, fn, 3, nxi);

    // Compute the transformation at the quadrature point
    TacsScalar T[9];
    transform->computeTransform(Xxi, T);

    // Assemble the terms Xd = [Xxi; n] and Xdz
    TacsScalar Xd[9], Xdz[9];
    assembleFrame(Xxi, n, Xd);
    assembleFrame(nxi, Xdz);

    // Compute the inverse of the 3x3 Jacobian transformation
    TacsScalar Xdinv[9];
    TacsScalar detXd = inv3x3(Xd, Xdinv);
    detXd *= weight;

    // Compute negXdinvXdz = -Xdinv*Xdz
    TacsScalar negXdinvXdz[9];
    mat3x3MatMult(Xdinv, Xdz, negXdinvXdz);
    for ( int i = 0; i < 9; i++ ){
      negXdinvXdz[i] *= -1.0;
    }

    // Compute XdinvTT = Xdinv*T^{T}
    TacsScalar XdinvTT[9];
    mat3x3MatTransMult(Xdinv, T, XdinvTT);

    // Compute the director field and the gradient of the director
    // field at the specified point
    TacsScalar d0[3], d0xi[6];
    basis::interpFields(pt, 3, d, 3, d0);
    basis::interpFieldsGrad(pt, 3, d, 3, d0xi);

    // Compute the gradient of the displacement solution at the quadrature points
    TacsScalar u0xi[6];
    basis::interpFieldsGrad(pt, vars_per_node, vars, 3, u0xi);

    // Input: u0xi, d0, d0xi, T, negXdinvXdz, XdinvTT
    // Output: u0x, u1x

    // Compute the transformation u0x = T*ueta*Xdinv*T^{T}
    // u0x = T*u0d*Xdinv*T^{T}
    TacsScalar u0x[9];
    assembleFrame(u0xi, d0, u0x); // Use u0x to store [u0,xi; d0]

    // u1x = T*(u0d*(-Xdinv*Xdz) + u1d)*Xdinv*T^{T}
    TacsScalar u1x[9], tmp[9];
    assembleFrame(d0xi, u1x); // Use u1x to store [d0,xi; 0]
    mat3x3MatMultAdd(u0x, negXdinvXdz, u1x);
    mat3x3MatMult(T, u1x, tmp);
    mat3x3MatMult(tmp, XdinvTT, u1x);

    // Compute u0x = T*[u0,xi ; d]*Xdinv*T^{T}
    mat3x3MatMult(u0x, XdinvTT, tmp);
    mat3x3MatMult(T, tmp, u0x);

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    model::template interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvTT^{T}*gty*XdinvTT
    mat3x3SymmTransformTranspose(XdinvTT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

    // Compute the derivative of the product of the stress and strain
    // with respect to u0x, u1x and e0ty
    TacsScalar du0x[9], du1x[9], de0ty[6];
    model::evalStrainSens(detXd, s, u0x, u1x, du0x, du1x, de0ty);

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgty[6];
    mat3x3SymmTransformTransSens(XdinvTT, de0ty, dgty);

    // Evaluate the tying strain
    model::template addInterpTyingStrainTranspose<basis>(pt, dgty, dety);

    // Compute du0d = T^{T}*(du0x*XdinvTT^{T} + du1x*XdinvTT^{T}*negXdinvXdz^{T})
    TacsScalar du0d[9];
    mat3x3MatTransMult(du1x, XdinvTT, du0d);
    mat3x3MatTransMult(du0d, negXdinvXdz, tmp);
    mat3x3MatTransMultAdd(du0x, XdinvTT, tmp);
    mat3x3TransMatMult(T, tmp, du0d);

    // Compute du1d = T^{T}*du1x*XdinvTT^{T}
    TacsScalar du1d[9];
    mat3x3MatTransMult(du1x, XdinvTT, tmp);
    mat3x3TransMatMult(T, tmp, du1d);

    // du0d = [du0xi; dd0]
    TacsScalar du0xi[6], dd0[3];
    extractFrame(du0d, du0xi, dd0);

    TacsScalar dd0xi[6];
    extractFrame(du1d, dd0xi);

    // Compute the director field and the gradient of the director
    // field at the specified point
    basis::addInterpFieldsTranspose(pt, 3, dd0, 3, dd);
    basis::addInterpFieldsGradTranspose(pt, 3, dd0xi, 3, dd);

    // Compute the gradient of the displacement solution at the quadrature points
    basis::addInterpFieldsGradTranspose(pt, 3, du0xi, vars_per_node, res);
  }

  // Set the total number of tying points needed for this element
  model::template addComputeTyingStrainTranspose<basis>(Xpts, fn, vars_per_node, vars, d, dety, res, dd);

  // Add the contributions to the director field
  for ( int i = 0, offset = 3; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::addDirectorResidual(&vars[offset], &dvars[offset], &ddvars[offset],
                                   &fn[3*i], &dd[3*i], &res[offset]);
  }

  for ( int i = 0; i < basis::NUM_NODES; i++ ){
    res[5 + 6*i] += 0.5*vars[5 + 6*i];
  }
}


/*
  Get the element data for the basis
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  getOutputData( int elemIndex,
                 ElementType etype,
                 int write_flag,
                 const TacsScalar Xpts[],
                 const TacsScalar vars[],
                 const TacsScalar dvars[],
                 const TacsScalar ddvars[],
                 int ld_data,
                 TacsScalar *data ){
  // Get the number of nodes associated with the visualization
  int num_vis_nodes = TacsGetNumVisNodes(basis::getLayoutType());

  // Compute the number of quadrature points
  const int vars_per_node = 3 + director::NUM_PARAMETERS;

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals(Xpts, fn);

  // Compute the frame normal and directors at each node
  TacsScalar d[3*basis::NUM_NODES];
  TacsScalar ddot[3*basis::NUM_NODES];
  for ( int i = 0, offset = 3; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::computeDirectorRates(&vars[offset], &dvars[offset],
                                   &fn[3*i], &d[3*i], &ddot[3*i]);
  }

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<basis>(Xpts, fn, vars_per_node, vars, d, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int index = 0; index < num_vis_nodes; index++ ){
    // Get the quadrature weight
    double pt[3];
    basis::getNodePoint(index, pt);

    // Compute X,xi = [dX/dxi1 ; dX/dxi2] and n,xi = [dn/dxi1; dn/dxi2]
    TacsScalar X[3], Xxi[6], n[3], nxi[6];
    basis::interpFields(pt, 3, Xpts, 3, X);
    basis::interpFields(pt, 3, fn, 3, n);
    basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
    basis::interpFieldsGrad(pt, 3, fn, 3, nxi);

    // Compute the transformation at the quadrature point
    TacsScalar T[9];
    transform->computeTransform(Xxi, T);

    // Assemble the terms Xd = [Xxi; n] and Xdz
    TacsScalar Xd[9], Xdz[9];
    assembleFrame(Xxi, n, Xd);
    assembleFrame(nxi, Xdz);

    // Compute the inverse of the 3x3 Jacobian transformation
    TacsScalar Xdinv[9];
    inv3x3(Xd, Xdinv);

    // Compute negXdinvXdz = -Xdinv*Xdz
    TacsScalar negXdinvXdz[9];
    mat3x3MatMult(Xdinv, Xdz, negXdinvXdz);
    for ( int i = 0; i < 9; i++ ){
      negXdinvXdz[i] *= -1.0;
    }

    // Compute XdinvTT = Xdinv*T^{T}
    TacsScalar XdinvTT[9];
    mat3x3MatTransMult(Xdinv, T, XdinvTT);

    // Compute the director field and the gradient of the director
    // field at the specified point
    TacsScalar d0[3], d0xi[6];
    basis::interpFields(pt, 3, d, 3, d0);
    basis::interpFieldsGrad(pt, 3, d, 3, d0xi);

    // Compute the gradient of the displacement solution at the quadrature points
    TacsScalar u0xi[6];
    basis::interpFieldsGrad(pt, vars_per_node, vars, 3, u0xi);

    // Input: u0xi, d0, d0xi, T, negXdinvXdz, XdinvTT
    // Output: u0x, u1x

    // Compute the transformation u0x = T*ueta*Xdinv*T^{T}
    // u0x = T*u0d*Xdinv*T^{T}
    TacsScalar u0x[9];
    assembleFrame(u0xi, d0, u0x); // Use u0x to store [u0,xi; d0]

    // u1x = T*(u0d*(-Xdinv*Xdz) + u1d)*Xdinv*T^{T}
    TacsScalar u1x[9], tmp[9];
    assembleFrame(d0xi, u1x); // Use u1x to store [d0,xi; 0]
    mat3x3MatMultAdd(u0x, negXdinvXdz, u1x);
    mat3x3MatMult(T, u1x, tmp);
    mat3x3MatMult(tmp, XdinvTT, u1x);

    // Compute u0x = T*[u0,xi ; d]*Xdinv*T^{T}
    mat3x3MatMult(u0x, XdinvTT, tmp);
    mat3x3MatMult(T, tmp, u0x);

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    model::template interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvTT^{T}*gty*XdinvTT
    mat3x3SymmTransformTranspose(XdinvTT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

    if (etype == TACS_BEAM_OR_SHELL_ELEMENT){
      if (write_flag & TACS_OUTPUT_NODES){
        data[0] = X[0];
        data[1] = X[1];
        data[2] = X[2];
        data += 3;
      }
      if (write_flag & TACS_OUTPUT_DISPLACEMENTS){
        int len = vars_per_node;
        if (len > 6){
          len = 6;
        }
        for ( int i = 0; i < len; i++ ){
          data[i] = vars[i + vars_per_node*index];
        }
        for ( int i = len; i < 6; i++ ){
          data[i] = 0.0;
        }
        data += 6;
      }
      if (write_flag & TACS_OUTPUT_STRAINS){
        for ( int i = 0; i < 9; i++ ){
          data[i] = e[i];
        }
        data += 9;
      }
      if (write_flag & TACS_OUTPUT_STRESSES){
        for ( int i = 0; i < 9; i++ ){
          data[i] = s[i];
        }
        data += 9;
      }
      if (write_flag & TACS_OUTPUT_EXTRAS){
        data[0] = con->evalFailure(elemIndex, pt, X, e);
        data[1] = con->evalDesignFieldValue(elemIndex, pt, X, 0);
        data[2] = con->evalDesignFieldValue(elemIndex, pt, X, 1);
        data[3] = con->evalDesignFieldValue(elemIndex, pt, X, 2);
        data += 4;
      }
    }
  }
}

#endif // TACS_SHELL_ELEMENT_H
