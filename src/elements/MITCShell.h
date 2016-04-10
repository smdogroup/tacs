#ifndef TACS_MITC_SHELL_H
#define TACS_MITC_SHELL_H

/*
  MITC (Mixed interpolation of tensorial components)-based shell
  element. 

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSElement.h"
#include "FSDTStiffness.h"
#include "FElibrary.h"
#include "TensorToolbox.h"
#include "TACSShell.h"
#include "ShellUtils.h"
#include "LargeRotUtils.h"

// Include all the functions from the namespace shellutils
using namespace shellutils;
using namespace largerot;

/*
  A Shell element for general linear and geometrically nonlinear
  analysis.

  This element employs a mixed interpolation of tensorial (strain)
  components (MITC) method to avoid shear locking problems. The flag
  MITCShellType signals whether to use the full nonlinear terms in the
  strain expressions. However, the element unknowns are in-plane
  displacements and small-angle rotations. This limits the element to
  relative small rotations. Stability and post-buckling can be explored
  as long as the rotations remain moderate. 

  In the MITC approach, the strain compoennts susceptible to locking
  are interpolated with appropriate low-order polynomails. In plance
  of interpolating the shear and in-plane components directly, the
  tensorial strains are interpolated and transformed back to the local
  coordinate strains at the interpolation points (also called tying
  points). This means that the evaluation of the residuals, stiffness
  matrix and design-type derivatives are more computationally expensive
  than an equivalent displacement-based shell. However, avoiding shear 
  locking is a key advantage and these elements should be used whenever
  possible!
*/
template <int order>
class MITCShell : public TACSShell {
 public:
  MITCShell( FSDTStiffness * _stiff,
	     TACSShellType _type = LINEAR, 
	     int _componentNum = 0 );
  ~MITCShell();
  
  // How many nodes/variables per element
  // ------------------------------------
  int numNodes();
  int numVariables();
  
  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies( TacsScalar *_Te, TacsScalar *_Pe,
			const TacsScalar Xpts[],
			const TacsScalar vars[],
			const TacsScalar dvars[] );

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void getResidual( TacsScalar res[], const TacsScalar Xpts[],
		    const TacsScalar vars[],
		    const TacsScalar dvars[],
		    const TacsScalar ddvars[] );

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void getJacobian( TacsScalar J[],
		    double alpha, double beta, double gamma,
		    const TacsScalar Xpts[],
		    const TacsScalar vars[],
		    const TacsScalar dvars[],
		    const TacsScalar ddvars[] );

  // Add the product of the adjoint with the derivative of the design variables
  // --------------------------------------------------------------------------
  void addAdjResProduct( double scale,
			 TacsScalar dvSens[], int dvLen,
			 const TacsScalar psi[],
			 const TacsScalar Xpts[],
			 const TacsScalar vars[],
			 const TacsScalar dvars[],
			 const TacsScalar ddvars[] );

  // Add the product of the adjoint with the derivative of the design variables
  // --------------------------------------------------------------------------  
  void getAdjResXptProduct( TacsScalar XptSens[],
			    const TacsScalar psi[],
			    const TacsScalar Xpts[],
			    const TacsScalar vars[],
			    const TacsScalar dvars[],
			    const TacsScalar ddvars[] );

  // Retrieve a specific time-independent matrix from the element
  // ------------------------------------------------------------
  void getMatType( ElementMatrixType matType, 
		   TacsScalar mat[], 
		   const TacsScalar Xpts[],
		   const TacsScalar vars[] );

  // Compute the derivative of the inner product w.r.t. design variables
  // -------------------------------------------------------------------
  void addMatDVSensInnerProduct( ElementMatrixType matType, 
				 double scale,
				 TacsScalar dvSens[], int dvLen,
				 const TacsScalar psi[], 
				 const TacsScalar phi[],
				 const TacsScalar Xpts[],
				 const TacsScalar vars[] );

  // Compute the derivative of the inner product w.r.t. vars[]
  // ---------------------------------------------------------
  void getMatSVSensInnerProduct( ElementMatrixType matType, 
				 TacsScalar res[],
				 const TacsScalar psi[], 
				 const TacsScalar phi[],
				 const TacsScalar Xpts[],
				 const TacsScalar vars[] );

  // Member functions for evaluating global functions of interest
  // ------------------------------------------------------------
  TACSConstitutive * getConstitutive();

  // Get the number of Gauss quadrature points
  // -----------------------------------------
  int getNumGaussPts();

  // Get the quadrature points and weights
  // -------------------------------------
  double getGaussWtsPts( const int num, double * pt );

  // Get the shape functions from the element
  // ----------------------------------------
  void getShapeFunctions( const double pt[], double N[] );
  
  // Return the determinant of the Jacobian at this point
  // ----------------------------------------------------
  TacsScalar getDetJacobian( const double * pt, 
                             const TacsScalar Xpts[] );

  // Return the determinant of the Jacobian and its sensitivity at this point
  // ------------------------------------------------------------------------
  TacsScalar getDetJacobianXptSens( TacsScalar * hXptSens, 
                                    const double * pt, 
                                    const TacsScalar Xpts[] );

  // This function returns the strain evaluated at pt
  // ------------------------------------------------
  void getStrain( TacsScalar strain[], 
		  const double pt[], 
		  const TacsScalar Xpts[],
		  const TacsScalar vars[] );

  // This function returns the sensitivity of the strain w.r.t. Xpts
  // ---------------------------------------------------------------
  void addStrainXptSens( TacsScalar strainXptSens[],
			 const double pt[], 
			 const TacsScalar scale,
			 const TacsScalar strainSens[], 
			 const TacsScalar Xpts[],
			 const TacsScalar vars[] );
 
  // This function adds the sensitivity of the strain to the state variables
  // -----------------------------------------------------------------------
  void addStrainSVSens( TacsScalar strainSVSens[], 
			const double pt[], 
			const TacsScalar scale,
			const TacsScalar strainSens[], 
			const TacsScalar Xpts[],
			const TacsScalar vars[] );

  // Functions for post-processing
  // -----------------------------
  void addOutputCount( int * nelems, int * nnodes, int * ncsr );
  void getOutputData( unsigned int out_type, 
		      double * data, int ld_data,
		      const TacsScalar Xpts[],
		      const TacsScalar vars[] );
  void getOutputConnectivity( int * con, int node );

  // Set the number of nodes/number of variables per element
  // -------------------------------------------------------
  static const int NUM_NODES = order*order;
  static const int NUM_VARIABLES = 6*order*order;

 private:
  static const int NUM_G11 = (order-1)*order;
  static const int NUM_G22 = (order-1)*order;
  static const int NUM_G12 = (order-1)*(order-1);
  static const int NUM_G13 = (order-1)*order;
  static const int NUM_G23 = (order-1)*order;

  inline static TacsScalar strain_product( const TacsScalar a[], 
					   const TacsScalar b[] ){
    return (a[0]*b[0] + a[1]*b[1] + 
            a[2]*b[2] + a[3]*b[3] + 
            a[4]*b[4] + a[5]*b[5] + 
            a[6]*b[6] + a[7]*b[7]);
  }

  // The type of strain expressions to use
  TACSShell::TACSShellType type;

  // The quadrature scheme for integrating the residual/stiffness
  // matrix -- this is Gauss quadrature
  int numGauss;
  const double *gaussWts, *gaussPts;

  // The knot locations
  const double * knots; // "order" Gauss points
  const double * pknots; // "order"-1 Gauss points
};

const double MITCShellFirstOrderKnots[2] = {-1.0, 1.0};

template <int order>
MITCShell<order>::MITCShell( FSDTStiffness * _stiff, 
			     TACSShell::TACSShellType _type, 
			     int _componentNum ):
TACSShell(_stiff, _componentNum){
  type = _type;
  
  unsigned int gaussOrder = order;
  numGauss = FElibrary::getGaussPtsWts(gaussOrder, &gaussPts, &gaussWts);

  // Get the knot points - the order and order-1-th Gauss points
  if (order == 2){
    knots = MITCShellFirstOrderKnots;
  }
  else { 
    FElibrary::getGaussPtsWts(order, &knots, NULL);
  }
  FElibrary::getGaussPtsWts(order-1, &pknots, NULL);
}

template <int order>
MITCShell<order>::~MITCShell(){}

template <int order>
int MITCShell<order>::numNodes(){ 
  return NUM_NODES; 
}

template <int order>
int MITCShell<order>::numVariables(){
  return NUM_VARIABLES; 
}

/*
  Get the element residuals corresponding to the strain energy
  contributions (e.g. not the work terms)

  output:
  res:     the element residual
  
  input:
  vars:    the element variables
  Xpts:    the element nodal locations in R^{3}
*/
template <int order>
void MITCShell<order>::getRes( TacsScalar * res, 
			       const TacsScalar vars[],
			       const TacsScalar Xpts[] ){
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9]; 

  TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the displacement strain
  TacsScalar b11[3*NUM_NODES*NUM_G11], b22[3*NUM_NODES*NUM_G22];
  TacsScalar b12[3*NUM_NODES*NUM_G12];
  TacsScalar b13[NUM_VARIABLES*NUM_G13], b23[NUM_VARIABLES*NUM_G23];

  double gpt[2] = {0.0, 0.0};

  TacsScalar strain[NUM_STRESSES], stress[NUM_STRESSES];
  TacsScalar drot[NUM_VARIABLES];
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  // Store the stiffness information for the element
  TacsScalar At[6], Bt[6], Dt[6], Ats[3];
  TacsScalar k_penalty = 0.0;

  // Zero the residual
  memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));

  // Evaluate the stiffness matrix
  if (!stiff->isVariableStiffness()){
    k_penalty = stiff->getStiffness(gpt, At, Bt, Dt, Ats);
  }
  
  if (type == LARGE_ROTATION){
    compute_lr_tying_bmat<order>(g11, g22, g12, g23, g13, 
				 b11, b22, b12, b23, b13, 
				 knots, pknots, vars, Xpts);
  }
  else {
    compute_tying_bmat<order>((type == LINEAR),
			      g11, g22, g12, g23, g13,
			      b11, b22, b12, b23, b13,
			      knots, pknots, vars, Xpts);
  }

  for ( int m = 0; m < numGauss; m++ ){
    for ( int n = 0; n < numGauss; n++ ){
      gpt[0] = gaussPts[n];
      gpt[1] = gaussPts[m];

      // Evaluate the stiffness matrix if required
      if (stiff->isVariableStiffness()){
        k_penalty = stiff->getStiffness(gpt, At, Bt, Dt, Ats);
      }
	
      // Calculate the shape functions
      shell_hessian(order, X, Xd, Xdd, 
                    N, Na, Nb, Naa, Nab, Nbb,
                    gpt, Xpts);
      compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

      TacsScalar h = 0.0;
      if (stiff->getTransformType() == FSDTStiffness::NATURAL){
	h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, 
                              Xd, Xdd);
      }
      else {
	const TacsScalar * axis = stiff->getRefAxis();
	h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, 
                                      normal_eta, axis, Xd, Xdd);
      }
      h = gaussWts[n]*gaussWts[m]*h;

      // Calculate the deformation at the current point... 
      TacsScalar rot = 0.0; // Difference between inplane/drilling rotation
      if (type == LINEAR){
	linear_bend_strain(strain, &rot, U, Ud, 
                           t, tx, ztx, 
                           normal, normal_xi, normal_eta);
	linear_bend_bmat(B, drot, NUM_NODES,
                         N, Na, Nb, t, tx, ztx,
                         normal, normal_xi, normal_eta);
      }
      else if (type == NONLINEAR){
	nonlinear_bend_strain(strain, &rot, U, Ud, 
                              t, tx, ztx, 
                              normal, normal_xi, normal_eta);
	nonlinear_bend_bmat(B, drot, NUM_NODES,
                            N, Na, Nb, U, Ud, t, tx, ztx,
                            normal, normal_xi, normal_eta);
      }
      else {
	// Rotation matrix data
	TacsScalar C[9], Ct[27], Ctt[54];

	// Compute the rotation matrices
	TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
	TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
	TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
	compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
	compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);

	// Evaluate the in-plane rotation term
	rot = compute_inplane_penalty(drot, NUM_NODES, Xd, Ud, 
				      C, Ct, N, Na, Nb);
	
	// Calculate the deformation at the current point... 
	large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, 
			      normal, normal_xi, normal_eta);
	large_rot_bend_bmat(B, NUM_NODES, N, Na, Nb, U, Ud, C, Ct, Ctt, 
			    t, tx, ztx, normal, normal_xi, normal_eta);
      }

      // Evaluate the strain interpolation at this point
      tying_interpolation<order>(gpt, N11, N22, N12, 
                                 knots, pknots);
      add_tying_strain<order>(strain, tx, 
                              g11, g22, g12, g23, g13,
                              N11, N22, N12);
      add_tying_bmat<order>(B, NUM_NODES, tx,
                            b11, b22, b12, b23, b13, 
                            N11, N22, N12);
      
      // Compute the stress at the current Gauss point
      stiff->calcStress(At, Bt, Dt, Ats, strain, stress);

      for ( int i = 0; i < NUM_NODES; i++ ){
	for ( int ii = 0; ii < NUM_DISPS; ii++ ){	    
	  int row = ii + NUM_DISPS*i;

	  res[row] += 
	    h*(strain_product(&B[NUM_STRESSES*row], stress) +
	       k_penalty*rot*drot[row]);
	}
      }
    }
  }
}

/*
  Get the element tangent stiffness matrix - the exact Jacobian of the
  residual expressions.

  output:
  mat:     the element tangent stiffness matrix
  res:     the element residual
  
  input:
  vars:    the element variables
  Xpts:    the element nodal locations in R^{3}
  matOr:   the matrix orientation (NORMAL or TRANSPOSE)
*/
template <int order>
void MITCShell<order>::getMat( TacsScalar * mat, TacsScalar * res, 
			       const TacsScalar vars[],
			       const TacsScalar Xpts[], 
			       MatrixOrientation matOr ){
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9]; 
  
  TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the displacement strain
  TacsScalar b11[3*NUM_NODES*NUM_G11], b22[3*NUM_NODES*NUM_G22];
  TacsScalar b12[3*NUM_NODES*NUM_G12];
  TacsScalar b13[NUM_VARIABLES*NUM_G13], b23[NUM_VARIABLES*NUM_G23];
  
  // The second derivatives of the tying strain
  TacsScalar n13[12*order*order*(order*order+1)*NUM_G13];
  TacsScalar n23[12*order*order*(order*order+1)*NUM_G23];

  // The stress and strain information
  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];
  TacsScalar BStress[NUM_STRESSES];
  TacsScalar drot[NUM_VARIABLES];

  // Store the stiffness information for the element
  TacsScalar At[6], Bt[6], Dt[6], Ats[3];
  TacsScalar k_penalty = 0.0;

  // Zero the element matrix and the residual
  memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));
  memset(mat, 0, NUM_VARIABLES*NUM_VARIABLES*sizeof(TacsScalar));

  // Evaluate the stiffness matrix
  double gpt[2] = {0.0, 0.0};
  if (!stiff->isVariableStiffness()){
    k_penalty = stiff->getStiffness(gpt, At, Bt, Dt, Ats);
  }
  
  if (type == LARGE_ROTATION){
    compute_lr_tying_bmat<order>(g11, g22, g12, g23, g13, 
				 b11, b22, b12, b23, b13, 
				 knots, pknots, vars, Xpts);
  }
  else {
    compute_tying_bmat<order>((type == LINEAR),
			      g11, g22, g12, g23, g13,
			      b11, b22, b12, b23, b13,
			      knots, pknots, vars, Xpts);
  }

  for ( int m = 0; m < numGauss; m++ ){
    for ( int n = 0; n < numGauss; n++ ){
      gpt[0] = gaussPts[n];
      gpt[1] = gaussPts[m];
      
      // Evaluate the stiffness matrix if required
      if (stiff->isVariableStiffness()){
        k_penalty = stiff->getStiffness(gpt, At, Bt, Dt, Ats);
      }

      // Calculate the shape functions	
      shell_hessian(order, X, Xd, Xdd, 
                    N, Na, Nb, Naa, Nab, Nbb,
                    gpt, Xpts);
      compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);
      
      TacsScalar h = 0.0;
      if (stiff->getTransformType() == FSDTStiffness::NATURAL){
	h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, 
                              Xd, Xdd);
      }
      else {
	const TacsScalar * axis = stiff->getRefAxis();
	h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, 
                                      normal_eta, axis, Xd, Xdd);
      }
      h = gaussWts[n]*gaussWts[m]*h;


      // Store the difference between the rotation variable
      // and the in-plane rotation
      TacsScalar rot = 0.0;

      // Rotation matrix data
      TacsScalar C[9], Ct[27], Ctt[54], Cttt[63];
      
      // Compute the strain at the current point
      if (type == LINEAR){
	linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, 
                           normal, normal_xi, normal_eta);
	linear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, t, tx, ztx,
                         normal, normal_xi, normal_eta);
      }
      else if (type == NONLINEAR){
	nonlinear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, 
                              normal, normal_xi, normal_eta);
	nonlinear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, U, Ud, t, tx, ztx,
                            normal, normal_xi, normal_eta);
      }
      else {
	// Compute the rotation matrices
	TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
	TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
	TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
	compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
	compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);
	compute_3rd_rate_matrix(Cttt, c1, s1, c2, s2, c3, s3);
	
	// Evaluate the in-plane rotation term
	rot = compute_inplane_penalty(drot, NUM_NODES, Xd, Ud, 
				      C, Ct, N, Na, Nb);

	// Calculate the strain/bmat at the current point 
	large_rot_bend_strain(strain, U, Ud, C, Ct,
			      t, tx, ztx, normal, normal_xi, normal_eta);
	large_rot_bend_bmat(B, NUM_NODES, N, Na, Nb, U, Ud, C, Ct, Ctt, 
			    t, tx, ztx, normal, normal_xi, normal_eta);
      }
      
      // Evaluate the strain interpolation at this point
      tying_interpolation<order>(gpt, N11, N22, N12, 
                                 knots, pknots);

      // Add the interpolated strain and the interpolated b-matrix to the
      // point-wise strain and strain-derivative (B)
      add_tying_strain<order>(strain, tx, g11, g22, g12, g23, g13,
                              N11, N22, N12);
      add_tying_bmat<order>(B, NUM_NODES, tx, b11, b22, b12, b23, b13, 
                            N11, N22, N12);

      stiff->calcStress(At, Bt, Dt, Ats, strain, stress);

      // Loop through all displacements except for the z-rotation
      for ( int i = 0; i < NUM_NODES; i++ ){
	for ( int ii = 0; ii < NUM_DISPS; ii++ ){	  
	  int row = ii + NUM_DISPS*i;

	  // Determine the contribution to the residual
	  res[row] += 
	    h*(strain_product(&B[row*NUM_STRESSES], stress) + 
               k_penalty*rot*drot[row]);	  
	}
      }

      for ( int i = 0; i < NUM_NODES; i++ ){
	for ( int ii = 0; ii < NUM_DISPS; ii++ ){
	  int row = ii + NUM_DISPS*i;
	  
	  // Calculate the stress associated with B
          stiff->calcStress(At, Bt, Dt, Ats, 
                                 &B[row*NUM_STRESSES], BStress);
	  
	  for ( int j = 0; j <= i; j++ ){
	    int end = ( j == i ? ii : NUM_DISPS-1);
	    for ( int jj = 0; jj <= end; jj++ ){
	      int col = jj + NUM_DISPS*j;
	      
	      // The regular element matrix
	      mat[col + row*NUM_VARIABLES] += 
		h*(strain_product(BStress, &B[col*NUM_STRESSES]) +    
		  k_penalty*drot[row]*drot[col]);	      
	    }
	  }
	}
      }

      if (type == NONLINEAR){
        nonlinear_bend_stress_bmat(mat, NUM_NODES, h, stress, 
                                   N, Na, Nb, t, tx, ztx,
                                   normal, normal_xi, normal_eta);
        add_nonlinear_tying_stress_nmat<order>(mat, h, stress,
                                               tx, N11, N22, N12,
                                               knots, pknots, Xpts);
      }
      else if (type == LARGE_ROTATION){
	// Add the second-derivative contributions from the bending strain
	add_large_rot_bend_stress_bmat(mat, NUM_NODES, h, stress, 
				       N, Na, Nb, U, Ud, C, Ct, Ctt, Cttt, 
				       t, tx, ztx, normal, normal_xi, normal_eta);
	
	// Add the contributions to the second derivative of the tying strain
	add_lr_tying_stress_nmat<order>(mat, h, stress, n13, n23, tx, 
					N11, N22, N12, knots, pknots); 
	
	// Add the second derivative of the in-plane penalty
	add_inplane_penalty(mat, NUM_NODES, h*k_penalty*rot, Xd, Ud, 
			    Ct, Ctt, N, Na, Nb);
      }
    }
  }

  // Copy over the matrix
  // Take the lower triangle and copy to the upper triangle
  for ( int row = 0; row < NUM_VARIABLES; row++ ){
    for ( int col = row+1; col < NUM_VARIABLES; col++ ){
      mat[col + row*NUM_VARIABLES] = mat[row + col*NUM_VARIABLES];
    }
  }
}

/*
  Evaluate the element matrix of a specified type

  output:
  mat:         the element matrix of the specified type 

  input:
  matType:     the matrix type (e.g. MASS_MATRIX)
  scaleFactor: scale factor such that mat = scaleFactor*M
  vars:        the element variables
  Xpts:        the nodal coordinates in R^{3}
  matOr:       the matrix orientation either NORMAL or TRANSPOSE
*/
template <int order>
void MITCShell<order>::getMatType( ElementMatrixTypes matType, 
				   TacsScalar scaleFactor, TacsScalar * mat, 
				   const TacsScalar vars[], 
				   const TacsScalar Xpts[], 
				   MatrixOrientation matOr ){  
  memset(mat, 0, NUM_VARIABLES*NUM_VARIABLES*sizeof(TacsScalar));
  
  if (matType == MASS_MATRIX){
    double gpt[2] = {0.0, 0.0};
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];

    for ( int m = 0; m < numGauss; m++ ){
      for ( int n = 0; n < numGauss; n++ ){
	gpt[0] = gaussPts[n];
	gpt[1] = gaussPts[m];

	TacsScalar X[3], Xd[9], normal[3];
	shell_jacobian(order, X, Xd, N, Na, Nb, gpt, Xpts);
	Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
	Tensor::normalize3D(normal);

	for ( int i = 0; i < 3; i++ ){
	  Xd[6+i] = normal[i]; 
	}

	TacsScalar h = FElibrary::jacobian3d(Xd);
	h = scaleFactor*gaussWts[n]*gaussWts[m]*h;

        TacsScalar mass[3];
	stiff->pointwiseMass(gpt, mass);

	for ( int i = 0; i < NUM_NODES; i++ ){
	  for ( int ii = 0; ii < 3; ii++ ){
	    int row = ii + NUM_DISPS*i;
	    for ( int j = 0; j < NUM_NODES; j++ ){
	      int col = ii + NUM_DISPS*j;	      
	      mat[col + row*NUM_VARIABLES] += mass[0]*h*N[i]*N[j];
	    }
	  }
	}

	for ( int ii = 3; ii < NUM_DISPS; ii++ ){	 
	  for ( int jj = 3; jj < NUM_DISPS; jj++ ){
	    TacsScalar a[3] = { 0.0, 0.0, 0.0 };
	    TacsScalar b[3] = { 0.0, 0.0, 0.0 };

	    TacsScalar phi_a[3], phi_b[3];
	    a[ii-3] = 1.0;
	    b[jj-3] = 1.0;

	    Tensor::crossProduct3D(phi_a, a, normal);
	    Tensor::crossProduct3D(phi_b, b, normal);
	    TacsScalar weight = Tensor::dot3D(phi_a, phi_b);

	    for ( int i = 0; i < NUM_NODES; i++ ){
	      int row = ii + NUM_DISPS*i;
	      for ( int j = 0; j < NUM_NODES; j++ ){
		int col = jj + NUM_DISPS*j;		
		mat[col + row*NUM_VARIABLES] += mass[2]*weight*h*N[i]*N[j];
	      }
	    }
	  }
	}
      }
    }    
  }
  else if (matType == GEOMETRIC_STIFFNESS_MATRIX &&
	   type == LINEAR){
    // Geometric data
    TacsScalar X[3], Xd[9], Xdd[9];
    TacsScalar normal[3], normal_xi[3], normal_eta[3];
    
    // Transformation and the transformation derivative w.r.t. zeta
    TacsScalar t[9], tx[9], ztx[9]; 
    
    // Displacements/rotations and shape functions
    TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];
    
    // Interpolations for the shear components
    double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];
    
    // The interpolated tensorial shear components
    TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
    TacsScalar g13[NUM_G13], g23[NUM_G23];

    // The stress/strain values
    TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];

    // Store the stiffness information for the element
    TacsScalar At[6], Bt[6], Dt[6], Ats[3];

    // Evaluate the stiffness matrix
    if (!stiff->isVariableStiffness()){
      double pt[2] = {0.0, 0.0};
      stiff->getStiffness(pt, At, Bt, Dt, Ats);
    }

    const int is_linear = 1;
    compute_tying_strain<order>(is_linear,
				g11, g22, g12, g23, g13, 
				knots, pknots, vars, Xpts);

    for ( int m = 0; m < numGauss; m++ ){
      for ( int n = 0; n < numGauss; n++ ){
	double gpt[2];
	gpt[0] = gaussPts[n];
	gpt[1] = gaussPts[m];

        // Evaluate the stiffness matrix if required
        if (stiff->isVariableStiffness()){
          stiff->getStiffness(gpt, At, Bt, Dt, Ats);
        }
		
	// Calculate the shape functions
	shell_hessian(order, X, Xd, Xdd, 
                      N, Na, Nb, Naa, Nab, Nbb,
                      gpt, Xpts);
	compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);
      
	TacsScalar h = 0.0;
	if (stiff->getTransformType() == FSDTStiffness::NATURAL){
	  h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, 
                                Xd, Xdd);
	}
	else {
	  const TacsScalar * axis = stiff->getRefAxis();
	  h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, 
                                        normal_eta, axis, Xd, Xdd);
	}
	h = gaussWts[n]*gaussWts[m]*h*scaleFactor;

	// Calculate the strain/bmat at the current point... 
	TacsScalar rot;
	linear_bend_strain(strain, &rot, U, Ud, 
                           t, tx, ztx, 
                           normal, normal_xi, normal_eta);
      
	// Evaluate the strain interpolation at this point
	tying_interpolation<order>(gpt, N11, N22, N12, 
                                   knots, pknots);
	add_tying_strain<order>(strain, tx, g11, g22, g12, g23, g13,
                                N11, N22, N12);

	// Compute the stress at this Gauss point
        stiff->calcStress(At, Bt, Dt, Ats, strain, stress);
      
	// Add the geometric stiffness terms
	nonlinear_bend_stress_bmat(mat, NUM_NODES, h, stress, 
				   N, Na, Nb, t, tx, ztx,
				   normal, normal_xi, normal_eta);
	add_nonlinear_tying_stress_nmat<order>(mat, h, stress,
					       tx, N11, N22, N12,
					       knots, pknots, Xpts);
      }
    }
  
    // Copy over the matrix
    // Take the lower triangle and copy to the upper triangle
    for ( int row = 0; row < NUM_VARIABLES; row++ ){
      for ( int col = row+1; col < NUM_VARIABLES; col++ ){
	mat[col + row*NUM_VARIABLES] = mat[row + col*NUM_VARIABLES];
      }
    }
  }
}

/*
  Evaluate the derivative of the element residuals with respect
  to the nodal coordinates e.g res = dR/dXpts

  output:
  res:  the derivative of the residuals w.r.t. the element nodes
  
  input:
  vars:    the element variables
  Xpts:    the element nodal locations
*/
template <int order>
void MITCShell<order>::getResXptSens( TacsScalar * res, 
				      const TacsScalar vars[], 
				      const TacsScalar Xpts[] ){
  memset(res, 0, 3*NUM_NODES*NUM_VARIABLES*sizeof(TacsScalar));
  
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];
  TacsScalar dnormal[9*NUM_NODES], dnormal_xi[9*NUM_NODES], 
    dnormal_eta[9*NUM_NODES];
    
  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9];
  TacsScalar dt[27*NUM_NODES], dtx[27*NUM_NODES], dztx[27*NUM_NODES];
  
  TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];
    
  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];
    
  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the tensorial tying strain components
  TacsScalar dg11[3*NUM_G11*NUM_NODES], dg22[3*NUM_G22*NUM_NODES], 
    dg12[3*NUM_G12*NUM_NODES];
  TacsScalar dg13[3*NUM_G13*NUM_NODES], dg23[3*NUM_G23*NUM_NODES];
    
  // The derivatives of the displacement strain
  TacsScalar b11[3*NUM_NODES*NUM_G11], b22[3*NUM_NODES*NUM_G22];
  TacsScalar b12[3*NUM_NODES*NUM_G12];
  TacsScalar b13[NUM_VARIABLES*NUM_G13], b23[NUM_VARIABLES*NUM_G23];
    
  double gpt[2] = {0.0, 0.0};
    
  // The stress and strain and their sensitivities
  TacsScalar strain[NUM_STRESSES], stress[NUM_STRESSES];
  TacsScalar dstrain[3*NUM_STRESSES*NUM_NODES], dstress[NUM_STRESSES];

  // The rotational contributions
  TacsScalar drot[NUM_VARIABLES], srot[NUM_VARIABLES];

  // The derivative of the strain w.r.t. the nodal displacements
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  // The sensitivity of the determinant of the Jacobian
  TacsScalar dh[3*NUM_NODES];

  // Store the stiffness information for the element
  TacsScalar At[6], Bt[6], Dt[6], Ats[3];
  TacsScalar k_penalty = 0.0;
  
  // Evaluate the stiffness matrix
  if (!stiff->isVariableStiffness()){
    k_penalty = stiff->getStiffness(gpt, At, Bt, Dt, Ats);
  }

  compute_tying_strain_sens<order>((type == LINEAR),
				   g11, g22, g12, g23, g13,
				   dg11, dg22, dg12, dg23, dg13,
				   knots, pknots, vars, Xpts);
  compute_tying_bmat<order>((type == LINEAR),
			    g11, g22, g12, g23, g13,
			    b11, b22, b12, b23, b13,
			    knots, pknots, vars, Xpts);
  
  for ( int m = 0; m < numGauss; m++ ){
    for ( int n = 0; n < numGauss; n++ ){
      gpt[0] = gaussPts[n];
      gpt[1] = gaussPts[m];

      // Evaluate the stiffness matrix if required
      if (stiff->isVariableStiffness()){
        k_penalty = stiff->getStiffness(gpt, At, Bt, Dt, Ats);
      }
		            
      // Calculate the shape functions and the surface derivatives
      shell_hessian(order, X, Xd, Xdd, 
                    N, Na, Nb, Naa, Nab, Nbb, 
                    gpt, Xpts);

      // Compute the values of U and Ud - the variables and their
      // derivatives w.r.t. the local coordinates
      compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

      // Evaluate the strain interpolation at this point
      tying_interpolation<order>(gpt, N11, N22, N12, knots, pknots);

      TacsScalar h = 0.0;  
      if (stiff->getTransformType() == FSDTStiffness::NATURAL){
        h = compute_transform_sens(dh, t, dt, tx, dtx, ztx, dztx, 
                                   normal, dnormal, 
                                   normal_xi, dnormal_xi, 
                                   normal_eta, dnormal_eta, 
                                   Xd, Xdd, Na, Nb, Naa, Nab, Nbb, NUM_NODES);
      }
      else {
        const TacsScalar * axis = stiff->getRefAxis();
        h = compute_transform_refaxis_sens(dh, t, dt, tx, dtx, ztx, dztx, 
                                           normal, dnormal, 
                                           normal_xi, dnormal_xi, 
                                           normal_eta, dnormal_eta, axis,
                                           Xd, Xdd, Na, Nb, Naa, Nab, Nbb, 
                                           NUM_NODES);
      }
      h = gaussWts[n]*gaussWts[m]*h;

      // Evaluate the strain and the derivative of the strain w.r.t. 
      // the displacements and the nodal coordinates
      TacsScalar rot;
      if (type == LINEAR){
	linear_bend_strain_sens(strain, dstrain, &rot, srot, U, Ud, 
                                t, dt, tx, dtx, ztx, dztx,
				normal, dnormal, normal_xi, 
				dnormal_xi, normal_eta, dnormal_eta, 
                                3*NUM_NODES);
	linear_bend_bmat(B, drot, NUM_NODES,
                         N, Na, Nb, t, tx, ztx,
                         normal, normal_xi, normal_eta);
      }
      else if (type == NONLINEAR){
	nonlinear_bend_strain_sens(strain, dstrain, &rot, srot, U, Ud, 
				   t, dt, tx, dtx, ztx, dztx,
				   normal, dnormal, normal_xi, 
				   dnormal_xi, normal_eta, dnormal_eta,
                                   3*NUM_NODES);
	nonlinear_bend_bmat(B, drot, NUM_NODES,
                            N, Na, Nb, U, Ud, t, tx, ztx,
                            normal, normal_xi, normal_eta);
      }
                        
      // Add the sensitivities of the tying strain points into the
      // strain senstivities
      add_tying_strain_sens<order>(strain, dstrain, tx, dtx,
                                   g11, g22, g12, g23, g13,
                                   dg11, dg22, dg12, dg23, dg13,
                                   N11, N22, N12);

      // Add the tying strain to the B-matrix
      add_tying_bmat<order>(B, NUM_NODES, tx, 
                            b11, b22, b12, b23, b13, 
                            N11, N22, N12);
            
      // Calculate the stress at the current point
      stiff->calcStress(At, Bt, Dt, Ats, strain, stress);

      for ( int k = 0; k < 3*NUM_NODES; k++ ){
        dh[k] = gaussWts[n]*gaussWts[m]*dh[k];
        stiff->calcStress(At, Bt, Dt, Ats, 
                               &dstrain[k*NUM_STRESSES], dstress);
            
        for ( int i = 0; i < NUM_NODES; i++ ){
          for ( int ii = 0; ii < NUM_DISPS; ii++ ){	    
            int row = NUM_DISPS*i + ii;

            res[NUM_VARIABLES*k + row] += 
              h*(strain_product(&B[NUM_STRESSES*row], dstress) +
                 k_penalty*(srot[k]*drot[row])) + 
              dh[k]*(strain_product(&B[NUM_STRESSES*row], stress) +
                     k_penalty*rot*drot[row]);
          }
        }
      }

      if (type == LINEAR){
        add_linear_bend_bmat_sens(res, NUM_NODES, h, h*k_penalty*rot,
                                  stress, N, Na, Nb, 
                                  t, dt, tx, dtx, ztx, dztx,
                                  normal, dnormal, normal_xi, dnormal_xi,
                                  normal_eta, dnormal_eta, 3*NUM_NODES);
      }
      else if (type == NONLINEAR){
        add_nonlinear_bend_bmat_sens(res, NUM_NODES, h, h*k_penalty*rot,
                                     stress, N, Na, Nb, U, Ud,
                                     t, dt, tx, dtx, ztx, dztx,
                                     normal, dnormal, normal_xi, dnormal_xi,
                                     normal_eta, dnormal_eta, 3*NUM_NODES);
      }

      add_tying_bmat_sens<order>((type == LINEAR), 
				 res, h, stress, tx, dtx,
				 knots, pknots, vars, Xpts, N11, N22, N12);
    }
  }
}

/*
  Compute the derivative of the product of the residual with the
  adjoint vector psi with respect to the design variables.  The result
  is added to a design-variable array.

  input:
  alpha:    the scaling factor applied to the derivative
  dvLen:    the length of the design variable vector
  psi:      the left-multiplying vector
  phi:      the right-multiplying vector
  vars:     the element state variables
  Xpts:     the nodal locations for this element

  output:
  dvSens:   the result is added to this vector
*/
template <int order>
void MITCShell<order>::addAdjResDVSensProduct( TacsScalar alpha,
					       TacsScalar dvSens[], int dvLen,
					       const TacsScalar psi[],
					       const TacsScalar vars[],
					       const TacsScalar Xpts[] ){
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9]; 

  TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the displacement strain
  TacsScalar b11[3*NUM_NODES*NUM_G11], b22[3*NUM_NODES*NUM_G22];
  TacsScalar b12[3*NUM_NODES*NUM_G12];
  TacsScalar b13[NUM_VARIABLES*NUM_G13], b23[NUM_VARIABLES*NUM_G23];

  TacsScalar strain[NUM_STRESSES], stress[NUM_STRESSES];
  TacsScalar drot[NUM_VARIABLES];
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  // Compute the tying terms in the matrix
  if (type == LARGE_ROTATION){
    compute_lr_tying_bmat<order>(g11, g22, g12, g23, g13, 
				 b11, b22, b12, b23, b13, 
				 knots, pknots, vars, Xpts);
  }
  else {
    compute_tying_bmat<order>((type == LINEAR),
			      g11, g22, g12, g23, g13,
			      b11, b22, b12, b23, b13,
			      knots, pknots, vars, Xpts);
  }

  for ( int m = 0; m < numGauss; m++ ){
    for ( int n = 0; n < numGauss; n++ ){
      double gpt[2];
      gpt[0] = gaussPts[n];
      gpt[1] = gaussPts[m];
	
      // Calculate the shape functions
      shell_hessian(order, X, Xd, Xdd, 
                    N, Na, Nb, Naa, Nab, Nbb,
                    gpt, Xpts);
      compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

      TacsScalar h = 0.0;
      if (stiff->getTransformType() == FSDTStiffness::NATURAL){
	h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, 
                              Xd, Xdd);
      }
      else {
	const TacsScalar * axis = stiff->getRefAxis();
	h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, 
                                      normal_eta, axis, Xd, Xdd);
      }
      h = gaussWts[n]*gaussWts[m]*h;

      // Calculate the deformation at the current point... 
      TacsScalar rot = 0.0; // Difference between inplane/drilling rotation
      if (type == LINEAR){
	linear_bend_strain(strain, &rot, U, Ud, 
                           t, tx, ztx, 
                           normal, normal_xi, normal_eta);
	linear_bend_bmat(B, drot, NUM_NODES,
                         N, Na, Nb, t, tx, ztx,
                         normal, normal_xi, normal_eta);
      }
      else if (type == NONLINEAR){
	nonlinear_bend_strain(strain, &rot, U, Ud, 
                              t, tx, ztx, 
                              normal, normal_xi, normal_eta);
	nonlinear_bend_bmat(B, drot, NUM_NODES,
                            N, Na, Nb, U, Ud, t, tx, ztx,
                            normal, normal_xi, normal_eta);
      }
      else {
	// Rotation matrix data
	TacsScalar C[9], Ct[27], Ctt[54];

	// Compute the rotation matrices
	TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
	TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
	TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
	compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
	compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);

	// Evaluate the in-plane rotation term
	rot = compute_inplane_penalty(drot, NUM_NODES, Xd, Ud, 
				      C, Ct, N, Na, Nb);
	
	// Calculate the deformation at the current point... 
	large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, 
			      normal, normal_xi, normal_eta);
	large_rot_bend_bmat(B, NUM_NODES, N, Na, Nb, U, Ud, C, Ct, Ctt, 
			    t, tx, ztx, normal, normal_xi, normal_eta);
      }

      // Evaluate the strain interpolation at this point
      tying_interpolation<order>(gpt, N11, N22, N12, 
                                 knots, pknots);
      add_tying_strain<order>(strain, tx, 
                              g11, g22, g12, g23, g13,
                              N11, N22, N12);
      add_tying_bmat<order>(B, NUM_NODES, tx,
                            b11, b22, b12, b23, b13, 
                            N11, N22, N12);

      // Compute the product of psi^{T}*B^{T}
      TacsScalar bpsi[NUM_STRESSES], brot = 0.0;
      memset(bpsi, 0, NUM_STRESSES*sizeof(TacsScalar));

      TacsScalar *b = B, *dr = drot;
      const TacsScalar *ps = psi;
      for ( int i = 0; i < NUM_VARIABLES; i++ ){
	bpsi[0] += ps[0]*b[0];
	bpsi[1] += ps[0]*b[1];
	bpsi[2] += ps[0]*b[2];
	bpsi[3] += ps[0]*b[3];
	bpsi[4] += ps[0]*b[4];
	bpsi[5] += ps[0]*b[5];
	bpsi[6] += ps[0]*b[6];	
	bpsi[7] += ps[0]*b[7];

	brot += ps[0]*dr[0];
	b += NUM_STRESSES;
	dr++;
	ps++;
      }

      // Add the term: alpha*psi^{T}*B^{T}*dC/dx*strain to the vector
      // dvSens - Note that this is much more efficient than computing
      // the terms component by component
      stiff->addStiffnessDVSens(gpt, alpha*h, bpsi, strain, brot, rot,
				dvSens, dvLen);
    }
  }
}

/*
  Compute the derivative of the inner product of the stiffness or mass
  matrices with the given vectors psi and phi with respect to the
  design variables. The result is a vector which is the length of the
  number of design variables.

  input:
  matType:  the type of matrix
  alpha:    the scaling factor applied to the derivative
  dvLen:    the length of the design variable vector
  psi:      the left-multiplying vector
  phi:      the right-multiplying vector
  vars:     the element state variables
  Xpts:     the nodal locations for this element

  output:
  dvSens:   the result is added to this vector
*/
template <int order>
void MITCShell<order>::addMatDVSensInnerProduct( ElementMatrixTypes matType,
						 TacsScalar alpha,
						 TacsScalar dvSens[], int dvLen,
						 const TacsScalar psi[], 
						 const TacsScalar phi[],
						 const TacsScalar vars[],
						 const TacsScalar Xpts[] ){
  if (matType == STIFFNESS_MATRIX){
    // Geometric data
    TacsScalar X[3], Xd[9], Xdd[9];
    TacsScalar normal[3], normal_xi[3], normal_eta[3];
    
    // Transformation and the transformation derivative w.r.t. zeta
    TacsScalar t[9], tx[9], ztx[9]; 
    
    TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];
    
    // Interpolations for the shear components
    double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];
    
    // The interpolated tensorial shear components
    TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
    TacsScalar g13[NUM_G13], g23[NUM_G23];
    
    // The derivatives of the displacement strain
    TacsScalar b11[3*NUM_NODES*NUM_G11], b22[3*NUM_NODES*NUM_G22];
    TacsScalar b12[3*NUM_NODES*NUM_G12];
    TacsScalar b13[NUM_VARIABLES*NUM_G13], b23[NUM_VARIABLES*NUM_G23];
    
    TacsScalar strain[NUM_STRESSES], stress[NUM_STRESSES];
    TacsScalar drot[NUM_VARIABLES];
    TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

    // Compute the tying terms in the matrix
    if (type == LARGE_ROTATION){
      compute_lr_tying_bmat<order>(g11, g22, g12, g23, g13, 
				   b11, b22, b12, b23, b13, 
				   knots, pknots, vars, Xpts);
    }
    else {
      compute_tying_bmat<order>((type == LINEAR),
				g11, g22, g12, g23, g13,
				b11, b22, b12, b23, b13,
				knots, pknots, vars, Xpts);
    }
    
    for ( int m = 0; m < numGauss; m++ ){
      for ( int n = 0; n < numGauss; n++ ){
	double gpt[2];
	gpt[0] = gaussPts[n];
	gpt[1] = gaussPts[m];
	
	// Calculate the shape functions
	shell_hessian(order, X, Xd, Xdd, 
		      N, Na, Nb, Naa, Nab, Nbb,
		      gpt, Xpts);
	compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);
	
	TacsScalar h = 0.0;
	if (stiff->getTransformType() == FSDTStiffness::NATURAL){
	  h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, 
				Xd, Xdd);
	}
	else {
	  const TacsScalar * axis = stiff->getRefAxis();
	  h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, 
					normal_eta, axis, Xd, Xdd);
	}
	h = gaussWts[n]*gaussWts[m]*h;
	
	// Calculate the deformation at the current point... 
	TacsScalar rot = 0.0; // Difference between inplane/drilling rotation
	if (type == LINEAR){
	  linear_bend_strain(strain, &rot, U, Ud, 
			     t, tx, ztx, 
			     normal, normal_xi, normal_eta);
	  linear_bend_bmat(B, drot, NUM_NODES,
			   N, Na, Nb, t, tx, ztx,
			   normal, normal_xi, normal_eta);
	}

	// Evaluate the strain interpolation at this point
	tying_interpolation<order>(gpt, N11, N22, N12, 
				   knots, pknots);
	add_tying_strain<order>(strain, tx, 
				g11, g22, g12, g23, g13,
				N11, N22, N12);
	add_tying_bmat<order>(B, NUM_NODES, tx,
			      b11, b22, b12, b23, b13, 
			      N11, N22, N12);

	// Compute the product of psi^{T}*B^{T}
	TacsScalar bpsi[NUM_STRESSES], bphi[NUM_STRESSES];
	memset(bpsi, 0, NUM_STRESSES*sizeof(TacsScalar));
	memset(bphi, 0, NUM_STRESSES*sizeof(TacsScalar));

	TacsScalar brpsi = 0.0, brphi = 0.0;
	TacsScalar *b = B, *dr = drot;
	const TacsScalar *ps = psi, *ph = phi;
	for ( int i = 0; i < NUM_VARIABLES; i++ ){
	  bpsi[0] += ps[0]*b[0];  bpsi[1] += ps[0]*b[1];
	  bpsi[2] += ps[0]*b[2];  bpsi[3] += ps[0]*b[3];
	  bpsi[4] += ps[0]*b[4];  bpsi[5] += ps[0]*b[5];
	  bpsi[6] += ps[0]*b[6];  bpsi[7] += ps[0]*b[7];

	  bphi[0] += ph[0]*b[0];  bphi[1] += ph[0]*b[1];
	  bphi[2] += ph[0]*b[2];  bphi[3] += ph[0]*b[3];
	  bphi[4] += ph[0]*b[4];  bphi[5] += ph[0]*b[5];
	  bphi[6] += ph[0]*b[6];  bphi[7] += ph[0]*b[7];
	  
	  brpsi += ps[0]*dr[0];
	  brphi += ph[0]*dr[0];
	  
	  // Increment the pointers
	  b += NUM_STRESSES;
	  dr++; ps++; ph++;
	}

	// Add the term: alpha*psi^{T}*B^{T}*dC/dx*strain to the vector
	// dvSens - Note that this is much more efficient than computing
	// the terms component by component
	stiff->addStiffnessDVSens(gpt, alpha*h, bpsi, bphi, brpsi, brphi,
				  dvSens, dvLen);
      }
    }
  }
  else if (matType == MASS_MATRIX){
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    
    for ( int m = 0; m < numGauss; m++ ){
      for ( int n = 0; n < numGauss; n++ ){
	double gpt[2];
	gpt[0] = gaussPts[n];
	gpt[1] = gaussPts[m];
	
	TacsScalar X[3], Xd[9], normal[3];
	shell_jacobian(order, X, Xd, N, Na, Nb, gpt, Xpts);
	Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
	Tensor::normalize3D(normal);
	
	for ( int i = 0; i < 3; i++ ){
	  Xd[6+i] = normal[i]; 
	}
	
	TacsScalar h = FElibrary::jacobian3d(Xd);
	h = gaussWts[n]*gaussWts[m]*h;

	// Compute the nodal accelerations at the quadrature point
	TacsScalar upsi[6], uphi[6];
	upsi[0] = upsi[1] = upsi[2] = 
	  upsi[3] = upsi[4] = upsi[5] = 0.0;
	uphi[0] = uphi[1] = uphi[2] = 
	  uphi[3] = uphi[4] = uphi[5] = 0.0;
	
	double *ns = N;
	const TacsScalar *ps = psi, *ph = phi;
	for ( int i = 0; i < NUM_NODES; i++ ){
	  upsi[0] += ns[0]*ps[0];  upsi[1] += ns[0]*ps[1];  
	  upsi[2] += ns[0]*ps[2];  upsi[3] += ns[0]*ps[3];
	  upsi[4] += ns[0]*ps[4];  upsi[5] += ns[0]*ps[5];

	  uphi[0] += ns[0]*ph[0];  uphi[1] += ns[0]*ph[1];  
	  uphi[2] += ns[0]*ph[2];  uphi[3] += ns[0]*ph[3];
	  uphi[4] += ns[0]*ph[4];  uphi[5] += ns[0]*ph[5];
	  
	  ps += 6; ph += 6; ns++;
	}
	
	// Compute the weights on each component of the mass moments
	TacsScalar rho_alpha[3];
	rho_alpha[0] = alpha*h*(upsi[0]*uphi[0] + upsi[1]*uphi[1] + 
				upsi[2]*uphi[2]);
	rho_alpha[1] = 0.0;

	TacsScalar tphi[3], tpsi[3];
	Tensor::crossProduct3D(tphi, &uphi[3], normal);
	Tensor::crossProduct3D(tpsi, &upsi[3], normal);
	rho_alpha[2] = alpha*h*(tphi[0]*tpsi[0] + tphi[1]*tpsi[1] + 
				tphi[2]*tpsi[2]);
	
	// Add the result to the design variable vector
	stiff->addPointwiseMassDVSens(gpt, rho_alpha, dvSens, dvLen);
      }
    }
  }
  else if (matType == GEOMETRIC_STIFFNESS_MATRIX && 
	   type == LINEAR){
    // Geometric data
    TacsScalar X[3], Xd[9], Xdd[9];
    TacsScalar normal[3], normal_xi[3], normal_eta[3];
    
    // Transformation and the transformation derivative w.r.t. zeta
    TacsScalar t[9], tx[9], ztx[9]; 
    
    // Displacements/rotations and shape functions
    TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];
    
    // Interpolations for the shear components
    double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];
    
    // The interpolated tensorial shear components
    TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
    TacsScalar g13[NUM_G13], g23[NUM_G23];

    // The strain values
    TacsScalar strain[NUM_STRESSES];
    TacsScalar bstrain[NUM_STRESSES];

    // Compute the tying strain values
    const int is_linear = 1;
    compute_tying_strain<order>(is_linear,
				g11, g22, g12, g23, g13, 
				knots, pknots, vars, Xpts);

    for ( int m = 0; m < numGauss; m++ ){
      for ( int n = 0; n < numGauss; n++ ){
	double gpt[2];
	gpt[0] = gaussPts[n];
	gpt[1] = gaussPts[m];
		
	// Calculate the shape functions
	shell_hessian(order, X, Xd, Xdd, 
                      N, Na, Nb, Naa, Nab, Nbb,
                      gpt, Xpts);
	compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);
      
	TacsScalar h = 0.0;
	if (stiff->getTransformType() == FSDTStiffness::NATURAL){
	  h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, 
                                Xd, Xdd);
	}
	else {
	  const TacsScalar * axis = stiff->getRefAxis();
	  h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, 
                                        normal_eta, axis, Xd, Xdd);
	}
	h = gaussWts[n]*gaussWts[m]*h;

	// Calculate the strain/bmat at the current point... 
	TacsScalar rot;
	linear_bend_strain(strain, &rot, U, Ud, 
                           t, tx, ztx, 
                           normal, normal_xi, normal_eta);
      
	// Evaluate the strain interpolation at this point
	tying_interpolation<order>(gpt, N11, N22, N12, 
                                   knots, pknots);
	add_tying_strain<order>(strain, tx, g11, g22, g12, g23, g13,
                                N11, N22, N12);

	// Compute the derivatives of the phi/psi vectors
	TacsScalar Upsi[NUM_DISPS], Udpsi[2*NUM_DISPS];
	TacsScalar Uphi[NUM_DISPS], Udphi[2*NUM_DISPS];
	compute_shell_Ud(NUM_NODES, Upsi, Udpsi, psi, N, Na, Nb);
	compute_shell_Ud(NUM_NODES, Uphi, Udphi, phi, N, Na, Nb);

	// Compute the inner product of the second derivatives of the
	// stiffness matrix with the psi and phi vectors
	inner_nonlinear_bend_bmat(bstrain, Upsi, Udpsi, Uphi, Udphi,
				  t, tx, ztx, normal, normal_xi, normal_eta);
	add_nonlinear_tying_inner_nmat<order>(bstrain, psi, phi,
					      tx, N11, N22, N12,
					      knots, pknots, Xpts);

	// Add the term: alpha*psi^{T}*B^{T}*dC/dx*strain to the vector
	// dvSens - Note that this is much more efficient than computing
	// the terms component by component
	stiff->addStiffnessDVSens(gpt, alpha*h, strain, bstrain, 0.0, 0.0,
				  dvSens, dvLen);
      }
    }
  }
}

/*
  Add the derivative of the inner product of the given matrix type
  with the given psi and phi vectors with respect to the state
  variables.  This only makes sense for nonlinear Jacobian matrices
  such as the geometric stiffness matrix.

  input:
  matType:  the type of matrix
  alpha:    the scaling factor applied to the derivative
  psi:      the left-multiplying vector
  phi:      the right-multiplying vector
  vars:     the element state variables
  Xpts:     the nodal locations for this element

  output:
  res:      the derivative of the inner product w.r.t. vars
*/
template <int order>
void MITCShell<order>::addMatSVSensInnerProduct( ElementMatrixTypes matType,
						 TacsScalar alpha, 
						 TacsScalar res[],
						 const TacsScalar psi[], 
						 const TacsScalar phi[],
						 const TacsScalar vars[], 
						 const TacsScalar Xpts[] ){
  if (matType == GEOMETRIC_STIFFNESS_MATRIX && 
      type == LINEAR){  
    // Geometric data
    TacsScalar X[3], Xd[9], Xdd[9];
    TacsScalar normal[3], normal_xi[3], normal_eta[3];
    
    // Transformation and the transformation derivative w.r.t. zeta
    TacsScalar t[9], tx[9], ztx[9]; 
    
    // Displacements/rotations and shape functions
    TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];
    
    // Interpolations for the shear components
    double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];
    
    // The interpolated tensorial shear components
    TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
    TacsScalar g13[NUM_G13], g23[NUM_G23];

    // The derivatives of the displacement strain
    TacsScalar b11[3*NUM_NODES*NUM_G11], b22[3*NUM_NODES*NUM_G22];
    TacsScalar b12[3*NUM_NODES*NUM_G12];
    TacsScalar b13[NUM_VARIABLES*NUM_G13], b23[NUM_VARIABLES*NUM_G23];

    // The strain values
    TacsScalar strain[NUM_STRESSES];
    TacsScalar bstrain[NUM_STRESSES], bstress[NUM_STRESSES];

    // The derivative of the strain w.r.t. the displacements
    TacsScalar B[NUM_STRESSES*NUM_VARIABLES];
    TacsScalar drot[NUM_VARIABLES];

    // Store the stiffness information for the element
    TacsScalar At[6], Bt[6], Dt[6], Ats[3];

    // Evaluate the stiffness matrix
    if (!stiff->isVariableStiffness()){
      double gpt[2] = {0.0, 0.0};
      stiff->getStiffness(gpt, At, Bt, Dt, Ats);
    }

    // Compute the tying strain values
    const int is_linear = 1;
    compute_tying_bmat<order>(is_linear,
			      g11, g22, g12, g23, g13,
			      b11, b22, b12, b23, b13,
			      knots, pknots, vars, Xpts);

    for ( int m = 0; m < numGauss; m++ ){
      for ( int n = 0; n < numGauss; n++ ){
	double gpt[2];
	gpt[0] = gaussPts[n];
	gpt[1] = gaussPts[m];
		
	// Evaluate the stiffness matrix if required
	if (stiff->isVariableStiffness()){
	  stiff->getStiffness(gpt, At, Bt, Dt, Ats);
	}

	// Calculate the shape functions
	shell_hessian(order, X, Xd, Xdd, 
                      N, Na, Nb, Naa, Nab, Nbb,
                      gpt, Xpts);
	compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);
      
	TacsScalar h = 0.0;
	if (stiff->getTransformType() == FSDTStiffness::NATURAL){
	  h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, 
                                Xd, Xdd);
	}
	else {
	  const TacsScalar * axis = stiff->getRefAxis();
	  h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, 
                                        normal_eta, axis, Xd, Xdd);
	}
	h = alpha*gaussWts[n]*gaussWts[m]*h;

	// Calculate the strain/bmat at the current point... 
	TacsScalar rot;
	linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, 
                           normal, normal_xi, normal_eta);
	linear_bend_bmat(B, drot, NUM_NODES,
                         N, Na, Nb, t, tx, ztx,
			 normal, normal_xi, normal_eta);
      
	// Evaluate the strain interpolation at this point
	tying_interpolation<order>(gpt, N11, N22, N12, 
                                   knots, pknots);
	add_tying_strain<order>(strain, tx, g11, g22, g12, g23, g13,
                                N11, N22, N12);
	add_tying_bmat<order>(B, NUM_NODES, tx,
			      b11, b22, b12, b23, b13, 
			      N11, N22, N12);

	// Compute the derivatives of the phi/psi vectors
	TacsScalar Upsi[NUM_DISPS], Udpsi[2*NUM_DISPS];
	TacsScalar Uphi[NUM_DISPS], Udphi[2*NUM_DISPS];
	compute_shell_Ud(NUM_NODES, Upsi, Udpsi, psi, N, Na, Nb);
	compute_shell_Ud(NUM_NODES, Uphi, Udphi, phi, N, Na, Nb);

	// Compute the inner product of the second derivatives of the
	// stiffness matrix with the psi and phi vectors
	inner_nonlinear_bend_bmat(bstrain, Upsi, Udpsi, Uphi, Udphi,
				  t, tx, ztx, normal, normal_xi, normal_eta);
	add_nonlinear_tying_inner_nmat<order>(bstrain, psi, phi,
					      tx, N11, N22, N12,
					      knots, pknots, Xpts);

	// Compute the stress at the current Gauss point
	stiff->calcStress(At, Bt, Dt, Ats, bstrain, bstress);

	for ( int i = 0; i < NUM_NODES; i++ ){
	  for ( int ii = 0; ii < NUM_DISPS; ii++ ){	    
	    int row = ii + NUM_DISPS*i;	    
	    res[row] += 
	      h*strain_product(&B[NUM_STRESSES*row], bstress);
	  }
	}
      }
    }
  }
}

/*
  Compute the derivative of the residuals with respect to one of the
  material design variables.

  output:
  res:    the derivative of the residual w.r.t. the material design var
  
  input:
  dvNum:   the design variable number 
  vars:    the variables
  Xpts:    the element nodal locations
*/
template <int order>
void MITCShell<order>::getResDVSens( int dvNum, TacsScalar * res, 
                                     const TacsScalar vars[], 
				     const TacsScalar Xpts[] ){
  memset(res, 0, NUM_VARIABLES*sizeof(TacsScalar));
  
  if (dvNum >= 0 && stiff->ownsDesignVar(dvNum)){
    // Geometric data
    TacsScalar X[3], Xd[9], Xdd[9];
    TacsScalar normal[3], normal_xi[3], normal_eta[3];
    
    // Transformation and the transformation derivative w.r.t. zeta
    TacsScalar t[9], tx[9], ztx[9]; 
    
    TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];
    
    // Interpolations for the shear components
    double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];
    
    // The interpolated tensorial shear components
    TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
    TacsScalar g13[NUM_G13], g23[NUM_G23];
    
    // The derivatives of the displacement strain
    TacsScalar b11[3*NUM_NODES*NUM_G11], b22[3*NUM_NODES*NUM_G22];
    TacsScalar b12[3*NUM_NODES*NUM_G12];
    TacsScalar b13[NUM_VARIABLES*NUM_G13], b23[NUM_VARIABLES*NUM_G23];
    
    double gpt[2] = {0.0, 0.0};
    
    TacsScalar strain[NUM_STRESSES], stress[NUM_STRESSES];
    TacsScalar drot[NUM_VARIABLES];
    TacsScalar B[NUM_STRESSES*NUM_VARIABLES];
    
    // Store the stiffness information for the element
    TacsScalar sAt[6], sBt[6], sDt[6], sAts[3];
    TacsScalar k_penalty_sens = 0.0;

    // Evaluate the stiffness matrix
    if (!stiff->isVariableStiffness()){
      k_penalty_sens = stiff->getStiffnessDVSens(dvNum, gpt, 
						 sAt, sBt, sDt, sAts);
    }
    
    if (type == LARGE_ROTATION){
      compute_lr_tying_bmat<order>(g11, g22, g12, g23, g13, 
				   b11, b22, b12, b23, b13, 
				   knots, pknots, vars, Xpts);
    }
    else {
      compute_tying_bmat<order>((type == LINEAR),
				g11, g22, g12, g23, g13,
				b11, b22, b12, b23, b13,
				knots, pknots, vars, Xpts);
    }

    for ( int m = 0; m < numGauss; m++ ){
      for ( int n = 0; n < numGauss; n++ ){
	gpt[0] = gaussPts[n];
	gpt[1] = gaussPts[m];

        if (stiff->isVariableStiffness()){
          k_penalty_sens = 
            stiff->getStiffnessDVSens(dvNum, gpt, 
				      sAt, sBt, sDt, sAts);
        }
	
	// Calculate the shape functions
	shell_hessian(order, X, Xd, Xdd, 
                      N, Na, Nb, Naa, Nab, Nbb,
                      gpt, Xpts);
	compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);
	
	TacsScalar h = 0.0;
	if (stiff->getTransformType() == FSDTStiffness::NATURAL){
	  h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, 
                                Xd, Xdd);
	}
	else {
	  const TacsScalar * axis = stiff->getRefAxis();
	  h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, 
                                        normal_eta, axis, Xd, Xdd);
	}
	h = gaussWts[n]*gaussWts[m]*h;
	
	// Calculate the deformation at the current point... 
	TacsScalar rot;
	if (type == LINEAR){
	  linear_bend_strain(strain, &rot, U, Ud, 
                             t, tx, ztx, 
                             normal, normal_xi, normal_eta);
	  linear_bend_bmat(B, drot, NUM_NODES,
                           N, Na, Nb, t, tx, ztx,
                           normal, normal_xi, normal_eta);
	}
	else if (type == NONLINEAR){
	  nonlinear_bend_strain(strain, &rot, U, Ud,  t, tx, ztx, 
                                normal, normal_xi, normal_eta);
	  nonlinear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, U, Ud, t, tx, ztx,
                              normal, normal_xi, normal_eta);
	}
	else {
	  // Rotation matrix data
	  TacsScalar C[9], Ct[27], Ctt[54];
	  
	  // Compute the rotation matrices
	  TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
	  TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
	  TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
	  compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
	  compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);
	  
	  // Evaluate the in-plane rotation term
	  rot = compute_inplane_penalty(drot, NUM_NODES, Xd, Ud, 
					C, Ct, N, Na, Nb);
	  
	  // Calculate the deformation at the current point... 
	  large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, 
				normal, normal_xi, normal_eta);
	  large_rot_bend_bmat(B, NUM_NODES, N, Na, Nb, U, Ud, C, Ct, Ctt, 
			      t, tx, ztx, normal, normal_xi, normal_eta);
	}

	// Evaluate the strain interpolation at this point
	tying_interpolation<order>(gpt, N11, N22, N12, 
                                   knots, pknots);
	add_tying_strain<order>(strain, tx, 
                                g11, g22, g12, g23, g13,
                                N11, N22, N12);
	add_tying_bmat<order>(B, NUM_NODES, tx,
                              b11, b22, b12, b23, b13, 
                              N11, N22, N12);
	
        // Calculate the sensitivity of the stress at the current point
	stiff->calcStress(sAt, sBt, sDt, sAts, strain, stress);

	for ( int i = 0; i < NUM_NODES; i++ ){
	  for ( int ii = 0; ii < NUM_DISPS; ii++ ){	    
	    int row = ii + NUM_DISPS*i;	    
	    res[row] += 
	      h*(strain_product(&B[NUM_STRESSES*row], stress) +
		 k_penalty_sens*rot*drot[row]);
	  }
	}
      }
    }
  }
}
  
/*
  Get the derivative of the specified matrix type with respect to the
  element nodal locations.
*/
template <int order>
void MITCShell<order>::getMatTypeXptSens( ElementMatrixTypes matType, 
					  TacsScalar scaleFactor, 
					  TacsScalar * mat, 
					  const TacsScalar vars[],
					  const TacsScalar Xpts[], 
					  const TacsScalar XptSens[], 
					  MatrixOrientation matOr ){
  memset(mat, 0, NUM_VARIABLES*NUM_VARIABLES*sizeof(TacsScalar));
}

/*
  Get the derivative of the specified matrix type with respect to the
  material design variables.
  
  output:
  mat:         the element matrix 

  input:
  dvNum:       the design variable number
  matType:     the element matrix type
  scaleFactor: the scale factor 
  vars:        the element variables
  Xpts:        the element nodes
  matOr:       the matrix orientation
*/
template <int order>
void MITCShell<order>::getMatTypeDVSens( int dvNum,
                                         ElementMatrixTypes matType, 
					 TacsScalar scaleFactor, 
					 TacsScalar * mat, 
					 const TacsScalar vars[], 
					 const TacsScalar Xpts[], 
					 MatrixOrientation matOr ){

  memset(mat, 0, NUM_VARIABLES*NUM_VARIABLES*sizeof(TacsScalar));
  
  if (dvNum >= 0 && stiff->ownsDesignVar(dvNum)){
    if (matType == MASS_MATRIX){
      double gpt[2] = {0.0, 0.0};
      double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];

      for ( int m = 0; m < numGauss; m++ ){
        for ( int n = 0; n < numGauss; n++ ){
          gpt[0] = gaussPts[n];
          gpt[1] = gaussPts[m];
	
          TacsScalar X[3], Xd[9], normal[3];
          shell_jacobian(order, X, Xd, N, Na, Nb, gpt, Xpts);
          Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
          Tensor::normalize3D(normal);
	
          for ( int i = 0; i < 3; i++ ){
            Xd[6+i] = normal[i]; 
          }
	  
          TacsScalar h = FElibrary::jacobian3d(Xd);
          h = scaleFactor*gaussWts[n]*gaussWts[m]*h;

          TacsScalar massDVSens[3];
          stiff->pointwiseMassDVSens(dvNum, gpt, massDVSens);
	  
          for ( int i = 0; i < NUM_NODES; i++ ){
            for ( int ii = 0; ii < 3; ii++ ){
              int row = ii + NUM_DISPS*i;	      
              for ( int j = 0; j < NUM_NODES; j++ ){
                int col = ii + NUM_DISPS*j;		
                mat[col + row*NUM_VARIABLES] += h*massDVSens[0]*N[i]*N[j];
              }
            }
          }
	  
          for ( int ii = 3; ii < NUM_DISPS; ii++ ){	 
            for ( int jj = 3; jj < NUM_DISPS; jj++ ){
              TacsScalar a[3] = { 0.0, 0.0, 0.0 };
              TacsScalar b[3] = { 0.0, 0.0, 0.0 };
	      
              TacsScalar phi_a[3], phi_b[3];
              a[ii-3] = 1.0;
              b[jj-3] = 1.0;
	      
              Tensor::crossProduct3D(phi_a, a, normal);
              Tensor::crossProduct3D(phi_b, b, normal);
              TacsScalar weight = Tensor::dot3D(phi_a, phi_b);
	      
              for ( int i = 0; i < NUM_NODES; i++ ){
                int row = ii + NUM_DISPS*i;

                for ( int j = 0; j < NUM_NODES; j++ ){
                  int col = jj + NUM_DISPS*j;
		  
                  mat[col + row*NUM_VARIABLES] += 
                    massDVSens[2]*weight*h*N[i]*N[j];
                }
              }
            }	  
          }
        }    
      }
    }    
    else if (matType == STIFFNESS_MATRIX &&
	     type != LARGE_ROTATION){
      TacsScalar X[3], Xd[9], Xdd[9];
      TacsScalar normal[3], normal_xi[3], normal_eta[3];
      // Transformation and the transformation derivative w.r.t. zeta
      TacsScalar t[9], tx[9], ztx[9]; 
  
      TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
      double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
      double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

      // Interpolations for the shear components
      double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

      // The interpolated tensorial shear components
      TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
      TacsScalar g13[NUM_G13], g23[NUM_G23];

      // The derivatives of the displacement strain
      TacsScalar b11[3*NUM_NODES*NUM_G11], b22[3*NUM_NODES*NUM_G22];
      TacsScalar b12[3*NUM_NODES*NUM_G12];
      TacsScalar b13[NUM_VARIABLES*NUM_G13], b23[NUM_VARIABLES*NUM_G23];
  
      double gpt[2] = {0.0, 0.0};

      TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];
      TacsScalar B[NUM_STRESSES*NUM_VARIABLES];
      TacsScalar BStress[NUM_STRESSES];

      TacsScalar drot[NUM_VARIABLES];

      // Store the stiffness information for the element
      TacsScalar sAt[6], sBt[6], sDt[6], sAts[3];
      TacsScalar k_penalty_sens = 0.0;

      // Evaluate the stiffness matrix
      if (!stiff->isVariableStiffness()){
        k_penalty_sens = stiff->getStiffnessDVSens(dvNum, gpt, 
                                                        sAt, sBt, sDt, sAts);
      }

      compute_tying_bmat<order>((type == LINEAR),
				g11, g22, g12, g23, g13,
				b11, b22, b12, b23, b13,
				knots, pknots, vars, Xpts);

      for ( int m = 0; m < numGauss; m++ ){
        for ( int n = 0; n < numGauss; n++ ){
          gpt[0] = gaussPts[n];
          gpt[1] = gaussPts[m];
	
          // Evaluate the stiffness matrix
          if (stiff->isVariableStiffness()){
            k_penalty_sens = 
              stiff->getStiffnessDVSens(dvNum, gpt, 
                                             sAt, sBt, sDt, sAts);
          }

          // Calculate the shape functions
          shell_hessian(order, X, Xd, Xdd, 
                        N, Na, Nb, Naa, Nab, Nbb,
                        gpt, Xpts);
          compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);
	     
          TacsScalar h = 0.0;
          if (stiff->getTransformType() == FSDTStiffness::NATURAL ){
            h = compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, 
                                  Xd, Xdd);
          }
          else {
            const TacsScalar * axis = stiff->getRefAxis();
            h = compute_transform_refaxis(t, tx, ztx, normal, normal_xi, 
                                          normal_eta, axis, Xd, Xdd);
          }
          h = gaussWts[n]*gaussWts[m]*h;
	
          // Calculate the strain at the current point... 
          TacsScalar rot;
          if (type == LINEAR){
            linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, 
                               normal, normal_xi, normal_eta);
            linear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, t, tx, ztx,
                             normal, normal_xi, normal_eta);
          }
          else {
            nonlinear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, 
                                  normal, normal_xi, normal_eta);
            nonlinear_bend_bmat(B, drot, NUM_NODES, N, Na, Nb, U, Ud, t, tx, ztx,
                                normal, normal_xi, normal_eta);
          }
	
          // Evaluate the strain interpolation at this point
          tying_interpolation<order>(gpt, N11, N22, N12, knots, pknots);
          add_tying_strain<order>(strain, tx, 
                                  g11, g22, g12, g23, g13,
                                  N11, N22, N12);
          add_tying_bmat<order>(B, NUM_NODES, tx,
                                b11, b22, b12, b23, b13, 
                                N11, N22, N12);
	
          stiff->calcStress(sAt, sBt, sDt, sAts, strain, stress);
	
          // Loop through all displacements except for the z-rotation
          for ( int i = 0; i < NUM_NODES; i++ ){
            for ( int ii = 0; ii < NUM_DISPS; ii++ ){
              int row = ii + NUM_DISPS*i;
	    
              // Calculate the stress associated with Bi
              stiff->calcStress(sAt, sBt, sDt, sAts, 
                                     &B[row*NUM_STRESSES], BStress);
	    
              for ( int j = 0; j <= i; j++ ){
                int end = ( j == i ? ii : NUM_DISPS-1 );
                for ( int jj = 0; jj <= end; jj++ ){
                  int col = jj + NUM_DISPS*j;
		
                  // The regular element matrix
                  mat[col + row*NUM_VARIABLES] += 
                    h*(strain_product(BStress, &B[col*NUM_STRESSES]) +    
                       k_penalty_sens*drot[row]*drot[col]);
                }
              }
            }
          }

          if (type == NONLINEAR){
            nonlinear_bend_stress_bmat(mat, NUM_NODES, h, stress, 
                                       N, Na, Nb, t, tx, ztx,
                                       normal, normal_xi, 
                                       normal_eta);
          }
        }
      }
    
      // Copy over the matrix
      for ( int row = 0; row < NUM_VARIABLES; row++ ){
        for ( int col = row+1; col < NUM_VARIABLES; col++ ){
          mat[col + row*NUM_VARIABLES] = mat[row + col*NUM_VARIABLES];
        }
      }
    }
  }
}

/*
  Return the number of points in the specified quadrature scheme
*/
template <int order>
int MITCShell<order>::getNumGaussPts(){
  return numGauss*numGauss;
}

/*
  Retrieve the Gauss points and weights for the given Gauss point
  number

  input:
  scheme: the Gauss quadrature scheme
  num:    the Gauss point index for the quadrature point
  
  returns: the Gauss quadrature weight for the given point

  output:
  pt:   the Gauss point for the given index
*/
template <int order>
TacsScalar MITCShell<order>::getGaussWtsPts( const int num, 
					     double * pt ){
  int m = (int)(num/numGauss);
  int n = num % numGauss;
  pt[0] = gaussPts[n];
  pt[1] = gaussPts[m];    
  
  return gaussWts[n]*gaussWts[m];
}

/*
  Evaluate the shape functions for this element at the specified point
  
  output:
  N:  the shape functions values evaluated at the parametric point pt

  input:
  pt: the parametric point within the element
*/
template <int order>
void MITCShell<order>::getShapeFunctions( const double pt[], double N[] ){
  double na[order], nb[order];
  FElibrary::lagrangeSF(na, pt[0], order);
  FElibrary::lagrangeSF(nb, pt[1], order);

  for ( int j = 0; j < order; j++ ){
    for ( int i = 0; i < order; i++ ){
      N[0] = na[i]*nb[j];
      N++;
    }
  }
}

/*
  Evaluate the determinant of the Jacobian for numerical integration

  returns: the determinant of the Jacobian
 
  input:
  pt: the parametric point within the element
  Xpts: the element nodes
*/
template <int order>
TacsScalar MITCShell<order>::getDetJacobian( const double * pt, 
                                             const TacsScalar Xpts[] ){
  TacsScalar X[3], Xd[9];
  TacsScalar normal[3];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];

  shell_jacobian(order, X, Xd, N, Na, Nb, pt, Xpts);
  Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
  Tensor::normalize3D(normal);

  for ( int i = 0; i < 3; i++ ){
    Xd[i+6] = normal[i];
  }

  return FElibrary::jacobian3d(Xd);
}

/*
  Evaluate the derivative of the determinant of the Jacobian with respect
  to the element nodal locations
  
  output:
  hXptSens:  the derivative of the determinant w.r.t. the nodal locations

  returns: the determinant of the Jacobian
 
  input:
  pt: the parametric point within the element
  Xpts: the element nodes
*/
template <int order>
TacsScalar MITCShell<order>::getDetJacobianXptSens( TacsScalar * hXptSens, 
                                                    const double * pt, 
                                                    const TacsScalar Xpts[] ){
  TacsScalar h = 0.0;
  TacsScalar X[3], Xd[9], XdSens[9];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];

  shell_jacobian(order, X, Xd, N, Na, Nb, pt, Xpts);

  for ( int i = 0; i < NUM_NODES; i++ ){
    for ( int k = 0; k < 3; k++ ){
      XdSens[0] = XdSens[1] = XdSens[2] = 0.0;
      XdSens[3] = XdSens[4] = XdSens[5] = 0.0;
      XdSens[k]   = Na[i];
      XdSens[3+k] = Nb[i];

      Tensor::crossProduct3DSens(&Xd[6], &XdSens[6],
				 &Xd[0], &Xd[3], &XdSens[0], &XdSens[3]);
      TacsScalar snrm;
      Tensor::normalize3DSens(&snrm, &Xd[6], &XdSens[6]);
      h = FElibrary::jacobian3dSens(Xd, XdSens, &hXptSens[0]);
      hXptSens++;
    }
  }

  return h;
}

/*
  Evaluate the strain at the specified point using the provided
  set of variables

  output:
  strain:   the strain evaluate at the specific parametric point
  
  input:
  vars:     the element variable values
  Xpts:     the element nodal locations
*/
template <int order>
void MITCShell<order>::getStrain( TacsScalar strain[], 
				  const double * pt,
				  const TacsScalar Xpts[],
				  const TacsScalar vars[] ){
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9]; 
  
  TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the displacement strain
  TacsScalar b11[3*NUM_NODES*NUM_G11], b22[3*NUM_NODES*NUM_G22];
  TacsScalar b12[3*NUM_NODES*NUM_G12];
  TacsScalar b13[NUM_VARIABLES*NUM_G13], b23[NUM_VARIABLES*NUM_G23];

  shell_hessian(order, X, Xd, Xdd, 
                N, Na, Nb, Naa, Nab, Nbb,
                pt, Xpts);
  compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

  if (stiff->getTransformType() == FSDTStiffness::NATURAL ){
    compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, 
                      Xd, Xdd);
  }
  else {
    const TacsScalar * axis = stiff->getRefAxis();
    compute_transform_refaxis(t, tx, ztx, normal, normal_xi, 
                              normal_eta, axis, Xd, Xdd);
  }

  if (type == LARGE_ROTATION){
    compute_lr_tying_bmat<order>(g11, g22, g12, g23, g13, 
				 b11, b22, b12, b23, b13, 
				 knots, pknots, vars, Xpts);
  }
  else {
    compute_tying_strain<order>((type == LINEAR),
				g11, g22, g12, g23, g13,
				knots, pknots, vars, Xpts);
  }

  TacsScalar rot;
  if (type == LINEAR){
    linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, 
                       normal, normal_xi, normal_eta);
  }
  else if (type == NONLINEAR){
    nonlinear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, 
                          normal, normal_xi, normal_eta);    
  }
  else {
    // Rotation matrix data
    TacsScalar C[9], Ct[27];

    // Compute the rotation matrices
    TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
    TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
    TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
    compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
    
    compute_lr_tying_strain<order>(g11, g22, g12, g23, g13, 
				   knots, pknots, vars, Xpts);
    large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, 
			  normal, normal_xi, normal_eta);
  }

  // Evaluate the strain interpolation at this point
  tying_interpolation<order>(pt, N11, N22, N12, 
                             knots, pknots);
  add_tying_strain<order>(strain, tx, 
                          g11, g22, g12, g23, g13,
                          N11, N22, N12);
}

/*
  Compute the strain and the derivative of the strain with respect to
  the nodal locations.

  output:
  strain:  the strain evaluate at the pamametric point pt
  strainXptSens: the derivative of the straint w.r.t. the nodal locations

  input:
  pt:            the parametric point within the element
  vars:          the element variables
  Xpts:          the nodal locations
*/
template <int order>
void MITCShell<order>::addStrainXptSens( TacsScalar strain[], 
					 TacsScalar strainXptSens[], 
					 const double * pt,
					 const TacsScalar vars[], 
					 const TacsScalar Xpts[] ){ 
  
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];
  TacsScalar dnormal[9*NUM_NODES], dnormal_xi[9*NUM_NODES], 
    dnormal_eta[9*NUM_NODES];

  // The derivative of the determinant of the Jacobian and the 
  // sensitivity of the drilling rotation
  TacsScalar dh[3*NUM_NODES], srot[3*NUM_NODES];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9]; 
  TacsScalar dt[27*NUM_NODES], dtx[27*NUM_NODES], dztx[27*NUM_NODES]; 
  
  TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial strain components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the tensorial tying strain components
  TacsScalar dg11[3*NUM_G11*NUM_NODES], dg22[3*NUM_G22*NUM_NODES], 
    dg12[3*NUM_G12*NUM_NODES];
  TacsScalar dg13[3*NUM_G13*NUM_NODES], dg23[3*NUM_G23*NUM_NODES];

  // Evaluate the tying interpolation
  tying_interpolation<order>(pt, N11, N22, N12, 
                             knots, pknots);

  // Calculate the shape functions
  shell_hessian(order, X, Xd, Xdd, 
                N, Na, Nb, Naa, Nab, Nbb,
                pt, Xpts);

  // Evaluate the displacements and their derivatives
  compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

  if (stiff->getTransformType() == FSDTStiffness::NATURAL){
    compute_transform_sens(dh, t, dt, tx, dtx, ztx, dztx, 
                           normal, dnormal, 
                           normal_xi, dnormal_xi, 
                           normal_eta, dnormal_eta, 
                           Xd, Xdd, Na, Nb, Naa, Nab, Nbb, NUM_NODES);
  }
  else {
    const TacsScalar * axis = stiff->getRefAxis();
    compute_transform_refaxis_sens(dh, t, dt, tx, dtx, ztx, dztx, 
                                   normal, dnormal, 
                                   normal_xi, dnormal_xi, 
                                   normal_eta, dnormal_eta, axis,
                                   Xd, Xdd, Na, Nb, Naa, Nab, Nbb, NUM_NODES);
  }

  if (type == LARGE_ROTATION){
    // Rotation matrix data
    TacsScalar C[9], Ct[27];
    
    // Evaluate the rotation matrices
    TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
    TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
    TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
    compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
    
    compute_lr_tying_strain_sens<order>(g11, g22, g12, g23, g13, 
					dg11, dg22, dg12, dg23, dg13,
					knots, pknots, vars, Xpts);
    large_rot_bend_strain_sens(strain, strainXptSens, 
			       U, Ud, C, Ct, t, dt, tx, dtx, ztx, dztx,
			       normal, dnormal, normal_xi, dnormal_xi, 
			       normal_eta, dnormal_eta, 3*NUM_NODES);
  }
  else {
    compute_tying_strain_sens<order>((type == LINEAR),
				     g11, g22, g12, g23, g13, 
				     dg11, dg22, dg12, dg23, dg13, 
				     knots, pknots, vars, Xpts);
    
    TacsScalar rot;
    if (type == LINEAR){
      linear_bend_strain_sens(strain, strainXptSens, 
			      &rot, srot, U, Ud, 
			      t, dt, tx, dtx, ztx, dztx,
			      normal, dnormal, normal_xi, dnormal_xi, 
			      normal_eta, dnormal_eta, 3*NUM_NODES);
    }
    else {
      nonlinear_bend_strain_sens(strain, strainXptSens, 
				 &rot, srot, U, Ud, 
				 t, dt, tx, dtx, ztx, dztx,
				 normal, dnormal, normal_xi, dnormal_xi, 
				 normal_eta, dnormal_eta, 3*NUM_NODES);
    }
  }
    
  // Evaluate the strain interpolation at this point
  add_tying_strain_sens<order>(strain, strainXptSens, tx, dtx, 
                               g11, g22, g12, g23, g13,
                               dg11, dg22, dg12, dg23, dg13,
                               N11, N22, N12);
}

/*
  Compute the derivative of the point-wise strain with respect to the
  element variables, multiply the derivative by a strain sensitivity
  vector and add the result, times a scalar multiple, to the outputt
  array. This can be used to evaluate the derivative of functions of
  the strain with respect to the state variables using the chain rule.

  output:
  elementSens: the output array - same length as the number of elem variables
 
  input:
  pt:          parametric point used to evaluate the derivative [-1, 1]^{2}
  scaleFactor: scale ther result by this scalar
  strainSens:  the sensitivity of each straint component 
  vars:        the element variables
  Xpts:        the element nodal locations
*/
template <int order>
void MITCShell<order>::addStrainSVSens( TacsScalar elementSens[], 
					const double * pt, 
					const TacsScalar scaleFactor, 
					const TacsScalar strainSens[],
					const TacsScalar vars[],  
					const TacsScalar Xpts[] ){

  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];
  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9]; 
  
  TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The tensorial components of the strain
  TacsScalar g11[NUM_G11], g22[NUM_G22];
  TacsScalar g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The derivatives of the displacement strain
  TacsScalar b11[3*NUM_NODES*NUM_G11], b22[3*NUM_NODES*NUM_G22];
  TacsScalar b12[3*NUM_NODES*NUM_G12];
  TacsScalar b13[NUM_VARIABLES*NUM_G13], b23[NUM_VARIABLES*NUM_G23];

  TacsScalar dinplane_rot[NUM_VARIABLES];
  TacsScalar B[NUM_STRESSES*NUM_VARIABLES];

  // Calculate the shape functions
  shell_hessian(order, X, Xd, Xdd, 
                N, Na, Nb, Naa, Nab, Nbb,
                pt, Xpts);
  compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);

  if (stiff->getTransformType() == FSDTStiffness::NATURAL ){
    compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, 
                      Xd, Xdd);
  }
  else {
    const TacsScalar * axis = stiff->getRefAxis();
    compute_transform_refaxis(t, tx, ztx, normal, normal_xi, 
                              normal_eta, axis, Xd, Xdd);
  }

  if (type == LARGE_ROTATION){
    // Rotational matrix data
    TacsScalar C[9], Ct[27], Ctt[54];
    
    // Compute the rotation matrices
    TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
    TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
    TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
    compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
    compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);
    
    compute_lr_tying_bmat<order>(g11, g22, g12, g23, g13,
				 b11, b22, b12, b23, b13, 
				 knots, pknots, vars, Xpts);
    large_rot_bend_bmat(B, NUM_NODES, N, Na, Nb, U, Ud, C, Ct, Ctt,
			t, tx, ztx, normal, normal_xi, normal_eta);
  }
  else {
    compute_tying_bmat<order>((type == LINEAR),
			      g11, g22, g12, g23, g13,
			      b11, b22, b12, b23, b13, 
			      knots, pknots, vars, Xpts);
    
    if (type == LINEAR){
      linear_bend_bmat(B, dinplane_rot, NUM_NODES,
		       N, Na, Nb, t, tx, ztx,
		       normal, normal_xi, normal_eta);
    }
    else {
      nonlinear_bend_bmat(B, dinplane_rot, NUM_NODES,
			  N, Na, Nb, U, Ud, t, tx, ztx,
			  normal, normal_xi, normal_eta);
    }
  }

  tying_interpolation<order>(pt, N11, N22, N12, 
                             knots, pknots);
  add_tying_bmat<order>(B, NUM_NODES, tx,
                        b11, b22, b12, b23, b13, 
                        N11, N22, N12);

  for ( int k = 0; k < NUM_VARIABLES; k++ ){
    elementSens[k] += scaleFactor*strain_product(strainSens, 
						 &B[k*NUM_STRESSES]);
  }
}

/*
  Determine the number of nodes and elements for visualization 
  generated by the data in this element. Note that higher-order
  elements are broken down into bi-linear elements for visualization.

  output:
  nelems:  the number of visualization elements
  nnodes:  the number of nodes used for visualization
  ncsr:    the number of entries in a CSR-type data structure used
  to store the connectivity
*/
template <int order>
void MITCShell<order>::addOutputCount( int * nelems, 
                                       int * nnodes, int * ncsr ){
  *nelems += (order-1)*(order-1);
  *nnodes += order*order;
  *ncsr += 4*(order-1)*(order-1);
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
void MITCShell<order>::getOutputData( unsigned int out_type,
				      double * data, int ld_data,
				      const TacsScalar Xpts[],
				      const TacsScalar vars[] ){
  // Geometric data
  TacsScalar X[3], Xd[9], Xdd[9];
  TacsScalar normal[3], normal_xi[3], normal_eta[3];

  // Transformation and the transformation derivative w.r.t. zeta
  TacsScalar t[9], tx[9], ztx[9]; 

  TacsScalar U[NUM_DISPS], Ud[2*NUM_DISPS];
  double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
  double Naa[NUM_NODES], Nab[NUM_NODES], Nbb[NUM_NODES];

  // Interpolations for the shear components
  double N11[NUM_G11], N22[NUM_G22], N12[NUM_G12];

  // The interpolated tensorial shear components
  TacsScalar g11[NUM_G11], g22[NUM_G22], g12[NUM_G12];
  TacsScalar g13[NUM_G13], g23[NUM_G23];

  // The stress and strain values
  TacsScalar stress[NUM_STRESSES], strain[NUM_STRESSES];

  if (type == LARGE_ROTATION){
    compute_lr_tying_strain<order>(g11, g22, g12, g23, g13, 
				   knots, pknots, vars, Xpts);
  }
  else {
    compute_tying_strain<order>((type == LINEAR),
				g11, g22, g12, g23, g13, 
				knots, pknots, vars, Xpts);
  }

  for ( int m = 0; m < order; m++ ){
    for ( int n = 0; n < order; n++ ){
      double pt[2];
      pt[0] = -1.0 + 2.0/(order - 1.0 )*n;
      pt[1] = -1.0 + 2.0/(order - 1.0 )*m;

      // Calculate the shape functions    
      shell_hessian(order, X, Xd, Xdd, 
                    N, Na, Nb, Naa, Nab, Nbb,
                    pt, Xpts);
      compute_shell_Ud(NUM_NODES, U, Ud, vars, N, Na, Nb);
	
      if (stiff->getTransformType() == FSDTStiffness::NATURAL){
        compute_transform(t, tx, ztx, normal, normal_xi, normal_eta, 
                          Xd, Xdd);
      }
      else {
        const TacsScalar * axis = stiff->getRefAxis();
        compute_transform_refaxis(t, tx, ztx, normal, normal_xi, 
                                  normal_eta, axis, Xd, Xdd);
      }
      
      // Evaluate the strain
      TacsScalar rot;
      if (type == LINEAR){
        linear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, 
                           normal, normal_xi, normal_eta);
      }
      else if (type == NONLINEAR){
        nonlinear_bend_strain(strain, &rot, U, Ud, t, tx, ztx, 
                              normal, normal_xi, normal_eta);
      }
      else {
	// Rotational matrix data
	TacsScalar C[9], Ct[27];
	
	// Compute the rotation matrices
	TacsScalar c1 = cos(U[3]), s1 = sin(U[3]);
	TacsScalar c2 = cos(U[4]), s2 = sin(U[4]);
	TacsScalar c3 = cos(U[5]), s3 = sin(U[5]);
	compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
	
	// Evaluate the strain
	large_rot_bend_strain(strain, U, Ud, C, Ct, t, tx, ztx, 
			      normal, normal_xi, normal_eta);
      }
      
      // Evaluate the strain interpolation at this point
      tying_interpolation<order>(pt, N11, N22, N12, 
                                 knots, pknots);
      add_tying_strain<order>(strain, tx, 
                              g11, g22, g12, g23, g13,
                              N11, N22, N12);

      int index = 0;
      int p = n + m*order;
      if (out_type & TACSElement::OUTPUT_NODES){
	for ( int k = 0; k < 3; k++ ){
	  data[index+k] = RealPart(Xpts[3*p+k]);
	}
        index += 3;
      }
      if (out_type & TACSElement::OUTPUT_DISPLACEMENTS){
	for ( int k = 0; k < NUM_DISPS; k++ ){
	  data[index+k] = RealPart(vars[NUM_DISPS*p+k]);
	}
        index += NUM_DISPS;
      }
      if (out_type & TACSElement::OUTPUT_STRAINS){
        for ( int k = 0; k < NUM_STRESSES; k++ ){
          data[index+k] = RealPart(strain[k]);
        }
        index += NUM_STRESSES;
      }
      if (out_type & TACSElement::OUTPUT_STRESSES){
        // Evaluate the stiffness at the current point 
        // and then calculate the stress
        stiff->calculateStress(pt, strain, stress);
        
        for ( int k = 0; k < NUM_STRESSES; k++ ){
          data[index+k] = RealPart(stress[k]);
	}
        index += NUM_STRESSES;
      }
      if (out_type & TACSElement::OUTPUT_EXTRAS){
	// Compute the failure value
	TacsScalar lambda;
	stiff->failure(pt, strain, &lambda);
	data[index] = RealPart(lambda);

	// Compute the buckling constraint value
	TacsScalar bval;
	stiff->buckling(strain, &bval);
	data[index+1] = RealPart(bval);

	data[index+2] = RealPart(stiff->getDVOutputValue(0, pt));
	data[index+3] = RealPart(stiff->getDVOutputValue(1, pt));

        index += NUM_EXTRAS;
      }
      if (out_type & TACSElement::OUTPUT_COORDINATES){
        // t is the transform from the global coordinates to the
        // local coordinates.
        for ( int i = 0; i < 3; i++ ){
          for ( int j = 0; j < 3; j++ ){
            data[index + 3*i + j] = RealPart(t[3*i + j]);
          }
        }
        index += 9;
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
void MITCShell<order>::getOutputConnectivity( int * con, int node ){
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
}

#endif
