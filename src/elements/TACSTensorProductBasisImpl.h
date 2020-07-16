

#ifndef TACS_TENSOR_PRODUCT_BASIS_IMPL_H
#define TACS_TENSOR_PRODUCT_BASIS_IMPL_H

#include "TACSObject.h"


/*
  3D Tensor product functions for p = 1
*/
void TACSInterpAllTensor3DInterp2( const int m,
                                   const TacsScalar N[],
                                   const TacsScalar Nx[],
                                   const TacsScalar values[],
                                   TacsScalar out[] );
// void TACSInterpAllTensor3DInterp2VarsPerNode1( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
// void TACSInterpAllTensor3DInterp2VarsPerNode3( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
// void TACSInterpAllTensor3DInterp2VarsPerNode4( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
void TacsAddAllTransTensor3DInterp2( const int m,
                                     const TacsScalar N[],
                                     const TacsScalar Nx[],
                                     const TacsScalar in[],
                                     TacsScalar values[] );
// void TACSInterpAllTensor3DInterp2VarsPerNode1( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
// void TACSInterpAllTensor3DInterp2VarsPerNode3( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
// void TACSInterpAllTensor3DInterp2VarsPerNode4( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );

/*
  3D Tensor product functions for p = 2
*/
void TACSInterpAllTensor3DInterp3( const int m,
                                   const TacsScalar N[],
                                   const TacsScalar Nx[],
                                   const TacsScalar values[],
                                   TacsScalar out[] );
// void TACSInterpAllTensor3DInterp5VarsPerNode1( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
// void TACSInterpAllTensor3DInterp5VarsPerNode3( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
// void TACSInterpAllTensor3DInter56VarsPerNode4( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
void TacsAddAllTransTensor3DInterp3( const int m,
                                     const TacsScalar N[],
                                     const TacsScalar Nx[],
                                     const TacsScalar in[],
                                     TacsScalar values[] );
// void TACSInterpAllTensor3DInterp3VarsPerNode1( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
// void TACSInterpAllTensor3DInterp3VarsPerNode3( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
// void TACSInterpAllTensor3DInterp3VarsPerNode4( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );


/*
  3D Tensor product functions for p = 3
*/
void TACSInterpAllTensor3DInterp4( const int m,
                                   const TacsScalar N[],
                                   const TacsScalar Nx[],
                                   const TacsScalar values[],
                                   TacsScalar out[] );
// void TACSInterpAllTensor3DInterp4VarsPerNode1( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
// void TACSInterpAllTensor3DInterp4VarsPerNode3( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
// void TACSInterpAllTensor3DInterp4VarsPerNode4( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
void TacsAddAllTransTensor3DInterp4( const int m,
                                     const TacsScalar N[],
                                     const TacsScalar Nx[],
                                     const TacsScalar in[],
                                     TacsScalar values[] );
// void TACSInterpAllTensor3DInterp4VarsPerNode1( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
// void TACSInterpAllTensor3DInterp4VarsPerNode3( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );
// void TACSInterpAllTensor3DInterp4VarsPerNode4( const TacsScalar N[],
//                                                const TacsScalar Nx[],
//                                                const TacsScalar values[],
//                                                TacsScalar out[] );

/*
  3D Tensor product functions for p = 4
*/
void TACSInterpAllTensor3DInterp5( const int m,
                                   const TacsScalar N[],
                                   const TacsScalar Nx[],
                                   const TacsScalar values[],
                                   TacsScalar out[] );
void TACSInterpAllTensor3DInterp5VarsPerNode1( const TacsScalar N[],
                                               const TacsScalar Nx[],
                                               const TacsScalar values[],
                                               TacsScalar out[] );
void TACSInterpAllTensor3DInterp5VarsPerNode3( const TacsScalar N[],
                                               const TacsScalar Nx[],
                                               const TacsScalar values[],
                                               TacsScalar out[] );
void TACSInterpAllTensor3DInterp5VarsPerNode4( const TacsScalar N[],
                                               const TacsScalar Nx[],
                                               const TacsScalar values[],
                                               TacsScalar out[] );

void TacsAddAllTransTensor3DInterp5( const int m,
                                     const TacsScalar N[],
                                     const TacsScalar Nx[],
                                     const TacsScalar in[],
                                     TacsScalar values[] );
void TacsAddAllTransTensor3DInterp5VarsPerNode1( const TacsScalar N[],
                                                 const TacsScalar Nx[],
                                                 const TacsScalar in[],
                                                 TacsScalar values[] );
void TacsAddAllTransTensor3DInterp5VarsPerNode3( const TacsScalar N[],
                                                 const TacsScalar Nx[],
                                                 const TacsScalar in[],
                                                 TacsScalar values[] );
void TacsAddAllTransTensor3DInterp5VarsPerNode4( const TacsScalar N[],
                                                 const TacsScalar Nx[],
                                                 const TacsScalar in[],
                                                 TacsScalar values[] );

/*
  3D Tensor product functions for p = 5
*/
void TACSInterpAllTensor3DInterp6( const int m,
                                   const TacsScalar N[],
                                   const TacsScalar Nx[],
                                   const TacsScalar values[],
                                   TacsScalar out[] );
void TACSInterpAllTensor3DInterp6VarsPerNode1( const TacsScalar N[],
                                               const TacsScalar Nx[],
                                               const TacsScalar values[],
                                               TacsScalar out[] );
void TACSInterpAllTensor3DInterp6VarsPerNode3( const TacsScalar N[],
                                               const TacsScalar Nx[],
                                               const TacsScalar values[],
                                               TacsScalar out[] );
void TACSInterpAllTensor3DInterp6VarsPerNode4( const TacsScalar N[],
                                               const TacsScalar Nx[],
                                               const TacsScalar values[],
                                               TacsScalar out[] );

void TacsAddAllTransTensor3DInterp6( const int m,
                                     const TacsScalar N[],
                                     const TacsScalar Nx[],
                                     const TacsScalar in[],
                                     TacsScalar values[] );
void TacsAddAllTransTensor3DInterp6VarsPerNode1( const TacsScalar N[],
                                                 const TacsScalar Nx[],
                                                 const TacsScalar in[],
                                                 TacsScalar values[] );
void TacsAddAllTransTensor3DInterp6VarsPerNode3( const TacsScalar N[],
                                                 const TacsScalar Nx[],
                                                 const TacsScalar in[],
                                                 TacsScalar values[] );
void TacsAddAllTransTensor3DInterp6VarsPerNode4( const TacsScalar N[],
                                                 const TacsScalar Nx[],
                                                 const TacsScalar in[],
                                                 TacsScalar values[] );

#endif // TACS_TENSOR_PRODUCT_BASIS_IMPL_H