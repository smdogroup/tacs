#ifndef TACS_AMD_INTERFACE_H
#define TACS_AMD_INTERFACE_H

/*
  An approximate minimum degree reordering scheme for the schur
  complement method. Complete documentation is provided in the .c
  file. Note that the contents of rowp/cols are overwritten during
  this call so be careful!

  Copyright (c) 2011 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

void amd_order_interface( int nvars, int * rowp, int * cols, int * perm,
			  int * interface_nodes, int ninterface_nodes,
			  int use_exact_degree );


#endif
