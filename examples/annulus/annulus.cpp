#include "PlaneStressQuad.h"
#include "TACSMeshLoader.h"

/*
  The following test illustrates the use of PlaneStressQuad elements
  and demonstrates how to load a BDF file using the TACSMeshLoader
  object. This example can be run in parallel.
*/
int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  // Create the mesh loader object on MPI_COMM_WORLD. The
  // TACSAssembler object will be created on the same comm
  TACSMeshLoader *mesh = new TACSMeshLoader(MPI_COMM_WORLD);
  mesh->incref();

  // Create the plane stress stiffness object
  PlaneStressStiffness *stiff = new PlaneStressStiffness(1.0, 70e3, 0.3);
  stiff->incref();

  // The TACSAssembler object - which should be allocated if the mesh
  // is loaded correctly
  TACSAssembler *tacs = NULL;

  // Try to load the input file as a BDF file through the
  // TACSMeshLoader class
  if (argc > 1){
    const char *filename = argv[1];
    FILE *fp = fopen(filename, "r");
    if (fp){
      fclose(fp);
      
      // Scan the BDF file
      int fail = mesh->scanBDFFile(filename);

      if (fail){
	fprintf(stderr, "Failed to read in the BDF file\n");
      }
      else {
	// Add the elements to the mesh loader class
	for ( int i = 0; i < mesh->getNumComponents(); i++ ){
	  TACSElement *elem = NULL;

	  // Get the BDF description of the element
	  const char *elem_descript = mesh->getElementDescript(i);
	  if (strcmp(elem_descript, "CQUAD4") == 0){
	    elem = new PlaneStressQuad<2>(stiff);
	  }
	  else if (strcmp(elem_descript, "CQUAD") == 0 ||
		   strcmp(elem_descript, "CQUAD9") == 0){
	    elem = new PlaneStressQuad<3>(stiff);
	  }
	  else if (strcmp(elem_descript, "CQUAD16") == 0){
	    elem = new PlaneStressQuad<4>(stiff);
	  }

	  // Set the element object into the mesh loader class
	  if (elem){
	    mesh->setElement(i, elem);
	  }
	}

	// Now, create the TACSAssembler object
        int vars_per_node = 2;
	tacs = mesh->createTACS(vars_per_node);
	tacs->incref();
      }
    }
    else {
      fprintf(stderr, "File %s does not exist\n", filename);
    }
  }
  else {
    fprintf(stderr, "No BDF file provided\n");
  }

  if (tacs){
    // Create the preconditioner
    TACSBVec *res = tacs->createVec();
    TACSBVec *ans = tacs->createVec();
    FEMat *mat = tacs->createFEMat();
  
    // Increment the reference count to the matrix/vectors
    res->incref();
    ans->incref();
    mat->incref();
    
    // Allocate the factorization
    int lev = 4500;
    double fill = 10.0;
    int reorder_schur = 1;
    PcScMat *pc = new PcScMat(mat, lev, fill, reorder_schur);
    pc->incref();

    // Allocate the GMRES object
    int gmres_iters = 80; 
    int nrestart = 2; // Number of allowed restarts
    int is_flexible = 0; // Is a flexible preconditioner?
    TACSKsm *ksm = new GMRES(mat, pc, gmres_iters, 
                             nrestart, is_flexible);
    ksm->incref();

    // Assemble and factor the stiffness/Jacobian matrix
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    tacs->assembleJacobian(alpha, beta, gamma, res, mat);
    pc->factor();
    
    res->set(1.0);
    tacs->applyBCs(res);
    ksm->solve(res, ans);
    tacs->setVariables(ans);

    // Create an TACSToFH5 object for writing output to files
    unsigned int write_flag = (TACSElement::OUTPUT_NODES |
			       TACSElement::OUTPUT_DISPLACEMENTS |
			       TACSElement::OUTPUT_STRAINS |
			       TACSElement::OUTPUT_STRESSES |
			       TACSElement::OUTPUT_EXTRAS);
    TACSToFH5 * f5 = new TACSToFH5(tacs, TACS_PLANE_STRESS, write_flag);
    f5->incref();
    f5->writeToFile("output.f5");
    
    // Free everything
    f5->decref();
    
    // Decrease the reference count to the linear algebra objects
    ksm->decref();
    pc->decref();
    mat->decref();
    ans->decref();
    res->decref();
  }

  // Deallocate the objects
  mesh->decref();
  stiff->decref();
  if (tacs){ tacs->decref(); }

  MPI_Finalize();
  return (0);
}
