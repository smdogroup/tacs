#include "TACSMeshLoader.h"
#include "MITCShell.h"
#include "isoFSDTStiffness.h"


int main( int argc, char **argv ){

  MPI_Init(&argc, &argv);

  const char *filename = "CRM_box_2nd.bdf";

  MPI_Comm comm = MPI_COMM_WORLD;

  TACSMeshLoader *mesh = new TACSMeshLoader(comm);
  mesh->incref();

  mesh->scanBDFFile(filename);

  int num_components = mesh->getNumComponents();

  double rho = 2500.0; // kg/m^3
  double E = 70e9;
  double nu = 0.3;
  double kcorr = 5.0/6.0;
  double ys = 350e6;

  for ( int i = 0; i < num_components; i++ ){
    const char *descriptor = mesh->getElementDescript(i);
    double min_thickness = 0.01;
    double max_thickness = 0.20;
    double thickness = 0.07;
    isoFSDTStiffness *stiff = new isoFSDTStiffness(rho, E, nu, kcorr, ys,
        thickness, i, min_thickness, max_thickness); 

    TACSElement *element = NULL;

    if ( strcmp(descriptor, "CQUAD") == 0 || 
         strcmp(descriptor, "CQUADR") == 0 ||
         strcmp(descriptor, "CQUAD4") == 0) {
      element = new MITCShell<2>(stiff, TACSShell::LINEAR, i);
    }
    mesh->setElement(i, element);
  }

  TACSAssembler *tacs = mesh->createTACS(6);
  tacs->incref();
  mesh->decref();

  FEMat *mat = tacs->createFEMat(); 
  mat->incref();

  tacs->assembleJacobian(NULL, mat, 1, 0, 0);
  mat->applyBCs();
  PcScMat *pc = new PcScMat(mat, 10000, 10, 1); 
  pc->incref();
  BVec *ans = tacs->createVec();
  ans->incref();
  BVec *f = tacs->createVec();
  f->incref();
  f->set(1.0);
  f->applyBCs();
  pc->factor();
  pc->applyFactor(f, ans);
  tacs->setVariables(ans);
  
  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 * f5 = new TACSToFH5(tacs, SHELL, write_flag);
  f5->incref();

  // Write the displacements
  f5->writeToFile("ucrm.f5");

  // Delete the viewer
  f5->decref();

  mat->decref();
  pc->decref();
  tacs->decref();

  MPI_Finalize();

  return 0;
}
