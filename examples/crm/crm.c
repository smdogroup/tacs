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

  pc->factor();

  mat->decref();
  pc->decref();
  tacs->decref();

  MPI_Finalize();

  return 0;
}
