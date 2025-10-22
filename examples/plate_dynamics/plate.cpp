#include "Compliance.h"
#include "InducedFailure.h"
#include "KSFailure.h"
#include "MITC9.h"
#include "StructuralMass.h"
#include "TACSIntegrator.h"
#include "TACSMeshLoader.h"
#include "TACSShellTraction.h"
#include "isoFSDTStiffness.h"

/*
  Code for testing adjoints with plate example. Use command line
  arguments as necessary.

  BDF1 BDF2 BDF3    : for BDF integrators
  DIRK2 DIRK3 DIRK4 : for DIRK integrators
  ABM1-6            : for ABM integrators
  NBG               : for Newmark integrator

  num_funcs     : 1 to 3 for the adjoint
  num_threads   : number of threads
  write_solution: write solution to f5 frequency
  print_level: 0, 1, 2
*/
const char *help_string[] = {
    "TACS time-dependent analysis of a plate located in plate.bdf",
    "num_funcs=1,2,3 and 12  : Number of functions for adjoint problem",
    "num_threads=1,2,3...    : Number of threads to use",
    "print_level=0,1,2       : Controls the amount of information to print",
    "write_solution=0,1,2... : Controls the frequency of f5 file output",
    "convert_mesh=0,1        : Converts the mesh to coordinate ordering"};

int main(int argc, char **argv) {
  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Parse command line arguments
  int num_funcs = 1;
  int num_threads = 1;
  int write_solution = 0;
  int print_level = 1;
  int convert_mesh = 0;
  double dh = 1e-7;

  for (int i = 0; i < argc; i++) {
    // Determine whether or not to test gradients with complex step
    if (strcmp("--help", argv[i]) == 0) {
      if (rank == 0) {
        for (int k = 0; k < 6; k++) {
          printf("%s\n", help_string[k]);
        }
      }
      MPI_Finalize();
      return 0;
    }

    // Determine the complex/finite-difference step interval
    if (sscanf(argv[i], "dh=%lf", &dh) == 1) {
      if (rank == 0) {
        printf("Difference interval : %g\n", dh);
      }
    }

    // Determine the number of functions for adjoint
    if (sscanf(argv[i], "num_funcs=%d", &num_funcs) == 1) {
      if (num_funcs < 0) {
        num_funcs = 1;
      }
      if (rank == 0) {
        printf("Number of functions : %d\n", num_funcs);
      }
    }

    // How frequent to write the f5 files
    if (sscanf(argv[i], "write_solution=%d", &write_solution) == 1) {
      if (write_solution < 0) {
        write_solution = 0;
      }
      if (rank == 0) {
        printf("Write solution freq : %d\n", write_solution);
      }
    }

    // Set the print level
    if (sscanf(argv[i], "print_level=%d", &print_level) == 1) {
      if (print_level < 0) {
        print_level = 1;
      }
      if (print_level > 3) {
        print_level = 3;
      }
      if (rank == 0) {
        printf("Print level : %d\n", print_level);
      }
    }

    // Determine the number of threads
    if (sscanf(argv[i], "num_threads=%d", &num_threads) == 1) {
      if (num_threads < 0) {
        num_threads = 1;
      }
      if (num_threads > 24) {
        num_threads = 24;
      }
      if (rank == 0) {
        printf("Number of threads : %d\n", num_threads);
      }
    }

    // Determine whether or not to convert the
    // connectivities to coordinate ordering
    if (sscanf(argv[i], "convert_mesh=%d", &convert_mesh) == 1) {
      if (convert_mesh != 1) {
        convert_mesh = 0;
      } else {
        convert_mesh = 1;
      }
      if (rank == 0) {
        printf("Convert mesh to coordinate order : %d\n", convert_mesh);
      }
    }
  }

  // Write name of BDF file to be load to char array
  const char *filename = "plate.bdf";

  // Create the mesh loader object and load file
  TACSMeshLoader *mesh = new TACSMeshLoader(comm);
  mesh->incref();
  mesh->setConvertToCoordinate(convert_mesh);
  mesh->scanBDFFile(filename);

  // Get number of components prescribed in BDF file
  int num_components = mesh->getNumComponents();

  // Set properties for structural elements
  double rho = 2500.0;       // density, kg/m^3
  double E = 70e9;           // elastic modulus, Pa
  double nu = 0.3;           // poisson's ratio
  double kcorr = 5.0 / 6.0;  // shear correction factor
  double ys = 350e6;         // yield stress, Pa

  // Set properties for dynamics
  TacsScalar g[] = {0.0, 0.0, -9.81};
  TacsScalar v_init[] = {0.0, 0.0, 0.0};
  TacsScalar omega_init[] = {0.0, 0.0, 0.0};

  /* TacsScalar v_init[] = {0.1, 0.1, 0.1};  */
  /* TacsScalar omega_init[] = {0.3, 0.1, 0.2}; */

  TACSGibbsVector *gravity = new TACSGibbsVector(g);
  gravity->incref();
  TACSGibbsVector *v0 = new TACSGibbsVector(v_init);
  v0->incref();
  TACSGibbsVector *omega0 = new TACSGibbsVector(omega_init);
  omega0->incref();

  // Variables per node
  int vars_per_node = 0;

  // Loop over components, creating constituitive object for each
  for (int i = 0; i < num_components; i++) {
    const char *descriptor = mesh->getElementDescript(i);
    double min_thickness = 5.0e-3;
    double max_thickness = 2.0e-2;
    double thickness = 5.0e-3;
    isoFSDTStiffness *stiff = new isoFSDTStiffness(
        rho, E, nu, kcorr, ys, thickness, i, min_thickness, max_thickness);
    stiff->incref();

    TacsScalar axis[] = {1.0, 0.0, 0.0};
    stiff->setRefAxis(axis);

    // Initialize element object
    TACSElement *element = NULL;

    // Create element object using constituitive information and type
    // defined in the descriptor
    if (strcmp(descriptor, "CQUAD9") == 0 || strcmp(descriptor, "CQUAD") == 0) {
      element = new MITC9(stiff, gravity, v0, omega0);
      element->incref();

      // Set the number of displacements
      vars_per_node = element->numDisplacements();
      mesh->setElement(i, element);

      stiff->decref();
      element->decref();
    } else {
      printf("[%d] TACS Warning: Unsupported element %s in BDF file\n", rank,
             descriptor);
    }
  }

  // Create tacs assembler from mesh loader object
  TACSAssembler *tacs = mesh->createTACS(vars_per_node);
  tacs->incref();

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag =
      (TACSElement::OUTPUT_NODES | TACSElement::OUTPUT_DISPLACEMENTS |
       TACSElement::OUTPUT_STRAINS | TACSElement::OUTPUT_STRESSES |
       TACSElement::OUTPUT_EXTRAS);

  TACSToFH5 *f5 = new TACSToFH5(tacs, TACS_SHELL, write_flag);
  f5->incref();

  /*-----------------------------------------------------------------*/
  /*------------------ Time Integration and Adjoint Solve -----------*/
  /*-----------------------------------------------------------------*/

  int num_dvs = num_components;

  // Create functions of interest
  TACSFunction *func[num_funcs];
  if (num_funcs == 1) {
    func[0] = new TACSCompliance(tacs);
  } else if (num_funcs == 2) {
    func[0] = new TACSKSFailure(tacs, 100.0);

    // Set the induced norm failure types
    TACSInducedFailure *ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::EXPONENTIAL);
    func[1] = ifunc;
  } else if (num_funcs == 3) {
    func[0] = new TACSKSFailure(tacs, 100.0);
    func[1] = new TACSCompliance(tacs);
    func[2] = new TACSStructuralMass(tacs);
  } else if (num_funcs == 12) {
    // Place functions into the func list
    func[0] = new TACSStructuralMass(tacs);
    func[1] = new TACSCompliance(tacs);

    // Set the discrete and continuous KS functions
    TACSKSFailure *ksfunc = new TACSKSFailure(tacs, 20.0);
    ksfunc->setKSFailureType(TACSKSFailure::DISCRETE);
    func[2] = ksfunc;

    ksfunc = new TACSKSFailure(tacs, 20.0);
    ksfunc->setKSFailureType(TACSKSFailure::CONTINUOUS);
    func[3] = ksfunc;

    // Set the induced norm failure types
    TACSInducedFailure *ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::EXPONENTIAL);
    func[4] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::DISCRETE_EXPONENTIAL);
    func[5] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::DISCRETE_EXPONENTIAL_SQUARED);
    func[6] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::EXPONENTIAL_SQUARED);
    func[7] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::POWER);
    func[8] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::DISCRETE_POWER);
    func[9] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::POWER_SQUARED);
    func[10] = ifunc;

    ifunc = new TACSInducedFailure(tacs, 20.0);
    ifunc->setInducedType(TACSInducedFailure::DISCRETE_POWER_SQUARED);
    func[11] = ifunc;
  }

  for (int i = 0; i < num_funcs; i++) {
    func[i]->incref();
  }

  TacsScalar *funcVals = new TacsScalar[num_funcs];     // adjoint
  TacsScalar *funcValsTmp = new TacsScalar[num_funcs];  // CSD
  TacsScalar *funcVals1 = new TacsScalar[num_funcs];    // forward/reverse

  TacsScalar *dfdx = new TacsScalar[num_funcs * num_dvs];     // adjoint
  TacsScalar *dfdx1 = new TacsScalar[num_funcs * num_dvs];    // CSD
  TacsScalar *dfdxTmp = new TacsScalar[num_funcs * num_dvs];  // forward/reverse

  // Create an array of design variables
  TacsScalar *x = new TacsScalar[num_dvs];
  x[0] = 0.02;

  // Set paramters for time marching
  double tinit = 0.0;
  double tfinal = 0.005;
  for (int k = 0; k < argc; k++) {
    if (sscanf(argv[k], "tfinal=%lf", &tfinal) == 1) {
    }
  }
  int num_steps_per_sec = 1000;

  int start_plane = 2;
  int end_plane = 5;

  // Check the BDF integrator
  int bdf_order = 3;
  TACSIntegrator *bdf =
      new TACSBDFIntegrator(tacs, tinit, tfinal, num_steps_per_sec, bdf_order);
  bdf->incref();

  // Set options
  bdf->setPrintLevel(print_level);
  bdf->setOutputFrequency(write_solution);

  // Set functions of interest for adjoint solve
  bdf->setFunctions(func, num_funcs, num_dvs, start_plane, end_plane);
  bdf->checkGradients(dh);
  bdf->decref();

  // Check the DIRK integrator
  int num_stages = 2;
  TACSIntegrator *dirk = new TACSDIRKIntegrator(tacs, tinit, tfinal,
                                                num_steps_per_sec, num_stages);
  dirk->incref();

  // Set options
  dirk->setPrintLevel(print_level);
  dirk->setOutputFrequency(write_solution);

  // Set functions of interest for adjoint solve
  dirk->setFunctions(func, num_funcs, num_dvs, start_plane, end_plane);
  dirk->checkGradients(dh);
  dirk->decref();

  mesh->decref();
  gravity->decref();
  v0->decref();
  omega0->decref();

  for (int i = 0; i < num_funcs; i++) {
    func[i]->decref();
  }

  tacs->decref();

  delete[] dfdx;
  delete[] dfdx1;
  delete[] dfdxTmp;

  delete[] x;

  delete[] funcVals;
  delete[] funcVals1;
  delete[] funcValsTmp;

  MPI_Finalize();
  return 0;
}
