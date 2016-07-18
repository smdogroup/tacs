# Import from TACS for definitions
from TACS cimport *

# Import numpy
import numpy as np
cimport numpy as np

# Ensure that numpy is initialized
np.import_array()

cdef extern from "mpi-compat.h":
   pass

cdef extern from "TACSIntegrator.h":
   # Declare the TACSIntegrator base class
   cdef cppclass TACSIntegrator(TACSObject):      
      # Integrate forward in time
      void integrate()
      
      # Returns the adjoint gradient for all functions that are set into TACS
      void getFuncGrad(int num_dv, TacsScalar *x, TacsScalar *fvals,
                       TacsScalar *dfdx)

      # Returns the finite-difference gradient for all the functions that are set into TACS
      void getFDFuncGrad(int num_dv, TacsScalar *x, TacsScalar *fvals,
                         TacsScalar *dfdx, double dh)

      # Setters for class variables
      void setFunction(TACSFunction **func, int num_funcs)
      void setPrintLevel(int print_level)
      void setRelTol(double rtol)
      void setAbsTol(double atol)
      void setMaxNewtonIters(int max_newton_iters)
      void setJacAssemblyFreq(int jac_comp_freq)
      void setUseLapack(int use_lapack)

      # Write F5 files
      void writeSolutionToF5()

   # BDF Implementation of the integrator
   cdef cppclass TACSBDFIntegrator(TACSIntegrator):
      TACSBDFIntegrator(TACSAssembler *tacs,
                        double tinit, double tfinal,
                        int num_steps_per_sec,
                        int max_bdf_order)

   # DIRK Implementation of the integrator
   cdef cppclass TACSDIRKIntegrator(TACSIntegrator):
      TACSDIRKIntegrator(TACSAssembler *tacs,
                        double tinit, double tfinal,
                        int num_steps_per_sec,
                        int num_stages)
