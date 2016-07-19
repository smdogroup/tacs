# Import the C++ definitions
from TACS cimport *
from solver cimport *

# Import numpy 
import numpy as np
cimport numpy as np

from libc.stdlib cimport malloc, free

# Ensure that numpy is initialized
np.import_array()

# A generic abstract class for all integrators implemented in TACS
cdef class Integrator:
    '''
    Class containing functions for solving the equations forward in
    time and adjoint.
    '''    
    cdef TACSIntegrator *ptr
    
    def __cinit__(self):
        self.ptr = NULL
        return
    
    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def integrate(self):
        '''
        Integrates the governing equations forward in time
        '''
        self.ptr.integrate()
        return

    def getFuncGrad(self,
                    int num_dv,
                    np.ndarray[TacsScalar, ndim=1, mode='c'] x,
                    np.ndarray[TacsScalar, ndim=1, mode='c'] fvals,
                    np.ndarray[TacsScalar, ndim=1, mode='c'] dfdx):
        '''
        Returns the function values and derivatives with respect to
        design variables for the given array of design variables

        input:
        x: array of design variables

        output:
        fvals: array of function values
        dfdx: gradient of function with respect to design variables
        '''
        self.ptr.getFuncGrad(num_dv,
                             <TacsScalar*>x.data,
                             <TacsScalar*>fvals.data,
                             <TacsScalar*>dfdx.data)
        return

    def getFDFuncGrad(self,
                      int num_dv,
                      np.ndarray[TacsScalar, ndim=1, mode='c'] x,
                      np.ndarray[TacsScalar, ndim=1, mode='c'] fvals,
                      np.ndarray[TacsScalar, ndim=1, mode='c'] dfdx,
                      double dh):
        '''
        Returns the function values and derivatives with respect to
        design variables for the given array of design variables using
        finite differences/complex step

        input:
        x: array of design variables
        dh: finite differences/complex step perturbation step size
        
        output:
        fvals: array of function values
        dfdx: gradient of function with respect to design variables
        '''
        self.ptr.getFDFuncGrad(num_dv,
                               <TacsScalar*>x.data,
                               <TacsScalar*>fvals.data,
                               <TacsScalar*>dfdx.data, dh)
        return

    def setFunction(self, funclist):
        '''
        Sets the functions for obtaining the derivatives. This
        function should be called before seeking the gradient from
        this class.
        '''
        # Allocate the array of TACSFunction pointers
        cdef TACSFunction **funcs
        funcs = <TACSFunction**>malloc(len(funclist)*sizeof(TACSFunction*))
        if funcs is NULL:
            raise MemoryError()
        num_funcs = len(funclist)
        for i in xrange(num_funcs):
            funcs[i] = (<Function>funclist[i]).ptr
        self.ptr.setFunction(funcs, num_funcs)
        # Need to free memory but integrator needs it ()
        # free(funcs)
        return

    def setPrintLevel(self, int print_level):
        self.ptr.setPrintLevel(print_level)
        return

    def setRelTol(self, double rtol):
        self.ptr.setRelTol(rtol)
        return

    def setAbsTol(self, double atol):
        self.ptr.setAbsTol(atol)
        return

    def setMaxNewtonIters(self, int max_newton_iters):
        self.ptr.setMaxNewtonIters(max_newton_iters)
        return

    def setJacAssemblyFreq(self, int freq):
        self.ptr.setJacAssemblyFreq(freq)
        return

    def setUseLapack(self, int use_lapack):
        self.ptr.setUseLapack(use_lapack)
        return

    ## def configureOutput(self,
    ##                     ToFH5 f5,
    ##                     int write_freq = 0,
    ##                     char* file_format   = ''):
    ##     if file_format is None:
    ##         file_format = 'solution_%4d.f5'            
    ##     self.ptr.configureOutput(f5, write_freq, &file_format[0])
    ##     return

cdef class BDFIntegrator(Integrator):
    '''
    Backward-Difference method for integration. This currently
    supports upto third order accuracy in time integration.
    '''    
    def __cinit__(self, Assembler tacs,
                  double tinit, double tfinal,
                  int num_steps_per_sec,
                  int max_bdf_order):
        '''
        Constructor for BDF Integrators of order 1, 2 and 3
        '''
        self.ptr = new TACSBDFIntegrator(tacs.ptr, tinit, tfinal, num_steps_per_sec, max_bdf_order)
        self.ptr.incref()
        return

cdef class DIRKIntegrator(Integrator):
    '''
    Diagonally-Implicit-Runge-Kutta integration class. This supports
    upto fourth order accuracy in time and domain. One stage DIRK is
    second order accurate, two stage DIRK is third order accurate and
    '''    
    def __cinit__(self, Assembler tacs,
                  double tinit, double tfinal,
                  int num_steps_per_sec,
                  int num_stages):
        self.ptr = new TACSDIRKIntegrator(tacs.ptr, tinit, tfinal, num_steps_per_sec, num_stages)
        self.ptr.incref()
        return
