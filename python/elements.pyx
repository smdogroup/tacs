# For the use of MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy 
import numpy as np
cimport numpy as np

# Ensure that numpy is initialized
np.import_array()

# Import the definition required for const strings
from libc.string cimport const_char

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

cdef extern from "mpi-compat.h":
   pass

# This class wraps a C++ array with a numpy array for later useage
cdef class NpArrayWrap:
   cdef int nptype
   cdef int dim1, dim2
   cdef void *data_ptr

   cdef set_data1d(self, int nptype, int dim1, void *data_ptr):
      '''Set data in the array'''
      self.nptype = nptype
      self.dim1 = dim1
      self.dim2 = -1
      self.data_ptr = data_ptr
      return

   cdef set_data2d(self, int nptype, int dim1, int dim2, void *data_ptr):
      '''Set data in the array'''
      self.nptype = nptype
      self.dim1 = dim1
      self.dim2 = dim2
      self.data_ptr = data_ptr
      return

   cdef as_ndarray(self):
      '''Return a numpy version of the array'''
      # Set the shape of the array
      cdef int size = 1
      cdef np.npy_intp shape[2]
      cdef np.ndarray ndarray

      shape[0] = <np.npy_intp> self.dim1
      if (self.dim2 > 0):
         size = 2
         shape[1] = <np.npy_intp> self.dim2
      
      # Create the array itself
      ndarray = np.PyArray_SimpleNewFromData(size, shape,
                                             self.nptype, self.data_ptr)
      
      # Set the base class who owns the memory
      ndarray.base = <PyObject*>self
      Py_INCREF(self)
      
      return ndarray

class cdef pyTestElement(pyTACSObject):
    cdef TestElement *this_ptr
    def __cinit__(self, pyTACSElement element,
                  np.ndarray[TacsScalar, ndim=1, mode='c']Xpts):
        '''
        This class is used to test the element implementation

        Each function tests a single function within the Element class for
        internal self-consistency. These functions test three types of 
        quantities:

        1. Consistency between the residual calculation and the calculation
        of the Jacobian of the residuals. 

        2. Consistency between the element calculations and derivatives
        with respect to the nodal coordinates (XptSens functions.)

        3. Consistency between element calculations and the derivatives
        with respect to material design variables (DVSens functions.)

        The error in each component is measured against the global relative
        tolerance fail_rtol and the global absolute tolerance fail_atol. If
        the component in any component of the derivatives is greater than
        either of these tolerances, the function returns a non-zero fail
        flag (indicating failure.) Otherwise the function returns 0 to
        indicate that no failure has occured. It is important to note that
        the absolute tolerance check is sensitive to the order of magnitude
        of the quantities. Setting the print level to greater than 0 prints
        all components.

        Note that if no variables are supplied to the class, it generates a
        random set of variables on the interval [-1, 1]. (Such large
        displacements often produce large residuals.)
        '''
        self.this_ptr = new TestElement(element.this_ptr, <TacsScalar*>Xpts.data)

    def __dealloc__(self):
        del self.this_ptr

    def setFailTolerance(self, double fail_rtol, double fail_atol):
        '''
        Set the fail tolerance
        '''
        self.this_ptr.setFailTolerance(fail_rtol, fail_atol)

    def setPrintLevel(self, int flag):
        '''
        Set the print level
        '''
        self.this_ptr.setPrintLevel(flag)

    def setStepSize(self, TacsScalar dh):
        '''
        Set the step size
        '''
        self.this_ptr.setStepSize(dh)

    def testStiffnessMat(self, int col = NULL):
        '''
        Test the stiffness matrix using residual function
        '''
        return self.this_ptr.testStiffnessMat(col)

    def testMatDVSens(self, ElementMatrixTypes mat_type):
        '''
        Test the sensitivity of the matrix w.r.t. the design variables
        '''
        return self.this_ptr.testMatDVSens(mat_type)

    def testStrainSVSens(self, np.ndarray[double, ndim=1, mode='c']pt):
        '''
        Test the derivative of the strain with respect to the state
        variables

        addPtwiseStrainSVSens adds the derivative of a function of the
        strain with respect to the element state variables.
        '''
        return self.this_ptr.testStrainSVSens(<double*>pt.data)

    def testJacobianXptSens(self, np.ndarray[double, ndim=1, mode='c']pt):
        '''
        Test the derivative of the Jacobian w.r.t. the nodal coordinates
        '''
        return self.this_ptr.testJacobianXptSens(<double*>pt.data)

    def testStrainXptSens(self, np.ndarray[double, ndim=1, mode='c']pt):
        '''
        Test the derivative of the strain w.r.t. the nodal coordinates
        '''
        return self.this_ptr.testStrainXptSens(<double*>pt.data)

# Wrap the abstract TACS2DElement class
cdef class pyTACS2DElement(pyTACSElement)[T]:
    cdef CyTACS2DElement[T] *this_ptr
    def __init__(self pyPlaneStressStiffness planeStiff,
                 int linear = 1, int compNum = 0):
        '''
        Wrap an abstract base TACS2DElement class
        '''
        self.this_ptr = new CyTACS2DElement[T](planeStiff.this_ptr, linear, compNum)
        self.this_ptr.setSelfPointer(<void*>self)
        
    def __dealloc__(self):
        del self.this_ptr
    
# Wrap the plane stress element
cdef class pyPlaneStress(pyTACS2DElement)[T]:
    cdef PlaneStress[T] *this_ptr
    def __cinit__(self, pyPlaneStressStiffness planeStiff,
                  int linear = 1, int compNum = 0):
        '''
        Plane stress element implementation
        '''
        self.this_ptr = new PlaneStress[T](planeStiff.this_ptr, linear, compNum)
        
    def __dealloc__(self):
        del self.this_ptr

    def elementName(self):
        '''
        Return the name of this element
        '''
        return self.this_ptr.elementName()
    
# Wrap the plane stress element
cdef class pyPlaneStressTri6(pyTACS2DElement)[T]:
    cdef PlaneStressTri6[T] *this_ptr
    def __cinit__(self, pyPlaneStressStiffness planeStiff,
                  int linear = 1, int compNum = 0):
        '''
        Plane stress Tri-6 element
        '''
        self.this_ptr = new PlaneStressTri6[T](planeStiff.this_ptr, linear, compNum)

    def __dealloc__(self):
        del self.this_ptr

    def elementName(self):
        '''
        Return the name of this element
        '''
        return self.this_ptr.elementName()

# Wrap the solid element
cdef class pySolid(pyTACS3DElement)[T]:
    cdef Solid[T] *this_ptr
    def __cinit__(self, pySolidStiffness stiff,
                    int linear = 1, int compNum = 0):
        '''
        Solid element implementation
        '''
        self.this_ptr = new Solid[T](stiff.this_ptr, linear, compNum)

    def __dealloc__(self):
        del self.this_ptr

# Wrap the TACSShell element
cdef class pyTACSShell(pyTACSElement):
    cdef CyTACSShell *this_ptr
    def __init__(self, pyFSDTStiffness stiff, int compNum):
        '''
        The following class defines the base class for all Shell elements
        used in TACS. This class defines some of the operations required 
        by the generic TACSElement base class.
        '''
        self.this_ptr = new CyTACSShell(stiff.this_ptr, compNum)
        self.this_ptr.setSelfPointer(<void*>self)

    def __dealloc__(self):
        del self.this_ptr
    
# Wrap the MITCShell element
cdef class pyMITCShell(pyTACSShell)[T]:
    cdef MITCShell[T] *this_ptr
    def __cinit__(self, pyFSDTStiffness stiff, int linear = 1,
                  int compNum = 0):
        '''
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
        '''
        self.this_ptr = new MITCShell[T](stiff.this_ptr, linear, compNum)

    def __dealloc__(self):
        del self.this_ptr

    def elementName(self):
        '''
        Return the name of this element
        '''
        return self.this_ptr.elementName()
