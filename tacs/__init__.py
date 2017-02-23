'''
TACS is a parallel finite-element package for analysis and
gradient-based design optimization.
'''

import os

def get_include():
	'''
	Get the include directory for the Cython .pxd files in TACS
	'''
	return os.path.abspath(os.path.dirname(__file__))
