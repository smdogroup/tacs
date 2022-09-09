from copy import deepcopy
from tacs import TACS
import unittest
from collections import namedtuple
from openmdao.utils.assert_utils import (
    assert_near_equal,
    assert_check_partials,
    assert_check_totals,
)

ErrorTuple = namedtuple("ErrorTuple", ["forward", "reverse", "forward_reverse"])

"""
This is a base class for running openmdao unit test cases.
This base class will test function evaluations, partial, and total 
sensitivities for the user-specified openmdao problems implemented by
the child test case. When the user creates a new test based 
on this class two methods are required to be defined in the child class. 
    1. setup_problem
    2. setup_funcs

See the virtual method implementations for each method 
below for more details.
NOTE: The child class must NOT implement its own setUp method 
for the unittest class. This is handled in the base class.
"""


class OpenMDAOTestCase:
    class OpenMDAOTest(unittest.TestCase):
        def setUp(self):
            self.dtype = TACS.dtype

            # Default fd/cs step size and tolerances
            # Can be overridden in child class
            if self.dtype == complex:
                self.rtol = 1e-8
                self.dh = 1e-50
                self.fd_method = "cs"
                self.fd_form = None
            else:
                self.rtol = 1e-2
                self.dh = 1e-6
                self.fd_method = "fd"
                self.fd_form = "central"

            # Basically only check rtol
            self.atol = 1e99

            # Setup user-specified openmdao problem for this test
            self.prob = self.setup_problem(self.dtype)
            self.prob.setup(mode="rev", force_alloc_complex=True)

            self.func_ref, self.wrt = self.setup_funcs()

            if self.wrt is None:
                self.wrt = []

        def setup_problem(self, dtype):
            """
            Setup openmdao problem object we will be testing.
            Must be defined in child class that inherits from this class.
            """
            raise NotImplementedError(
                "Child class must implement a 'setup_problem' method"
            )
            return

        def setup_funcs(self):
            """
            Create a dictionary of functions and their reference values to be tested for the problem.
            Also provides a list of variable names to test total sensitivity wrt.
            Must be defined in child class that inherits from this class.
            """
            raise NotImplementedError(
                "Child class must implement a 'setup_funcs' method"
            )
            return

        def test_solve(self):
            """
            Test linear solve and function evaluations
            """
            # solve
            self.prob.run_model()

            # Test functions values against historical values
            for func_name in self.func_ref:
                with self.subTest(function=func_name):
                    assert_near_equal(
                        self.prob[func_name], self.func_ref[func_name], self.rtol
                    )

        def test_partials(self):
            """
            Test partial sensitivities using fd/cs
            """
            # solve
            self.prob.run_model()

            # Test functions values against historical values
            data = self.prob.check_partials(
                compact_print=True,
                out_stream=None,
                method=self.fd_method,
                form=self.fd_form,
                step=self.dh,
            )
            # Remove forward checks from data, TACS only works in rev anyways
            clean_data = self.cleanup_fwd_data(data)
            assert_check_partials(clean_data, atol=self.atol, rtol=self.rtol)

        def test_totals(self):
            """
            Test total sensitivities using fd/cs
            """
            # solve
            self.prob.run_model()

            # Test functions total sensitivities
            of = list(self.func_ref.keys())
            for var_wrt in self.wrt:
                with self.subTest(wrt=var_wrt):
                    for var_of in of:
                        with self.subTest(of=var_of):
                            data = self.prob.check_totals(
                                of=var_of,
                                wrt=var_wrt,
                                compact_print=True,
                                out_stream=None,
                                method=self.fd_method,
                                form=self.fd_form,
                                step=self.dh,
                            )
                            assert_check_totals(data, atol=self.atol, rtol=self.rtol)

        def cleanup_fwd_data(self, data):
            """
            Replace the forward sensitivity error with 0.0.
            We know this is always going to be wrong because TACS only works in reverse.
            """
            clean_data = deepcopy(data)
            for component in clean_data:
                for in_out_tuple in clean_data[component]:
                    for error_type in ["abs error", "rel error"]:
                        rev_error = clean_data[component][in_out_tuple][
                            error_type
                        ].reverse
                        clean_data[component][in_out_tuple][error_type] = ErrorTuple(
                            0.0, rev_error, 0.0
                        )
            return clean_data
