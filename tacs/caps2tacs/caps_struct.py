
"""
Sean Engelstad
Georgia Tech SMDO
08/18/2022
"""

__all__ = ["CapsProblem", "CapsStruct", "CapsFluid"]

import pyCAPS
from typing import TYPE_CHECKING

# import each of the aim modules here
from .tacs_aim import TacsAim
from .egads_aim import EgadsAim

class CapsStruct:
    """
    Base class for Structure problems with ESP/CAPS
    Often uses TACS for structure solver
    """
    def __init__(self, problem:pyCAPS.Problem, comm=None):
        self.comm = comm
        self._problem = problem


    @classmethod
    def build(cls, csm_file:str, problem_name:str="capsStruct"):
        """
        auto build a caps struct problem
        syntax: CapsStruct.build(csm)
        """
        problem = pyCAPS.Problem(problemName=problem_name, capsFile=csm_file, outLevel=1)
        return cls(problem)

    @property
    def geometry(self):
        return self._problem.geometry

    @property
    def view(self):
        self.geometry.view()

    @property
    def design_parameters(self):
        return self.geometry.despmtr.keys()

    @property
    def tacsAim(self) -> TacsAim:
        return TacsAim(caps_problem=self._problem)

    @property
    def egadsAim(self) -> EgadsAim:
        return EgadsAim(caps_problem=self._problem)

    @property
    def design_parameters(self):
        return self.geometry.despmtr.keys()
