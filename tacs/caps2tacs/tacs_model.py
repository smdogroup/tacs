__all__ = ["TacsModel"]

import pyCAPS
from .tacs_aim import TacsAim
from .egads_aim import EgadsAim
from .analysis_function import AnalysisFunction, Derivative
from .materials import Material
from .constraints import Constraint
from .property import ShellProperty
from .loads import Load
from .variables import ShapeVariable, ThicknessVariable
from .egads_aim import EgadsAim
from typing import List
from tacs.pytacs import pyTACS


class TacsModel:
    def __init__(self, tacs_aim: TacsAim, egads_aim: EgadsAim, comm=None):
        self._tacs_aim = tacs_aim
        self._egads_aim = egads_aim
        self.comm = comm

        self._analysis_functions = []
        self.SPs = None
        self._setup = False

    @property
    def tacs_aim(self) -> TacsAim:
        return self._tacs_aim

    @property
    def egads_aim(self) -> EgadsAim:
        return self._egads_aim

    @classmethod
    def build(cls, csm_file, comm=None, problem_name: str = "capsStruct"):
        """
        make a pyCAPS problem with the tacsAIM and egadsAIM on serial / root proc

        Parameters
        ---------------------------------
        csm_file : filepath
            filename / full path of ESP/CAPS Constructive Solid Model or .CSM file
        comm : MPI.COMM
            MPI communicator
        """
        caps_problem = None
        if comm is None or comm.rank == 0:
            caps_problem = pyCAPS.Problem(
                problemName=problem_name, capsFile=csm_file, outLevel=1
            )
        tacs_aim = TacsAim(caps_problem, comm)
        egads_aim = EgadsAim(caps_problem)
        return cls(tacs_aim, egads_aim, comm)

    def get_config_parameter(self, param_name: str):
        return self.tacs_aim.get_config_parameter(param_name=param_name)

    def register(self, obj):
        """
        register each of the objects to the tacs model
        can also register any object for the tacs aim to the tacs model which passes it on
        """

        if isinstance(obj, AnalysisFunction):
            self._analysis_functions.append(obj)

        tacs_aim_objects = [
            Material,
            ThicknessVariable,
            ShapeVariable,
            ShellProperty,
            Constraint,
            Load,
            EgadsAim,
        ]
        for tacs_aim_obj in tacs_aim_objects:
            if isinstance(obj, tacs_aim_obj):
                self.tacs_aim.register(obj)
                return

        return

    def setup(self, include_aim: bool = True):
        """
        setup the analysis functions to store derivatives

        Parameters
        --------------------------------------------
        auto_tacs_aim : bool
            automatically setup the tacs aim too
        """
        # add each variable as a derivative object for each analysis function
        for func in self.analysis_functions:
            for var in self.variables:
                func._derivatives.append(Derivative(name=var.name, value=0.0))

        if include_aim:
            self.tacs_aim.setup_aim()

            # go ahead and generate the first input files and mesh for TACS
            if not self.tacs_aim.change_shape:
                self.tacs_aim.pre_analysis()

        return self

    @property
    def analysis_functions(self) -> List[AnalysisFunction]:
        """
        return the list of analysis function objects registered to the tacs aim wrapper class
        to add more functions use Function.(...).register_to(tacs_aim) or tacs_aim.register(my_analysis_function)
        """
        return self._analysis_functions

    @property
    def function_names(self) -> List[str]:
        """
        list of names of each analysis function
        """
        return [func.name for func in self.analysis_functions]

    @property
    def analysis_dir(self) -> str:
        return self.tacs_aim.analysis_dir

    @property
    def geometry(self):
        """
        link to pyCAPS geometry object to enable shape change in tacsAIM
        """
        return self.tacs_aim.geometry

    @property
    def variables(self) -> List[ShapeVariable or ThicknessVariable]:
        return self.tacs_aim.variables

    @property
    def shape_variables(self) -> List[ShapeVariable]:
        return self.tacs_aim.shape_variables

    @property
    def thickness_variables(self) -> List[ThicknessVariable]:
        return self.tacs_aim.thickness_variables

    @property
    def root_proc(self) -> bool:
        return self.comm is None or self.comm.rank == 0

    def update_design(self, input_dict: dict):
        """
        method to change the values of each design variable in tacsAim wrapper and ESP/CAPS
        """

        # track any design change to monitor capsDirty
        changed_design = False

        # change all shape variables in TacsAim and update CAD geometry
        for shape_var in self.shape_variables:
            if shape_var.name in input_dict:
                shape_var.value = float(input_dict[shape_var.name])

                # update the CAD geometry on root proc / serial since ESP/CAPS doesn't handle MPI directly
                if self.root_proc:
                    if self.geometry.despmtr[shape_var.name].value != shape_var.value:
                        self.geometry.despmtr[shape_var.name].value = shape_var.value
                        changed_design = True

        # change all thickness variables in TacsAim
        for thick_var in self.thickness_variables:
            if thick_var.name in input_dict:
                thick_var.value = float(input_dict[thick_var.name])

        # record whether the design has changed & first analysis flag as well
        if self._first_analysis:
            self._first_analysis = False
            changed_design = True

        return changed_design

    @property
    def fea_solver(self) -> pyTACS:
        """
        build pyTACS from nastran dat file and comm
        """
        return pyTACS(self.tacs_aim.dat_file_path, self.comm)

    def createTACSProbs(self, add_functions: bool = True):
        """
        creates TACS list of static, transient, or modal analysis TACS problems from the TacsAim class
        most important call method from the tacsAim class: SPs = tacs_aim.createTACSProbs
        """
        fea_solver = self.fea_solver
        fea_solver.initialize()
        SPs = fea_solver.createTACSProbsFromBDF()
        self.SPs = SPs  # store the static problems as well

        # add the analysis functions of the model into the static problems
        # add each analysis function into the static problems
        if add_functions:
            for caseID in self.SPs:
                for analysis_function in self.analysis_functions:
                    self.SPs[caseID].addFunction(
                        name=analysis_function.name,
                        handle=analysis_function.handle,
                        compIDs=analysis_function.compIDs,
                        **(analysis_function.kwargs),
                    )
        return self.SPs

    def pre_analysis(self):
        """
        call tacs aim pre_analysis to build TACS input files and mesh
        only regenerate the mesh each time if there are shape variables
        """
        if self.tacs_aim.change_shape:
            self.tacs_aim.pre_analysis()

    def run_analysis(self, write_f5: bool = True):
        """
        run the static problem analysis
        """

        # create a new set of static problems for w/ or w/o shape change
        self.SPs = self.tacs_aim.createTACSProbs(add_functions=True)

        # solve the forward and adjoint analysis for each struct problem
        self._tacs_funcs = {}
        self._tacs_sens = {}
        for caseID in self.SPs:
            # write in the new thickness variables
            xarray = self.SPs[caseID].x.getArray()
            for ithick, thick_var in enumerate(self.thickness_variables):
                xarray[ithick] = float(thick_var.value)

            self.SPs[caseID].solve()
            self.SPs[caseID].evalFunctions(
                self._tacs_funcs, evalFuncs=self.function_names
            )
            self.SPs[caseID].evalFunctionsSens(
                self._tacs_sens, evalFuncs=self.function_names
            )

            if (
                self.tacs_aim.change_shape
            ):  # if the shape changes write a sensitivity file to the tacsAim directory
                self.SPs[caseID].writeSensFile(
                    evalFuncs=self.function_names,
                    tacsAim=self.tacs_aim,
                )

            if write_f5:
                self.SPs[caseID].writeSolution(
                    baseName="tacs_output", outputDir=self.tacs_aim.analysis_dir
                )

        # return this object for method cascading
        return self

    def post_analysis(self):
        """
        call tacs aim wrapper postAnalysis and update analysis functions and gradients
        """

        if self.tacs_aim.change_shape:
            # call serial tacsAim postAnalysis if shape changes
            functions_dict = None
            gradients_dict = None

            if self.root_proc:
                self.tacs_aim.post_analysis()
                functions_dict = {}
                gradients_dict = {}
                # update functions and gradients on root proc from tacsAim dynout of ESP/CAPS serial
                for func in self.analysis_functions:
                    functions_dict[func.name] = self.tacs_aim.aim.dynout[
                        func.name
                    ].value
                    gradients_dict[func.name] = {}
                    for var in self.variables:
                        gradients_dict[func.name][var.name] = self.tacs_aim.aim.dynout[
                            func.name
                        ].deriv(var.name)

            # broadcast functions and gradients dict to all other processors from root proc
            if self.comm is not None:
                functions_dict = self.comm.bcast(functions_dict, root=0)
                gradients_dict = self.comm.bcast(gradients_dict, root=0)

            # update functions and gradients into the tacsAim analysis_functions
            for func in self.analysis_functions:
                func.value = functions_dict[func.name]
                for var in self.variables:
                    func.set_derivative(var, gradients_dict[func.name][var.name])

        # otherwise use struct problems to read in function values and gradients
        else:  # just thickness variable case
            for func in self.analysis_functions:
                # corresponding tacs key for single loadset (key="loadset#" + "func_name")
                for tacs_key in self._tacs_funcs:
                    if func.name in tacs_key:
                        break

                # add the function and gradients to each analysis function
                func.value = self._tacs_funcs[tacs_key].real

                struct_derivs = self._tacs_sens[tacs_key]["struct"]
                for ithick, thick_var in enumerate(self.thickness_variables):
                    func.set_derivative(thick_var, struct_derivs[ithick].real)

        return self