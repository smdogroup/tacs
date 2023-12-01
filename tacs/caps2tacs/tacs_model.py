"""
Written by Sean Engelstad, GT SMDO Lab, 2022-2023
"""
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
from .aflr_aim import AflrAim
from typing import List
from tacs.pytacs import pyTACS

# optional funtofem dependency for shape variable animations
# will still work without funtofem
import importlib

f2f_loader = importlib.util.find_spec("funtofem")


class TacsModel:
    MESH_AIMS = ["egads", "aflr"]

    def __init__(self, tacs_aim: TacsAim, mesh_aim, comm=None):
        self._tacs_aim = tacs_aim
        self._mesh_aim = mesh_aim
        self.comm = comm

        self._analysis_functions = []
        self.SPs = None
        self._setup = False
        self._first_analysis = True

    @property
    def tacs_aim(self) -> TacsAim:
        return self._tacs_aim

    @property
    def mesh_aim(self) -> AflrAim:
        return self._mesh_aim

    @property
    def uses_egads(self):
        return isinstance(self.mesh_aim, EgadsAim)

    @property
    def uses_aflr(self):
        return isinstance(self.mesh_aim, AflrAim)

    @classmethod
    def build(
        cls,
        csm_file,
        comm=None,
        mesh="egads",
        tacs_project="tacs",
        active_procs: list = [0],
        problem_name: str = "capsStruct",
        mesh_morph: bool = False,
        verbosity: int = 0,
    ):
        """
        make a pyCAPS problem with the tacsAIM and egadsAIM on serial / root proc

        Parameters
        ---------------------------------
        csm_file : filepath
            filename / full path of ESP/CAPS Constructive Solid Model or .CSM file
        comm : MPI.COMM
            MPI communicator
        mesh: str
            whether to use "egads" or "aflr" meshing tool
        tacs_project: str
            project name for all ESP/CAPS work directory files (affects filenames)
        active_procs: list
            list of active procs for parallel tacsAIM instances. Ex: if you want to use 5
            parallel tacsAIMs you can choose active_procs=[_ for _ in range(5)] or
            equivalently active_procs=[0,1,2,3,4]
        problem_name: str
            the folder name for the ESP/CAPS work directory (rename this to avoid conflicts)
        mesh_morph: bool
            whether to engage structural mesh morphing (usually not necessary)
        verbosity: int
            0 to suppress print statements to terminal, 1 to allow them.
        """

        assert 0 in active_procs
        caps_problem = None
        assert mesh in cls.MESH_AIMS
        assert max(active_procs) + 1 <= comm.Get_size()

        for iproc in active_procs:
            if comm.rank == iproc:
                caps_problem = pyCAPS.Problem(
                    problemName=problem_name + "_" + str(iproc),
                    capsFile=csm_file,
                    outLevel=verbosity,
                )

        tacs_aim = TacsAim(
            caps_problem,
            comm,
            project_name=tacs_project,
            mesh_morph=mesh_morph,
            active_procs=active_procs,
        )
        mesh_aim = None
        if mesh == "egads":
            mesh_aim = EgadsAim(caps_problem, comm, active_procs=active_procs)
        elif mesh == "aflr":
            mesh_aim = AflrAim(caps_problem, comm, active_procs=active_procs)
        return cls(tacs_aim, mesh_aim, comm)

    @property
    def active_procs(self) -> list:
        return self.tacs_aim.active_procs

    def get_config_parameter(self, param_name: str):
        return self.tacs_aim.get_config_parameter(param_name=param_name)

    def get_output_parameter(self, param_name: str):
        return self.tacs_aim.get_output_parameter(param_name=param_name)

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
            AflrAim,
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

            # Set additional options for meshing AIM through dictionaries
            if self.mesh_aim._dictOptions is not None:
                self.mesh_aim._set_dict_options()

            # go ahead and generate the first input files and mesh for TACS
            if not self.tacs_aim.change_shape:
                self.tacs_aim.pre_analysis()

        self._setup = True

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

    def analysis_dir(self, proc: int = 0) -> str:
        return self.tacs_aim.analysis_dir(proc)

    @property
    def mesh_morph(self) -> bool:
        return self.tacs_aim.mesh_morph

    @mesh_morph.setter
    def mesh_morph(self, new_bool):
        self.tacs_aim.mesh_morph = new_bool

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
    def variable_dict(self) -> dict:
        return {var.name: var.value for var in self.variables}

    @property
    def shape_variables(self) -> List[ShapeVariable]:
        return self.tacs_aim.shape_variables

    @property
    def thickness_variables(self) -> List[ThicknessVariable]:
        return self.tacs_aim.thickness_variables

    @property
    def root_proc_ind(self) -> int:
        return self.active_procs[0]

    @property
    def root_proc(self) -> bool:
        return self.comm.rank == self.root_proc_ind

    def update_design(self, input_dict: dict = None):
        """
        method to change the values of each design variable in tacsAim wrapper and ESP/CAPS
        """

        input_dict = input_dict if input_dict is not None else self.variable_dict

        # track any design change to monitor capsDirty
        changed_design = False

        # change all shape variables in TacsAim and update CAD geometry
        for shape_var in self.shape_variables:
            if shape_var.name in input_dict:
                if input_dict[shape_var.name] is not None:
                    shape_var.value = float(input_dict[shape_var.name])

                # update the CAD geometry on root proc / serial since ESP/CAPS doesn't handle MPI directly
                for proc in self.active_procs:
                    # apply geometry changes to all procs even if that shape variable is differentiated on that processor
                    if self.comm.rank == proc:
                        if (
                            self.geometry.despmtr[shape_var.name].value
                            != shape_var.value
                        ):
                            changed_design = True
                            if shape_var.value is not None:
                                self.geometry.despmtr[
                                    shape_var.name
                                ].value = shape_var.value
                            else:
                                shape_var.value = self.geometry.despmtr[
                                    shape_var.name
                                ].value

        # change all thickness variables in TacsAim
        for thick_var in self.thickness_variables:
            if thick_var.name in input_dict:
                if thick_var.value != float(input_dict[thick_var.name]):
                    thick_var.value = float(input_dict[thick_var.name])
                    changed_design = True

        # update thickness prop cards in t
        if self.tacs_aim.change_shape:
            self.tacs_aim.update_properties()

        # record whether the design has changed & first analysis flag as well
        if self._first_analysis:
            self._first_analysis = False
            changed_design = True

        return changed_design

    def fea_solver(self, root=0) -> pyTACS:
        """
        build pyTACS from nastran dat file on the root proc
        """
        return pyTACS(self.tacs_aim.dat_file_path(root), self.comm)

    def createTACSProbs(self, root=0, addFunctions: bool = True):
        """
        creates TACS list of static, transient, or modal analysis TACS problems from the TacsAim class
        most important call method from the tacsAim class: SPs = tacs_aim.createTACSProbs
        """
        fea_solver = self.fea_solver(root)
        fea_solver.initialize()
        SPs = fea_solver.createTACSProbsFromBDF()
        self.SPs = SPs  # store the static problems as well

        # add the analysis functions of the model into the static problems
        # add each analysis function into the static problems
        if addFunctions:
            for caseID in self.SPs:
                for analysis_function in self.analysis_functions:
                    self.SPs[caseID].addFunction(
                        funcName=analysis_function.name,
                        funcHandle=analysis_function.handle,
                        compIDs=analysis_function.compIDs,
                        **(analysis_function.kwargs),
                    )
        return self.SPs

    def _convert_shape_vars_dict(self, shape_vars_dict: dict):
        """convert the shape variable dict keys from str or funtofem variables to Tacs Shape Variables if necessary"""
        # convert shape variable dict to caps2tacs ShapeVariable keys if necessary
        shape_vars_dict2 = {}
        first_key = list(shape_vars_dict.keys())[0]
        if isinstance(first_key, ShapeVariable):
            shape_vars_dict2 = shape_vars_dict
        elif isinstance(first_key, str):
            # get matching caps2tacs shape variable for each string
            for key in shape_vars_dict:
                for var in self.tacs_aim.shape_variables:
                    if var.name == key:
                        shape_vars_dict2[var] = shape_vars_dict[key]
                        break  # go to next key
        elif f2f_loader is not None:  # optional funtofem dependency
            # import has to go here to prevent circular import
            from funtofem import Variable

            if isinstance(first_key, Variable):
                # get matching caps2tacs shape variable for each F2F shape variable
                for f2f_var in shape_vars_dict:
                    for tacs_var in self.tacs_aim.shape_variables:
                        if f2f_var.name == tacs_var.name:
                            # copy the value list to the new dictionary
                            shape_vars_dict2[tacs_var] = shape_vars_dict[f2f_var]
                            break  # go to next variable

            else:
                raise AssertionError(
                    f"The datatype {type(first_key)} is not allowed in caps2tacs shape variable dictionaries."
                )
        else:
            raise AssertionError(
                f"The datatype {type(first_key)} is not allowed in caps2tacs shape variable dictionaries."
            )

        return shape_vars_dict2

    def animate_shape_vars(self, shape_vars_dict: dict):
        """
        Animate the shape variables in ESP/CAPS.

        Parameters
        ----------
        shape_vars_dict: dict[ShapeVariable : list[float]]
            dict of the list of values for each shape variable to take

        e.g. if you wish to animate over the ShapeVariable objects var1, and var2
        create the following shape_vars_dict
        NOTE : the keys can be str, ShapeVariable or funtofem variable objects objects
        shape_vars_dict = {
          var1 : [_*0.1 for _ in range(1,4)],
          var2 : [_*0.05 for _ in range(-3,4)],
        }
        """
        # make an analysis function if none have been made
        if len(self.analysis_functions) == 0:
            if self.comm.rank == 0:
                print(
                    "Adding mass analysis function to enable caps2tacs postAnalysis for animating shape variables..."
                )
            AnalysisFunction.mass().register_to(self)

        # set all shape variables inactive
        for shape_var in self.tacs_aim.shape_variables:
            shape_var._active = False

        # convert the shape variables dict keys from str or F2F variable to Tacs Shape Variables
        shape_vars_dict = self._convert_shape_vars_dict(shape_vars_dict)

        for shape_var in shape_vars_dict:
            value_list = shape_vars_dict[shape_var]

            # want only this shape variable to be active so that it doesn't
            # try to do extra mesh sensitivity chain rule products in tacsAIM
            shape_var._active = True
            self.tacs_aim.setup_aim()

            # save the original value
            orig_value = shape_var.value * 1.0

            for i, value in enumerate(value_list):
                shape_var.value = value
                # update geometry on all active procs
                for proc in self.active_procs:
                    if self.comm.rank == proc:
                        self.geometry.despmtr[shape_var.name].value = value

                self.tacs_aim.pre_analysis()
                self.SPs = self.createTACSProbs(addFunctions=True)
                for caseID in self.SPs:
                    # note we don't solve the forward / adjoint analysis here
                    # as we only care about visualizing the change in the structural shape / mesh
                    # not the actual structural analysis results

                    self.SPs[caseID].writeSolution(
                        baseName=shape_var.name,
                        outputDir=self.tacs_aim.analysis_dir(self.root_proc_ind),
                        number=i,
                    )
                    self.SPs[caseID].writeSensFile(
                        evalFuncs=None, tacsAim=self.tacs_aim, root=0, proc=0
                    )

                del self.SPs
                self.tacs_aim.post_analysis()

            # return to the original value
            for proc in self.active_procs:
                if self.comm.rank == proc:
                    self.geometry.despmtr[shape_var.name].value = orig_value
            shape_var.value = orig_value

            # make the shape variable inactive again
            shape_var._active = False
        print(
            f"Done animating caps2tacs shape variables.. the f5 files are found in {self.tacs_aim.analysis_dir}."
        )
        print("Use the GifWriter script in examples/caps_wing/archive to ")

    def pre_analysis(self, mpi_barrier: bool = True):
        """
        call tacs aim pre_analysis to build TACS input files and mesh
        only regenerate the mesh each time if there are shape variables
        """
        if self.tacs_aim.change_shape:
            self.tacs_aim.pre_analysis()

        # other procs wait on the active procs to finish
        if mpi_barrier:
            self.comm.Barrier()
        return self

    def run_analysis(
        self, write_f5: bool = True, iteration: float = 0, mpi_barrier: bool = True
    ):
        """
        run the static problem analysis
        """

        assert self._setup

        # create a new set of static problems for w/ or w/o shape change
        self.SPs = self.createTACSProbs(addFunctions=True)

        # solve the forward and adjoint analysis for each struct problem
        self._tacs_funcs = {}
        self._tacs_sens = {}
        for caseID in self.SPs:
            # write in the new thickness variables for sizing only case
            if not self.tacs_aim.change_shape:
                if self.comm.rank == 0:  # only change TACS design vars on root proc
                    xarray = self.SPs[caseID].x.getArray()
                    for ithick, thick_var in enumerate(self.thickness_variables):
                        xarray[ithick] = float(thick_var.value)

            self.SPs[caseID].solve()

            if (
                self.tacs_aim.change_shape
            ):  # if the shape changes write a sensitivity file to the tacsAim directory
                for proc in self.active_procs:
                    # write a sens file to each of the capsStruct directories
                    self.SPs[caseID].writeSensFile(
                        evalFuncs=self.function_names,
                        tacsAim=self.tacs_aim,
                        root=0,
                        proc=proc,
                    )
            else:  # only call evalFunctions and evalFunctionSens if not shape change else redundant
                self.SPs[caseID].evalFunctions(
                    self._tacs_funcs, evalFuncs=self.function_names
                )
                self.SPs[caseID].evalFunctionsSens(
                    self._tacs_sens, evalFuncs=self.function_names
                )

                if self.comm.rank != 0:
                    self._tacs_funcs = None
                    self._tacs_sens = None

                # broadcast from root proc to other procs the sens
                #    since some derivatives are missing here
                # don't care about the Xpts sens which is distributed since
                #    this is the case with no shape variables
                self._tacs_funcs = self.comm.bcast(self._tacs_funcs, root=0)
                self._tacs_sens = self.comm.bcast(self._tacs_sens, root=0)

            if write_f5:
                self.SPs[caseID].writeSolution(
                    baseName="tacs_output",
                    outputDir=self.tacs_aim.analysis_dir(self.root_proc_ind),
                    number=iteration,
                )

        # other procs wait on the active procs to finish
        if mpi_barrier:
            self.comm.Barrier()

        # return this object for method cascading
        return self

    def post_analysis(self, mpi_barrier=True):
        """
        call tacs aim wrapper postAnalysis and update analysis functions and gradients
        """

        if self.tacs_aim.change_shape:
            self.tacs_aim.post_analysis()

            functions_dict = {}
            gradients_dict = {}
            # get shape variable derivatives from each processor
            for func in self.analysis_functions:
                func_value = None
                if self.comm.rank == self.root_proc_ind:
                    func_value = self.tacs_aim.aim.dynout[func.name].value

                # broadcast and save function value on other procs
                func_value = self.comm.bcast(func_value, root=self.root_proc_ind)
                functions_dict[func.name] = func_value
                gradients_dict[func.name] = {}

            for var in self.variables:
                # get derivative from one of the processors that has it and broadcast
                c_deriv_dict = (
                    {}
                )  # key is func name, holds derivatives for current variable
                if isinstance(var, ThicknessVariable):
                    c_root = self.root_proc_ind
                elif isinstance(var, ShapeVariable):
                    c_root = self.tacs_aim.get_proc_with_shape_var(var)
                else:
                    raise AssertionError(
                        f"caps2tacs - var {type(var)} is not a thickness or shape variable.."
                    )

                for func in self.analysis_functions:
                    if self.comm.rank == c_root:
                        c_deriv_dict[func.name] = self.tacs_aim.aim.dynout[
                            func.name
                        ].deriv(var.name)

                # broadcast c_deriv_dict to the other processors
                c_deriv_dict = self.comm.bcast(c_deriv_dict, root=c_root)

                # update the gradients dict
                for func in self.analysis_functions:
                    gradients_dict[func.name][var.name] = c_deriv_dict[func.name]

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

        if self.mesh_morph:
            self.tacs_aim.unlink()

        # other procs wait on the active procs to finish
        if mpi_barrier:
            self.comm.Barrier()

        return self
