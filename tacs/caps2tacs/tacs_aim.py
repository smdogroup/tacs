__all__ = ["TacsAim"]

from typing import TYPE_CHECKING, List
import os, numpy as np
from .proc_decorator import root_proc, root_broadcast
from .materials import Material
from .constraints import Constraint
from .property import ShellProperty, Property
from .loads import Load
from .variables import ShapeVariable, ThicknessVariable
from .egads_aim import EgadsAim
from .aflr_aim import AflrAim


class TacsAimMetadata:
    def __init__(self, analysis_dir, project_name):
        self.analysis_dir = analysis_dir
        self.project_name = project_name


class TacsAim:
    """
    Wrapper class for TacsAim with default build setting in different scenarios
    applies default settings and spits it back out at the end
    only supports shell properties at the moment
    """

    def __init__(self, caps_problem, comm=None):
        self.comm = comm

        # geometry and design parameters to change the design of the CSM file during an optimization
        self._aim = None
        self._geometry = None
        self._build_aim(caps_problem)

        self._analysis_functions = []
        self._materials = []
        self._loads = []
        self._properties = []
        self._constraints = []
        self._design_variables = []
        self._mesh_aim = None

        # build flags
        self._setup = False
        self._first_setup = True

        # broadcast TacsAimMetadata from root proc to other processors
        self._metadata = None
        self._broadcast_metadata()

    @root_proc
    def _build_aim(self, caps_problem):
        """
        build the TacsAim pyCAPS object inside our wrapper class on root proc
        """
        self._aim = caps_problem.analysis.create(aim="tacsAIM", name="tacs")
        self._geometry = caps_problem.geometry

    def _broadcast_metadata(self):
        """
        broadcast any tacs aim metadata needed for this class from root proc to other processors
        """
        if self.comm is None:
            self._metadata = TacsAimMetadata(
                analysis_dir=self._aim.analysisDir,
                project_name=self._aim.input.Proj_Name,
            )
        else:
            if self.comm.rank == 0:
                self._metadata = TacsAimMetadata(
                    analysis_dir=self._aim.analysisDir,
                    project_name=self._aim.input.Proj_Name,
                )
            self._metadata = self.comm.bcast(self._metadata, root=0)

    def register(self, obj):
        """
        register any one of the available objects: Materials, Properties, Variables, etc.
        """
        if isinstance(obj, Material):
            self._materials.append(obj)
        elif isinstance(obj, ThicknessVariable):
            self._design_variables.append(obj)
            if obj.can_make_shell:
                self._properties.append(obj.shell_property)
        elif isinstance(obj, ShapeVariable):
            self._design_variables.append(obj)
        elif isinstance(obj, Property):
            self._properties.append(obj)
        elif isinstance(obj, Constraint):
            self._constraints.append(obj)
        elif isinstance(obj, Load):
            self._loads.append(obj)
        elif isinstance(obj, EgadsAim):
            self._mesh_aim = obj
        elif isinstance(obj, AflrAim):
            self._mesh_aim = obj
        else:
            raise AssertionError(
                "Object could not be registered to TacsAim as it is not an appropriate type."
            )

    def setup_aim(
        self,
        large_format: bool = True,
        static: bool = True,
        auto_shape_variables: bool = False,
    ):
        # make sure there is at least one material, property, constraint, etc.
        assert len(self._materials) > 0
        assert len(self._properties) > 0
        assert len(self._constraints) > 0
        assert self._mesh_aim is not None

        # this part runs on serial
        if self.comm is None or self.comm.rank == 0:
            # write in the original shape variable values into the pyCAPS geometry
            for shape_var in self.shape_variables:
                if shape_var.value is not None:
                    self.geometry.despmtr[shape_var.name].value = shape_var.value

            # increase the precision in the BDF file
            self.aim.input.File_Format = "Large" if large_format else "Small"
            self.aim.input.Mesh_File_Format = "Large" if large_format else "Small"

            # set the analysis type
            if static:
                self.aim.input.Analysis_Type = "Static"
            else:
                raise AssertionError(
                    "Analysis types other than static analyses for tacsAim are not supported yet."
                )

            # add materials to tacsAim
            self.aim.input.Material = {
                material.name: material.dictionary for material in self._materials
            }

            # add properties to tacsAim
            self.aim.input.Property = {
                prop.caps_group: prop.dictionary for prop in self._properties
            }

            # add constraints to tacsAim
            self.aim.input.Constraint = {
                con.name: con.dictionary for con in self._constraints
            }

            # add loads to tacsAim
            if len(self._loads) > 0:
                self.aim.input.Load = {
                    load.name: load.dictionary for load in self._loads
                }

            if auto_shape_variables and self._first_setup:
                for despmtr in self._metadata.design_parameters:
                    # TODO : setup for dv arrays too but not yet
                    new_value = self._geometry.despmtr[despmtr].value
                    if isinstance(
                        new_value, float
                    ):  # make sure not a list despmtr, not supported yet
                        shape_var = ShapeVariable(name=despmtr, value=new_value)
                        self.add_variable(variable=shape_var)
                self._first_setup = False

            # link the egads aim to the tacs aim
            self.aim.input["Mesh"].link(self._mesh_aim.aim.output["Surface_Mesh"])

            # add the design variables to the DesignVariable and DesignVariableRelation properties
            if len(self.thickness_variables) > 0:
                self.aim.input.Design_Variable_Relation = {
                    dv.name: dv.DVR_dictionary
                    for dv in self._design_variables
                    if isinstance(dv, ThicknessVariable)
                }
            if len(self.variables) > 0:
                self.aim.input.Design_Variable = {
                    dv.name: dv.DV_dictionary for dv in self._design_variables
                }

        # end of serial or root proc section

        # note that setup is finished now
        self._setup = True
        return self  # return object for method cascading

    @root_proc
    def set_config_parameter(self, param_name: str, value: float):
        self.geometry.cfgpmtr[param_name].value = value
        return

    @root_broadcast
    def get_config_parameter(self, param_name: str):
        return self.geometry.cfgpmtr[param_name].value

    @root_broadcast
    def get_output_parameter(self, out_name: str):
        return self.geometry.outpmtr[out_name].value

    @property
    def geometry(self):
        """
        caps problem geometry object pyCAPS.problem.geometry
        """
        return self._geometry

    @property
    def variables(self) -> List[ShapeVariable or ThicknessVariable]:
        return self._design_variables

    @property
    def shape_variables(self) -> List[ShapeVariable]:
        return [dv for dv in self.variables if isinstance(dv, ShapeVariable)]

    @property
    def thickness_variables(self) -> List[ThicknessVariable]:
        """
        return sorted thickness vars so that the TACS derivatives can be appropriately obtained
        """
        thick_var_names = [
            dv.name for dv in self.variables if isinstance(dv, ThicknessVariable)
        ]
        thick_sorted_names = np.sort(np.array(thick_var_names))
        sorted_dvs = []
        for sort_name in thick_sorted_names:
            for var in self.variables:
                if isinstance(var, ThicknessVariable) and var.name == sort_name:
                    sorted_dvs.append(var)
                    break
        return sorted_dvs

    @property
    def analysis_dir(self) -> str:
        return self._metadata.analysis_dir

    @property
    def dat_file(self) -> str:
        return self.project_name + ".dat"

    @property
    def dat_file_path(self) -> str:
        return os.path.join(self.analysis_dir, self.dat_file)

    @property
    def sens_file(self) -> str:
        return self.project_name + ".sens"

    @property
    def sens_file_path(self) -> str:
        return os.path.join(self.analysis_dir, self.sens_file)

    @property
    def project_name(self) -> str:
        return self._metadata.project_name

    @property
    def is_setup(self) -> bool:
        return self._setup

    @property
    def aim(self):
        """
        returns the auto-built tacsAim object
        """
        return self._aim

    @property
    def change_shape(self) -> bool:
        """
        whether the aim will change shape (only if shape variables provided)
        """
        return len(self.shape_variables) > 0

    def update_properties(self):
        """
        update thickness properties and design variables in ESP/CAPS inputs
        if shape change w/ thickness variables
        """
        # exit if no thickness variables
        if len(self.thickness_variables) == 0:
            return

        # update property thicknesses by the modified thickness variables
        for property in self._properties:
            for dv in self._design_variables:
                if isinstance(property, ShellProperty) and isinstance(
                    dv, ThicknessVariable
                ):
                    if property.caps_group == dv.caps_group:
                        property.membrane_thickness == dv.value
                        break

        if self.comm is None or self.comm.rank == 0:
            # input new design var and property cards
            self.aim.input.Design_Variable = {
                dv.name: dv.DV_dictionary for dv in self._design_variables
            }
            self.aim.input.Property = {
                prop.caps_group: prop.dictionary for prop in self._properties
            }

        return

    @root_proc
    def pre_analysis(self):
        """
        provide access to the tacs aim preAnalysis for generating TACS input files and mesh
        """
        assert self._setup
        self.aim.preAnalysis()

        # return this object for method cascading
        return self

    @root_proc
    def post_analysis(self):
        """
        provide access to the tacs aim postAnalysis for collecting analysis outputs - functions and derivatives
        """
        self.aim.postAnalysis()
        return self
