__all__ = ["TacsAim"]

from typing import TYPE_CHECKING, List
import pyCAPS
import os
from .materials import Material
from .constraints import Constraint
from .property import ShellProperty
from .loads import *
from .variables import CapsShapeVariable, CapsThicknessVariable


class TacsAim:
    """
    Wrapper class for TacsAim with default build setting in different scenarios
    applies default settings and spits it back out at the end
    only supports shell properties at the moment
    """

    def __init__(self, caps_problem: pyCAPS.Problem):
        self._aim = caps_problem.analysis.create(aim="tacsAIM", name="tacs")

        # geometry and design parameters to change the design of the CSM file during an optimization
        self._design_parameters = caps_problem.geometry.despmtr.keys()
        self._geometry = caps_problem.geometry

        self._materials = []
        self._loads = []
        self._properties = []
        self._constraints = []
        self._design_variables = []

        # build flags
        self._setup = False
        self._first_setup = True
        self._first_analysis = True

    def update_design(self, design_dict: dict):
        """
        method to change the values of each design variable in tacs, caps
        input x is a dictionary of values for each variable {"name" : value}
        """

        changed_design: bool = False

        # change setup to False since we have a new design we have to set up the aim again

        for idx, design_variable in enumerate(self._design_variables):
            dv_name = design_variable.name
            if dv_name in design_dict:

                new_value = float(design_dict[dv_name])

                if isinstance(design_variable, CapsShapeVariable):
                    same_value: bool = (
                        self._geometry.despmtr[dv_name].value == new_value
                    )
                    if not same_value:
                        print(f"changing dv '{dv_name}' to {new_value}")
                        design_variable.value = new_value
                        self._geometry.despmtr[dv_name].value = new_value
                        changed_design = True
                        # self.aim.geometry.despmtr[dv_name].value = new_value

                elif isinstance(design_variable, CapsThicknessVariable):

                    # change the thickness design variable
                    same_value = design_variable.value == new_value
                    if not same_value:
                        self._design_variables[idx].value = new_value

                    # change that property
                    matching_property = False

                    for idx, property in enumerate(self._properties):

                        if isinstance(property, ShellProperty):
                            matching_property = (
                                property.caps_group == design_variable.caps_group
                            )

                            if matching_property:
                                same_value: bool = (
                                    property.membrane_thickness == new_value
                                )

                                if not same_value:
                                    self._properties[idx].membrane_thickness = new_value
                                    changed_design = True
                                break

                    if not matching_property:
                        raise AssertionError(
                            f"Couldn't find matching property '{dv_name}'..."
                        )

        self._setup = not (changed_design)
        if self._first_analysis:
            self._first_analysis = False
            return True
        else:
            return changed_design

    def register(self, obj):
        """
        register any one of the available objects: Materials, Properties, Variables, etc.
        """
        if isinstance(obj, Material):
            self._materials.append(obj)
        elif isinstance(obj, CapsThicknessVariable):
            self._design_variables.append(obj)
            if obj.can_make_shell:
                self._properties.append(obj.shell_property)
        elif isinstance(obj, CapsShapeVariable):
            self._design_variables.append(obj)
        elif isinstance(obj, ShellProperty):
            self._properties.append(obj)
        elif isinstance(obj, Constraint):
            self._constraints.append(obj)
        elif isinstance(obj, Load):
            self._loads.append(obj)
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

        # increase the precision in the BDF file
        self._aim.input.File_Format = "Large" if large_format else "Small"
        self._aim.input.Mesh_File_Format = "Large" if large_format else "Small"

        # set the analysis type
        if static:
            self._aim.input.Analysis_Type = "Static"
        else:
            raise AssertionError(
                "Analysis types other than static analyses for tacsAim are not supported yet."
            )

        # add materials to tacsAim
        self._aim.input.Material = {
            material.name: material.dictionary for material in self._materials
        }

        # add properties to tacsAim
        self._aim.input.Property = {
            prop.caps_group: prop.dictionary for prop in self._properties
        }

        # add constraints to tacsAim
        self._aim.input.Constraint = {
            con.name: con.dictionary for con in self._constraints
        }

        # add loads to tacsAim
        if len(self._loads) > 0:
            self._aim.input.Load = {load.name: load.dictionary for load in self._loads}

        if auto_shape_variables and self._first_setup:
            for despmtr in self._design_parameters:
                # TODO : setup for dv arrays too but not yet
                new_value = self._geometry.despmtr[despmtr].value
                if isinstance(
                    new_value, float
                ):  # make sure not a list despmtr, not supported yet
                    shape_var = CapsShapeVariable(name=despmtr, value=new_value)
                    self.add_variable(variable=shape_var)
            self._first_setup = False

        # add the design variables to the DesignVariable and DesignVariableRelation properties
        self._aim.input.Design_Variable_Relation = {
            dv.name: dv.DVR_dictionary
            for dv in self._design_variables
            if isinstance(dv, ThicknessVariable)
        }
        self._aim.input.Design_Variable = {
            dv.name: dv.DV_dictionary for dv in self._design_variables
        }

        # note that setup is finished now
        self._setup = True

    @property
    def variables(self) -> List[CapsShapeVariable or CapsThicknessVariable]:
        return self._design_variables

    @property
    def shape_variables(self) -> List[CapsShapeVariable]:
        return [dv for dv in self.variables if isinstance(dv, CapsShapeVariable)]

    @property
    def thickness_variables(self) -> List[CapsThicknessVariable]:
        return [dv for dv in self.variables if isinstance(dv, CapsThicknessVariable)]

    def pre_analysis(self):
        """
        provide access to the tacs aim preAnalysis for running
        """
        self.aim.preAnalysis()

    @property
    def analysis_dir(self) -> str:
        return self.aim.analysisDir

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
        return self.aim.input.Proj_Name

    def post_analysis(self):
        self.aim.postAnalysis()

    @property
    def is_setup(self) -> bool:
        return self._setup

    @property
    def aim(self):
        """
        returns the auto-built tacsAim object
        """
        return self._aim
