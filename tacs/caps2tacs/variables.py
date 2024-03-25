"""
Written by Sean Engelstad, GT SMDO Lab, 2022-2023
"""

__all__ = ["ShapeVariable", "ThicknessVariable"]

from .materials import *
from .property import *


class ShapeVariable:
    """
    shape variables in ESP/CAPS are design parameters that affect the structural geometry
    """

    def __init__(self, name: str, value=None, active: bool = True):
        """
        ESP/CAPS shape variable controls a design parameter in the CSM file
            name: corresponds to the design parameter in the CSM file
            value: can be used to modify the design parameter
        """
        self.name = name
        self._value = value
        self._active = active

    @property
    def DV_dictionary(self) -> dict:
        return {}

    def register_to(self, tacs_aim):
        """
        cascaded method to register this ShapeVariable to TacsAim
        """
        tacs_aim.register(self)
        return self

    @property
    def value(self) -> float:
        return self._value

    @value.setter
    def value(self, new_value: float):
        self._value = new_value


class ThicknessVariable:
    """
    Caps Thickness Variables control membrane thickness of shell elements
    """

    def __init__(
        self,
        caps_group: str,
        value: float = 1.0,
        name: str = None,
        lower_bound: float = None,
        upper_bound: float = None,
        max_delta: float = None,
        material: Material = None,
        ply_angle: float = None,
        active: bool = True,
    ):
        """
        ESP/CAPS Thickness variable sets the thickness over a portion of the geometry in the CSM file
            caps_group: corresponds to the geometry attribute CapsGroup in the CSM file
        Optional arguments
            value: sets the membrane thickness of that section of the geometry's shell elements
            material: corresponding material, can be used to auto-create a shell property
            name: name of the variable, can be set to the same as caps_group
            lower_bound, upper_bound, and max_delta are bounds for optimization
        """
        self.caps_group = caps_group
        self._value = value
        if name is not None:  # copy caps_group to name if not specified
            self.name = name
        else:
            self.name = caps_group
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.max_delta = max_delta
        self.ply_angle = ply_angle
        self._active = active

        # private variables used to create shell property
        self._material = material

    def material(self, material: Material):
        """
        method cascading setter method
        """
        self._material = material
        return self

    def set_bounds(self, lower_bound: float, value: float, upper_bound: float):
        self.lower_bound = lower_bound
        self.value = value
        self.upper_bound = upper_bound
        return self

    @property
    def value(self) -> float:
        return self._value

    @value.setter
    def value(self, new_value: float):
        self.lower_bound = 0.5 * new_value
        self._value = new_value
        self.upper_bound = 2.0 * new_value
        return

    @property
    def can_make_shell(self):
        return self._material is not None and self.value is not None

    @property
    def DV_dictionary(self) -> dict:
        """
        ESP/CAPS design variable dictionary
        """
        return {
            "groupName": self.caps_group,
            "initialValue": self.value,
            "lowerBound": (
                self.lower_bound if self.lower_bound is not None else self.value * 0.5
            ),
            "upperBound": (
                self.upper_bound if self.upper_bound is not None else self.value * 2.0
            ),
            "maxDelta": (
                self.max_delta if self.max_delta is not None else self.value * 0.1
            ),
        }

    @property
    def DVR_dictionary(self) -> dict:
        """
        ESP/CAPS design variable relations dictionary
        """
        return {
            "componentType": "Property",
            "fieldName": "T",
            "componentName": self.caps_group,
            "variableName": self.name,
        }

    @property
    def auto_property(self) -> BaseProperty:
        # automatic property is made if a material is provided
        assert self._material is not None
        if isinstance(self._material, Isotropic):
            return ShellProperty(
                caps_group=self.caps_group,
                material=self._material,
                membrane_thickness=self.value,
            )
        elif isinstance(self._material, Orthotropic):
            # auto-making a unidirectional laminate CompositeProperty (for ease of use in large cases with metal + smeared stringer models)
            return CompositeProperty(
                caps_group=self.caps_group,
                ply_materials=[self._material],
                ply_thicknesses=[self.value],
                ply_angles=[self.ply_angle],
            )
        else:
            raise AssertionError(
                f"Can't directly make a thickness variable {self.name} from a composite material."
                + "Don't provide a material => just give the caps_group for the thickness variable"
                + " and separately make the composite property (registering each to the tacs model or tacs aim)."
            )

    def register_to(self, tacs_aim):
        """
        cascaded method to register this ThicknessVariable to TacsAim
        """
        tacs_aim.register(self)
        return self
