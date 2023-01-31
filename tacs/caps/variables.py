__all__ = ["CapsShapeVariable", "CapsThicknessVariable"]


class CapsShapeVariable:
    """
    shape variables in ESP/CAPS are design parameters that affect the structural geometry
    """

    def __init__(self, name: str, value: float):
        """
        ESP/CAPS shape variable controls a design parameter in the CSM file
            name: corresponds to the design parameter in the CSM file
            value: can be used to modify the design parameter
        """
        self.name = name
        self.value = value

    @property
    def DV_dictionary(self) -> dict:
        return {}


class CapsThicknessVariable:
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
        self.value = value
        if name is not None:  # copy caps_group to name if not specified
            self.name = name
        else:
            self.name = caps_group
        self.lower_bound = lower_bound if lower_bound is not None else value * 0.5
        self.upper_bound = upper_bound if upper_bound is not None else value * 1.5
        self.max_delta = max_delta if max_delta is not None else value * 0.1

        # private variables used to create shell property
        self._material = None

    def material(self, material: Material):
        """
        method cascading setter method
        """
        self._material = material
        return self

    def set_bounds(self, lower_bound: float, value: float, upper_bound: float):

        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        return self

    @property
    def DV_dictionary(self) -> dict:
        """
        ESP/CAPS design variable dictionary
        """
        return {
            "groupName": self.caps_group,
            "initialValue": self.value,
            "lowerBound": self.lower_bound,
            "upperBound": self.upper_bound,
            "maxDelta": self.max_delta,
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
    def shell_property(self) -> ShellProperty:
        assert self._material is not None
        return ShellProperty(
            caps_group=self.caps_group,
            material=self._material,
            membrane_thickness=self.value,
        )
