__all__ = ["Load", "Pressure", "GridForce"]

from typing import TYPE_CHECKING, List


class Load:
    """
    Base class for FEA loads into the tacsAIM in ESP/CAPS
    load = {"groupName" : "plate",
            "loadType" : "Pressure",
            "pressureForce" : pressload}
    caps_load
        matches caps load attribute in CSM file
    """

    def __init__(self, name: str, caps_load: str, load_type: str):
        assert load_type in [
            "GridForce",
            "GridMoment",
            "Rotational",
            "Thermal",
            "Pressure",
            "PressureDistribute",
            "PressureExternal",
            "Gravity",
        ]
        self._name = name
        self._caps_load = caps_load
        self._load_type = load_type

    @property
    def name(self) -> str:
        return self._name

    def register_to(self, tacs_aim):
        """
        cascaded method to register this load to TacsAim
        """
        tacs_aim.register(self)
        return self


class Pressure(Load):
    """
    Apply pressure loads to the FEA structure
    """

    def __init__(self, caps_load: str, name: str = None, force: float = 2.0e6):
        if name is None:
            name = f"Pressure_{caps_load}"
        super(Pressure, self).__init__(
            name=name, caps_load=caps_load, load_type="Pressure"
        )
        self._pressure_force = force

    @property
    def dictionary(self) -> dict:
        return {
            "groupName": self._caps_load,
            "loadType": self._load_type,
            "pressureForce": self._pressure_force,
        }


class GridForce(Load):
    """
    Apply body forces to the FEA problem in tacsAIM
    """

    def __init__(
        self,
        caps_load: str,
        name: str = None,
        direction: List[float] = [0.0, 0.0, 1.0],
        magnitude: float = 1.0e3,
    ):
        if name is None:
            name = f"GridForce_{caps_load}"
        super(GridForce, self).__init__(
            name=name, caps_load=caps_load, load_type="GridForce"
        )
        self._direction_vector = direction
        self._force_scale_factor = magnitude

    @property
    def dictionary(self) -> dict:
        return {
            "groupName": self._caps_load,
            "loadType": self._load_type,
            "directionVector": self._direction_vector,
            "forceScaleFactor": self._force_scale_factor,
        }
