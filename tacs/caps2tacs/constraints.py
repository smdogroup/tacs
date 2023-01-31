__all__ = ["Constraint", "ZeroConstraint", "TemperatureConstraint"]


class Constraint:
    def __init__(
        self,
        name: str,
        caps_constraint: str,
        dof_constraint: int,
        grid_displacement: float = 0.0,
    ):
        """
        Base class for ESP/CAPS constraint wrapper class
        """
        dof_str = str(dof_constraint)
        for char in dof_str:
            assert int(char) in range(
                0, 7
            )  # means only allow dof 0-6, see pyMeshLoader _isDOFinString for more info
        self._name = name
        self._caps_constraint = caps_constraint
        self._dof_constraint = dof_constraint
        self._grid_displacement = grid_displacement

    @property
    def name(self) -> str:
        return self._name

    @property
    def dictionary(self) -> dict:
        return {
            "groupName": self._caps_constraint,
            "constraintType": "Displacement",
            "dofConstraint": self._dof_constraint,
            "gridDisplacement": self._grid_displacement,
        }

    def register_to(self, tacs_aim):
        """
        cascaded method to register this constraint to TacsAim
        """
        tacs_aim.register(self)
        return self


class ZeroConstraint(Constraint):
    def __init__(self, name: str, caps_constraint: str, dof_constraint: int = 123):
        """
        Elastic pin ESP/CAPS Constraint by default
        Can change the dof to other than u=v=w=0 aka 123
        """
        super(ZeroConstraint, self).__init__(
            name=name,
            caps_constraint=caps_constraint,
            constraint_type="Displacement",
            dof_constraint=dof_constraint,
        )

    @property
    def dictionary(self) -> dict:
        return {
            "groupName": self._caps_constraint,
            "constraintType": "ZeroDisplacement",
            "dofConstraint": self._dof_constraint,
        }


class TemperatureConstraint(Constraint):
    def __init__(self, name: str, caps_constraint: str, temperature: float = 0.0):
        """
        Isothermal constraints in ESP/CAPS for the TacsAim
        """
        super(TemperatureConstraint, self).__init__(
            name=name,
            caps_constraint=caps_constraint,
            constraint_type="Thermal",
            dof_constraint=0,
            grid_displacement=temperature,
        )

    @property
    def dictionary(self) -> dict:
        return super(TemperatureConstraint, self).dictionary
