__all__ = ["BaseProperty", "ShellProperty", "CompositeProperty"]
"""
Written by Sean Engelstad, GT SMDO Lab, 2022-2023
"""

from .materials import Material
from typing import TYPE_CHECKING


class BaseProperty:
    def __init__(self, caps_group: str, material: Material, property_type: str):
        self._caps_group = caps_group
        self._material = material
        self._property_type = property_type

    @property
    def caps_group(self) -> str:
        """
        return capsGroup attribute associated with this property
        """
        return self._caps_group

    @property
    def dictionary(self) -> dict:
        """
        return property dictionary, however this is only fully defined in subclasses
        """
        return {}

    def register_to(self, tacs_aim):
        """
        cascaded method to register this property to TacsAim
        """
        tacs_aim.register(self)
        return self


class ShellProperty(BaseProperty):
    """
    Example ShellProperty Dictionary
    shell  = {"propertyType" : "Shell",
                    "membraneThickness" : 0.006,
                    "material"        : "madeupium",
                    "bendingInertiaRatio" : 1.0, # Default
                    "shearMembraneRatio"  : 5.0/6.0} # Default
    self._aim.input.Property = {"plate": shell}
    """

    # TODO : add other available settings for shell properties -> mass, inertias, etc.
    def __init__(
        self,
        caps_group: str,
        material: Material,
        membrane_thickness: float,
        bending_inertia: float = 1.0,
        shear_membrane_ratio: float = 5.0 / 6.0,
    ):
        super(ShellProperty, self).__init__(
            caps_group=caps_group, material=material, property_type="Shell"
        )
        self._membrane_thickness = membrane_thickness
        self._bending_inertia = bending_inertia
        self._shear_membrane_ratio = shear_membrane_ratio

    @property
    def membrane_thickness(self) -> float:
        return self._membrane_thickness

    @membrane_thickness.setter
    def membrane_thickness(self, new_thickness: float):
        self._membrane_thickness = new_thickness

    @property
    def dictionary(self) -> dict:
        """
        return property dictionary to pass into tacsAim
        """
        return {
            "propertyType": self._property_type,
            "membraneThickness": self._membrane_thickness,
            "material": self._material.name,
            "bendingInertiaRatio": self._bending_inertia,  # Default
            "shearMembraneRatio": self._shear_membrane_ratio,
        }  # Default-

    def register_to(self, tacs_aim):
        """
        cascaded method to register this ShellProperty to TacsAim
        """
        tacs_aim.register(self)
        return self


class CompositeProperty(BaseProperty):
    """
    Define a composite material property from ESP/CAPS
    Tip: you need to use the csystem command and a 1x9 list to define
    "Ox,Oy,Oz, d1x, d1y, d1z, d2x, d2y, d2z"
    which is the origin, 1-direction, and 2-direction used for shells
    """

    def __init__(
        self,
        caps_group: str,
        ply_materials: list,
        ply_thicknesses: list,
        ply_angles: list,
        sym_laminate: bool = True,
        composite_failure_theory: str = "STRN",
        shear_bond_allowable: float = 1.0e9,  # high in case just stringer-one ply
        bending_inertia: float = 1.0,
        shear_membrane_ratio: float = 0.0,
    ):
        self._caps_group = caps_group
        self._ply_materials = ply_materials
        self._ply_thicknesses = ply_thicknesses
        self._ply_angles = ply_angles  # in degrees
        self._sym_laminate = sym_laminate
        self._composite_failure_theory = composite_failure_theory
        self._shear_bond_allowable = shear_bond_allowable
        self._bending_inertia = bending_inertia
        self._shear_membrane_ratio = shear_membrane_ratio

    @classmethod
    def one_ply(
        cls,
        caps_group: str,
        material: Material,
        thickness: float,
        ply_angle: float,
        sym_laminate: bool = True,
        composite_failure_theory: str = "STRN",
        shear_bond_allowable: float = 1.0e9,
        bending_inertia: float = 1.0,
        shear_membrane_ratio: float = 0.0,
    ) -> BaseProperty:
        return cls(
            caps_group=caps_group,
            ply_materials=[material],
            ply_thicknesses=[thickness],
            ply_angles=[ply_angle],
            sym_laminate=sym_laminate,
            composite_failure_theory=composite_failure_theory,
            shear_bond_allowable=shear_bond_allowable,
            bending_inertia=bending_inertia,
            shear_membrane_ratio=shear_membrane_ratio,
        )

    @property
    def ply_materials(self) -> list:
        if isinstance(self._ply_materials[0], str):
            return self._ply_materials
        elif isinstance(self._ply_materials[0], Material):
            return [_.name for _ in self._ply_materials]
        else:
            raise AssertionError(
                "caps2tacs error: Ply material objects need to be list of material name strings or material objects.."
            )

    @property
    def dictionary(self) -> dict:
        """
        return property dictionary to pass into tacsAim
        """
        return {
            "propertyType": "Composite",
            "shearBondAllowable": self._shear_bond_allowable,
            "bendingInertiaRatio": self._bending_inertia,
            "shearMembraneRatio": self._shear_membrane_ratio,
            "compositeMaterial": self.ply_materials,
            "compositeThickness": self._ply_thicknesses,
            "compositeOrientation": self._ply_angles,  # in degrees
            "symmetricLaminate": self._sym_laminate,
            "compositeFailureTheory": self._composite_failure_theory,
        }

    def register_to(self, tacs_aim):
        """
        cascaded method to register this ShellProperty to TacsAim
        """
        tacs_aim.register(self)
        return self
