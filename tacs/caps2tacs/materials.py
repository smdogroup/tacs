__all__ = ["Material", "Isotropic", "Orthotropic"]

from typing import TYPE_CHECKING


class Material:
    def __init__(
        self,
        name: str,
        material_type: str,
        young_modulus: float,
        poisson_ratio: float,
        density: float,
        tension_allow: float,
        compression_allow: float = None,
        shear_allow: float = None,
        yield_allow: float = None,
        thermExpCoeff: float = None,
    ):
        """
        Material base class to wrap ESP/CAPS material inputs to TACS AIM
        """
        assert material_type in [
            "Isotropic",
            "Anisothotropic",
            "Orthotropic",
            "Anisotropic",
        ]
        self._name = name
        self._material_type = material_type
        self._young_modulus = young_modulus
        self._poisson_ratio = poisson_ratio
        self._density = density
        self._tension_allow = tension_allow
        self._compression_allow = compression_allow
        self._shear_allow = shear_allow
        self._yield_allow = yield_allow
        self._thermExpCoeff = thermExpCoeff

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, new_name: str):
        self._name = new_name

    @property
    def dictionary(self) -> dict:
        """
        return dictionary of material settings
        """
        m_dict = {}
        m_dict["materialType"] = self._material_type
        m_dict["youngModulus"] = self._young_modulus
        m_dict["poissonRatio"] = self._poisson_ratio
        m_dict["density"] = self._density
        m_dict["thermalExpCoeff"] = self._thermExpCoeff
        m_dict["tensionAllow"] = self._tension_allow
        m_dict["compressionAllow"] = self._compression_allow
        m_dict["shearAllow"] = self._shear_allow
        m_dict["yieldAllow"] = self._yield_allow

        # return all items that are not None
        return {k: v for k, v in m_dict.items() if v is not None}

    @property
    def young_modulus(self) -> float:
        return self._young_modulus

    @young_modulus.setter
    def young_modulus(self, value: float):
        self._young_modulus = value

    def register_to(self, tacs_aim):
        """
        cascaded method to register this constraint to TacsAim
        """
        tacs_aim.register(self)
        return self


class Isotropic(Material):
    def __init__(
        self,
        name: str,
        young_modulus: float,
        poisson_ratio: float,
        density: float,
        tension_allow: float,
        compression_allow: float = None,
        shear_allow: float = None,
        yield_allow: float = None,
        thermExpCoeff: float = None,
    ):
        """
        wrapper class for ESP/CAPS isotropic materials
        """
        super(Isotropic, self).__init__(
            name=name,
            material_type="Isotropic",
            young_modulus=young_modulus,
            poisson_ratio=poisson_ratio,
            density=density,
            tension_allow=tension_allow,
            compression_allow=compression_allow,
            shear_allow=shear_allow,
            yield_allow=yield_allow,
            thermExpCoeff=thermExpCoeff,
        )

    @classmethod
    def madeupium(
        cls,
        young_modulus=72.0e9,
        poisson_ratio=0.33,
        density=2.8e3,
        tension_allow=20.0e7,
    ):
        return cls(
            name="Madeupium",
            young_modulus=young_modulus,
            poisson_ratio=poisson_ratio,
            density=density,
            tension_allow=tension_allow,
        )

    @classmethod
    def aluminum(cls):
        return cls(
            name="aluminum",
            young_modulus=70.0e9,
            poisson_ratio=0.35,
            density=2.7e3,
            tension_allow=20.0e7,
            compression_allow=20.0e7,
            yield_allow=20.0e7,
            thermExpCoeff=23.1e-6,
        )

    @classmethod
    def steel(cls):
        return cls(
            name="steel",
            young_modulus=200.0e9,
            poisson_ratio=0.30,
            density=7.8e3,
            tension_allow=1.0e9,
            compression_allow=1.7e9,
            yield_allow=0.9e9,
            thermExpCoeff=11.5e-6,
        )


class Orthotropic(Material):
    # TBD
    pass
