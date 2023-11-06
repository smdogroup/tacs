"""
Written by Sean Engelstad, GT SMDO Lab, 2022-2023
"""
__all__ = ["Material", "Isotropic", "Orthotropic"]

from typing import TYPE_CHECKING


class Material:
    def __init__(
        self,
        name: str,
        material_type: str,
        E1=None,
        E2=None,
        E3=None,
        nu12=None,
        nu13=None,
        nu23=None,
        rho=None,
        cp=None,
        kappa1=None,
        kappa2=None,
        kappa3=None,
        alpha1=None,
        alpha2=None,
        alpha3=None,
        G12=None,
        G13=None,
        G23=None,
        T1=None,
        T2=None,
        C1=None,
        C2=None,
        S1=None,
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
        self._E1 = E1
        self._E2 = E2
        self._E3 = E3
        self._nu12 = nu12
        self._nu13 = nu13
        self._nu23 = nu23
        self._rho = rho
        self._cp = cp
        self._kappa1 = kappa1
        self._kappa2 = kappa2
        self._kappa3 = kappa3
        self._alpha1 = alpha1
        self._alpha2 = alpha2
        self._alpha3 = alpha3
        self._G12 = G12
        self._G13 = G13
        self._G23 = G23
        self._T1 = T1
        self._T2 = T2
        self._C1 = C1
        self._C2 = C2
        self._S1 = S1

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
        m_dict["density"] = self._rho
        m_dict["specificHeat"] = self._cp
        if self._kappa2 is None:
            m_dict["kappa"] = self._kappa1
        else:
            # [KXX, KXY, KXZ, KYY, KYZ, KZZ] = [k1, 0, 0, k2, 0, k3]
            m_dict["K"] = [self._kappa1, 0.0, 0.0, self._kappa2, 0.0, self._kappa3]
        m_dict["thermalExpCoeff"] = self._alpha1
        m_dict["thermalExpCoeffLateral"] = self._alpha2
        m_dict["youngModulus"] = self._E1
        m_dict["youngModulusLateral"] = self._E2
        m_dict["poissonRatio"] = self._nu12
        m_dict["poissonRatio23"] = self._nu23
        m_dict["shearModulus"] = self._G12
        m_dict["shearModulusTrans1Z"] = self._G13
        m_dict["shearModulusTrans2Z"] = self._G23
        m_dict["tensionAllow"] = self._T1
        m_dict["tensionAllowLateral"] = self._T2
        m_dict["compressionAllow"] = self._C1
        m_dict["compressionAllowLateral"] = self._C2
        m_dict["shearAllow"] = self._S1
        # m_dict["yieldAllow"] = self._yield_allow

        # return all items that are not None
        return {k: v for k, v in m_dict.items() if v is not None}

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
        E,
        nu,
        rho,
        T1,
        C1=None,
        S1=None,
        alpha=None,
        kappa=None,
        cp=None,
        G=None,
    ):
        """
        wrapper class for ESP/CAPS isotropic materials
        """
        if G is None:
            G = E / 2.0 / (1 + nu)
        if C1 is None:
            C1 = T1
        super(Isotropic, self).__init__(
            name=name,
            material_type="Isotropic",
            E1=E,
            nu12=nu,
            rho=rho,
            T1=T1,
            C1=C1,
            S1=S1,
            alpha1=alpha,
            kappa1=kappa,
            cp=cp,
            G12=G,
        )

    @classmethod
    def madeupium(
        cls,
        E=72.0e9,
        nu=0.33,
        rho=2.8e3,
        tension_allow=20.0e7,
    ):
        return cls(
            name="Madeupium",
            E=E,
            nu=nu,
            rho=rho,
            T1=tension_allow,
        )

    @classmethod
    def aluminum(cls):
        return cls(
            name="aluminum",
            E=70.0e9,
            nu=0.35,
            rho=2.7e3,
            T1=20.0e7,
            C1=20.0e7,
            S1=270e6,
            alpha=23.1e-6,
            cp=903,
            kappa=237,
        )

    @classmethod
    def titanium(cls):
        return cls(
            name="titanium",
            E=120e9,
            nu=0.361,
            rho=4.51e3,
            T1=2.4e7,
            C1=2.4e7,
            S1=0.6 * 2.4e7,  # estimated shear strength
            alpha=8.41e-6,
            cp=522.3,
            kappa=11.4,
        )

    @classmethod
    def titanium_alloy(cls):
        # Ti 6AL-4V (grade 5)
        return cls(
            name="titanium-alloy",
            E=114e9,
            nu=0.361,
            rho=4.43e3,
            T1=880e6,
            C1=970e6,
            S1=0.6 * 880e6,  # estimated shear strength
            alpha=9.2e-6,
            cp=526.3,
            kappa=6.7,
        )

    @classmethod
    def aluminum_alloy(cls):
        # Aluminum alloy Al-MS89
        return cls(
            name="aluminum-alloy",
            E=90e9,
            rho=2.92e3,
            nu=0.3,
            T1=420e6,
            C1=420e6,
            S1=0.6 * 420e6,  # estimated
            alpha=19.0e-6,
            kappa=115.0,
            cp=903,  # guessed the cp (not provided in data sheet)
        )

    @classmethod
    def steel(cls):
        return cls(
            name="steel",
            E=200.0e9,
            nu=0.30,
            rho=7.8e3,
            T1=1.0e9,
            C1=1.7e9,
            S1=0.6 * 420e6,  # estimated
            alpha=11.5e-6,
            kappa=45,
            cp=420,
        )


class Orthotropic(Material):
    def __init__(
        self,
        name: str,
        rho,
        E1,
        E2,
        nu12,
        nu13=None,
        nu23=None,
        E3=None,
        cp=None,
        kappa1=None,
        kappa2=None,
        kappa3=None,
        alpha1=None,
        alpha2=None,
        alpha3=None,
        G12=None,
        G13=None,
        G23=None,
        T1=None,
        T2=None,
        C1=None,
        C2=None,
        S1=None,
    ):
        super(Orthotropic, self).__init__(
            name=name,
            material_type="Orthotropic",
            E1=E1,
            E2=E2,
            nu12=nu12,
            nu23=nu23,
            rho=rho,
            T1=T1,
            G12=G12,
            G13=G13,
            G23=G23,
            T2=T2,
            C1=C1,
            C2=C2,
            S1=S1,
            E3=E3,
            nu13=nu13,
            cp=cp,
            kappa1=kappa1,
            kappa2=kappa2,
            kappa3=kappa3,
            alpha1=alpha1,
            alpha2=alpha2,
            alpha3=alpha3,
        )

    @classmethod
    def carbon_fiber(cls):
        # STD CF UD (carbon-fiber fiber/epoxy resin)
        # TODO : add more kinds here
        return cls(
            name="carbon_fiber_UD",
            E1=135e9,
            E2=10e9,
            G12=5e9,
            nu12=0.3,
            T1=1.5e9,
            C1=1.2e9,
            T2=50e6,
            C2=250e6,
            S1=70e9,
            alpha1=-0.3e-6,
            alpha2=28e-6,
            rho=1.6e3,
            kappa1=14.5,  # W/m-K
            kappa2=4.8,  # W/m-K
            kappa3=4.8,  # W/m-K
            cp=1130.0,  # J / kg-K
        )
