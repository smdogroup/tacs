__all__ = ["AnalysisFunction"]

from tacs.functions import (
    KSFailure,
    StructuralMass,
    KSTemperature,
    AverageTemperature,
    Compliance,
)


class AnalysisFunction:
    def __init__(self, name: str, handle, compIDs=None, scale: float = 1.0, **kwargs):
        """
        Analysis functions/functionals for structural optimization
            Every function is added to new structural analysis problems in each iteration if the design changes
                with the TacsAim shape variables for instance.
            Calls addFunction in each new staticProblem for now

        Parameters
        -----------------------------------------------------
        name : str
            the name of the function i.e. mass, ks_vmfailure, etc.
        handle : TACS.function
            the function handle of the cython TACS function, comes from tacs.functions module
        compIDs : list
            list of component IDs to select
        scale : float
            the scale used for optimization
        **kwargs:
            any keyword arguments to pass to the
        """

        self.name = name
        self.handle = handle
        self.compIDs = compIDs
        self.scale = scale
        self.kwargs = kwargs

    def register_to(self, tacs_aim):
        # register the analysis function to the TacsAim wrapper class
        tacs_aim.register(self)
        return

    # class methods to add certain functions
    @classmethod
    def mass(cls, compIDs=None, scale: float = 1.0):
        return cls(name="mass", handle=StructuralMass, compIDs=compIDs, scale=scale)

    @classmethod
    def ksfailure(
        cls,
        compIDs=None,
        safetyFactor: float = 1.0,
        ksWeight: float = 50.0,
        scale: float = 1.0,
    ):
        return cls(
            name="ksfailure",
            handle=KSFailure,
            compIDs=compIDs,
            scale=scale,
            safetyFactor=safetyFactor,
            ksWeight=ksWeight,
        )

    @classmethod
    def ks_temperature(
        cls,
        compIDs=None,
        alpha: float = 1.0,
        ksWeight: float = 50.0,
        scale: float = 1.0,
    ):
        return cls(
            name="ks_temperature",
            handle=KSTemperature,
            compIDs=compIDs,
            scale=scale,
            alpha=alpha,
            ksWeight=ksWeight,
        )

    @classmethod
    def avg_temperature(cls, compIDs=None, volume: float = 1.0, scale: float = 1.0):
        return cls(
            name="avg_temperature",
            handle=AverageTemperature,
            compIDs=compIDs,
            scale=scale,
            volume=volume,
        )

    @classmethod
    def compliance(cls, compIDs=None, scale: float = 1.0):
        return cls(name="compliance", handle=Compliance, compIDs=compIDs, scale=scale)
