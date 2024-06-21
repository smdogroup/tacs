#  This file is part of TACS: The Toolkit for the Analysis of Composite
#  Structures, a parallel finite-element code for structural and
#  multidisciplinary design optimization.
#
#  Copyright (C) 2014 Georgia Tech Research Corporation
#
#  TACS is licensed under the Apache License, Version 2.0 (the
#  "License"); you may not use this software except in compliance with
#  the License.  You may obtain a copy of the License at
#
#  http://www.apache.org/licenses/LICENSE-2.0

#distutils: language=c++

# Import from TACS for definitions
from tacs.cpp_headers.TACS cimport *

cdef extern from "TACSMaterialProperties.h":
    cdef cppclass TACSMaterialProperties(TACSObject):
        TACSMaterialProperties(TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar)
        TACSMaterialProperties(TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar,
                               TacsScalar, TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar, TacsScalar)
        void setDensity(TacsScalar)
        TacsScalar getDensity();
        void setSpecificHeat(TacsScalar)
        void getIsotropicProperties(TacsScalar*, TacsScalar*)
        void getOrthotropicProperties(TacsScalar*, TacsScalar*, TacsScalar*,
                                      TacsScalar*, TacsScalar*, TacsScalar*,
                                      TacsScalar*, TacsScalar*, TacsScalar*)
        void getStrengthProperties(TacsScalar*, TacsScalar*, TacsScalar*,
                                   TacsScalar*, TacsScalar*, TacsScalar*,
                                   TacsScalar*, TacsScalar*, TacsScalar*)
        void getCoefThermalExpansion(TacsScalar*, TacsScalar*, TacsScalar*)
        void getThermalConductivity(TacsScalar*, TacsScalar*, TacsScalar*)

        MaterialType getMaterialType();

    cdef cppclass TACSOrthotropicPly(TACSObject):
        TACSOrthotropicPly(TacsScalar, TACSMaterialProperties*)
        void setKSWeight(TacsScalar)
        void setUseMaxStrainCriterion()
        void setUseTsaiWuCriterion()

cdef extern from "TACSPlaneStressConstitutive.h":
    cdef cppclass TACSPlaneStressConstitutive(TACSConstitutive):
        TACSPlaneStressConstitutive(TACSMaterialProperties*,
                                    TacsScalar, int, TacsScalar, TacsScalar)

cdef extern from "TACSPhaseChangeMaterialConstitutive.h":
    cdef cppclass TACSPhaseChangeMaterialConstitutive(TACSConstitutive):
        TACSPhaseChangeMaterialConstitutive(TACSMaterialProperties*,
                                            TACSMaterialProperties*,
                                            TacsScalar, TacsScalar,
                                            TacsScalar, TacsScalar,
                                            int, TacsScalar, TacsScalar)

cdef extern from "TACSSolidConstitutive.h":
    cdef cppclass TACSSolidConstitutive(TACSConstitutive):
        TACSSolidConstitutive(TACSMaterialProperties*,
                              TacsScalar, int, TacsScalar, TacsScalar)
        TACSMaterialProperties* getMaterialProperties()

cdef extern from "TACSShellConstitutive.h":
    cdef cppclass TACSShellConstitutive(TACSConstitutive):
        void setDrillingRegularization(double)

cdef extern from "TACSIsoShellConstitutive.h":
    cdef cppclass TACSIsoShellConstitutive(TACSShellConstitutive):
        TACSIsoShellConstitutive(TACSMaterialProperties*, TacsScalar, int,
                                 TacsScalar, TacsScalar, TacsScalar)

cdef extern from "TACSCompositeShellConstitutive.h":
    cdef cppclass TACSCompositeShellConstitutive(TACSShellConstitutive):
        TACSCompositeShellConstitutive(int, TACSOrthotropicPly**, const TacsScalar*,
                                       const TacsScalar*, TacsScalar, TacsScalar)
        void getPlyThicknesses(TacsScalar*);
        void getPlyAngles(TacsScalar*);
        TacsScalar getThicknessOffset();

cdef extern from "TACSLamParamShellConstitutive.h":
    cdef cppclass TACSLamParamShellConstitutive(TACSShellConstitutive):
        TACSLamParamShellConstitutive(TACSOrthotropicPly*, TacsScalar, int, TacsScalar, TacsScalar,
                                      TacsScalar, TacsScalar, TacsScalar, int, int, int,
                                      TacsScalar, TacsScalar, TacsScalar,
                                      TacsScalar, TacsScalar, int, int, TacsScalar, TacsScalar)

cdef extern from "TACSIsoShellConstitutive.h":
    cdef cppclass TACSIsoShellConstitutive(TACSShellConstitutive):
        TACSIsoShellConstitutive(TACSMaterialProperties*, TacsScalar, int,
                                 TacsScalar, TacsScalar)

cdef extern from "TACSCompositeShellConstitutive.h":
    cdef cppclass TACSCompositeShellConstitutive(TACSShellConstitutive):
        TACSCompositeShellConstitutive(int, TACSOrthotropicPly**, const TacsScalar*,
                                       const TacsScalar*, TacsScalar)
        void getPlyThicknesses(TacsScalar*);
        void getPlyAngles(TacsScalar*);

cdef extern from "TACSSmearedCompositeShellConstitutive.h":
    cdef cppclass TACSSmearedCompositeShellConstitutive(TACSShellConstitutive):
        TACSSmearedCompositeShellConstitutive(int, TACSOrthotropicPly**, TacsScalar,
                                       const TacsScalar*, const TacsScalar*,
                                       int, const int*,
                                       TacsScalar, TacsScalar,
                                       const TacsScalar*, const TacsScalar*,
                                       TacsScalar)
        TacsScalar getLaminateThickness();
        void getPlyAngles(TacsScalar*);
        void getPlyFractions(TacsScalar*);
        TacsScalar getThicknessOffset();

cdef extern from "TACSLamParamShellConstitutive.h":
    cdef cppclass TACSLamParamShellConstitutive(TACSShellConstitutive):
        TACSLamParamShellConstitutive(TACSOrthotropicPly*, TacsScalar, int, TacsScalar, TacsScalar,
                                      TacsScalar, TacsScalar, TacsScalar, int, int, int,
                                      TacsScalar, TacsScalar, TacsScalar,
                                      TacsScalar, TacsScalar, int, int, TacsScalar, TacsScalar)

cdef extern from "TACSBladeStiffenedShellConstitutive.h":
    cdef cppclass TACSBladeStiffenedShellConstitutive(TACSShellConstitutive):
        TACSBladeStiffenedShellConstitutive(
            TACSOrthotropicPly*, # panelPly
            TACSOrthotropicPly*, # stiffenerPly
            TacsScalar, # kcorr
            TacsScalar, # panelLength
            int, # panelLengthNum
            TacsScalar, # stiffenerPitch
            int, # stiffenerPitchNum
            TacsScalar, # panelThick
            int, # panelThickNum
            int, # numPanelPlies
            TacsScalar[], # panelPlyAngles
            TacsScalar[], # panelPlyFracs
            int[], # panelPlyFracNums
            TacsScalar, # stiffenerHeight
            int, # stiffenerHeightNum
            TacsScalar, # stiffenerThick
            int, # stiffenerThickNum
            int, # numStiffenerPlies
            TacsScalar[], # stiffenerPlyAngles
            TacsScalar[], # stiffenerPlyFracs
            int[], # stiffenerPlyFracNums
            TacsScalar, # flangeFraction
        )
        int getNumPanelPlies()
        int getNumStiffenerPlies()
        void setIncludePanelMaterialFailure(bool _includePanelMaterialFailure)
        void setIncludeStiffenerMaterialFailure(bool _includeStiffenerMaterialFailure)
        void setIncludeGlobalBuckling(bool _includeGlobalBuckling)
        void setIncludeLocalBuckling(bool _includeLocalBuckling)
        void setIncludeStiffenerColumnBuckling(bool _includeStiffenerColumnBuckling)
        void setIncludeStiffenerCrippling(bool _includeStiffenerCrippling)
        void setKSWeight(double ksWeight)
        void setStiffenerPitchBounds(TacsScalar lowerBound, TacsScalar upperBound)
        void setStiffenerHeightBounds(TacsScalar lowerBound, TacsScalar upperBound)
        void setStiffenerThicknessBounds(TacsScalar lowerBound, TacsScalar upperBound)
        void setPanelThicknessBounds(TacsScalar lowerBound, TacsScalar upperBound)
        void setStiffenerPlyFractionBounds(TacsScalar[] lowerBound, TacsScalar[] upperBound)
        void setPanelPlyFractionBounds(TacsScalar[] lowerBound, TacsScalar[] upperBound)

cdef extern from "TACSBeamConstitutive.h":
    cdef cppclass TACSBeamConstitutive(TACSConstitutive):
        pass

cdef extern from "TACSBasicBeamConstitutive.h":
    cdef cppclass TACSBasicBeamConstitutive(TACSBeamConstitutive):
        TACSBasicBeamConstitutive(TACSMaterialProperties*, TacsScalar, TacsScalar, TacsScalar,
                                  TacsScalar, TacsScalar, TacsScalar, TacsScalar)

cdef extern from "TACSIsoTubeBeamConstitutive.h":
    cdef cppclass TACSIsoTubeBeamConstitutive(TACSBeamConstitutive):
        TACSIsoTubeBeamConstitutive(TACSMaterialProperties*, TacsScalar, TacsScalar,
                                    int, int, TacsScalar, TacsScalar, TacsScalar, TacsScalar)

cdef extern from "TACSIsoRectangleBeamConstitutive.h":
    cdef cppclass TACSIsoRectangleBeamConstitutive(TACSBeamConstitutive):
        TACSIsoRectangleBeamConstitutive(TACSMaterialProperties*, TacsScalar, TacsScalar,
                                         int, int, TacsScalar, TacsScalar, TacsScalar, TacsScalar,
                                         TacsScalar, TacsScalar)

cdef extern from "TACSGeneralMassConstitutive.h":
    cdef cppclass TACSGeneralMassConstitutive(TACSConstitutive):
        TACSGeneralMassConstitutive(const TacsScalar*)
        void evalMassMatrix(int, const double*, const TacsScalar*, TacsScalar*)

cdef extern from "TACSPointMassConstitutive.h":
    cdef cppclass TACSPointMassConstitutive(TACSGeneralMassConstitutive):
        TACSPointMassConstitutive(TacsScalar, TacsScalar, TacsScalar, TacsScalar, TacsScalar, TacsScalar, TacsScalar,
                                  int, TacsScalar, TacsScalar, int, TacsScalar, TacsScalar, int, TacsScalar, TacsScalar,
                                  int, TacsScalar, TacsScalar, int, TacsScalar, TacsScalar, int, TacsScalar, TacsScalar,
                                  int, TacsScalar, TacsScalar)

cdef extern from "TACSGeneralSpringConstitutive.h":
    cdef cppclass TACSGeneralSpringConstitutive(TACSConstitutive):
        TACSGeneralSpringConstitutive(TacsScalar*)
        void evalMassMatrix(int, const double *, const TacsScalar *, TacsScalar *)

cdef extern from "TACSDOFSpringConstitutive.h":
    cdef cppclass TACSDOFSpringConstitutive(TACSGeneralSpringConstitutive):
        TACSDOFSpringConstitutive(TacsScalar*)

# Special functions required for converting pointers
cdef extern from "":
    TACSBeamConstitutive* _dynamicBeamConstitutive"dynamic_cast<TACSBeamConstitutive*>"(TACSConstitutive*)

cdef extern from "TACSConstitutiveVerification.h":
    int TacsTestConstitutiveDensity(TACSConstitutive*, int, const double*, const TacsScalar*, int,
                                    const TacsScalar*, double, int, double, double)
    int TacsTestConstitutiveSpecificHeat(TACSConstitutive*, int, const double*, const TacsScalar*, int,
                                         const TacsScalar*, double, int, double, double)
    int TacsTestConstitutiveHeatFlux(TACSConstitutive*, int, const double*, const TacsScalar*, int,
                                     const TacsScalar*, double, int, double, double)
    int TacsTestConstitutiveStress(TACSConstitutive*, int, const double*, const TacsScalar*, int,
                                   const TacsScalar*, double, int, double, double)
    int TacsTestConstitutiveThermalStrain(TACSConstitutive*, int, const double*, const TacsScalar*, int,
                                          const TacsScalar*, double, int, double, double)
    int TacsTestConstitutiveFailure(TACSConstitutive*, int, const double*, const TacsScalar*, int,
                                    const TacsScalar*, double, int, double, double)
    int TacsTestConstitutiveFailureStrainSens(TACSConstitutive*, int, const double*, const TacsScalar*,
                                              double, int, double, double)
    int TacsTestConstitutive(TACSConstitutive*, int, double, int, double, double)
