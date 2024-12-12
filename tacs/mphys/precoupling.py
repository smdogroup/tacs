import numpy as np

import openmdao.api as om
from mphys.core import UnmaskedConverter, MaskedVariableDescription

from .dv import TacsDVComp


class TacsPrecouplingGroup(om.Group):
    def initialize(self):
        self.options.declare(
            "fea_assembler",
            default=None,
            desc="the pytacs object itself",
            recordable=False,
        )
        self.options.declare(
            "discipline_vars",
            default=None,
            desc="the pytacs object itself",
            recordable=False,
        )
        self.options.declare(
            "initial_dv_vals",
            default=None,
            desc="initial values for global design variable vector",
        )
        self.options.declare(
            "separate_mass_dvs",
            default=False,
            desc="Flag for whether or not to separate out point mass dvs using user-defined names",
        )
        self.options.declare("discipline_vars")

    def setup(self):
        # Promote state variables/rhs with physics-specific tag that MPhys expects
        promotes_inputs = ["*"]

        fea_assembler = self.options["fea_assembler"]
        initial_dv_vals = self.options["initial_dv_vals"]
        separate_mass_dvs = self.options["separate_mass_dvs"]
        discipline_vars = self.options["discipline_vars"]

        coords_name = discipline_vars.COORDINATES

        self.add_subsystem(
            "distributor",
            TacsDVComp(
                fea_assembler=fea_assembler,
                initial_dv_vals=initial_dv_vals,
                separate_mass_dvs=separate_mass_dvs,
            ),
            promotes_inputs=promotes_inputs,
        )

        nnodes = fea_assembler.getNumOwnedNodes()
        nmult = fea_assembler.getNumOwnedMultiplierNodes()
        unmask_output = MaskedVariableDescription(
            coords_name, shape=nnodes * 3, tags=["mphys_coordinates"]
        )
        unmask_input = MaskedVariableDescription(
            f"{coords_name}_masked",
            shape=(nnodes - nmult) * 3,
            tags=["mphys_coordinates"],
        )
        mult_ids = fea_assembler.getLocalMultiplierNodeIDs()
        mask = np.zeros([nnodes, 3], dtype=bool)
        mask[:, :] = True
        mask[mult_ids, :] = False
        vals = fea_assembler.getOrigNodes()
        unmasker = UnmaskedConverter(
            input=unmask_input,
            output=unmask_output,
            mask=mask.flatten(),
            default_values=vals,
            distributed=True,
        )
        self.add_subsystem(
            "unmasker",
            unmasker,
            promotes_inputs=[(f"{coords_name}_masked", coords_name)],
        )
