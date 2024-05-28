import numpy as np

import openmdao.api as om

from mphys.core import MaskedConverter, MaskedVariableDescription


class TacsMesh(om.IndepVarComp):
    """
    Component to read the initial mesh coordinates with TACS
    """

    def initialize(self):
        self.options.declare(
            "fea_assembler",
            default=None,
            desc="the pytacs object itself",
            recordable=False,
        )
        self.options.declare(
            "discipline_vars",
            desc="the MPhys disciplinary variable enum structure object",
            default=False,
        )

    def setup(self):
        fea_assembler = self.options["fea_assembler"]
        self.discipline_vars = self.options["discipline_vars"]

        self.mesh_name = self.discipline_vars.Mesh.COORDINATES

        xpts = fea_assembler.getOrigNodes()
        self.add_output(
            self.mesh_name,
            distributed=True,
            val=xpts,
            shape=xpts.size,
            desc="fem node coordinates",
            tags=["mphys_coordinates"],
        )


class TacsMeshGroup(om.Group):
    def initialize(self):
        self.options.declare(
            "fea_assembler",
            default=None,
            desc="the pytacs object itself",
            recordable=False,
        )
        self.options.declare(
            "discipline_vars",
            desc="the MPhys disciplinary variable enum structure object",
            recordable=False,
        )

    def setup(self):
        fea_assembler = self.options["fea_assembler"]
        self.discipline_vars = self.options["discipline_vars"]

        self.mesh_name = self.discipline_vars.Mesh.COORDINATES
        self.add_subsystem(
            "fea_mesh",
            TacsMesh(fea_assembler=fea_assembler, discipline_vars=self.discipline_vars),
        )

        # Identify tacs nodes corresponding to lagrange multipliers. These are nodes that are typically added in tacs
        # whenever an element using a lagrange multiplier formulation is used (such as an RBE). It is important that
        # these nodes not be included in the aerostructural coupling procedure, as they a purely mathematical constructs.
        # We'll use this information later to create a mask for filtering out these nodes in the coupling procedure.
        nnodes = fea_assembler.getNumOwnedNodes()
        nmult = fea_assembler.getNumOwnedMultiplierNodes()
        mask_input = MaskedVariableDescription(
            self.mesh_name, shape=nnodes * 3, tags=["mphys_coordinates"]
        )
        mask_output = MaskedVariableDescription(
            f"{self.mesh_name}_masked",
            shape=(nnodes - nmult) * 3,
            tags=["mphys_coordinates"],
        )
        mult_ids = fea_assembler.getLocalMultiplierNodeIDs()
        mask = np.zeros([nnodes, 3], dtype=bool)
        mask[:, :] = True
        mask[mult_ids, :] = False
        x_orig = fea_assembler.getOrigNodes()
        x_masked = x_orig[mask.flatten()]
        masker = MaskedConverter(
            input=mask_input,
            output=mask_output,
            mask=mask.flatten(),
            init_output=x_masked,
            distributed=True,
        )
        self.add_subsystem(
            "masker",
            masker,
            promotes_outputs=[(f"{self.mesh_name}_masked", self.mesh_name)],
        )

        self.connect(f"fea_mesh.{self.mesh_name}", f"masker.{self.mesh_name}")
