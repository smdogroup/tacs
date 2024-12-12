import numpy as np

import openmdao.api as om

from mphys.core import (
    MaskedConverter,
    UnmaskedConverter,
    MaskedVariableDescription,
    DistributedSummer,
    DistributedVariableDescription,
)

from .solver import TacsSolver


class TacsCouplingGroup(om.Group):
    def initialize(self):
        self.options.declare("fea_assembler", recordable=False)
        self.options.declare("check_partials")
        self.options.declare("discipline_vars")
        self.options.declare("coupling_loads", default=[])
        self.options.declare("scenario_name", default=None)
        self.options.declare("problem_setup", default=None)
        self.options.declare("res_ref", default=None)

    def setup(self):
        self.fea_assembler = self.options["fea_assembler"]
        self.check_partials = self.options["check_partials"]
        self.coupling_loads = self.options["coupling_loads"]
        self.coupled = len(self.coupling_loads) > 0
        self.discipline_vars = self.options["discipline_vars"]
        self.res_ref = self.options["res_ref"]

        self.coords_name = self.discipline_vars.COORDINATES
        self.states_name = self.discipline_vars.STATES
        self.rhs_name = "rhs"

        # Promote state variables/rhs with physics-specific tag that MPhys expects
        promotes_inputs = [
            (self.coords_name, f"unmasker.{self.coords_name}"),
            ("tacs_dvs", "distributor.tacs_dvs"),
        ]

        # Identify tacs nodes corresponding to lagrange multipliers. These are nodes that are typically added in tacs
        # whenever an element using a lagrange multiplier formulation is used (such as an RBE). It is important that
        # these nodes not be included included in the aerostructural coupling procedure, as they a purely mathematical constructs.
        # We'll use this information later to create a mask for filtering out these nodes in the coupling procedure.
        nnodes = self.fea_assembler.getNumOwnedNodes()
        nmult = self.fea_assembler.getNumOwnedMultiplierNodes()
        vpn = self.fea_assembler.getVarsPerNode()
        mult_ids = self.fea_assembler.getLocalMultiplierNodeIDs()
        mask = np.zeros([nnodes, vpn], dtype=bool)
        mask[:, :] = True
        mask[mult_ids, :] = False

        # Create an unmasking component to process the coupled force vector from the aerostructural
        # load transfer component. This component takes the masked load vector (length = vpn * (nnodes - nmult))
        # and creates an unmasked load vector (length = vpn * nnodes), by inserting 0.0 for the forces on
        # the multiplier nodes
        if self.coupled:
            self.rhs_name_masked = f"{self.rhs_name}_masked"
            input_metadata = []
            for i, load_name in enumerate(self.coupling_loads):
                var_desc = DistributedVariableDescription(
                    name=load_name,
                    shape=(nnodes - nmult) * vpn,
                    tags=["mphys_coupling"],
                )
                input_metadata.append(var_desc)
            output_metadata = DistributedVariableDescription(
                name=self.rhs_name_masked, shape=(nnodes - nmult) * vpn
            )
            rhs_sum = DistributedSummer(inputs=input_metadata, output=output_metadata)
            self.add_subsystem("rhs_sum", rhs_sum, promotes_inputs=self.coupling_loads)

            unmask_output = MaskedVariableDescription(self.rhs_name, shape=nnodes * vpn)
            unmask_input = MaskedVariableDescription(
                self.rhs_name + "_masked",
                shape=(nnodes - nmult) * vpn,
            )
            unmasker = UnmaskedConverter(
                input=unmask_input,
                output=unmask_output,
                mask=mask.flatten(),
                distributed=True,
            )
            self.add_subsystem("unmasker", unmasker)
            self.connect(
                f"rhs_sum.{self.rhs_name}_masked", f"unmasker.{self.rhs_name}_masked"
            )

        self.add_subsystem(
            "solver",
            TacsSolver(
                fea_assembler=self.fea_assembler,
                check_partials=self.check_partials,
                coupling_loads=self.coupling_loads,
                discipline_vars=self.discipline_vars,
                res_ref=self.res_ref,
            ),
            promotes_inputs=promotes_inputs,
        )

        # Create a masking component to process the full structural displacement vector (length = vpn * nnodes)
        # and remove indices corresponding to multiplier nodes (length = vpn*(nnodes - nmult)).
        mask_input = MaskedVariableDescription(
            self.states_name, shape=nnodes * vpn, tags=["mphys_coupling"]
        )
        mask_output = MaskedVariableDescription(
            self.states_name + "_masked",
            shape=(nnodes - nmult) * vpn,
            tags=["mphys_coupling"],
        )
        masker = MaskedConverter(
            input=mask_input,
            output=mask_output,
            mask=mask.flatten(),
            distributed=True,
            init_output=0.0,
        )
        self.add_subsystem(
            "masker",
            masker,
            promotes_outputs=[(self.states_name + "_masked", self.states_name)],
        )

        self.connect("solver." + self.states_name, "masker." + self.states_name)

        if self.coupled:
            self.connect("unmasker." + self.rhs_name, "solver." + self.rhs_name)

        # Name problem based on scenario that's calling builder
        scenario_name = self.options["scenario_name"]
        if scenario_name is None:
            # Default structural problem
            name = "default"
        else:
            name = scenario_name
        sp = self.fea_assembler.createStaticProblem(name=name)

        # Setup TACS problem with user-defined structural loads
        problem_setup = self.options["problem_setup"]
        if problem_setup is not None:
            new_problem = problem_setup(scenario_name, self.fea_assembler, sp)
            # Check if the user provided back a new problem to overwrite the default
            if new_problem is not None:
                sp = new_problem

        # Set problem
        self.solver.set_sp(sp)

        self.sp = sp

    def write_bdf(self, file_name):
        """
        Write optimized structure and loads to BDF file.
        """
        self.fea_assembler.writeBDF(file_name, self.sp)
