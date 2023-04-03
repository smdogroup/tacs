import numpy as np
import copy

import openmdao.api as om
from openmdao.utils.mpi import MPI

from mphys.builder import Builder
from mphys import MaskedConverter, UnmaskedConverter, MaskedVariableDescription

from .. import pyTACS, functions


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

    def setup(self):
        fea_assembler = self.options["fea_assembler"]
        xpts = fea_assembler.getOrigNodes()
        self.add_output(
            "x_struct0",
            distributed=True,
            val=xpts,
            shape=xpts.size,
            desc="structural node coordinates",
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

    def setup(self):
        fea_assembler = self.options["fea_assembler"]
        self.add_subsystem("fea_mesh", TacsMesh(fea_assembler=fea_assembler))

        # Identify tacs nodes corresponding to lagrange multipliers. These are nodes that are typically added in tacs
        # whenever an element using a lagrange multiplier formulation is used (such as an RBE). It is important that
        # these nodes not be included included in the aerostructural coupling procedure, as they a purely mathematical constructs.
        # We'll use this information later to create a mask for filtering out these nodes in the coupling procedure.
        nnodes = fea_assembler.getNumOwnedNodes()
        nmult = fea_assembler.getNumOwnedMultiplierNodes()
        mask_input = MaskedVariableDescription(
            "x_struct0", shape=nnodes * 3, tags=["mphys_coordinates"]
        )
        mask_output = MaskedVariableDescription(
            "x_struct0_masked", shape=(nnodes - nmult) * 3, tags=["mphys_coordinates"]
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
            "masker", masker, promotes_outputs=[("x_struct0_masked", "x_struct0")]
        )

        self.connect("fea_mesh.x_struct0", "masker.x_struct0")


class TacsDVComp(om.ExplicitComponent):
    """
    Component for splitting serial tacs design variable from top level
    into distributed vector used by tacs.
    """

    def initialize(self):
        self.options.declare(
            "fea_assembler",
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

    def setup(self):
        self.fea_assembler = self.options["fea_assembler"]
        self.src_indices = self.get_dv_src_indices()
        vals = self.options["initial_dv_vals"]
        ndv = self.fea_assembler.getNumDesignVars()
        # Keep a list of dvs corresponding to regular struct dvs and mass dvs
        self.struct_dvs = list(range(len(vals)))
        self.mass_dvs = {}
        # Check if user wants to separate out point mass dvs by user-defined names
        if self.options["separate_mass_dvs"]:
            g_dvs = self.fea_assembler.getGlobalDVs()
            for dv_name in g_dvs:
                dv_dict = g_dvs[dv_name]
                if dv_dict["isMassDV"]:
                    dv_num = dv_dict["num"]
                    mass_val = vals[dv_num]
                    # Store mass dv num with user-defined dv name
                    self.mass_dvs[f"dv_mass_{dv_name}"] = dv_num
                    # Remove mass dv from struct list
                    self.struct_dvs.remove(dv_num)
                    # Add user-defined dv name as input
                    self.add_input(
                        f"dv_mass_{dv_name}",
                        desc="serial mass design variable holding one mass design variable instance for tacs",
                        val=mass_val,
                        distributed=False,
                        tags=["mphys_input"],
                    )
            # Remove all mass dvs from vals
            vals = vals[self.struct_dvs]

        self.add_input(
            "dv_struct",
            desc="serial vector holding all structural tacs design variable values",
            val=vals,
            distributed=False,
            tags=["mphys_input"],
        )
        self.add_output(
            "tacs_dvs",
            desc="distributed vector holding tacs design variable values\
                        for this proc",
            shape=ndv,
            distributed=True,
            tags=["mphys_coupling"],
        )

    def compute(self, inputs, outputs):
        # Create serial array to holding all dv vals
        tot_ndv = len(self.struct_dvs) + len(self.mass_dvs)
        full_dv_array = np.zeros(tot_ndv, dtype=inputs["dv_struct"].dtype)
        # Place struct dvs in full array
        full_dv_array[self.struct_dvs] = inputs["dv_struct"]
        # Place mass dvs (if they were defined) in full array
        for dv_name in self.mass_dvs:
            full_dv_array[self.mass_dvs[dv_name]] = inputs[dv_name]
        # Slice full array with src_indices to get distributed dv array
        outputs["tacs_dvs"] = full_dv_array[self.src_indices]

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        tot_ndv = len(self.struct_dvs) + len(self.mass_dvs)
        dfull_dv_array = np.zeros(tot_ndv, dtype=inputs["dv_struct"].dtype)
        if mode == "fwd":
            if "tacs_dvs" in d_outputs:
                if "dv_struct" in d_inputs:
                    dfull_dv_array[self.struct_dvs] += d_inputs["dv_struct"]
                for dv_name in self.mass_dvs:
                    if dv_name in d_inputs:
                        dfull_dv_array[self.mass_dvs[dv_name]] += d_inputs[dv_name]
                d_outputs["tacs_dvs"] += dfull_dv_array[self.src_indices]
        else:  # mode == 'rev'
            if "tacs_dvs" in d_outputs:
                dfull_dv_array[self.src_indices] += d_outputs["tacs_dvs"]
                if "dv_struct" in d_inputs:
                    d_inputs["dv_struct"] += self.comm.allreduce(
                        dfull_dv_array[self.struct_dvs]
                    )
                for dv_name in self.mass_dvs:
                    if dv_name in d_inputs:
                        d_inputs[dv_name] += self.comm.allreduce(
                            dfull_dv_array[self.mass_dvs[dv_name]]
                        )

    def get_dv_src_indices(self):
        """
        Method to get src_indices on each processor
        for tacs distributed design variable vec
        """
        if MPI is not None and self.comm.size > 1:
            local_ndvs = self.fea_assembler.getNumDesignVars()
            all_proc_ndvs = self.comm.gather(local_ndvs, root=0)
            all_proc_indices = []
            if self.comm.rank == 0:
                tot_ndvs = 0
                for proc_i in range(self.comm.size):
                    local_ndvs = all_proc_ndvs[proc_i]
                    proc_indices = np.arange(tot_ndvs, tot_ndvs + local_ndvs)
                    all_proc_indices.append(proc_indices)
                    tot_ndvs += local_ndvs
            local_dv_indices = self.comm.scatter(all_proc_indices, root=0)
            return local_dv_indices
        else:
            ndvs = len(self.options["initial_dv_vals"])
            all_dv_indices = np.arange(ndvs)
            return all_dv_indices


class TacsPrecouplingGroup(om.Group):
    def initialize(self):
        self.options.declare(
            "fea_assembler",
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

    def setup(self):
        # Promote state variables/rhs with physics-specific tag that MPhys expects
        promotes_inputs = ["*"]

        fea_assembler = self.options["fea_assembler"]
        initial_dv_vals = self.options["initial_dv_vals"]
        separate_mass_dvs = self.options["separate_mass_dvs"]

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
            "x_struct0", shape=nnodes * 3, tags=["mphys_coordinates"]
        )
        unmask_input = MaskedVariableDescription(
            "x_struct0_masked", shape=(nnodes - nmult) * 3, tags=["mphys_coordinates"]
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
            "unmasker", unmasker, promotes_inputs=[("x_struct0_masked", "x_struct0")]
        )


class TacsSolver(om.ImplicitComponent):
    """
    Component to perform TACS steady analysis

    Assumptions:
        - The TACS steady residual is R = K * u_s - f_s = 0

    """

    def initialize(self):
        self.options.declare("fea_assembler", recordable=False)
        self.options.declare("conduction", default=False)
        self.options.declare("check_partials")
        self.options.declare("coupled", default=False)

        self.fea_assembler = None

        self.transposed = False
        self.check_partials = False

        self.old_dvs = None
        self.old_xs = None

    def setup(self):
        self.check_partials = self.options["check_partials"]
        self.fea_assembler = self.options["fea_assembler"]
        self.conduction = self.options["conduction"]
        self.coupled = self.options["coupled"]

        if self.conduction:
            self.states_name = "T_conduct"
            self.rhs_name = "q_conduct"
        else:
            self.states_name = "u_struct"
            self.rhs_name = "f_struct"

        # OpenMDAO setup
        local_ndvs = self.fea_assembler.getNumDesignVars()
        self.ndof = self.fea_assembler.getVarsPerNode()
        state_size = self.fea_assembler.getNumOwnedNodes() * self.ndof

        # inputs
        self.add_input(
            "tacs_dvs",
            distributed=True,
            shape=local_ndvs,
            desc="tacs distributed design variables",
            tags=["mphys_coupling"],
        )
        self.add_input(
            "x_struct0",
            distributed=True,
            shape_by_conn=True,
            desc="distributed structural node coordinates",
            tags=["mphys_coordinates"],
        )
        if self.coupled:
            self.add_input(
                self.rhs_name,
                distributed=True,
                shape=state_size,
                val=0.0,
                desc="coupling load vector",
                tags=["mphys_coupling"],
            )

        # outputs
        # its important that we set this to zero since this displacement value is used for the first iteration of the aero
        self.add_output(
            self.states_name,
            distributed=True,
            shape=state_size,
            val=np.zeros(state_size),
            desc="structural state vector",
            tags=["mphys_coupling"],
        )

    def _need_update(self, inputs):
        update = False

        dvs = inputs["tacs_dvs"]
        xs = inputs["x_struct0"]

        if self.old_dvs is None:
            self.old_dvs = inputs["tacs_dvs"].copy()
            update = True

        elif len(dvs) > 0:
            if max(np.abs(dvs - self.old_dvs)) > 0.0:  # 1e-7:
                self.old_dvs = inputs["tacs_dvs"].copy()
                update = True

        if self.old_xs is None:
            self.old_xs = inputs["x_struct0"].copy()
            update = True

        elif len(xs) > 0:
            if max(np.abs(xs - self.old_xs)) > 0.0:  # 1e-7:
                self.old_xs = inputs["x_struct0"].copy()
                update = True

        tmp = update
        # Perform all reduce to check if any other procs came back True
        update = self.comm.allreduce(tmp)
        return update

    def _update_internal(self, inputs, outputs=None):
        if self._need_update(inputs):
            self.sp.setDesignVars(inputs["tacs_dvs"])
            self.sp.setNodes(inputs["x_struct0"])
        if outputs is not None:
            self.sp.setVariables(outputs[self.states_name])
        self.sp._updateAssemblerVars()

    def apply_nonlinear(self, inputs, outputs, residuals):
        self._update_internal(inputs, outputs)

        if self.coupled:
            Fext = inputs[self.rhs_name]
        else:
            Fext = None

        self.sp.getResidual(res=residuals[self.states_name], Fext=Fext)

    def solve_nonlinear(self, inputs, outputs):
        self._update_internal(inputs)

        if self.coupled:
            Fext = inputs[self.rhs_name]
        else:
            Fext = None

        hasConverged = self.sp.solve(Fext=Fext)
        if not hasConverged:
            # TODO: In future we could add something here to distinguish between fatal failures and those that could be recovered from
            self.sp.zeroVariables()
            raise om.AnalysisError("TACS solver did not converge")
        self.sp.getVariables(states=outputs[self.states_name])

    def solve_linear(self, d_outputs, d_residuals, mode):
        if mode == "fwd":
            if self.check_partials:
                print("solver fwd")
            else:
                raise ValueError("forward mode requested but not implemented")

        if mode == "rev":
            self.sp.solveAdjoint(
                d_outputs[self.states_name], d_residuals[self.states_name]
            )

    def apply_linear(self, inputs, outputs, d_inputs, d_outputs, d_residuals, mode):
        self._update_internal(inputs, outputs)
        if mode == "fwd":
            if not self.check_partials:
                raise ValueError("TACS forward mode requested but not implemented")

        if mode == "rev":
            if self.states_name in d_residuals:
                if self.states_name in d_outputs:
                    self.sp.addTransposeJacVecProduct(
                        d_residuals[self.states_name], d_outputs[self.states_name]
                    )

                if self.rhs_name in d_inputs:
                    array_w_bcs = d_residuals[self.states_name].copy()
                    self.fea_assembler.applyBCsToVec(array_w_bcs)
                    d_inputs[self.rhs_name] -= array_w_bcs

                if "x_struct0" in d_inputs:
                    self.sp.addAdjointResXptSensProducts(
                        [d_residuals[self.states_name]],
                        [d_inputs["x_struct0"]],
                        scale=1.0,
                    )

                if "tacs_dvs" in d_inputs:
                    self.sp.addAdjointResProducts(
                        [d_residuals[self.states_name]],
                        [d_inputs["tacs_dvs"]],
                        scale=1.0,
                    )

    def _design_vector_changed(self, x):
        if self.x_save is None:
            self.x_save = x.copy()
            return True
        elif not np.allclose(x, self.x_save, rtol=1e-10, atol=1e-10):
            self.x_save = x.copy()
            return True
        else:
            return False

    def set_sp(self, sp):
        self.sp = sp


# All TACS function types that should be included under mass funcs group
MASS_FUNCS_CLASSES = [
    functions.StructuralMass,
    functions.CenterOfMass,
    functions.MomentOfInertia,
]


class TacsFunctions(om.ExplicitComponent):
    """
    Component to compute non-mass TACS functions
    """

    def initialize(self):
        self.options.declare("fea_assembler", recordable=False)
        self.options.declare("conduction", default=False)
        self.options.declare("check_partials")
        self.options.declare("write_solution")

        self.fea_assembler = None

        self.check_partials = False

    def setup(self):
        self.fea_assembler = self.options["fea_assembler"]
        self.check_partials = self.options["check_partials"]
        self.write_solution = self.options["write_solution"]
        self.conduction = self.options["conduction"]
        self.solution_counter = 0

        if self.conduction:
            self.states_name = "T_conduct"
        else:
            self.states_name = "u_struct"

        # TACS part of setup
        local_ndvs = self.fea_assembler.getNumDesignVars()

        # OpenMDAO part of setup
        self.add_input(
            "tacs_dvs",
            distributed=True,
            shape=local_ndvs,
            desc="tacs design variables",
            tags=["mphys_coupling"],
        )
        self.add_input(
            "x_struct0",
            distributed=True,
            shape_by_conn=True,
            desc="structural node coordinates",
            tags=["mphys_coordinates"],
        )
        self.add_input(
            self.states_name,
            distributed=True,
            shape_by_conn=True,
            desc="structural state vector",
            tags=["mphys_coupling"],
        )

    def mphys_set_sp(self, sp):
        # this is the external function to set the sp to this component
        self.sp = sp

        # Add eval funcs as outputs
        for func_name in self.sp.functionList:
            func_handle = self.sp.functionList[func_name]
            # Skip any mass functions
            if type(func_handle) not in MASS_FUNCS_CLASSES:
                self.add_output(
                    func_name, distributed=False, shape=1, tags=["mphys_result"]
                )

    def _update_internal(self, inputs):
        self.sp.setDesignVars(inputs["tacs_dvs"])
        self.sp.setNodes(inputs["x_struct0"])
        self.sp.setVariables(inputs[self.states_name])

    def compute(self, inputs, outputs):
        self._update_internal(inputs)

        # Evaluate functions
        funcs = {}
        self.sp.evalFunctions(funcs, evalFuncs=outputs.keys())
        for func_name in outputs:
            # Add struct problem name from key
            key = self.sp.name + "_" + func_name
            outputs[func_name] = funcs[key]

        if self.write_solution:
            # write the solution files.
            self.sp.writeSolution(number=self.solution_counter)
            self.solution_counter += 1

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        if mode == "fwd":
            if not self.check_partials:
                raise ValueError("TACS forward mode requested but not implemented")
        if mode == "rev":
            # always update internal because same tacs object could be used by multiple scenarios
            # and we need to load this scenario's state back into TACS before doing derivatives
            self._update_internal(inputs)

            for func_name in d_outputs:
                d_func = d_outputs[func_name]

                if "tacs_dvs" in d_inputs:
                    self.sp.addDVSens([func_name], [d_inputs["tacs_dvs"]], scale=d_func)

                if "x_struct0" in d_inputs:
                    self.sp.addXptSens(
                        [func_name], [d_inputs["x_struct0"]], scale=d_func
                    )

                if self.states_name in d_inputs:
                    sv_sens = np.zeros_like(d_inputs[self.states_name])
                    self.sp.addSVSens([func_name], [sv_sens])
                    d_inputs[self.states_name][:] += sv_sens * d_func


class MassFunctions(om.ExplicitComponent):
    """
    Component to compute TACS mass-specific functions
    """

    def initialize(self):
        self.options.declare("fea_assembler", recordable=False)
        self.options.declare("check_partials")

        self.fea_assembler = None
        self.check_partials = False

    def setup(self):
        self.fea_assembler = self.options["fea_assembler"]
        self.check_partials = self.options["check_partials"]

        # TACS part of setup
        local_ndvs = self.fea_assembler.getNumDesignVars()

        # OpenMDAO part of setup
        self.add_input(
            "tacs_dvs",
            distributed=True,
            shape=local_ndvs,
            desc="tacs design variables",
            tags=["mphys_coupling"],
        )
        self.add_input(
            "x_struct0",
            distributed=True,
            shape_by_conn=True,
            desc="structural node coordinates",
            tags=["mphys_coordinates"],
        )

    def mphys_set_sp(self, sp):
        # this is the external function to set the sp to this component
        self.sp = sp

        # Add eval funcs as outputs
        for func_name in self.sp.functionList:
            func_handle = self.sp.functionList[func_name]
            # Only include mass functions
            if type(func_handle) in MASS_FUNCS_CLASSES:
                self.add_output(
                    func_name, distributed=False, shape=1, tags=["mphys_result"]
                )

    def _update_internal(self, inputs):
        self.sp.setDesignVars(inputs["tacs_dvs"])
        self.sp.setNodes(inputs["x_struct0"])

    def compute(self, inputs, outputs):
        self._update_internal(inputs)

        # Evaluate functions
        funcs = {}
        self.sp.evalFunctions(funcs, evalFuncs=outputs.keys())
        for func_name in outputs:
            # Add struct problem name from key
            key = self.sp.name + "_" + func_name
            outputs[func_name] = funcs[key]

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        if mode == "fwd":
            if not self.check_partials:
                raise ValueError("TACS forward mode requested but not implemented")
        if mode == "rev":
            # always update internal because same tacs object could be used by multiple scenarios
            # and we need to load this scenario's state back into TACS before doing derivatives
            self._update_internal(inputs)

            for func_name in d_outputs:
                d_func = d_outputs[func_name]

                if "tacs_dvs" in d_inputs:
                    self.sp.addDVSens([func_name], [d_inputs["tacs_dvs"]], scale=d_func)

                if "x_struct0" in d_inputs:
                    self.sp.addXptSens(
                        [func_name], [d_inputs["x_struct0"]], scale=d_func
                    )


class TacsCouplingGroup(om.Group):
    def initialize(self):
        self.options.declare("fea_assembler", recordable=False)
        self.options.declare("check_partials")
        self.options.declare("conduction", default=False)
        self.options.declare("coupled", default=False)
        self.options.declare("scenario_name", default=None)
        self.options.declare("problem_setup", default=None)

    def setup(self):
        self.fea_assembler = self.options["fea_assembler"]
        self.check_partials = self.options["check_partials"]
        self.coupled = self.options["coupled"]
        self.conduction = self.options["conduction"]

        # Promote state variables/rhs with physics-specific tag that MPhys expects
        promotes_inputs = [
            ("x_struct0", "unmasker.x_struct0"),
            ("tacs_dvs", "distributor.tacs_dvs"),
        ]
        if self.conduction:
            self.states_name = "T_conduct"
            self.rhs_name = "q_conduct"
        else:
            self.states_name = "u_struct"
            self.rhs_name = "f_struct"

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
            unmask_output = MaskedVariableDescription(
                self.rhs_name, shape=nnodes * vpn, tags=["mphys_coupling"]
            )
            unmask_input = MaskedVariableDescription(
                self.rhs_name + "_masked",
                shape=(nnodes - nmult) * vpn,
                tags=["mphys_coupling"],
            )
            unmasker = UnmaskedConverter(
                input=unmask_input,
                output=unmask_output,
                mask=mask.flatten(),
                distributed=True,
            )
            self.add_subsystem(
                "unmasker",
                unmasker,
                promotes_inputs=[(self.rhs_name + "_masked", self.rhs_name)],
            )

        self.add_subsystem(
            "solver",
            TacsSolver(
                fea_assembler=self.fea_assembler,
                check_partials=self.check_partials,
                coupled=self.coupled,
                conduction=self.conduction,
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


class TacsFuncsGroup(om.Group):
    def initialize(self):
        self.options.declare("fea_assembler", recordable=False)
        self.options.declare("check_partials")
        self.options.declare("conduction", default=False)
        self.options.declare("scenario_name", default=None)
        self.options.declare("problem_setup", default=None)
        self.options.declare("write_solution")

    def setup(self):
        self.fea_assembler = self.options["fea_assembler"]
        self.check_partials = self.options["check_partials"]
        self.conduction = self.options["conduction"]
        self.write_solution = self.options["write_solution"]

        # Setup problem based on scenario that's calling builder
        scenario_name = self.options["scenario_name"]
        if scenario_name is None:
            # Default structural problem
            name = "default"
        else:
            name = scenario_name
        sp = self.fea_assembler.createStaticProblem(name=name)

        # Setup TACS problem with user-defined output functions
        problem_setup = self.options["problem_setup"]
        if problem_setup is not None:
            new_problem = problem_setup(scenario_name, self.fea_assembler, sp)
            # Check if the user provided back a new problem to overwrite the default
            if new_problem is not None:
                sp = new_problem

        # Promote state variables with physics-specific tag that MPhys expects
        promotes_inputs = [
            ("x_struct0", "unmasker.x_struct0"),
            ("tacs_dvs", "distributor.tacs_dvs"),
        ]
        if self.conduction:
            promotes_states = [("T_conduct", "solver.T_conduct")]
        else:
            promotes_states = [("u_struct", "solver.u_struct")]

        # Add function evaluation component for non-mass outputs
        self.add_subsystem(
            "eval_funcs",
            TacsFunctions(
                fea_assembler=self.fea_assembler,
                check_partials=self.check_partials,
                conduction=self.conduction,
                write_solution=self.write_solution,
            ),
            promotes_inputs=promotes_inputs + promotes_states,
            promotes_outputs=["*"],
        )
        self.eval_funcs.mphys_set_sp(sp)

        # Check if there are any mass functions added by user
        mass_funcs = False
        for func_handle in sp.functionList.values():
            if type(func_handle) in MASS_FUNCS_CLASSES:
                mass_funcs = True

        # Mass functions are handled in a seperate component to prevent useless adjoint solves
        if mass_funcs:
            # Note: these functions do not depend on the states
            self.add_subsystem(
                "mass_funcs",
                MassFunctions(
                    fea_assembler=self.fea_assembler, check_partials=self.check_partials
                ),
                promotes_inputs=promotes_inputs,
                promotes_outputs=["*"],
            )
            self.mass_funcs.mphys_set_sp(sp)


class TacsBuilder(Builder):
    def __init__(
        self,
        options,
        check_partials=False,
        conduction=False,
        coupled=True,
        write_solution=True,
        separate_mass_dvs=False,
    ):
        self.options = copy.deepcopy(options)
        self.check_partials = check_partials
        # Flag to switch to tacs conduction solver (False->structural)
        self.conduction = conduction
        # Flag to turn on f5 file writer
        self.write_solution = write_solution
        # Flag to turn on coupling variables
        self.coupled = coupled
        # Flag to separate point mass dvs from struct dvs in openmdao input array
        self.separate_mass_dvs = separate_mass_dvs

    def initialize(self, comm):
        pytacs_options = copy.deepcopy(self.options)
        bdf_file = pytacs_options.pop("mesh_file")

        # Load optional user-defined callback function for setting up tacs elements
        if "assembler_setup" in pytacs_options:
            assembler_setup = pytacs_options.pop("assembler_setup")
        else:
            assembler_setup = None

        # Load optional user-defined callback function for setting up tacs elements
        if "element_callback" in pytacs_options:
            element_callback = pytacs_options.pop("element_callback")
        else:
            element_callback = None

        # Load optional user-defined callback function for setting up tacs elements
        if "problem_setup" in pytacs_options:
            self.problem_setup = pytacs_options.pop("problem_setup")
        else:
            self.problem_setup = None

        # Create pytacs instance
        self.fea_assembler = pyTACS(bdf_file, options=pytacs_options, comm=comm)
        self.comm = comm

        # Do any pre-initialize setup requested by user
        if assembler_setup is not None:
            assembler_setup(self.fea_assembler)

        # Set up elements and TACS assembler
        self.fea_assembler.initialize(element_callback)

    def get_coupling_group_subsystem(self, scenario_name=None):
        return TacsCouplingGroup(
            fea_assembler=self.fea_assembler,
            conduction=self.conduction,
            check_partials=self.check_partials,
            coupled=self.coupled,
            scenario_name=scenario_name,
            problem_setup=self.problem_setup,
        )

    def get_mesh_coordinate_subsystem(self, scenario_name=None):
        return TacsMeshGroup(fea_assembler=self.fea_assembler)

    def get_pre_coupling_subsystem(self, scenario_name=None):
        initial_dvs = self.get_initial_dvs()
        return TacsPrecouplingGroup(
            fea_assembler=self.fea_assembler,
            initial_dv_vals=initial_dvs,
            separate_mass_dvs=self.separate_mass_dvs,
        )

    def get_post_coupling_subsystem(self, scenario_name=None):
        return TacsFuncsGroup(
            fea_assembler=self.fea_assembler,
            check_partials=self.check_partials,
            conduction=self.conduction,
            write_solution=self.write_solution,
            scenario_name=scenario_name,
            problem_setup=self.problem_setup,
        )

    def get_ndof(self):
        return self.fea_assembler.getVarsPerNode()

    def get_number_of_nodes(self):
        """
        Get the number of nodes on this processor,
        not including lagrange multiplier nodes
        """
        nnodes = self.fea_assembler.getNumOwnedNodes()
        nmult = self.fea_assembler.getNumOwnedMultiplierNodes()
        return nnodes - nmult

    def get_initial_dvs(self):
        """
        Get an array holding all dvs values that have been added to TACS
        """
        if MPI is not None and self.comm.size > 1:
            # Get DVs locally owned by this processor
            local_dvs = self.fea_assembler.getOrigDesignVars()
            local_dvs = local_dvs.astype(float)
            # Size of design variable on this processor
            local_ndvs = self.fea_assembler.getNumDesignVars()
            # Size of design variable vector on each processor
            dv_sizes = self.comm.allgather(local_ndvs)
            # Offsets for global design variable vector
            offsets = np.zeros(self.comm.size, dtype=int)
            offsets[1:] = np.cumsum(dv_sizes)[:-1]
            # Gather the portions of the design variable array distributed across each processor
            tot_ndvs = sum(dv_sizes)
            global_dvs = np.zeros(tot_ndvs, dtype=local_dvs.dtype)
            self.comm.Allgatherv(local_dvs, [global_dvs, dv_sizes, offsets, MPI.DOUBLE])
            # return the global dv array
            return global_dvs
        else:
            return self.fea_assembler.getOrigDesignVars()

    def get_ndv(self):
        """
        Get total number of structural design variables across all procs
        """
        return self.fea_assembler.getTotalNumDesignVars()

    def get_solver(self):
        # this method is only used by the RLT transfer scheme
        return self.fea_assembler.assembler

    def get_fea_assembler(self):
        """
        Returns underlying pytacs object.
        """
        return self.fea_assembler
