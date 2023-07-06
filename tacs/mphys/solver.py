import numpy as np

import openmdao.api as om


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

        self.sp.solve(Fext=Fext)
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

    def set_sp(self, sp):
        self.sp = sp
