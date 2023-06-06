import openmdao.api as om


class ConstraintComponent(om.ExplicitComponent):
    """
    Component to compute TACS constraint functions
    """

    def initialize(self):
        self.options.declare("fea_assembler", recordable=False)
        self.options.declare("constraint_object")

        self.fea_assembler = None
        self.constr = None

    def setup(self):
        self.fea_assembler = self.options["fea_assembler"]
        self.constr = self.options["constraint_object"]

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

        # Add eval funcs as outputs
        con_names = self.constr.getConstraintKeys()
        con_sizes = {}
        self.constr.getConstraintSizes(con_sizes)
        for con_name in con_names:
            con_key = f"{self.constr.name}_{con_name}"
            ncon = con_sizes[con_key]
            self.add_output(
                con_name, distributed=False, shape=ncon, tags=["mphys_result"]
            )

    def _update_internal(self, inputs):
        self.constr.setDesignVars(inputs["tacs_dvs"])
        self.constr.setNodes(inputs["x_struct0"])

    def compute(self, inputs, outputs):
        self._update_internal(inputs)

        # Evaluate functions
        funcs = {}
        self.constr.evalConstraints(funcs, evalCons=outputs.keys())
        for con_name in outputs:
            # Add struct problem name from key
            key = self.constr.name + "_" + con_name
            outputs[con_name] = funcs[key]

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        # always update internal because same tacs object could be used by multiple scenarios
        # and we need to load this scenario's state back into TACS before doing derivatives
        self._update_internal(inputs)
        funcs_sens = {}
        self.constr.evalConstraintsSens(funcs_sens)

        for out_name in d_outputs:
            output_key = f"{self.constr.name}_{out_name}"
            Jdv = funcs_sens[output_key]["struct"].astype(float)
            Jxpt = funcs_sens[output_key]["Xpts"].astype(float)

            if mode == "fwd":
                if "tacs_dvs" in d_inputs:
                    out = Jdv.dot(d_inputs["tacs_dvs"])
                    d_outputs[out_name] += self.comm.allreduce(out)

                if "x_struct0" in d_inputs:
                    out = Jxpt.dot(d_inputs["x_struct0"])
                    d_outputs[out_name] += self.comm.allreduce(out)

            elif mode == "rev":
                if "tacs_dvs" in d_inputs:
                    d_inputs["tacs_dvs"] += Jdv.T.dot(d_outputs[out_name])

                if "x_struct0" in d_inputs:
                    d_inputs["x_struct0"] += Jxpt.T.dot(d_outputs[out_name])
