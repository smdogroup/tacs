import numpy as np
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

        self.old_dvs = None
        self.old_xs = None

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
        dvsNeedUpdate, xsNeedUpdate = self._need_update(inputs)
        if dvsNeedUpdate:
            self.constr.setDesignVars(inputs["tacs_dvs"])
        if xsNeedUpdate:
            self.constr.setNodes(inputs["x_struct0"])

    def _need_update(self, inputs):
        """Checks whether the design variables or coordinates being passed
        in by OpenMDAO are different from those currently stored in TACS

        Parameters
        ----------
        inputs : OpenMDAO input vector
            _description_

        Returns
        -------
        (bool, bool)
            Whether the design variables or coordinates need to be updated
            respectively
        """
        dvsNeedUpdate = False
        xsNeedUpdate = False

        dvs = inputs["tacs_dvs"]
        xs = inputs["x_struct0"]

        if self.old_dvs is None:
            self.old_dvs = inputs["tacs_dvs"].copy()
            dvsNeedUpdate = True

        elif len(dvs) > 0:
            if max(np.abs(dvs - self.old_dvs)) > 0.0:  # 1e-7:
                self.old_dvs = inputs["tacs_dvs"].copy()
                dvsNeedUpdate = True

        if self.old_xs is None:
            self.old_xs = inputs["x_struct0"].copy()
            xsNeedUpdate = True

        elif len(xs) > 0:
            if max(np.abs(xs - self.old_xs)) > 0.0:  # 1e-7:
                self.old_xs = inputs["x_struct0"].copy()
                xsNeedUpdate = True

        tmp1 = dvsNeedUpdate
        tmp2 = xsNeedUpdate
        # Perform all reduce to check if any other procs came back True
        dvsNeedUpdate = self.comm.allreduce(tmp1)
        xsNeedUpdate = self.comm.allreduce(tmp2)
        return dvsNeedUpdate, xsNeedUpdate

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
