import numpy as np

import openmdao.api as om


class TacsBuckling(om.ExplicitComponent):
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

    def mphys_set_bp(self, bp):
        # this is the external function to set the bp to this component
        self.bp = bp

        # Add eval funcs as outputs
        for func_name in self.bp.functionList:
            func_name = func_name.replace(".", "_")
            self.add_output(
                func_name, distributed=False, shape=1, tags=["mphys_result"]
            )

    def _update_internal(self, inputs):
        self.bp.setDesignVars(inputs["tacs_dvs"])
        self.bp.setNodes(inputs["x_struct0"])

    def compute(self, inputs, outputs):
        self._update_internal(inputs)

        # Solve
        self.bp.solve(u0=inputs[self.states_name])

        # Evaluate functions
        funcs = {}
        self.bp.evalFunctions(funcs, evalFuncs=outputs.keys())
        for out_name in outputs:
            eig_name = out_name.replace("_", ".")
            func_key = f"{self.bp.name}_{eig_name}"
            self.bp.evalFunctions(funcs, evalFuncs=[eig_name])
            outputs[out_name] = funcs[func_key]

        if self.write_solution:
            # write the solution files.
            self.bp.writeSolution(number=self.solution_counter)
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
                eig_name = func_name.replace("_", ".")
                mode_i = self.bp.functionList[eig_name]

                if d_func[0] != 0.0:
                    if "tacs_dvs" in d_inputs:
                        self.bp.addDVSens(
                            [mode_i], [d_inputs["tacs_dvs"]], scale=d_func
                        )

                    if "x_struct0" in d_inputs:
                        self.bp.addXptSens(
                            [mode_i], [d_inputs["x_struct0"]], scale=d_func
                        )

                    if self.states_name in d_inputs:
                        sv_sens = np.zeros_like(d_inputs[self.states_name])
                        self.bp.evalSVSens([mode_i], [sv_sens])
                        d_inputs[self.states_name][:] += sv_sens * d_func
