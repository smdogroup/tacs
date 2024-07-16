import numpy as np

import openmdao.api as om

from tacs import functions


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
        self.auto_write_solution = self.options["write_solution"]
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

    def write_solution(self):
        # write the solution files.
        self.sp.writeSolution(number=self.solution_counter)
        self.solution_counter += 1

    def compute(self, inputs, outputs):
        self._update_internal(inputs)

        # Evaluate functions
        funcs = {}
        self.sp.evalFunctions(funcs, evalFuncs=outputs.keys())
        for func_name in outputs:
            # Add struct problem name from key
            key = self.sp.name + "_" + func_name
            outputs[func_name] = funcs[key]

        if self.auto_write_solution:
            self.write_solution()

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
