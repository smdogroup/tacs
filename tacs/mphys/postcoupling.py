import openmdao.api as om

import tacs.problems
from .functions import TacsFunctions, MassFunctions, MASS_FUNCS_CLASSES
from .buckling import TacsBuckling
from .constraints import ConstraintComponent


class TacsPostcouplingGroup(om.Group):
    def initialize(self):
        self.options.declare("fea_assembler", recordable=False)
        self.options.declare("check_partials")
        self.options.declare("conduction", default=False)
        self.options.declare("scenario_name", default=None)
        self.options.declare("problem_setup", default=None)
        self.options.declare("write_solution")
        self.options.declare("constraint_setup", default=None)
        self.options.declare("buckling_setup", default=None)

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

        # Mass functions are handled in a separate component to prevent useless adjoint solves
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

        # Setup TACS problem with user-defined output functions
        buckling_setup = self.options["buckling_setup"]
        if buckling_setup is not None:
            new_problem = buckling_setup(scenario_name, self.fea_assembler)
            # Check if the user provided back a new problem to overwrite the default
            if isinstance(new_problem, tacs.problems.BucklingProblem):
                bp = new_problem

                # Add buckling evaluation component for eigenvalue outputs
                self.add_subsystem(
                    "buckling",
                    TacsBuckling(
                        fea_assembler=self.fea_assembler,
                        check_partials=self.check_partials,
                        conduction=self.conduction,
                        write_solution=self.write_solution,
                    ),
                    promotes_inputs=promotes_inputs + promotes_states,
                    promotes_outputs=["*"],
                )
                self.buckling.mphys_set_bp(bp)

        # Check if there are any user-defined TACS constraints
        # Constraints behave similar to "mass" functions (i.e. they don't depend on the solution state)
        constraint_setup = self.options["constraint_setup"]
        tacs_constraints = []
        if constraint_setup is not None:
            new_constraints = constraint_setup(
                scenario_name, self.fea_assembler, tacs_constraints
            )
            # Check if the user provided back new constraints to overwrite the default
            if new_constraints is not None:
                tacs_constraints = new_constraints

        # Only add constraint group if there are constraints to add
        if len(tacs_constraints) > 0:
            con_group = self.add_subsystem("constraints", om.Group(), promotes=["*"])
            # Loop through each constraint in lista and add to group
            for constraint in tacs_constraints:
                con_comp = ConstraintComponent(
                    fea_assembler=self.fea_assembler,
                    constraint_object=constraint,
                )
                con_group.add_subsystem(
                    constraint.name, con_comp, promotes_inputs=promotes_inputs
                )
