"""
Sean Engelstad, January 2023
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example
"""

# also requires pySNOPT
import os, numpy as np, matplotlib.pyplot as plt
import openmdao, openmdao.api as om
from tacs import functions, caps2tacs
from mpi4py import MPI


class AnalysisManager(om.ExplicitComponent):
    def initialize(self):
        """
        Declare the capsProblem to the openMdao component
        """
        # CAPS problem input
        self.options.declare("tacs_aim", types=object)
        self.options.declare("struct_problems", types=object)
        self.options.declare("write_f5", types=bool)

    def setup(self):
        """
        Declare the inputs and outputs of the capsProblem to openmdao
        """
        tacs_aim = self.options["tacs_aim"]

        # add the thickness variables as openmdao inputs, with starting uniform thickness
        for thick_var in tacs_aim.thickness_variables:
            self.add_input(thick_var.name, val=thick_var.value)

        # add output variables and functions
        self._func_names = ["mass", "ks_vmfailure"]
        for func_name in self._func_names:
            self.add_output(func_name)

        # first design boolean
        self._first_design = True

        # function histories
        self._func_history = {func_name: [] for func_name in self._func_names}

        # design history file
        self._design_hdl = open(
            os.path.join(tacs_aim.analysis_dir, "design_hist.txt"), "w"
        )

    def setup_partials(self):
        """
        declare partial derivatives
        """
        tacs_aim = self.options["tacs_aim"]

        for func_name in self._func_names:
            self.declare_partials(
                func_name,
                [thick_var.name for thick_var in tacs_aim.thickness_variables],
            )

    def compute(self, inputs, outputs):
        """
        compute the objective functions
        """
        # obtain the aim from openmdao storage
        tacs_aim = self.options["tacs_aim"]

        # update the input values
        design_change = False
        if not self._first_design:
            for caseID in self.SPs:
                xarray = self.SPs[
                    caseID
                ].x.getArray()  # gets the xarray of last Struct Problem
                for ithick, thick_var in enumerate(tacs_aim.thickness_variables):
                    if xarray[ithick].real != inputs[thick_var.name]:
                        design_change = True

        # count a first design as a design change
        if not (design_change) and self._first_design:
            design_change = True
            self._first_design = False

        # print new design
        if design_change:
            self._print_design(inputs)

            # run a forward + adjoint analysis
            self._run_analysis(inputs)

            # update func history and report to design file
            self._update_history()
            self._function_report()

        # Grab objectives and attach as openmdao outputs
        for func_name in self._func_names:
            outputs[func_name] = self._functions[func_name]

        return

    def _print_design(self, inputs):
        tacs_aim = self.options["tacs_aim"]
        self._design_hdl.write("New Design...\n")
        self._design_hdl.write(
            f"\tthick dvs = {[_.name for _ in tacs_aim.thickness_variables]}\n"
        )
        real_xarray = [float(inputs[key]) for key in inputs]
        self._design_hdl.write(f"\tvalues = {real_xarray}\n")
        self._design_hdl.flush()

    def _update_history(self):
        tacs_aim = self.options["tacs_aim"]
        for func_name in self._func_names:
            self._func_history[func_name].append(self._functions[func_name])
        self._plot_history(directory=tacs_aim.analysis_dir, filename="opt_history.png")

    def _function_report(self):
        self._design_hdl.write("Analysis result:\n")
        for func_name in self._func_names:
            self._design_hdl.write(
                f"\tfunc {func_name} = {self._func_history[func_name][-1]}\n"
            )
        self._design_hdl.write("\n")
        self._design_hdl.flush()

    def _plot_history(self, directory, filename):
        num_iterations = len(self._func_history[self._func_names[0]])
        iterations = [_ for _ in range(num_iterations)]
        styles = ["k-", "b-"]
        plt.figure()
        for ifunc, func_name in enumerate(self._func_names):
            yvec = self._func_history[func_name]
            if func_name == "mass":
                yvec /= max(np.array(yvec))
            plt.plot(iterations, yvec, styles[ifunc], linewidth=2, label=func_name)
        plt.legend()
        plt.xlabel("iterations")
        plt.ylabel("func values")
        plt.yscale("log")
        plot_filepath = os.path.join(directory, filename)
        plt.savefig(plot_filepath)
        plt.close("all")

    def compute_partials(self, inputs, partials):
        """
        the actual value of partial derivatives assigned here
        """
        tacs_aim = self.options["tacs_aim"]

        # Update input values
        design_change = False
        if not self._first_design:
            for caseID in self.SPs:  # get tacs DVs for last caseID Struct Problem
                xarray = self.SPs[caseID].x.getArray()
            for ithick, thick_var in enumerate(tacs_aim.thickness_variables):
                if xarray[ithick].real != inputs[thick_var.name]:
                    design_change = True

        if not (design_change) and self._first_design:
            design_change = True
            self._first_design = False

        if design_change:
            self._run_analysis(inputs)

        # Get derivatives and partials
        for func_name in self._func_names:
            for thick_var in tacs_aim.thickness_variables:
                partials[func_name, thick_var.name] = self._gradients[func_name][
                    thick_var.name
                ]

    def _run_analysis(self, inputs):
        """
        run forward and adjoint TACS analysis on fixed structural mesh geometry from tacsAIM
        """
        tacs_aim = self.options["tacs_aim"]
        write_f5 = self.options["write_f5"]

        print("running analysis...")

        self.SPs = tacs_aim.createTACSProbs()
        for caseID in self.SPs:
            self.SPs[caseID].addFunction("mass", functions.StructuralMass)
            self.SPs[caseID].addFunction(
                "ks_vmfailure", functions.KSFailure, safetyFactor=1.5, ksWeight=50
            )

        # solve the forward and adjoint analysis for each struct problem
        tacs_funcs = {}
        tacs_sens = {}
        for caseID in self.SPs:
            xarray = self.SPs[caseID].x.getArray()
            for ithick, thick_var in enumerate(tacs_aim.thickness_variables):
                xarray[ithick] = float(inputs[thick_var.name])

            self.SPs[caseID].solve()
            self.SPs[caseID].evalFunctions(tacs_funcs, evalFuncs=self._func_names)
            self.SPs[caseID].evalFunctionsSens(tacs_sens, evalFuncs=self._func_names)
            if write_f5:
                self.SPs[caseID].writeSolution(
                    baseName="tacs_output", outputDir=tacs_aim.analysis_dir
                )

        self._functions = {}
        self._gradients = {}
        for func_name in self._func_names:
            # corresponding tacs key for single loadset (key="loadset#" + "func_name")
            for tacs_key in tacs_funcs:
                if func_name in tacs_key:
                    break

            self._functions[func_name] = tacs_funcs[tacs_key].real
            self._gradients[func_name] = {}

            # grab structDV derivatives, not coordinate derivatives from tacs_sens
            struct_derivs = tacs_sens[tacs_key]["struct"]
            for ithick, thick_var in enumerate(tacs_aim.thickness_variables):
                self._gradients[func_name][thick_var.name] = struct_derivs[ithick].real

        print("done running analysis...")
        return


# --------------------------------------------------------------#
# Setup CAPS Problem
# --------------------------------------------------------------#
comm = MPI.COMM_WORLD
# can also switch to large_naca_wing.csm file here if you want and it will automatically update the DVs
caps_struct = caps2tacs.CapsStruct.build(csm_file="large_naca_wing.csm", comm=comm)
tacs_aim = caps_struct.tacs_aim
caps_struct.egads_aim.set_mesh(  # need a refined-enough mesh for the derivative test to pass
    edge_pt_min=15,
    edge_pt_max=20,
    global_mesh_size=0.01,
    max_surf_offset=0.01,
    max_dihedral_angle=5,
).register_to(
    tacs_aim
)

aluminum = caps2tacs.Isotropic.aluminum().register_to(tacs_aim)

# setup the thickness design variables + automatic shell properties
nribs = int(tacs_aim.get_config_parameter("nribs"))
nspars = int(tacs_aim.get_config_parameter("nspars"))
nOML = nribs - 1
for irib in range(1, nribs + 1):
    caps2tacs.ThicknessVariable(
        caps_group=f"rib{irib}", value=0.03, material=aluminum
    ).register_to(tacs_aim)
for ispar in range(1, nspars + 1):
    caps2tacs.ThicknessVariable(
        caps_group=f"spar{ispar}", value=0.03, material=aluminum
    ).register_to(tacs_aim)
for iOML in range(1, nOML + 1):
    caps2tacs.ThicknessVariable(
        caps_group=f"OML{iOML}", value=0.1, material=aluminum
    ).register_to(tacs_aim)

# add constraints and loads
caps2tacs.PinConstraint("root").register_to(tacs_aim)
caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=10).register_to(tacs_aim)

# run the pre analysis to build tacs input files
tacs_aim.setup_aim().pre_analysis()

# --------------------------------------------------------------------------#
# Setup OpenMDAO Problem
# --------------------------------------------------------------------------#

# setup the OpenMDAO Problem object
prob = om.Problem()

# Create the OpenMDAO component
tacs_system = AnalysisManager(tacs_aim=tacs_aim, write_f5=False)
prob.model.add_subsystem("tacsSystem", tacs_system)


# setup the optimizer settings # COBYLA for auto-FDing
prob.driver = om.ScipyOptimizeDriver(optimizer="SLSQP", tol=1.0e-9, disp=True)

# add design variables to the model
for thick_var in tacs_aim.thickness_variables:
    prob.model.add_design_var(
        f"tacsSystem.{thick_var.name}", lower=0.001, upper=0.5, scaler=100.0
    )

# add objectives & constraints to the model
prob.model.add_objective("tacsSystem.mass", scaler=1.0e-2)
prob.model.add_constraint("tacsSystem.ks_vmfailure", upper=0.267)

# Start the optimization
print("\n==> Starting optimization...")
prob.setup()

debug = False
if debug:
    print("Checking partials...", flush=True)
    prob.check_partials(compact_print=True)

else:
    prob.run_driver()

    # report the final optimal design
    design_hdl = tacs_system._design_hdl
    for thick_var in tacs_aim.thickness_variables:
        opt_value = prob.get_val(f"tacsSystem.{thick_var.name}")
        design_hdl.write(f"\t{thick_var.name} = {opt_value}")

    prob.cleanup()
    # close the design hdl file
    design_hdl.close()
