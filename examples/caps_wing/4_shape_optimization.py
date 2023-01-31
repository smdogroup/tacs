"""
Sean Engelstad, October 2022
GT SMDO Lab, Dr. Graeme Kennedy
Caps to TACS example
"""


import os, numpy as np, matplotlib.pyplot as plt
from tacs import functions, caps2tacs
from mpi4py import MPI
import openmdao, openmdao.api as om

class Caps2TacsComponent(om.ExplicitComponent):
    def initialize(self):
        """
        Declare the capsProblem to the openMdao component
        """
        # CAPS Problem input
        self.options.declare("tacs_aim", types=object)

    def setup(self):
        """
        Declare the inputs and outputs of the caps problem to openmdao
        """
        tacs_aim_wrapper = self.options["tacs_aim"]
        tacs_aim = tacs_aim_wrapper.aim

        # attach parameters to the openMdao object
        for shape_var in tacs_aim_wrapper.shape_variables:
            self.add_input(shape_var.name, val=tacs_aim.geometry.despmtr[shape_var.name].value)

        # add output variables
        self._func_names = ["mass"]
        for func_name in self._func_names:
            self.add_output(func_name)

        # private metadata variables
        self._first_design = True
        self._func_history = {func_name: [] for func_name in self._func_names}
        self._design_hdl = open(
            os.path.join(tacs_aim.analysisDir, "design_hist.txt"), "w"
        )

    def compute(self, inputs, outputs):
        """
        Compute the objective function
        """

        # obtain the aim from openmdao storage
        tacs_aim_wrapper = self.options["tacs_aim"]
        tacs_aim = tacs_aim_wrapper.aim

        # Update input values
        design_change = False
        for shape_var in tacs_aim_wrapper.shape_variables:
            if tacs_aim.geometry.despmtr[shape_var.name].value != inputs[shape_var.name]:
                design_change = True
                tacs_aim.geometry.despmtr[shape_var.name].value = inputs[shape_var.name]

        if not (design_change) and self._first_design:
            design_change = True
            self._first_design = False

        if design_change:
            self._design_hdl.write("New Design...\n")
            for shape_var in tacs_aim_wrapper.shape_variables:
                self._design_hdl.write(f"\tDV {shape_var.name} = {inputs[shape_var.name][0]}\n")
            self._design_hdl.flush()

        # forward analysis here
        if design_change:
            self._run_analysis(tacs_aim_wrapper)
            self._plot_history(
                tacs_aim, filename="opt_history.png"
            )
            self._write_history()

        # Grab objectives and attach as openmdao outputs
        for func_name in self._func_names:
            outputs[func_name] = self._functions[func_name]

    def _plot_history(self, tacs_aim, filename):

        directory = tacs_aim.analysisDir
        for func_name in self._func_names:
            self._func_history[func_name].append(tacs_aim.dynout[func_name].value)
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

        plot_filepath = os.path.join(directory, filename)
        plt.savefig(plot_filepath)
        plt.close("all")

    def _write_history(self):
        # write new function values to the design file
        self._design_hdl.write(f"Analysis result:\n")
        for func_name in self._func_names:
            self._design_hdl.write(
                f"\tfunc {func_name} = {self._func_history[func_name][-1]}\n"
            )
        self._design_hdl.write("\n")
        self._design_hdl.flush()

    def setup_partials(self):
        """
        declare partial derivatives
        """
        tacs_aim_wrapper = self.options["tacs_aim"]
        for func_name in self._func_names:
            for shape_var in tacs_aim_wrapper.shape_variables:
                self.declare_partials(func_name, shape_var.name)

    def compute_partials(self, inputs, partials):
        """
        the actual value of the partial derivatives are assigned here
        """

        tacs_aim_wrapper = self.options["tacs_aim"]
        tacs_aim = tacs_aim_wrapper.aim

        # Update input values
        design_change = False
        for shape_var in tacs_aim_wrapper.shape_variables:
            if tacs_aim.geometry.despmtr[shape_var.name].value != inputs[shape_var.name]:
                design_change = True
                tacs_aim.geometry.despmtr[shape_var.name].value = inputs[shape_var.name]

        if not (design_change) and self._first_design:
            design_change = True
            self._first_design = False

        if design_change:
            self._design_hdl.write("New Design...\n")
            for shape_var in tacs_aim_wrapper.shape_variables:
                self._design_hdl.write(f"\tDV {shape_var.name} = {inputs[shape_var.name][0]}\n")
            self._design_hdl.flush()

        if design_change:
            self._run_analysis(tacs_aim_wrapper)

        # Get derivatives and partials
        for func_name in self._func_names:
            for shape_var in tacs_aim_wrapper.shape_variables:
                partials[func_name, shape_var.name] = self._gradients[func_name][shape_var.name]

    def _run_analysis(self, tacs_aim_wrapper):
        """
        perform the TACS + ESP/CAPS analysis
        """

        tacs_aim_wrapper.pre_analysis()
        fea_solver = tacs_aim_wrapper.fea_solver.initialize()

        # create the SPs from the BDF file (Struct Problem)
        SPs = fea_solver.createTACSProbsFromBDF()
        for caseID in SPs:
            SPs[caseID].addFunction("mass", functions.StructuralMass)
            # SPs[caseID].addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5, ksWeight=50)

        # solve the forward and adjoint analysis for each struct problem
        for caseID in SPs:
            SPs[caseID].solve()
            SPs[caseID].write_sensitivity_file(
                evalFuncs=self._func_names,
                sens_file=tacs_aim_wrapper.sens_file_path,
            )

        # done running the write sens file operation
        # run the tacs AIM postAnalysis which reads the coordinate derivatives & computes geom DV derivatives
        tacs_aim_wrapper.post_analysis()

        # functions and gradients are stored in the tacs AIM dynout methods
        self._functions = tacs_aim_wrapper.get_functions(self._func_names)
        self._gradients = tacs_aim_wrapper.get_gradients(self._func_names, tacs_aim_wrapper.shape_variables)
        return


# --------------------------------------------------------------#
# Setup CAPS Problem
# --------------------------------------------------------------#
comm = MPI.COMM_WORLD
caps_struct = caps2tacs.CapsStruct.build(csm_file="simple_naca_wing.csm", comm=comm)
tacs_aim = caps_struct.tacs_aim
caps_struct.egads_aim.set_mesh(  # need a refined-enough mesh for the derivative test to pass
    edge_pt_min=15,
    edge_pt_max=20,
    global_mesh_size=0.1,
    max_surf_offset=0.01,
    max_dihedral_angle=5,
).register_to(
    tacs_aim
)

aluminum = caps2tacs.Isotropic.aluminum().register_to(tacs_aim)

# setup the thickness design variables + automatic shell properties
nribs = int(tacs_aim.get_config_parameter("nribs"))
nspars = int(tacs_aim.get_config_parameter("nspars"))
for irib in range(1, nribs + 1):
    caps2tacs.ShellProperty(
        caps_group=f"rib{irib}", material=aluminum, membrane_thickness=0.05
    ).register_to(tacs_aim)
for ispar in range(1, nspars + 1):
    caps2tacs.ShellProperty(
        caps_group=f"spar{ispar}", material=aluminum, membrane_thickness=0.05
    ).register_to(tacs_aim)
caps2tacs.ShellProperty(
    caps_group="OML", material=aluminum, membrane_thickness=0.03
).register_to(tacs_aim)

# register one shape variable rib_a1
caps2tacs.ShapeVariable("rib_a1").register_to(tacs_aim)
caps2tacs.ShapeVariable("rib_a2").register_to(tacs_aim)

# add constraints and loads
caps2tacs.PinConstraint("root").register_to(tacs_aim)
caps2tacs.GridForce("OML", direction=[0, 0, 1.0], magnitude=100).register_to(
    tacs_aim
)

# run the pre analysis to build tacs input files
tacs_aim.setup_aim()

# ---------------------------------------------------------------------------------#
# Setup OpenMDAO Problem
# ---------------------------------------------------------------------------------#

# setup the OpenMDAO problem object
om_problem = om.Problem()

# Create the OpenMDAO component
tacs_system = Caps2TacsComponent(tacs_aim=tacs_aim)
om_problem.model.add_subsystem("tacsSystem", tacs_system)

# setup the optimization
om_problem.driver = om.ScipyOptimizeDriver()
om_problem.driver.options["optimizer"] = "SLSQP"
om_problem.driver.options["tol"] = 1.0e-9
om_problem.driver.options["disp"] = True

# add design variables to the model
om_problem.model.add_design_var("tacsSystem.rib_a1", lower=0.6, upper=1.4)
om_problem.model.add_design_var("tacsSystem.rib_a2", lower=-0.3, upper=0.3)

# add objectives to the model
om_problem.model.add_objective("tacsSystem.mass")
# om_problem.model.add_objective('tacsSystem.ks_vmfailure')

# Start the optimization
print("\n==> Starting Optimization...")
om_problem.setup()
om_problem.run_driver()

tacs_system._design_hdl.write("--> Optimized values:\n")
rib_a1_value = om_problem.get_val("tacsSystem.rib_a1")
rib_a2_value = om_problem.get_val("tacsSystem.rib_a2")
tacs_system._design_hdl.write(f"\trib_a1 = {rib_a1_value}\n")
tacs_system._design_hdl.write(f"\trib_a2 = {rib_a2_value}\n")
# close the design hdl file and cleanup
tacs_system._design_hdl.close()
om_problem.cleanup()

