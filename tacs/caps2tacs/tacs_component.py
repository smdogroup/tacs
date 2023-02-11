__all__ = ["TacsStaticComponent"]

import os, numpy as np, matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import openmdao, openmdao.api as om
from .tacs_aim import TacsAim
from tacs import functions
from mpi4py import MPI


class TacsStaticComponent(om.ExplicitComponent):
    """
    OpenMDAO component for sizing and/or shape optimization of TACS static analysis problems
    TODO : future work could equivalent unsteady / modal analysis component
    """

    def initialize(self):
        """
        Declare the capsProblem to the openMdao component
        """
        self.options.declare(
            "tacs_aim", types=object
        )  # takes in the TacsAim wrapper class

    def setup(self):
        tacs_aim = self.options["tacs_aim"]
        assert isinstance(tacs_aim, TacsAim)

        # add the thickness variables as openmdao inputs, with starting uniform thickness
        for thick_var in tacs_aim.thickness_variables:
            self.add_input(thick_var.name, val=thick_var.value)

        # add the shape variables as openmdao inputs
        for shape_var in tacs_aim.shape_variables:
            self.add_input(shape_var.name, val=tacs_aim.get_shape_var_value(shape_var))

        # add output analysis functions
        assert (
            len(tacs_aim.analysis_functions) > 0
        )  # makes sure we have some analysis functions before running an analysis)
        for func in tacs_aim.analysis_functions:
            self.add_output(func.name)

        # first design boolean
        self._first_design = True

        # function histories
        self._func_history = {func_name: [] for func_name in tacs_aim.function_names}

        # design history file
        self._design_hdl = open(
            os.path.join(tacs_aim.analysis_dir, "design_hist.txt"), "w"
        )

    def setup_partials(self):
        tacs_aim = self.options["tacs_aim"]

        for func in tacs_aim.analysis_functions:
            self.declare_partials(func.name, tacs_aim.variables)

    def compute(self, inputs, outputs):
        """
        compute the objective functions
        """
        # obtain the aim from openmdao storage
        tacs_aim = self.options["tacs_aim"]

        # update the input values
        design_change = False
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

        # create the list of steady structural analysis problems in TACS
        # this one is to check the design variables at the start, recreated multiple times
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
        return
