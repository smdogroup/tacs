"""
Written by Sean Engelstad, GT SMDO Lab, 2022-2023
"""
__all__ = ["TacsStaticComponent"]

import os, numpy as np, matplotlib.pyplot as plt
import openmdao.api as om


class TacsStaticComponent(om.ExplicitComponent):
    """
    OpenMDAO component for sizing and/or shape optimization of TACS static analysis problems
    TODO : future work could equivalent unsteady / modal analysis component
    """

    def initialize(self):
        """
        Declare the capsProblem to the openMdao component
        Makes the __init__ construct : TacsStaticComponent(tacs_model, write_f5, track_history)
        """
        self.options.declare(
            "tacs_model", types=object
        )  # takes in the TacsAim wrapper class

        # whether to write f5 files of each iteration
        self.options.declare("write_f5", types=bool, default=True)

        # whether to write history files and plot history
        self.options.declare("track_history", types=bool, default=True)

        self._iteration = 0  # track iteration number

    def setup(self):
        tacs_model = self.options["tacs_model"]
        track_history = self.options["track_history"]

        # make sure we have more than zero tacs aim variables
        assert len(tacs_model.variables) > 0

        # add the thickness variables as openmdao inputs, with starting uniform thickness
        for thick_var in tacs_model.thickness_variables:
            self.add_input(thick_var.name, val=thick_var.value)

        # add the shape variables as openmdao inputs
        for shape_var in tacs_model.shape_variables:
            self.add_input(shape_var.name, val=shape_var.value)

        # add output analysis functions
        assert (
            len(tacs_model.analysis_functions) > 0
        )  # makes sure we have some analysis functions before running an analysis)
        for func in tacs_model.analysis_functions:
            self.add_output(func.name)

        if track_history:
            # function histories
            self._func_history = {
                func_name: [] for func_name in tacs_model.function_names
            }

            # design history file
            if tacs_model.root_proc:
                self._design_hdl = open(
                    os.path.join(
                        tacs_model.analysis_dir(tacs_model.root_proc_ind),
                        "design_hist.txt",
                    ),
                    "w",
                )

    def setup_partials(self):
        tacs_model = self.options["tacs_model"]

        for func in tacs_model.analysis_functions:
            for var in tacs_model.variables:
                self.declare_partials(func.name, var.name)

    def compute(self, inputs, outputs):
        """
        compute the objective functions
        """
        # obtain the aim from openmdao storage
        tacs_model = self.options["tacs_model"]
        track_history = self.options["track_history"]
        write_f5 = self.options["write_f5"]

        # update the design
        new_design = tacs_model.update_design(inputs)

        if new_design:
            if track_history:
                self._print_design(inputs)

            # run a forward + adjoint analysis and apply any shape changes if necessary
            tacs_model.pre_analysis()
            tacs_model.run_analysis(write_f5=write_f5, iteration=self._iteration)
            tacs_model.post_analysis()

            self._iteration += 1

            # update func history and report to design file
            if track_history:
                self._update_history()
                self._function_report()

        # Grab the function values and attach as openmdao outputs
        for func in tacs_model.analysis_functions:
            outputs[func.name] = func.value

        return

    def compute_partials(self, inputs, partials):
        """
        the actual value of partial derivatives assigned here
        """
        # obtain the aim from openmdao storage
        tacs_model = self.options["tacs_model"]
        track_history = self.options["track_history"]
        write_f5 = self.options["write_f5"]

        # update the design
        new_design = tacs_model.update_design(inputs)

        if new_design:
            if track_history:
                self._print_design(inputs)

            # run a forward + adjoint analysis and apply any shape changes if necessary
            tacs_model.pre_analysis()
            tacs_model.run_analysis(write_f5=write_f5, iteration=self._iteration)
            tacs_model.post_analysis()

            self._iteration += 1

            # update func history and report to design file
            if track_history:
                self._update_history()
                self._function_report()

        # Grab the function values and attach as openmdao outputs
        for func in tacs_model.analysis_functions:
            for var in tacs_model.variables:
                partials[func.name, var.name] = func.get_derivative(var)

        return

    # helper methods for writing history, plotting history, etc.
    def _update_history(self):
        tacs_model = self.options["tacs_model"]
        for func in tacs_model.analysis_functions:
            self._func_history[func.name].append(func.value)

        if tacs_model.root_proc:
            self._plot_history(
                directory=tacs_model.analysis_dir(tacs_model.root_proc_ind),
                filename="opt_history.png",
            )

    def _function_report(self):
        tacs_model = self.options["tacs_model"]

        if tacs_model.root_proc:
            self._design_hdl.write("Analysis result:\n")
            for func_name in tacs_model.function_names:
                self._design_hdl.write(
                    f"\tfunc {func_name} = {self._func_history[func_name][-1]}\n"
                )
            self._design_hdl.write("\n")
            self._design_hdl.flush()

    def _plot_history(self, directory, filename):
        tacs_model = self.options["tacs_model"]

        if tacs_model.root_proc:
            num_iterations = len(self._func_history[tacs_model.function_names[0]])
            iterations = [_ for _ in range(num_iterations)]
            plt.figure()
            for func_name in tacs_model.function_names:
                yvec = self._func_history[func_name]
                yvec /= max(np.array(yvec))
                plt.plot(iterations, yvec, linewidth=2, label=func_name)
            plt.legend()
            plt.xlabel("iterations")
            plt.ylabel("func values")
            plt.yscale("log")
            plot_filepath = os.path.join(directory, filename)
            plt.savefig(plot_filepath)
            plt.close("all")

    def _print_design(self, inputs):
        tacs_model = self.options["tacs_model"]

        if tacs_model.root_proc:
            self._design_hdl.write("New Design...\n")
            self._design_hdl.write(
                f"\tthick dvs = {[_.name for _ in tacs_model.variables]}\n"
            )
            real_xarray = [float(inputs[key]) for key in inputs]
            self._design_hdl.write(f"\tvalues = {real_xarray}\n")
            self._design_hdl.flush()
