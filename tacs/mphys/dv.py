import numpy as np

import openmdao.api as om
from openmdao.utils.mpi import MPI


class TacsDVComp(om.ExplicitComponent):
    """
    Component for splitting serial tacs design variable from top level
    into distributed vector used by tacs.
    """

    def initialize(self):
        self.options.declare(
            "fea_assembler",
            default=None,
            desc="the pytacs object itself",
            recordable=False,
        )
        self.options.declare(
            "initial_dv_vals",
            default=None,
            desc="initial values for global design variable vector",
        )
        self.options.declare(
            "separate_mass_dvs",
            default=False,
            desc="Flag for whether or not to separate out point mass dvs using user-defined names",
        )

    def setup(self):
        self.fea_assembler = self.options["fea_assembler"]
        self.src_indices = self.get_dv_src_indices()
        vals = self.options["initial_dv_vals"]
        ndv = self.fea_assembler.getNumDesignVars()
        # Keep a list of dvs corresponding to regular struct dvs and mass dvs
        self.struct_dvs = list(range(len(vals)))
        self.mass_dvs = {}
        # Check if user wants to separate out point mass dvs by user-defined names
        if self.options["separate_mass_dvs"]:
            g_dvs = self.fea_assembler.getGlobalDVs()
            for dv_name in g_dvs:
                dv_dict = g_dvs[dv_name]
                # For array-valued global DVs, "isMassDV" is a per-entry bool mask;
                # only the entries assigned through pyTACS.assignMassDV are mass dvs
                if np.any(dv_dict["isMassDV"]):
                    mass_mask = np.atleast_1d(dv_dict["isMassDV"])
                    dv_nums = np.atleast_1d(dv_dict["num"])[mass_mask]
                    mass_vals = vals[dv_nums]
                    # Store mass dv nums with user-defined dv name
                    self.mass_dvs[f"dv_mass_{dv_name}"] = dv_nums
                    # Remove mass dvs from struct list
                    for dv_num in dv_nums:
                        self.struct_dvs.remove(dv_num)
                    # Add user-defined dv name as input
                    self.add_input(
                        f"dv_mass_{dv_name}",
                        desc="serial mass design variable holding mass design variable instances for tacs",
                        val=mass_vals,
                        distributed=False,
                        tags=["mphys_input"],
                    )
            # Remove all mass dvs from vals
            vals = vals[self.struct_dvs]

        self.add_input(
            "dv_struct",
            desc="serial vector holding all structural tacs design variable values",
            val=vals,
            distributed=False,
            tags=["mphys_input"],
        )
        self.add_output(
            "tacs_dvs",
            desc="distributed vector holding tacs design variable values\
                        for this proc",
            shape=ndv,
            distributed=True,
            tags=["mphys_coupling"],
        )

    def get_tot_ndv(self):
        """
        Get total number of serial design variables (struct + mass)
        """
        num_mass_dvs = sum(len(dv_nums) for dv_nums in self.mass_dvs.values())
        return len(self.struct_dvs) + num_mass_dvs

    def compute(self, inputs, outputs):
        # Create serial array to holding all dv vals
        tot_ndv = self.get_tot_ndv()
        full_dv_array = np.zeros(tot_ndv, dtype=inputs["dv_struct"].dtype)
        # Place struct dvs in full array
        full_dv_array[self.struct_dvs] = inputs["dv_struct"]
        # Place mass dvs (if they were defined) in full array
        for dv_name in self.mass_dvs:
            full_dv_array[self.mass_dvs[dv_name]] = inputs[dv_name]
        # Slice full array with src_indices to get distributed dv array
        outputs["tacs_dvs"] = full_dv_array[self.src_indices]

    def compute_jacvec_product(self, inputs, d_inputs, d_outputs, mode):
        tot_ndv = self.get_tot_ndv()
        dfull_dv_array = np.zeros(tot_ndv, dtype=inputs["dv_struct"].dtype)
        if mode == "fwd":
            if "tacs_dvs" in d_outputs:
                if "dv_struct" in d_inputs:
                    dfull_dv_array[self.struct_dvs] += d_inputs["dv_struct"]
                for dv_name in self.mass_dvs:
                    if dv_name in d_inputs:
                        dfull_dv_array[self.mass_dvs[dv_name]] += d_inputs[dv_name]
                d_outputs["tacs_dvs"] += dfull_dv_array[self.src_indices]
        else:  # mode == 'rev'
            if "tacs_dvs" in d_outputs:
                dfull_dv_array[self.src_indices] += d_outputs["tacs_dvs"]
                if "dv_struct" in d_inputs:
                    d_inputs["dv_struct"] += self.comm.allreduce(
                        dfull_dv_array[self.struct_dvs]
                    )
                for dv_name in self.mass_dvs:
                    if dv_name in d_inputs:
                        d_inputs[dv_name] += self.comm.allreduce(
                            dfull_dv_array[self.mass_dvs[dv_name]]
                        )

    def get_dv_src_indices(self):
        """
        Method to get src_indices on each processor
        for tacs distributed design variable vec
        """
        if MPI is not None and self.comm.size > 1:
            local_ndvs = self.fea_assembler.getNumDesignVars()
            all_proc_ndvs = self.comm.gather(local_ndvs, root=0)
            all_proc_indices = []
            if self.comm.rank == 0:
                tot_ndvs = 0
                for proc_i in range(self.comm.size):
                    local_ndvs = all_proc_ndvs[proc_i]
                    proc_indices = np.arange(tot_ndvs, tot_ndvs + local_ndvs)
                    all_proc_indices.append(proc_indices)
                    tot_ndvs += local_ndvs
            local_dv_indices = self.comm.scatter(all_proc_indices, root=0)
            return local_dv_indices
        else:
            ndvs = len(self.options["initial_dv_vals"])
            all_dv_indices = np.arange(ndvs)
            return all_dv_indices
