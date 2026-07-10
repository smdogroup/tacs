Extracting and Restoring Component Design Variables
***************************************************
.. note:: The script for this example can be found under the ``examples/component_dv_extraction/component_dv_extraction.py`` file.

This example demonstrates how to extract a complete, named description of a model's design variables with :meth:`StaticProblem.getComponentDesignVars <tacs.problems.StaticProblem.getComponentDesignVars>` and use it to recreate the same sizing state in a separate TACS execution.

The model combines four shell components with two tube-beam components that cross at the plate center.

The returned dictionary is keyed by component description, and each entry maps design variable group names to their current values.
Two properties make it suitable as a save/restore format:

#. Group names and value types match the keyword arguments of the corresponding constitutive class constructor, so values can be passed straight back to the constructor inside an ``elemCallBack`` function.
#. Every design variable group is always included, whether or not its entries are active design variables, so the receiving model is free to make a different choice of active design variables.

The dictionary is identical on every processor, and active entries always reflect the calling problem's or constraint's current design variable values at the time :meth:`StaticProblem.getComponentDesignVars <tacs.problems.StaticProblem.getComponentDesignVars>` is called.
A companion method, :meth:`StaticProblem.getComponentDesignVarNums <tacs.problems.StaticProblem.getComponentDesignVarNums>`, returns the global design variable numbers in the same structure (with -1 marking inactive entries).

First, import the required libraries:

.. literalinclude:: ../../../examples/component_dv_extraction/component_dv_extraction.py
   :language: python
   :start-after: # [docs:imports-start]
   :end-before: # [docs:imports-end]

Define the component thicknesses and choose which components will have active design variables:

.. literalinclude:: ../../../examples/component_dv_extraction/component_dv_extraction.py
   :language: python
   :start-after: # [docs:parameters-start]
   :end-before: # [docs:parameters-end]

The element callback assigns each plate component its own thickness and maps the beam components to tube-beam constitutive objects:

.. literalinclude:: ../../../examples/component_dv_extraction/component_dv_extraction.py
   :language: python
   :start-after: # [docs:element-callback-start]
   :end-before: # [docs:element-callback-end]

Set up the assembler and a static problem:

.. literalinclude:: ../../../examples/component_dv_extraction/component_dv_extraction.py
   :language: python
   :start-after: # [docs:setup-start]
   :end-before: # [docs:setup-end]

Perturb the active design variables, and solve the problem to produce a new set of function values:

.. literalinclude:: ../../../examples/component_dv_extraction/component_dv_extraction.py
   :language: python
   :start-after: # [docs:perturb-start]
   :end-before: # [docs:perturb-end]

Extract the component design variable dictionary and save it to disk.
The dictionary contains plain Python scalars and numpy arrays, so any serialization format that handles those will work; here we use pickle:

.. literalinclude:: ../../../examples/component_dv_extraction/component_dv_extraction.py
   :language: python
   :start-after: # [docs:extract-start]
   :end-before: # [docs:extract-end]

.. code-block:: text

   {'BEAM_X': {'d': 0.015, 't': 0.0015},
   'BEAM_Y': {'d': 0.015, 't': 0.0015},
   'PLATE.00': {'t': 0.015},
   'PLATE.01': {'t': 0.012},
   'PLATE.02': {'t': 0.021},
   'PLATE.03': {'t': 0.016}}


Finally, in a separate TACS execution, load the file and use it inside the element callback to rebuild the model with the saved sizing.
The restored model produces the same function values as the original:

.. literalinclude:: ../../../examples/component_dv_extraction/component_dv_extraction.py
   :language: python
   :start-after: # [docs:restore-start]
   :end-before: # [docs:restore-end]

.. code-block:: text

   Original model functions:
   {'gravity_ks_vmfailure': np.float64(0.01738452431014043),
   'gravity_mass': np.float64(11.336157282530241)}

   Restored model functions:
   {'gravity_ks_vmfailure': np.float64(0.01738452431014043),
   'gravity_mass': np.float64(11.336157282530241)}
