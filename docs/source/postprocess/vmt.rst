Shear, Moment and Torque Diagrams
=================================

.. automodule:: tacs.postprocess.vmt

Command-line usage
------------------

.. program-output:: python -m tacs.postprocess.vmt --help

API Reference
-------------

.. autofunction:: tacs.postprocess.vmt.computeVMTFromF5

.. autofunction:: tacs.postprocess.vmt.computeVMT

.. autoclass:: tacs.postprocess.vmt.VMTResult
  :members:

.. autofunction:: tacs.postprocess.vmt.plotVMT

.. autofunction:: tacs.postprocess.vmt.scaleVMTResult

.. autofunction:: tacs.postprocess.vmt.writeVMTCsv

Several load cases and envelopes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: tacs.postprocess.vmt.computeVMTCases

.. autofunction:: tacs.postprocess.vmt.computeEnvelope

.. autoclass:: tacs.postprocess.vmt.VMTEnvelope
  :members:

.. autofunction:: tacs.postprocess.vmt.plotVMTEnvelope

.. autofunction:: tacs.postprocess.vmt.plotVMTReport

.. autofunction:: tacs.postprocess.vmt.writeVMTEnvelopeCsv

Reading f5 files
^^^^^^^^^^^^^^^^

.. autofunction:: tacs.postprocess.vmt.loadF5

.. autoclass:: tacs.postprocess.vmt.F5Data
  :members:

.. autofunction:: tacs.postprocess.vmt.nodalLoads

.. autofunction:: tacs.postprocess.vmt.selectComponents

.. autofunction:: tacs.postprocess.vmt.listComponents

.. autofunction:: tacs.postprocess.vmt.buildSilhouette

.. autoclass:: tacs.postprocess.vmt.Silhouette
  :members:
