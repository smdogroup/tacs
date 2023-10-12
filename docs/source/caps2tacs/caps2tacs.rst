caps2tacs
*********
The caps2tacs module is a python interface for running TACS analysis on geometries built in ESP/CAPS (Engineering Sketch Pad).
Caps2tacs is built on top of the pytacs interface and tacsAIM (Analysis Interface Module) from ESP/CAPS. 
Caps2tacs is used in FUNtoFEM for sizing and shape optimization of an aerodynmamic structure under oneway-coupled and fully-coupled
aerostructural optimizations. Caps2tacs also provides thermoelastic properties through ESP/CAPS and TACS. For each structural
design, caps2tacs + the tacsAIM build nastran files for the mesh \*.bdf, \*.dat which are located in the ESP/CAPS work directory
usually under ``capsStruct/Scratch/tacs``. Output files can also be stored there as well such as \*.f5 files.

Developer:  
* Sean Engelstad

For more examples using caps2tacs for thermoelastic analysis and involving CFD, please see the 
`FUNtoFEM github <https://github.com/smdogroup/funtofem/>`_.

.. toctree::
  :maxdepth: 1

  main