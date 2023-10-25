# Nonlinear cantilever beam example

This example runs an analysis of a cantilever beam modeled with shell elements subject to a vertical tip force.
The problem is taken from section 3.1 of ["Popular benchmark problems for geometric nonlinear analysis of shells" by Sze et al](https://doi.org/10.1016/j.finel.2003.11.001).

The code in this directory can be used to run analyses using the different linear and nonlinear shell formulations available in TACS and then to generate plots comparing TACS with both Abaqus and analytic results.

## How to run

### Run analyses

```bash
python analysis.py --strainType nonlinear --rotationType linear
python analysis.py --strainType nonlinear --rotationType quadratic
```

### Generate validation plot

```bash
python ValidationPlots.py
```
