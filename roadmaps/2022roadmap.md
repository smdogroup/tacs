TACS 2022 Development Roadmap
=================================

Date: 2022-03-28

## What is this roadmap for?
This document represents the perspective of the TACS development team on what we feel will be the most important areas to focus on for framework development in 2022.
It is intended to communicate this perspective to the TACS user community, and to provide a basis for feedback.
It is not intended to be all-encompassing or a binding  development plan, but should be regarded as an indication of the current direction of the core development team.
Also, it should be clearly understood that there is very probably more listed here than can actually be accomplished in a single year!

# 2022 main focus areas
- Adding Automatic Differentiation (AD) tools to TACS source
- Allow for use of external numerical solver libraries with TACS
- Simplify TACS installation procedure
- Add more documentation

----------------------------------------------
----------------------------------------------
----------------------------------------------

# Adding Automatic Differentiation (AD) tools to TACS source

One of the painful processes in creating new element implementations in TACS has always been providing the element derivatives necessary for computing the system-level adjoint.
Up until now, these derivatives had derived analytically by the user, which often a labor-intensive process.
We feel that the addition of a common AD library for the TACS C++ source code will make future development of new elements easier.

## 1) Develop AD library for common numerical operations in TACS source code

### Goal:
Provide a common C++ library within TACS from which developers include in their element source code that will simplify the procedure for implementing derivative-based methods.
This library should be simple, fast, and memory-efficient.

### Potential Challenges:
Currently the scope of this AD toolbox is only to be used within TACS element classes. Other class types, like constitutive classes,
will still need to implement their derivatives by hand.

## 2) Use AD to implement missing elements sensitivities

### Goal:
There are still a number of elements in the TACS library that do not full implementation of all sensitivity methods required for adjoint computations.
By default, if these are not provided, then they are estimated using a FD or CS approximation.
We want to re-implement these missing methods using the AD tools, leading to more accurate and efficient computations.

# Allow for use of external numerical solver libraries with TACS
There are many options for numerical solver libraries available to users (ex. PETSc, scipy, etc.).
However, currently most typical TACS problems (linear systems, eigenvalue problems, etc.),
can only be solved with native TACS solver implementation (TACS.KSM, TACS.FrequencyAnalysis, etc.).
We hope to relax this restriction and allow for use of these various libraries in conjunction with TACS numerical objects (vectors, matrices, etc.)

# Simplify TACS installation procedure
Currently the procedure for installing TACS through Python doesn't follow a lot of the Python conventions.
This can lead to confusion for new users when setting up TACS for the first time.

## 1) Clean up setup.py

### Goal:
Properly compile and expose C++ source code dependencies through setup.py such that all dependencies are properly
compiled and installed as expected using a `pip install .` command.

## 2) Package TACS into a PyPI and/or conda package

### Goal:
Once a more standard Python installation process is defined, we'd like to package and deploy the TACS library into a PyPI or conda-forge.
Users should then be able to install TACS, and all of it's dependencies, with python using one command:
```
$pip install tacs
```
or
```
$conda install tacs
```

# Add more documentation
Currently the TACS library documentation is still limited. Effort will be put forth to expand this coverage out further.

## 1) Add more examples to docs

### Goal:
Expand out the docs examples so we demonstrate a broader variety of TACS methods and use cases..

## 2) Document MPHYS interface and add examples/tests

### Goal:
The recently implemented MPHYS wrapper for TACS is not documented and no tests or examples include it in the repo.
These three things will need to be added to ensure stable development of the module in the future.
