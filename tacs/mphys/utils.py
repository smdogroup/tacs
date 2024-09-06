import warnings
import tacs


def add_tacs_constraints(scenario):
    """Call this in the configure method to add all TACS constraints functions
    in a TacsPostcouplingGroup as constraints in an OpenMDAO model. This saves
    you having to call OpenMDAO's `add_constraint` method for each constraint
    function yourself.

    This function relies on you having defined the lower and/or upper bounds
    of the constraints when you created them in your `constraint_setup`
    function.

    Parameters
    ----------
    scenario : MPhys scenario
        Scenario containing a TACS postcoupling group with constraints
    """
    # First, verify that the model has a struct_post component that is a
    # TacsPostcouplingGroup
    if not hasattr(scenario, "struct_post"):
        raise ValueError(
            f"Scenario {scenario.name} does not have a struct_post component"
        )
    if not isinstance(
        scenario.struct_post, tacs.mphys.postcoupling.TacsPostcouplingGroup
    ):
        raise ValueError(
            f"{scenario.name}.struct_post is not a TacsPostcouplingGroup, cannot add TACS constraints"
        )
    # If the provided scenario has no constraints component, warn the user and return
    if not hasattr(scenario.struct_post, "constraints"):
        warnings.warn(
            f"Scenario {scenario.name} contains no TACS constraints to add",
            stacklevel=2,
        )
        return
    else:
        for system in scenario.struct_post.constraints.system_iter():
            constraint = system.constr
            constraintFuncNames = constraint.getConstraintKeys()
            bounds = {}
            constraint.getConstraintBounds(bounds)
            for conName in constraintFuncNames:
                if scenario.comm.rank == 0:
                    print("Adding TACS constraint: ", conName)
                name = f"{system.name}.{conName}"
                # Panel length constraints are nonlinear, all other constrain
                # types (DV and adjacency) are linear
                if isinstance(
                    constraint,
                    tacs.constraints.panel_length.PanelLengthConstraint,
                ):
                    scenario.add_constraint(name, equals=0.0, scaler=1.0)
                else:
                    lb = bounds[f"{system.name}_{conName}"][0]
                    ub = bounds[f"{system.name}_{conName}"][1]
                    if all(lb == ub):
                        scenario.add_constraint(name, equals=lb, linear=True)
                    else:
                        scenario.add_constraint(name, lower=lb, upper=ub, linear=True)
