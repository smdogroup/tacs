"""
==============================================================================

==============================================================================
@File    :   GenerateBDF.py
@Date    :   2025/09/29
@Author  :   Alasdair Christison Gray
@Description : This script generates a NASTRAN input file (.bdf) for a tapered
cantilever beam analysis using pyNastran.

The file models a 1m long beam made of Aluminium 7075. The
rectangular cross-section tapers linearly from the fixed end to the
free end.

Then we define a modal analysis (SOL 103) to find the first 10 natural frequencies.

Command-Line Usage:
    # Generate a model with 30 elements using PBEAM cards
    python generate_tapered_beam.py --num-elements 10 --prop-type PBEAM --output beam_pbeam.bdf

    # Generate a model with 30 elements using PBEAML cards
    python generate_tapered_beam.py --num-elements 20 --prop-type PBEAML --output beam_pbeaml.bdf
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
from typing import List, Tuple

# ==============================================================================
# External Python modules
# ==============================================================================
from pyNastran.bdf.bdf import BDF, CaseControlDeck
import numpy as np

# ==============================================================================
# Extension modules
# ==============================================================================


def getBeamSectionProperties(
    width: float, depth: float
) -> Tuple[float, float, float, float]:
    """
    Calculates the cross-sectional properties for a rectangular section.

    Args:
        width (float): The width of the rectangle.
        depth (float): The depth (height) of the rectangle.

    Returns:
        tuple: A tuple containing:
            - area (float): Cross-sectional area.
            - i1 (float): Moment of inertia about the local y-axis.
            - i2 (float): Moment of inertia about the local z-axis.
            - j (float): Torsional constant.
    """
    area = width * depth
    i1 = width * depth**3 / 12.0
    i2 = depth * width**3 / 12.0

    # Torsional constant J using Roark's formula for a solid rectangular bar
    a = 0.5 * max(width, depth)
    b = 0.5 * min(width, depth)
    j = (a * b**3) * (16.0 / 3.0 - 3.36 * (b / a) * (1.0 - (b**4) / (12.0 * a**4)))
    return area, i1, i2, j


def addBeamPropertyCard(
    model: BDF,
    propType: str,
    propertyId: int,
    materialId: int,
    width1: float,
    depth1: float,
    width2: float,
    depth2: float,
) -> None:
    """Create the property card for a tapered beam element

    Parameters
    ----------
    model : BDF
        pyNastran BDF model object
    propType : str
        Property card type ('PBEAM' or 'PBEAML')
    propertyId : int
        ID for the property card
    materialId : int
        Material ID to link to the property
    width1 : float
        Rectangular section width at end 1
    depth1 : float
        Rectangular section depth at end 1
    width2 : float
        Rectangular section width at end 2
    depth2 : float
        Rectangular section depth at end 2
    """
    stations = [0.0, 1.0]  # Stations along the beam element length
    if propType.upper() == "PBEAML":
        # PBEAML is ideal for linearly tapered beams.
        # We define the dimensions at the start (A) and end (B) of the beam.
        beamType = "BAR"
        dimsA = [width1, depth1]  # Dimensions at End A
        dimsB = [width2, depth2]  # Dimensions at End B
        model.add_pbeaml(
            pid=propertyId,
            mid=materialId,
            beam_type=beamType,
            xxb=stations,
            dims=[dimsA, dimsB],
        )

    elif propType.upper() == "PBEAM":
        # PBEAM requires explicit calculation of section properties at stations.
        # We will define properties at each end of the single property region.
        areas, i1s, i2s, js = [], [], [], []

        # Calculate properties at the two ends of the full beam length
        for width, depth in [(width1, depth1), (width2, depth2)]:
            area, i1, i2, j = getBeamSectionProperties(width, depth)
            areas.append(area)
            i1s.append(i1)
            i2s.append(i2)
            js.append(j)

        model.add_pbeam(
            propertyId,
            materialId,
            so="NO",
            area=areas,
            i1=i1s,
            i2=i2s,
            i12=[0.0] * len(areas),
            j=js,
            xxb=stations,
        )


def generateTaperedBeamBdf(numElements: int, propType: str, output: str) -> None:
    """
    Generates the NASTRAN BDF file for the tapered beam analysis.

    Args:
        numElements (int): The number of CBEAM elements to create.
        propType (str): The property card type to use ('PBEAM' or 'PBEAML').
        output (str): The path to write the output BDF file.
    """
    if numElements <= 0:
        raise ValueError("Number of elements must be a positive integer.")
    if propType.upper() not in ["PBEAM", "PBEAML"]:
        raise ValueError("Property type must be either 'PBEAM' or 'PBEAML'.")

    # --- Model Definition ---
    model = BDF(debug=False)

    # --- Define offsets for ID numbers ---
    gridOffset = 1000
    materialOffset = 2000
    propertyOffset = 3000
    elementOffset = 4000
    SPCOffset = 5000
    loadOffset = 6000

    # --- Initial IDs ---
    nodeId = 1 + gridOffset
    elementId = 1 + elementOffset
    materialId = 1 + materialOffset
    propertyId = 1 + propertyOffset
    spcId = 1 + SPCOffset
    gravLoadId = 1 + loadOffset
    eigrlId = 1
    coordId = 0  # Basic coordinate system

    # --- Beam Geometry ---
    beamLength = 1.0  # meters
    # Dimensions at the fixed end (x=0)
    width1, depth1 = 0.20, 0.05
    # Dimensions at the free end (x=L)
    width2, depth2 = 0.10, 0.02

    widths = np.linspace(width1, width2, numElements + 1)
    depths = np.linspace(depth1, depth2, numElements + 1)

    # --- Material Properties (Aluminium 7075) ---
    E = 71.7e9  # Young's Modulus in N/m^2
    nu = 0.33
    G = E / (2 * (1 + nu))  # Shear Modulus in N/m^2
    rho = 2810.0  # Density in kg/m^3
    model.add_mat1(mid=materialId, E=E, G=G, nu=nu, rho=rho, comment="Aluminium 7075")

    # --- Nodes and Elements ---
    xNodes = np.linspace(0.0, beamLength, numElements + 1)
    nodeIds = []
    # Create GRID cards along the beam's length (x-axis)
    for xNode in xNodes:
        model.add_grid(nodeId, [xNode, 0.0, 0.0])
        nodeIds.append(nodeId)
        nodeId += 1

    # Create CBEAM elements connecting the nodes and corresponding properties
    # The z axis is vertical (depth direction), y axis is width direction
    orientationVector = [0.0, 0.0, 1.0]
    componentPrefix = "Elements and Element Properties for region :"
    for ii in range(numElements):
        elementId
    for i in range(numElements):
        # Create property card
        addBeamPropertyCard(
            model=model,
            propType=propType,
            propertyId=propertyId,
            materialId=materialId,
            width1=widths[i],
            depth1=depths[i],
            width2=widths[i + 1],
            depth2=depths[i + 1],
        )

        # Create element card
        nodeA = nodeIds[i]
        nodeB = nodeIds[i + 1]
        model.add_cbeam(
            elementId,
            propertyId,
            [nodeA, nodeB],
            orientationVector,
            comment=f"{componentPrefix} Element {elementId-elementOffset}",
        )
        elementId += 1
        propertyId += 1

    # --- Boundary Conditions ---
    # Fix the root of the cantilever beam (node 1) in all 6 DOFs
    fixedNodeId = nodeIds[0]
    constrainedDofs = "123456"
    model.add_spc1(spcId, constrainedDofs, [fixedNodeId])

    # --- Loads ---
    # Apply standard gravity (9.81 m/s^2) in the negative Z direction
    g = 9.81
    gravityVector = [0.0, 0.0, -1.0]
    model.add_grav(gravLoadId, g, gravityVector)

    # --- Analysis Case Control ---
    # Just do the modal analysis for now
    eigrl = model.add_eigrl(42, nd=10)

    model.sol = 103  # start=103
    cc = CaseControlDeck(
        [
            "SUBCASE 1",
            "  LABEL = Modal Analysis",
            "  METHOD = 42",  # TODO: remove
            "  SPC = %s" % spcId,
            "  DISP(PLOT) = ALL",
        ]
    )
    # print(cc)
    model.case_control_deck = cc
    model.validate()

    # --- Write BDF File ---
    # Set PARAM POST 1 to output an op2 file
    model.add_param("POST", 1)
    model.write_bdf(output, enddata=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate a NASTRAN BDF file for a tapered cantilever beam analysis.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-n",
        "--num-elements",
        type=int,
        default=10,
        help="Number of CBEAM elements to use along the beam length.",
    )
    parser.add_argument(
        "-p",
        "--prop-type",
        type=str,
        default="PBEAM",
        choices=["PBEAM", "PBEAML"],
        help="The type of beam property card to use.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="tapered_beam.bdf",
        help="Path for the output BDF file.",
    )
    args = parser.parse_args()
    generateTaperedBeamBdf(
        numElements=args.num_elements, propType=args.prop_type, output=args.output
    )
