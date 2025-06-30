"""
==============================================================================
Plotting stiffened plate optimization results
==============================================================================
@File    :   CompareCrossSections.py
@Date    :   2025/04/27
@Author  :   Alasdair Christison Gray
@Description :
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import os

# ==============================================================================
# External Python modules
# ==============================================================================
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import niceplots
from pyoptsparse import History
import tabulate
import numpy as np

# ==============================================================================
# Extension modules
# ==============================================================================


def plotDesign(
    ax,
    stiffPitch,
    skinThickness,
    stiffenerHeight,
    stiffThickness,
    flangeFraction,
    colors=None,
    label=None,
    **kwargs,
):
    """Plot a stiffened panel cross section

    Parameters
    ----------
    ax : matplotlib axis
        Axis to plot on
    stiffPitch : float
        Stiffener pitch
    skinThickness : float
        Panel skin thickness
    stiffenerHeight : float
        Stiffener height
    stiffThickness : float
        Stiffener thickness
    """

    if colors is None:
        colors = ["blue", "orange"]

    totalWidth = stiffPitch + stiffenerHeight

    # Plot skin
    skin = mpatches.Rectangle(
        (-stiffenerHeight / 2, 0),
        totalWidth,
        skinThickness,
        facecolor=colors[0],
        **kwargs,
    )
    ax.add_artist(skin)

    # Stiffeners
    for xCentre in [0, stiffPitch]:
        flange = mpatches.Rectangle(
            (xCentre - flangeFraction * stiffenerHeight / 2, skinThickness),
            flangeFraction * stiffenerHeight,
            stiffThickness,
            facecolor=colors[1],
            **kwargs,
        )
        web = mpatches.Rectangle(
            (xCentre - stiffThickness / 2, skinThickness + stiffThickness),
            stiffThickness,
            stiffenerHeight,
            facecolor=colors[1],
            **kwargs,
        )
        ax.add_artist(flange)
        ax.add_artist(web)

    # Place a label to the right of the top of the stiffener
    if label is not None:
        ax.annotate(
            label,
            (2 * stiffThickness / 2, stiffenerHeight + skinThickness),
            verticalalignment="center",
            horizontalalignment="left",
            color=colors[1],
        )

    xMargin = 0.05 * totalWidth
    yMargin = 0.05 * stiffenerHeight
    # Compute the axis limits
    xMin = -stiffenerHeight / 2 - xMargin
    xMax = -stiffenerHeight / 2 + totalWidth + xMargin
    yMin = -yMargin
    yMax = skinThickness + stiffThickness + stiffenerHeight + yMargin
    return xMin, xMax, yMin, yMax


def getValuesFromHistory(values, iteration):
    # Struct DVs
    structDVs = values["dv_struct"][iteration]
    currentDVNum = 0
    optStiffPitch = structDVs[currentDVNum] * 1e3 / 10
    currentDVNum += 1
    optSkinThickness = structDVs[currentDVNum] * 1e3 / 100
    currentDVNum += 1
    optSkinPlyFrac = structDVs[currentDVNum : currentDVNum + 4]
    currentDVNum += 4
    optStiffHeight = structDVs[currentDVNum] * 1e3 / 10
    currentDVNum += 1
    optStiffThickness = structDVs[currentDVNum] * 1e3 / 100

    # Strength ratios
    inPlaneFail = 1.5 / values["CompAndShear.KSFailure"][iteration][0]
    pressureFail = 1.5 / values["Pressure.KSFailure"][iteration][0]

    # Max displacement
    maxDisp = values["Pressure.MaxDispZ"][iteration][0]

    # Mass
    mass = values["CompAndShear.Mass"][iteration][0]

    return {
        "Skin thickness": optSkinThickness,
        "Skin ply fractions": list(optSkinPlyFrac * 100),
        "Stiffener thickness": optStiffThickness,
        "Stiffener height": optStiffHeight,
        "Stiffener pitch": optStiffPitch,
        "Stiffener web aspect ratio": optStiffHeight / optStiffThickness,
        "In-plane FOS": inPlaneFail,
        "Pressure FOS": pressureFail,
        "Max displacement": maxDisp,
        "Mass": mass,
    }


# ==============================================================================
# Load the history files
# ==============================================================================

caseDirs = [
    "WithoutStiffenerBuckling",
    "WithStiffenerBuckling",
]
caseNames = [
    "Without stiffener buckling",
    "With stiffener buckling",
]
results = {}
for ii in range(len(caseDirs)):
    histFile = os.path.join(caseDirs[ii], "structOpt.hst")
    hist = History(histFile)
    values = hist.getValues()
    # Get the optimal design
    results[caseNames[ii]] = getValuesFromHistory(values, -1)

# Get the baseline design
results["Baseline"] = getValuesFromHistory(values, 0)
# ==============================================================================
# Convert the results to a table that tabulate can use, one column per case
# ==============================================================================
tableData = []
tableData.append(["Quantity"] + list(results.keys()))
for quantity in results["Baseline"].keys():
    row = [quantity]
    for case in results.keys():
        if quantity == "Skin ply fractions":
            fractions = "[ "
            for frac in results[case][quantity]:
                fractions += f"{frac:.2g}% "
            fractions += "]"
            row.append(fractions)
        else:
            row.append(f"{results[case][quantity]:.4g}")
    tableData.append(row)
print(tabulate.tabulate(tableData, headers="firstrow"))

plt.style.use(niceplots.get_style())
niceColours = niceplots.get_colors()
niceColoursList = niceplots.get_colors_list()


# Create two figures, one where each design is in it's own subplot, and another where they're plotted on top of one another
axWidth = 10
axHeight = 4
combinedFig, combinedAx = plt.subplots(figsize=(axWidth, axHeight))
sepFig, sepAxes = plt.subplots(
    3, 1, sharex=True, sharey=True, figsize=(axWidth, len(results) * axHeight)
)

xLim = [np.inf, -np.inf]
yLim = [np.inf, -np.inf]
for ii, (case, values) in enumerate(results.items()):
    colour = niceColoursList[ii]

    xMin, xMax, yMin, yMax = plotDesign(
        combinedAx,
        values["Stiffener pitch"],
        values["Skin thickness"],
        values["Stiffener height"],
        values["Stiffener thickness"],
        1.0,
        colors=[colour] * 2,
        label=case,
        # alpha=0.5,
        edgecolor=niceColours["Axis"],
        linewidth=0.75,
    )
    xLim[0] = min(xLim[0], xMin)
    xLim[1] = max(xLim[1], xMax)
    yLim[0] = min(yLim[0], yMin)
    yLim[1] = max(yLim[1], yMax)

    plotDesign(
        sepAxes[ii],
        values["Stiffener pitch"],
        values["Skin thickness"],
        values["Stiffener height"],
        values["Stiffener thickness"],
        1.0,
        colors=[colour] * 2,
        # alpha=0.5,
        edgecolor=niceColours["Axis"],
        linewidth=0.75,
    )
    sepAxes[ii].set_aspect("equal")
    niceplots.adjust_spines(sepAxes[ii])


combinedAx.set_aspect("equal")
niceplots.adjust_spines(combinedAx)
combinedAx.set_xlim(xLim[0], xLim[1])
combinedAx.set_ylim(yLim[0], yLim[1])
combinedAx.set_xlabel("$x_2$ (mm)")
combinedAx.set_ylabel("$-x_3$ (mm)")
niceplots.save_figs(combinedFig, "CrossSectionCombined", formats=["pdf", "png", "svg"])

sepAxes[-1].set_xlim(xLim[0], xLim[1])
sepAxes[-1].set_ylim(yLim[0], yLim[1])
sepAxes[-1].set_xlabel("$x_2$ (mm)")
sepAxes[-1].set_ylabel("$-x_3$ (mm)")
niceplots.save_figs(sepFig, "CrossSectionSeparate", formats=["pdf", "png", "svg"])

# plt.show()
