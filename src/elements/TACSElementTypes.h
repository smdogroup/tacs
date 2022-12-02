/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_ELEMENT_TYPES_H
#define TACS_ELEMENT_TYPES_H

/*
  Quantity output types defined globally
*/
const int TACS_ELEMENT_DENSITY = 1;
const int TACS_STRAIN_ENERGY_DENSITY = 2;
const int TACS_FAILURE_INDEX = 3;
const int TACS_HEAT_FLUX = 4;
const int TACS_TEMPERATURE = 5;
const int TACS_TOTAL_STRAIN_ENERGY_DENSITY = 6;
const int TACS_ELEMENT_DISPLACEMENT = 7;
const int TACS_ELEMENT_STRAIN = 8;
const int TACS_ELEMENT_STRESS = 9;
const int TACS_ELEMENT_DENSITY_MOMENT = 10;
const int TACS_ELEMENT_MOMENT_OF_INERTIA = 11;

/*
  The output flags from each element type
*/
const int TACS_OUTPUT_CONNECTIVITY = 1;
const int TACS_OUTPUT_NODES = 2;
const int TACS_OUTPUT_DISPLACEMENTS = 4;
const int TACS_OUTPUT_STRAINS = 8;
const int TACS_OUTPUT_STRESSES = 16;
const int TACS_OUTPUT_EXTRAS = 32;

/**
  The Element type defines how many displacement, stress, and strain
  components are written to an output file.
*/
enum ElementType {
  TACS_ELEMENT_NONE = 0,
  TACS_SCALAR_2D_ELEMENT = 1,
  TACS_SCALAR_3D_ELEMENT = 2,
  TACS_BEAM_OR_SHELL_ELEMENT = 4,
  TACS_PLANE_STRESS_ELEMENT = 5,
  TACS_SOLID_ELEMENT = 6,
  TACS_RIGID_ELEMENT = 7,
  TACS_MASS_ELEMENT = 8,
  TACS_SPRING_ELEMENT = 9,
  TACS_PCM_ELEMENT = 10
};

/**
  Defines the element layout used for writing to an output file
*/
enum ElementLayout {
  TACS_LAYOUT_NONE = 0,
  TACS_POINT_ELEMENT = 1,

  TACS_LINE_ELEMENT = 2,
  TACS_LINE_QUADRATIC_ELEMENT = 3,
  TACS_LINE_CUBIC_ELEMENT = 4,

  TACS_TRI_ELEMENT = 5,
  TACS_TRI_QUADRATIC_ELEMENT = 6,
  TACS_TRI_CUBIC_ELEMENT = 7,

  TACS_QUAD_ELEMENT = 8,
  TACS_QUAD_QUADRATIC_ELEMENT = 9,
  TACS_QUAD_CUBIC_ELEMENT = 10,
  TACS_QUAD_QUARTIC_ELEMENT = 11,
  TACS_QUAD_QUINTIC_ELEMENT = 12,

  TACS_TETRA_ELEMENT = 13,
  TACS_TETRA_QUADRATIC_ELEMENT = 14,
  TACS_TETRA_CUBIC_ELEMENT = 15,

  TACS_HEXA_ELEMENT = 16,
  TACS_HEXA_QUADRATIC_ELEMENT = 17,
  TACS_HEXA_CUBIC_ELEMENT = 18,
  TACS_HEXA_QUARTIC_ELEMENT = 19,
  TACS_HEXA_QUINTIC_ELEMENT = 20,

  TACS_PENTA_ELEMENT = 21,
  TACS_PENTA_QUADRATIC_ELEMENT = 22,
  TACS_PENTA_CUBIC_ELEMENT = 23
};

/**
   The different element matrix types
*/
enum ElementMatrixType {
  TACS_JACOBIAN_MATRIX,
  TACS_STIFFNESS_MATRIX,
  TACS_MASS_MATRIX,
  TACS_GEOMETRIC_STIFFNESS_MATRIX,
  TACS_STIFFNESS_PRODUCT_DERIVATIVE
};

/**
  Get the total number of components associated with the element type
  and flag

  @param etype Element type
  @param flag A binary flag indicating all the categories of component
  @return The total number of components
*/
int TacsGetTotalOutputCount(ElementType etype, int flag);

/**
  Get the number of components associated with the output

  @param etype Element type
  @param comp The category of component
  @return The number of components for that element type
*/
int TacsGetOutputComponentCount(ElementType etype, int comp);

/**
  Get the number of components associated with the output

  @param etype Element type
  @param comp The category of component
  @param index The index of the component
  @return The name of the component
*/
const char *TacsGetOutputComponentName(ElementType etype, int comp, int index);

/**
   Get the number of visualization nodes for the given element layout

   @param ltype The element layout type for the visualization
   @return The number of nodes in the visual representation of the element
*/
int TacsGetNumVisNodes(ElementLayout ltype);

/**
  Get the element layout count and number of entries in the new
  connectivity.

  Each higher-order element type can be converted to a series of
  basic element types for visualiztaion. This function returns the
  number of element types that are needed and the number of new
  entries required in the connectivity array. The basic element types
  are:

  TACS_POINT_ELEMENT
  TACS_LINE_ELEMENT
  TACS_TRI_ELEMENT
  TACS_QUAD_ELEMENT
  TACS_TETRA_ELEMENT
  TACS_HEXA_ELEMENT

  @param ltype The element layout type
  @param ntypes The number of basic element types.
*/
void TacsConvertVisLayoutToBasicCount(ElementLayout ltype, int *ntypes,
                                      int *nconn);

/**
  Retrieve the new element types and new element connectivity for the
  basic elements used for visualization

  @param ltype The element layout type
  @param conn The original element connectivity
  @param basic_ltypes The basic element types (cast of ElementLayout to int)
  @param basic_conn The basic element connectivity
*/
void TacsConvertVisLayoutToBasic(ElementLayout ltype, const int conn[],
                                 int basic_ltypes[], int basic_conn[]);

#endif  // TACS_ELEMENT_TYPES_H
