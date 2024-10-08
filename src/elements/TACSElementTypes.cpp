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

#include "TACSElementTypes.h"

#include "TACSObject.h"

/*
  Count up the total number of values associate with all outputs
*/
int TacsGetTotalOutputCount(ElementType etype, int flag) {
  int nvals = 0;
  if (flag & TACS_OUTPUT_NODES) {
    nvals += TacsGetOutputComponentCount(etype, TACS_OUTPUT_NODES);
  }
  if (flag & TACS_OUTPUT_DISPLACEMENTS) {
    nvals += TacsGetOutputComponentCount(etype, TACS_OUTPUT_DISPLACEMENTS);
  }
  if (flag & TACS_OUTPUT_STRAINS) {
    nvals += TacsGetOutputComponentCount(etype, TACS_OUTPUT_STRESSES);
  }
  if (flag & TACS_OUTPUT_STRESSES) {
    nvals += TacsGetOutputComponentCount(etype, TACS_OUTPUT_STRAINS);
  }
  if (flag & TACS_OUTPUT_EXTRAS) {
    nvals += TacsGetOutputComponentCount(etype, TACS_OUTPUT_EXTRAS);
  }
  if (flag & TACS_OUTPUT_LOADS) {
    nvals += TacsGetOutputComponentCount(etype, TACS_OUTPUT_LOADS);
  }

  return nvals;
}

/*
  Get the number of components associated with the output
*/
int TacsGetOutputComponentCount(ElementType etype, int comp) {
  if (comp == TACS_OUTPUT_NODES) {
    return 3;
  }

  if (etype == TACS_ELEMENT_NONE) {
    return 0;
  } else if (etype == TACS_SCALAR_2D_ELEMENT) {
    if (comp == TACS_OUTPUT_DISPLACEMENTS) {
      return 1;
    } else if (comp == TACS_OUTPUT_STRAINS) {
      return 2;
    } else if (comp == TACS_OUTPUT_STRESSES) {
      return 2;
    } else if (comp == TACS_OUTPUT_EXTRAS) {
      return 4;
    } else if (comp == TACS_OUTPUT_LOADS) {
      return 1;
    }
  } else if (etype == TACS_SCALAR_3D_ELEMENT) {
    if (comp == TACS_OUTPUT_DISPLACEMENTS) {
      return 1;
    } else if (comp == TACS_OUTPUT_STRAINS) {
      return 3;
    } else if (comp == TACS_OUTPUT_STRESSES) {
      return 3;
    } else if (comp == TACS_OUTPUT_EXTRAS) {
      return 4;
    } else if (comp == TACS_OUTPUT_LOADS) {
      return 1;
    }
  } else if (etype == TACS_BEAM_OR_SHELL_ELEMENT) {
    if (comp == TACS_OUTPUT_DISPLACEMENTS) {
      return 6;
    } else if (comp == TACS_OUTPUT_STRAINS) {
      return 9;
    } else if (comp == TACS_OUTPUT_STRESSES) {
      return 9;
    } else if (comp == TACS_OUTPUT_EXTRAS) {
      return 14;
    } else if (comp == TACS_OUTPUT_LOADS) {
      return 6;
    }
  } else if (etype == TACS_PLANE_STRESS_ELEMENT) {
    if (comp == TACS_OUTPUT_DISPLACEMENTS) {
      return 2;
    } else if (comp == TACS_OUTPUT_STRAINS) {
      return 3;
    } else if (comp == TACS_OUTPUT_STRESSES) {
      return 3;
    } else if (comp == TACS_OUTPUT_EXTRAS) {
      return 4;
    } else if (comp == TACS_OUTPUT_LOADS) {
      return 2;
    }
  } else if (etype == TACS_SOLID_ELEMENT) {
    if (comp == TACS_OUTPUT_DISPLACEMENTS) {
      return 3;
    } else if (comp == TACS_OUTPUT_STRAINS) {
      return 6;
    } else if (comp == TACS_OUTPUT_STRESSES) {
      return 6;
    } else if (comp == TACS_OUTPUT_EXTRAS) {
      return 4;
    } else if (comp == TACS_OUTPUT_LOADS) {
      return 3;
    }
  } else if (etype == TACS_PCM_ELEMENT) {
    if (comp == TACS_OUTPUT_DISPLACEMENTS) {
      return 1;
    } else if (comp == TACS_OUTPUT_STRAINS) {
      return 2;
    } else if (comp == TACS_OUTPUT_STRESSES) {
      return 2;
    } else if (comp == TACS_OUTPUT_EXTRAS) {
      return 5;
    } else if (comp == TACS_OUTPUT_LOADS) {
      return 1;
    }
  }

  return 0;
}

/*
  Get the output string associated with the component and index
*/
const char *TacsGetOutputComponentName(ElementType etype, int comp, int index) {
  if (index < 0 || index > TacsGetOutputComponentCount(etype, comp)) {
    return NULL;
  }

  if (comp == TACS_OUTPUT_NODES) {
    switch (index) {
      case 0:
        return "X";
      case 1:
        return "Y";
      case 2:
        return "Z";
      default:
        return NULL;
    }
  }

  if (etype == TACS_ELEMENT_NONE) {
    return NULL;
  } else if (etype == TACS_SCALAR_2D_ELEMENT) {
    if (comp == TACS_OUTPUT_DISPLACEMENTS) {
      return "u";
    } else if (comp == TACS_OUTPUT_STRAINS) {
      switch (index) {
        case 0:
          return "ux";
        case 1:
          return "uy";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_STRESSES) {
      switch (index) {
        case 0:
          return "sx";
        case 1:
          return "sy";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_EXTRAS) {
      switch (index) {
        case 0:
          return "failure";
        case 1:
          return "dv1";
        case 2:
          return "dv2";
        case 3:
          return "dv3";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_LOADS) {
      return "f";
    }
  } else if (etype == TACS_SCALAR_3D_ELEMENT) {
    if (comp == TACS_OUTPUT_DISPLACEMENTS) {
      return "u";
    } else if (comp == TACS_OUTPUT_STRAINS) {
      switch (index) {
        case 0:
          return "ux";
        case 1:
          return "uy";
        case 2:
          return "uz";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_STRESSES) {
      switch (index) {
        case 0:
          return "sx";
        case 1:
          return "sy";
        case 2:
          return "sz";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_EXTRAS) {
      switch (index) {
        case 0:
          return "failure";
        case 1:
          return "dv1";
        case 2:
          return "dv2";
        case 3:
          return "dv3";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_LOADS) {
      return "f";
    }
  } else if (etype == TACS_BEAM_OR_SHELL_ELEMENT) {
    if (comp == TACS_OUTPUT_DISPLACEMENTS) {
      switch (index) {
        case 0:
          return "u";
        case 1:
          return "v";
        case 2:
          return "w";
        case 3:
          return "rotx";
        case 4:
          return "roty";
        case 5:
          return "rotz";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_STRAINS) {
      switch (index) {
        case 0:
          return "ex0";
        case 1:
          return "ey0";
        case 2:
          return "exy0";
        case 3:
          return "ex1";
        case 4:
          return "ey1";
        case 5:
          return "exy1";
        case 6:
          return "eyz0";
        case 7:
          return "exz0";
        case 8:
          return "erot";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_STRESSES) {
      switch (index) {
        case 0:
          return "sx0";
        case 1:
          return "sy0";
        case 2:
          return "sxy0";
        case 3:
          return "sx1";
        case 4:
          return "sy1";
        case 5:
          return "sxy1";
        case 6:
          return "syz0";
        case 7:
          return "sxz0";
        case 8:
          return "srot";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_EXTRAS) {
      switch (index) {
        case 0:
          return "failure0";
        case 1:
          return "failure1";
        case 2:
          return "failure2";
        case 3:
          return "failure3";
        case 4:
          return "failure4";
        case 5:
          return "failure5";
        case 6:
          return "failure6";
        case 7:
          return "dv1";
        case 8:
          return "dv2";
        case 9:
          return "dv3";
        case 10:
          return "dv4";
        case 11:
          return "dv5";
        case 12:
          return "dv6";
        case 13:
          return "dv7";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_LOADS) {
      switch (index) {
        case 0:
          return "fx";
        case 1:
          return "fy";
        case 2:
          return "fz";
        case 3:
          return "mx";
        case 4:
          return "my";
        case 5:
          return "mz";
        default:
          return NULL;
      }
    }
  } else if (etype == TACS_PLANE_STRESS_ELEMENT) {
    if (comp == TACS_OUTPUT_DISPLACEMENTS) {
      switch (index) {
        case 0:
          return "u";
        case 1:
          return "v";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_STRAINS) {
      switch (index) {
        case 0:
          return "exx";
        case 1:
          return "eyy";
        case 2:
          return "gxy";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_STRESSES) {
      switch (index) {
        case 0:
          return "sxx";
        case 1:
          return "syy";
        case 2:
          return "sxy";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_EXTRAS) {
      switch (index) {
        case 0:
          return "failure";
        case 1:
          return "dv1";
        case 2:
          return "dv2";
        case 3:
          return "dv3";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_LOADS) {
      switch (index) {
        case 0:
          return "fx";
        case 1:
          return "fy";
        default:
          return NULL;
      }
    }
  } else if (etype == TACS_SOLID_ELEMENT) {
    if (comp == TACS_OUTPUT_DISPLACEMENTS) {
      switch (index) {
        case 0:
          return "u";
        case 1:
          return "v";
        case 2:
          return "w";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_STRAINS) {
      switch (index) {
        case 0:
          return "exx";
        case 1:
          return "eyy";
        case 2:
          return "ezz";
        case 3:
          return "gyz";
        case 4:
          return "gxz";
        case 5:
          return "gxy";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_STRESSES) {
      switch (index) {
        case 0:
          return "sxx";
        case 1:
          return "syy";
        case 2:
          return "szz";
        case 3:
          return "syz";
        case 4:
          return "sxz";
        case 5:
          return "sxy";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_EXTRAS) {
      switch (index) {
        case 0:
          return "failure";
        case 1:
          return "dv1";
        case 2:
          return "dv2";
        case 3:
          return "dv3";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_LOADS) {
      switch (index) {
        case 0:
          return "fx";
        case 1:
          return "fy";
        case 2:
          return "fz";
        default:
          return NULL;
      }
    }
  } else if (etype == TACS_PCM_ELEMENT) {
    if (comp == TACS_OUTPUT_DISPLACEMENTS) {
      switch (index) {
        case 0:
          return "dT";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_STRAINS) {
      switch (index) {
        case 0:
          return "gradx";
        case 1:
          return "grady";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_STRESSES) {
      switch (index) {
        case 0:
          return "fluxx";
        case 1:
          return "fluxy";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_EXTRAS) {
      switch (index) {
        case 0:
          return "rho";
        case 1:
          return "dv1";
        case 2:
          return "dv2";
        case 3:
          return "dv3";
        case 4:
          return "phase";
        default:
          return NULL;
      }
    } else if (comp == TACS_OUTPUT_LOADS) {
      switch (index) {
        case 0:
          return "Q";
        default:
          return NULL;
      }
    }
  }

  return NULL;
}

/*
  Get the number of visualization nodes for the given element layout
*/
int TacsGetNumVisNodes(ElementLayout ltype) {
  switch (ltype) {
    case TACS_LAYOUT_NONE:
      return 0;
    case TACS_POINT_ELEMENT:
      return 1;

    case TACS_LINE_ELEMENT:
      return 2;
    case TACS_LINE_QUADRATIC_ELEMENT:
      return 3;
    case TACS_LINE_CUBIC_ELEMENT:
      return 4;

    case TACS_TRI_ELEMENT:
      return 3;
    case TACS_TRI_QUADRATIC_ELEMENT:
      return 6;
    case TACS_TRI_CUBIC_ELEMENT:
      return 10;

    case TACS_QUAD_ELEMENT:
      return 4;
    case TACS_QUAD_QUADRATIC_ELEMENT:
      return 9;
    case TACS_QUAD_CUBIC_ELEMENT:
      return 16;
    case TACS_QUAD_QUARTIC_ELEMENT:
      return 25;
    case TACS_QUAD_QUINTIC_ELEMENT:
      return 36;

    case TACS_TETRA_ELEMENT:
      return 4;
    case TACS_TETRA_QUADRATIC_ELEMENT:
      return 10;
    case TACS_TETRA_CUBIC_ELEMENT:
      return 17;

    case TACS_HEXA_ELEMENT:
      return 8;
    case TACS_HEXA_QUADRATIC_ELEMENT:
      return 27;
    case TACS_HEXA_CUBIC_ELEMENT:
      return 64;
    case TACS_HEXA_QUARTIC_ELEMENT:
      return 125;
    case TACS_HEXA_QUINTIC_ELEMENT:
      return 216;

    case TACS_PENTA_ELEMENT:
      return 6;
    case TACS_PENTA_QUADRATIC_ELEMENT:
      return 0;  /// Fix me later
    case TACS_PENTA_CUBIC_ELEMENT:
      return 0;  /// Fix me later
    default:
      break;
  }

  return -1;
}

/*
  Get the element layout count and number of entries in the new
  connectivity.
*/
void TacsConvertVisLayoutToBasicCount(ElementLayout ltype, int *ntypes,
                                      int *nconn) {
  *ntypes = 0;
  *nconn = 0;

  switch (ltype) {
    case TACS_POINT_ELEMENT:
      break;

    case TACS_LINE_ELEMENT:
      *ntypes = 1;
      *nconn = 2;
      break;
    case TACS_LINE_QUADRATIC_ELEMENT:
      *ntypes = 2;
      *nconn = 4;
      break;
    case TACS_LINE_CUBIC_ELEMENT:
      *ntypes = 3;
      *nconn = 6;
      break;

    case TACS_TRI_ELEMENT:
      *ntypes = 1;
      *nconn = 3;
      break;
    case TACS_TRI_QUADRATIC_ELEMENT:
      *ntypes = 4;
      *nconn = 12;
      break;
    case TACS_TRI_CUBIC_ELEMENT:
      *ntypes = 9;
      *nconn = 27;
      break;

    case TACS_QUAD_ELEMENT:
      *ntypes = 1;
      *nconn = 4;
      break;
    case TACS_QUAD_QUADRATIC_ELEMENT:
      *ntypes = 4;
      *nconn = 4 * 4;
      break;
    case TACS_QUAD_CUBIC_ELEMENT:
      *ntypes = 9;
      *nconn = 4 * 9;
      break;
    case TACS_QUAD_QUARTIC_ELEMENT:
      *ntypes = 16;
      *nconn = 4 * 16;
      break;
    case TACS_QUAD_QUINTIC_ELEMENT:
      *ntypes = 25;
      *nconn = 4 * 25;
      break;

    case TACS_TETRA_ELEMENT:
      *ntypes = 1;
      *nconn = 4;
      break;
    case TACS_TETRA_QUADRATIC_ELEMENT:
    case TACS_TETRA_CUBIC_ELEMENT:
      break;

    case TACS_HEXA_ELEMENT:
      *ntypes = 1;
      *nconn = 8;
      break;
    case TACS_HEXA_QUADRATIC_ELEMENT:
      *ntypes = 8;
      *nconn = 8 * 8;
      break;
    case TACS_HEXA_CUBIC_ELEMENT:
      *ntypes = 27;
      *nconn = 8 * 27;
      break;
    case TACS_HEXA_QUARTIC_ELEMENT:
      *ntypes = 64;
      *nconn = 8 * 64;
      break;
    case TACS_HEXA_QUINTIC_ELEMENT:
      *ntypes = 125;
      *nconn = 8 * 125;
      break;

    case TACS_PENTA_ELEMENT:
    case TACS_PENTA_QUADRATIC_ELEMENT:
    case TACS_PENTA_CUBIC_ELEMENT:
      break;

    default:
      break;
  }
}

/*
  Retrieve the new element types and new element connectivity for the
  basic elements used for visualization
*/
void TacsConvertVisLayoutToBasic(ElementLayout ltype, const int conn[],
                                 int basic_ltypes[], int basic_conn[]) {
  if (ltype == TACS_LINE_ELEMENT || ltype == TACS_LINE_QUADRATIC_ELEMENT ||
      ltype == TACS_LINE_CUBIC_ELEMENT) {
    int order = 2;
    if (ltype == TACS_QUAD_QUADRATIC_ELEMENT) {
      order = 3;
    } else if (ltype == TACS_QUAD_CUBIC_ELEMENT) {
      order = 4;
    }

    for (int i = 0; i < order - 1; i++) {
      basic_ltypes[i] = TACS_LINE_ELEMENT;

      basic_conn[2 * i] = conn[i];
      basic_conn[2 * i + 1] = conn[i + 1];
    }
  } else if (ltype == TACS_TRI_ELEMENT) {
    basic_ltypes[0] = TACS_TRI_ELEMENT;
    basic_conn[0] = conn[0];
    basic_conn[1] = conn[1];
    basic_conn[2] = conn[2];
  } else if (ltype == TACS_TRI_QUADRATIC_ELEMENT) {
    for (int i = 0; i < 4; i++) {
      basic_ltypes[i] = TACS_TRI_ELEMENT;
    }
    const int tris[] = {0, 3, 5, 3, 1, 4, 5, 3, 4, 5, 4, 2};
    for (int i = 0; i < 12; i++) {
      basic_conn[i] = conn[tris[i]];
    }
  } else if (ltype == TACS_TRI_CUBIC_ELEMENT) {
  } else if (ltype == TACS_QUAD_ELEMENT ||
             ltype == TACS_QUAD_QUADRATIC_ELEMENT ||
             ltype == TACS_QUAD_CUBIC_ELEMENT ||
             ltype == TACS_QUAD_QUARTIC_ELEMENT ||
             ltype == TACS_QUAD_QUINTIC_ELEMENT) {
    int order = 2;
    if (ltype == TACS_QUAD_QUADRATIC_ELEMENT) {
      order = 3;
    } else if (ltype == TACS_QUAD_CUBIC_ELEMENT) {
      order = 4;
    } else if (ltype == TACS_QUAD_QUARTIC_ELEMENT) {
      order = 5;
    } else if (ltype == TACS_QUAD_QUINTIC_ELEMENT) {
      order = 6;
    }

    for (int j = 0, index = 0; j < order - 1; j++) {
      for (int i = 0; i < order - 1; i++, index++) {
        basic_ltypes[index] = TACS_QUAD_ELEMENT;

        basic_conn[4 * index] = conn[i + j * order];
        basic_conn[4 * index + 1] = conn[i + 1 + j * order];
        basic_conn[4 * index + 2] = conn[i + (j + 1) * order];
        basic_conn[4 * index + 3] = conn[i + 1 + (j + 1) * order];
      }
    }
  } else if (ltype == TACS_TETRA_ELEMENT) {
    basic_ltypes[0] = TACS_TETRA_ELEMENT;
    basic_conn[0] = conn[0];
    basic_conn[1] = conn[1];
    basic_conn[2] = conn[2];
    basic_conn[3] = conn[3];
  } else if (ltype == TACS_TETRA_QUADRATIC_ELEMENT) {
  } else if (ltype == TACS_TETRA_CUBIC_ELEMENT) {
  } else if (ltype == TACS_HEXA_ELEMENT ||
             ltype == TACS_HEXA_QUADRATIC_ELEMENT ||
             ltype == TACS_HEXA_CUBIC_ELEMENT ||
             ltype == TACS_HEXA_QUARTIC_ELEMENT ||
             ltype == TACS_HEXA_QUINTIC_ELEMENT) {
    int order = 2;
    if (ltype == TACS_HEXA_QUADRATIC_ELEMENT) {
      order = 3;
    } else if (ltype == TACS_HEXA_CUBIC_ELEMENT) {
      order = 4;
    } else if (ltype == TACS_HEXA_QUARTIC_ELEMENT) {
      order = 5;
    } else if (ltype == TACS_HEXA_QUINTIC_ELEMENT) {
      order = 6;
    }

    for (int k = 0, index = 0; k < order - 1; k++) {
      for (int j = 0; j < order - 1; j++) {
        for (int i = 0; i < order - 1; i++, index++) {
          basic_ltypes[index] = TACS_HEXA_ELEMENT;

          basic_conn[8 * index] = conn[i + j * order + k * order * order];
          basic_conn[8 * index + 1] =
              conn[i + 1 + j * order + k * order * order];
          basic_conn[8 * index + 2] =
              conn[i + (j + 1) * order + k * order * order];
          basic_conn[8 * index + 3] =
              conn[i + 1 + (j + 1) * order + k * order * order];
          basic_conn[8 * index + 4] =
              conn[i + j * order + (k + 1) * order * order];
          basic_conn[8 * index + 5] =
              conn[i + 1 + j * order + (k + 1) * order * order];
          basic_conn[8 * index + 6] =
              conn[i + (j + 1) * order + (k + 1) * order * order];
          basic_conn[8 * index + 7] =
              conn[i + 1 + (j + 1) * order + (k + 1) * order * order];
        }
      }
    }
  } else if (ltype == TACS_PENTA_ELEMENT ||
             ltype == TACS_PENTA_QUADRATIC_ELEMENT ||
             ltype == TACS_PENTA_CUBIC_ELEMENT) {
  }
}
