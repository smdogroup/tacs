#ifndef TACS_HEAT_FLUX_H
#define TACS_HEAT_FLUX_H

#include "TACSFunction.h"

/*
  Compute the KS functional of the heat flux on a given face or edge
*/
class HeatFluxIntegral : public TACSFunction {
 public:
  HeatFluxIntegral( TACSAssembler *_tacs, int *_elem_index,
                    int *_surface_index, int _num_elems );
  ~HeatFluxIntegral();
  // Retrieve the name of the function
  // ---------------------------------
  const char *functionName();

  // Create the function context for evaluation
  // ------------------------------------------
  TACSFunctionCtx *createFunctionCtx();

  // Collective calls on the TACS MPI Comm
  // -------------------------------------
  void initEvaluation( EvaluationType ftype );
  void finalEvaluation( EvaluationType ftype );

  // Functions for integration over the structural domain on each thread
  // -------------------------------------------------------------------
  void initThread( double tcoef,
                   EvaluationType ftype,
                   TACSFunctionCtx *ctx );
  void elementWiseEval( EvaluationType ftype,
                        TACSElement *element, int elemNum,
                        const TacsScalar Xpts[], const TacsScalar vars[],
                        const TacsScalar dvars[], const TacsScalar ddvars[],
                        TACSFunctionCtx *ctx );
  void finalThread( double tcoef,
                    EvaluationType ftype,
                    TACSFunctionCtx *ctx );

  // Return the value of the function
  // --------------------------------
  TacsScalar getFunctionValue();

  // State variable sensitivities
  // ----------------------------
  void getElementSVSens( double alpha, double beta, double gamma,
                         TacsScalar *elemSVSens,
                         TACSElement *element, int elemNum,
                         const TacsScalar Xpts[], const TacsScalar vars[],
                         const TacsScalar dvars[], const TacsScalar ddvars[],
                         TACSFunctionCtx *ctx );

  // Design variable sensitivity evaluation
  // --------------------------------------
  void addElementDVSens( double tcoef, TacsScalar *fdvSens, int numDVs,
                         TACSElement *element, int elemNum,
                         const TacsScalar Xpts[], const TacsScalar vars[],
                         const TacsScalar dvars[], const TacsScalar ddvars[],
                         TACSFunctionCtx *ctx );

  // Nodal sensitivities
  // -------------------
  void getElementXptSens( double tcoef, TacsScalar fXptSens[],
                          TACSElement *element, int elemNum,
                          const TacsScalar Xpts[], const TacsScalar vars[],
                          const TacsScalar dvars[], const TacsScalar ddvars[],
                          TACSFunctionCtx *ctx );

 private:
  // The name of the function
  static const char *funcName;

  // The value of the KS weight
  TacsScalar value;

  // The max number of nodes
  int maxNumNodes;

  // List of elements and its associated surfaces
  int *elem_index;
  int *surface_index;
  int num_elems;
};
#endif // TACS_HEAT_FLUX_H
