#include "TACSToFH5.h"

/*
  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/  

/*
  Create the TACSToFH5 file creation object

  input:
  tacs:        the instance of the TACSAssembler object
  elem_type:   the type of element to be used
  write_flag:  a variable indicating the type of output to write

  For example if:

  write_flag = 
  (Element::OUTPUT_NODES | 
   Element::OUTPUT_DISPLACEMENTS | 
   Element::OUTPUT_STRESSES)

  then the TACSToFH5 object will output the nodes, displacements and
  stresses for each element of type 'elem_type' in the TACSAssembler
  object. Note that this is a bit-wise OR operation.  
*/
TACSToFH5::TACSToFH5( TACSAssembler * _tacs,
                      enum ElementType _elem_type,
                      unsigned int _write_flag ){
  tacs = _tacs;
  tacs->incref();
  elem_type = _elem_type;
  write_flag = _write_flag;

  // 3 coordinates for each of 3 coordinate axes 
  ncoordinates = 9; 

  // Get the number of displacements/stresses
  ndisplacements = 0;
  nstresses = 0;
  nextras = 0;

  max_design_vars = 0;
  int numElements = tacs->getNumElements();
  TACSElement ** elements = tacs->getElements();
  for ( int i = 0; i < numElements; i++ ){
    if (elements[i] && elements[i]->getElementType() == elem_type){
      ndisplacements = elements[i]->numDisplacements();
      nstresses = elements[i]->numStresses();
      nextras = elements[i]->numExtras();
      int ndvs = elements[i]->getNumDesignVars();
      if (ndvs > max_design_vars){
	max_design_vars = ndvs;
      }
    }
  }

  int max_dvs = max_design_vars;
  MPI_Allreduce(&max_dvs, &max_design_vars, 1, MPI_INT, 
		MPI_MAX, tacs->getMPIComm());

  num_design_vars = tacs->getNumDesignVars();

  // Count up the number of values that will be output for each point
  // in the mesh
  nvals = 0;
  if (write_flag & TACSElement::OUTPUT_NODES){
    nvals += 3; // Always 3 coordinates - even for planar problems
  }
  if (write_flag & TACSElement::OUTPUT_DISPLACEMENTS){
    nvals += ndisplacements;
  }
  if (write_flag & TACSElement::OUTPUT_STRAINS){
    nvals += nstresses; // Number of stresses == number of strains 
  }
  if (write_flag & TACSElement::OUTPUT_STRESSES){
    nvals += nstresses;
  }
  if (write_flag & TACSElement::OUTPUT_EXTRAS){
    nvals += nextras;
  }
  if (write_flag & TACSElement::OUTPUT_COORDINATES){
    nvals += ncoordinates;
  }
  if (write_flag & TACSElement::OUTPUT_DESIGN_VARIABLES){
    nvals += max_design_vars;
  }

  // Get a comma separated list of the variable names
  variable_names = getElementVarNames();
  
  // Retrieve the number of components
  num_components = tacs->getNumComponents();
  
  // Allocate space for the component names
  component_names = new char*[ num_components ];
  memset(component_names, 0, num_components*sizeof(char*));

  for ( int k = 0; k < num_components; k++ ){
    char comp_name[128];
    sprintf(comp_name, "Component %d", k);
    setComponentName(k, comp_name);
  }
}

TACSToFH5::~TACSToFH5(){
  tacs->decref();
  
  // Deallocate the comma separated list of variable names
  delete [] variable_names;

  // Deallocate the component names
  for ( int k = 0; k < num_components; k++ ){
    if (component_names[k]){
      delete [] component_names[k];
    }
  }
  delete [] component_names;
}

/*
  Set the component names for the elements
*/
void TACSToFH5::setComponentName( int comp_num, const char * group_name ){
  if (comp_num >= 0 && comp_num < num_components){
    // If the name already exists, over-write it
    if (component_names[comp_num]){
      delete [] component_names[comp_num];
    }

    // Allocate space for the name
    size_t len = strlen(group_name)+1;
    component_names[comp_num] = new char[ len ];

    // Copy the component name
    strcpy(component_names[comp_num], group_name);
  }
}

/*
  Write the data stored in the TACSAssembler object to a file
  
  input:
  load_case: the load case number corresponding to the output
  filename:  the name of the file to create
*/
void TACSToFH5::writeToFile( int load_case, const char * filename ){
  int rank, size;
  MPI_Comm_rank(tacs->getMPIComm(), &rank);
  MPI_Comm_size(tacs->getMPIComm(), &size);

  // Create the FH5 file object for writting
  FH5File * file = new FH5File(tacs->getMPIComm());
  file->incref();

  // Open the file - if possible for writing
  int write_err = file->createFile(filename, component_names,
                                   num_components);

  if (write_err){
    file->decref();
    if (rank == 0){
      fprintf(stderr, "[%d] TACSToFH5 error: Could not create file\n",
              rank);
    }
    return;
  }

  // Write the connectivity
  int *comp_nums = NULL; // The component numbers of each element
  int *csr = NULL;       // Local segment of the global CSR data structure
  int *csr_range = NULL; // node/csr range over the number of processes
  int *node_range = NULL;
  tacs->getOutputConnectivity(elem_type, &comp_nums,
                              &csr, &csr_range, &node_range);

  int con_size = 4;
  if (elem_type == EULER_BEAM || 
      elem_type == TIMOSHENKO_BEAM){
    con_size = 2;
  }
  else if (elem_type == SOLID){
    con_size = 8;
  }

  // Write the component numbers to a zone
  int dim1 = (csr_range[rank+1] - csr_range[rank])/con_size;
  int dim2 = 1;
  char comp_name[] = "components";
  file->writeZoneData(comp_name, FH5File::FH5_INT,
                      comp_name, comp_nums, dim1, dim2);
  if (comp_nums){ delete [] comp_nums; }

  // Write the data to a zone
  dim1 = (csr_range[rank+1] - csr_range[rank])/con_size;
  dim2 = con_size;
  char conn_name[] = "connectivity";
  file->writeZoneData(conn_name, FH5File::FH5_INT,
                      conn_name, csr, dim1, dim2);  
  if (csr){ delete [] csr; }
  if (csr_range){ delete [] csr_range; }

  // Allocate space for the output data - the nodes, displacements, stresses etc.
  int len = nvals*(node_range[rank+1] - node_range[rank]);
  double * data = new double[ len ];
  memset(data, 0, len*sizeof(double));

  // Get the output data from TACS
  tacs->getOutputData(load_case, elem_type,
		      write_flag, data, nvals);
  
  if (write_flag & TACSElement::OUTPUT_DESIGN_VARIABLES){
    // Get the output data for the design variables
    TacsScalar * x = new TacsScalar[ num_design_vars ];
    tacs->getDesignVars(x, num_design_vars);
    tacs->getOutputDesignVarData(x, num_design_vars, elem_type, 
				 &data[nvals-max_design_vars],
				 nvals);
    delete [] x;
  }

  dim1 = node_range[rank+1] - node_range[rank];
  dim2 = nvals;

  char data_name[] = "data"; 
  file->writeZoneData(data_name, FH5File::FH5_DOUBLE, variable_names,
                      data, dim1, dim2);
  delete [] data;

  file->close();
  file->decref();
}

/*
  Create a comma-separated list of the element variable names
*/
char * TACSToFH5::getElementVarNames(){
  // Find the first variable name
  int numElements = tacs->getNumElements();
  TACSElement ** elements = tacs->getElements();
  TACSElement * elem_match = NULL;
  char * elem_vars = NULL;

  for ( int i = 0; i < numElements; i++ ){
    if (elements[i] && elements[i]->getElementType() == elem_type){
      elem_match = elements[i];
      break;
    }
  }

  if (!elem_match){
    int rank;
    MPI_Comm_rank(tacs->getMPIComm(), &rank);
    fprintf(stderr, "[%d] TACSToFH5: Could not find an element match\n",
            rank);
    return elem_vars;
  }

  char * output_names[7] = {NULL, NULL, NULL, 
                            NULL, NULL, NULL, NULL};

  if (write_flag & TACSElement::OUTPUT_NODES){ 
    output_names[0] = new char[ 6 ];
    sprintf(output_names[0], "X,Y,Z");
  }
  if (write_flag & TACSElement::OUTPUT_DISPLACEMENTS){
    int nd = elem_match->numDisplacements();
    size_t str_len = 2;
    for ( int i = 0; i < nd; i++ ){
      str_len += strlen(elem_match->displacementName(i))+1;
    }
    char * temp = new char[ str_len ];
    if (nd > 0){
      strcpy(temp, elem_match->displacementName(0));
      for ( int i = 1; i < nd; i++ ){
	size_t len = strlen(temp);
	sprintf(&(temp[len]), ",%s", elem_match->displacementName(i));
      }
    }
    output_names[1] = temp;
  }
  if (write_flag & TACSElement::OUTPUT_STRAINS){
    int ns = elem_match->numStresses();
    size_t str_len = 2;
    for ( int i = 0; i < ns; i++ ){
      str_len += strlen(elem_match->strainName(i))+1;
    }
    char * temp = new char[ str_len ];
    if (ns > 0){
      strcpy(temp, elem_match->strainName(0));
      for ( int i = 1; i < ns; i++ ){
	size_t len = strlen(temp);
	sprintf(&temp[len], ",%s", elem_match->strainName(i));
      }
    }
    output_names[2] = temp;
  }
  if (write_flag & TACSElement::OUTPUT_STRESSES){
    int ns = elem_match->numStresses();
    size_t str_len = 2;
    for ( int i = 0; i < ns; i++ ){
      str_len += strlen(elem_match->stressName(i))+1;
    }
    char * temp = new char[ str_len ];
    if (ns > 0){
      strcpy(temp, elem_match->stressName(0));
      for ( int i = 1; i < ns; i++ ){
	size_t len = strlen(temp);
	sprintf(&temp[len], ",%s", elem_match->stressName(i));
      }
    }
    output_names[3] = temp;
  }
  if (write_flag & TACSElement::OUTPUT_EXTRAS){
    int ne = elem_match->numExtras();
    size_t str_len = 2;
    for ( int i = 0; i < ne; i++ ){
      str_len += strlen(elem_match->extraName(i))+1;
    }
    char * temp = new char[ str_len ];
    if (ne > 0){
      strcpy(temp, elem_match->extraName(0));
      for ( int i = 1; i < ne; i++ ){
	size_t len = strlen(temp);
	sprintf(&temp[len], ",%s", elem_match->extraName(i));
      }
    }
    output_names[4] = temp;
  }
  if (write_flag & TACSElement::OUTPUT_COORDINATES){
    int ne = 9;
    size_t str_len = 4*ne;
    output_names[5] = new char[ str_len ];
    strcpy(output_names[5], "1x,1y,1z,2x,2y,2z,3x,3y,3z");
  }
  if (write_flag & TACSElement::OUTPUT_DESIGN_VARIABLES){
    size_t str_len = 1 + 5*max_design_vars;   
    char * temp = new char[ str_len ];
    sprintf(temp, "x%d", 0);
    for ( int i = 1; i < max_design_vars; i++ ){
      size_t len = strlen(temp);
      sprintf(&temp[len], ",x%d", i);
    }
    
    output_names[6] = temp;
  }

  // Count up the size of the elem_vars string
  int elem_size = 14; // Extra space for either a comma or \0
  for ( int k = 0; k < 7; k++ ){
    if (output_names[k]){
      elem_size += strlen(output_names[k]);
    }
  }

  elem_vars = new char[ elem_size ];

  // Copy the first zone into the list directly
  int k = 0; 
  for ( ; k < 7; k++ ){
    if (output_names[k]){
      strcpy(elem_vars, output_names[k]);
      k++;
      break;
    }
  }

  // For subsequent non-zero zones - add a comma before adding
  // the remainder of the list
  for ( ; k < 7; k++ ){
    if (output_names[k]){
      int len = strlen(elem_vars);
      sprintf(&elem_vars[len], ",%s", output_names[k]);
    }
  }

  for ( int k = 0; k < 7; k++ ){
    if (output_names[k]){ delete [] output_names[k]; }
  }

  return elem_vars;
}
