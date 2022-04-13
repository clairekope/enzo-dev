/***********************************************************************
/
/  HLLD-PLM MHD SOLVER
/
/  written by: J. S. Oishi
/  date:       April 2011
/  modified1:  
/              
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "EOS.h"
#include "../hydro_rk/ReconstructionRoutines.h"

int plm(float **prim, float **priml, float **primr, int ActiveSize, int Neq);
int cons_plm(float **prim, float **priml, float **primr, int ActiveSize, int Neq, char direc);
inline void plm_point(float &vm1, float &v, float &vp1, float &vl_plm);
int plm_species(float **prim, int is, float **species, float *flux0, int ActiveSize);
int plm_color(float **prim, int is, float **color, float *flux0, int ActiveSize);
int hlld_mhd(float **FluxLine, float **priml, float **primr, float **prim, int ActiveSize);
int hlle_cr(float **FluxLine, float **priml, float **primr, int ActiveSize, float *v_cr1);

int HLLD_PLM_MHD(float **prim, float **priml, float **primr,
		float **species, float **colors,  float **FluxLine, int ActiveSize,
		char direc, int jj, int kk, float *v_cr1)
{
  int idual = (DualEnergyFormalism) ? 1 : 0;

  // compute priml and primr
  if (ConservativeReconstruction == 1)
    cons_plm(prim, priml, primr, ActiveSize, NEQ_MHD-idual, direc);
  else
    plm(prim, priml, primr, ActiveSize, NEQ_MHD-idual);

  // compute FluxLine
  if (hlld_mhd(FluxLine, priml, primr, prim, ActiveSize)==FAIL) {
    return FAIL;
  }

  // compute CR flux, if applicable, with it's dedicated Riemann solver
  if (CRModel > 1) {
    if (hlle_cr(FluxLine, priml, primr, ActiveSize, v_cr1)==FAIL) {
      return FAIL;
    }
  }

  if (NSpecies > 0) {
    plm_species(prim, NEQ_MHD-idual, species, FluxLine[iD], ActiveSize);
    for (int field = NEQ_MHD; field < NEQ_MHD+NSpecies; field++) {
      for (int i = 0; i < ActiveSize+1; i++) {
	      FluxLine[field][i] = FluxLine[iD][i]*species[field-NEQ_MHD][i];
      }
    }
  }

  if (NColor > 0) {
    plm_color(prim, NEQ_MHD-idual, colors, FluxLine[iD], ActiveSize);
    for (int field = NEQ_MHD+NSpecies; field < NEQ_MHD+NSpecies+NColor; field++) {
      for (int i = 0; i < ActiveSize+1; i++) {
      	FluxLine[field][i] = FluxLine[iD][i]*colors[field-NEQ_MHD-NSpecies][i];
      }
    }
  }
  
  return SUCCESS;
}
