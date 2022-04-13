/***********************************************************************
/
/  HLLC RIEMANN SOLVER
/
/  written by: Peng Wang
/  date:       May, 2007
/  modified1:
/
/
************************************************************************/

#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ReconstructionRoutines.h"
#include "EOS.h"  

int hlle_cr(float **FluxLine, float **priml, float **primr, int ActiveSize, float *v_cr1)
{
  float Ul[4], Ur[4], Fl[4], Fr[4];
  float rho_l, rho_r, vx_l, vx_r, cs_l, cs_r, vcr_l, vcr_r, 
        ecr_l, ecr_r, fcrx_l, fcrx_r, fcry_l, fcry_r, fcrz_l, fcrz_r;
  float Sp_l, Sm_l, S_l, Sp_r, Sm_r, S_r;

  
  float vmax2 = CRMaxVelocity*CRMaxVelocity;

  for (int n = 0; n < ActiveSize+1; n++) {

    // CR transport velocity computed elsewhere
    // Consistent with zero reconstruction
    vcr_l = v_cr1[n];
    vcr_r = v_cr1[n+1];

    // First, compute Fl and Ul
    rho_l  = priml[ 0][n];
    vx_l   = priml[ 2][n]; // current sweep direction is always set to x
    ecr_l  = priml[ 9][n];
    fcrx_l = priml[10][n];
    fcry_l = priml[11][n];
    fcrz_l = priml[12][n];

    Ul[0] = ecr_l;
    Ul[1] = fcrx_l / vmax2;
    Ul[2] = fcry_l / vmax2;
    Ul[3] = fcrz_l / vmax2;

    Fl[0] = fcrx_l;
    Fl[1] = 1.0/3.0;
    Fl[2] = 0.0; // if you want a non-isotropic pressure tensor,
    Fl[3] = 0.0; // change these!

    cs_l = sqrt(CRgamma * ecr_l/3.0 / rho_l);

    // Full CR transport velocity includes advection speed of the gas
    Sp_l = vx_l + vcr_l + cs_l;
    Sm_l = vx_l + vcr_l - cs_l;
   

    // Then, Fr and Ur
    rho_r  = priml[ 0][n];
    vx_r   = priml[ 2][n];
    ecr_r  = priml[ 9][n];
    fcrx_r = priml[10][n];
    fcry_r = priml[11][n];
    fcrz_r = priml[12][n];

    Ul[0] = ecr_r;
    Ul[1] = fcrx_r / vmax2;
    Ul[2] = fcry_r / vmax2;
    Ul[3] = fcrz_r / vmax2;

    Fl[0] = fcrx_r;
    Fl[1] = 1.0/3.0;
    Fl[2] = 0.0; // if you want a non-isotropic pressure tensor,
    Fl[3] = 0.0; // change these!

    cs_r = sqrt( CRgamma * ecr_r/3.0 / rho_r);

    Sp_r = vx_r + vcr_r + cs_r;
    Sm_r = vx_r + vcr_r - cs_r;


    // Determine region based on overall fastest signal speeds
    S_l = Min(0.0, Sm_l, Sm_r);
    S_r = Max(0.0, Sp_l, Sp_r);

    // Limit signal to Vmax/sqrt(3), where sqrt(3) comes from P = E/3
    S_l = max(S_l, -CRMaxVelocity/sqrt(3.0));
    S_r = min(S_r,  CRMaxVelocity/sqrt(3.0));
    
    for (int field=iCRE; field < NEQ_MHD; ++field)
      FluxLine[field][n] = (S_r*Fl[field]-S_l*Fr[field]+S_l*S_r*(Ur[field]-Ul[field]))/(S_r-S_l);
  }

  return SUCCESS;
}
