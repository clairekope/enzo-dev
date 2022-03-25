/***********************************************************************
/
/  GRID CLASS (Compute and apply cosmic ray diffusion)
/
/  written by:  Munier A. Salem
/  date:        January, 2011
/
/  PURPOSE:  Calculates and applies cosmic ray diffusion 
/  
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <cmath>
#include <math.h> 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "phys_constants.h"
#include "CosmologyParameters.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
            float *TemperatureUnits, float *TimeUnits,
            float *VelocityUnits, double *MassUnits, FLOAT Time);

void RotateVec(const float sint, const float cost, 
               const float sinp, const float cosp, 
               float &v1, float &v2, float &v3);

void InvRotateVec(const float sint, const float cost, 
                  const float sinp, const float cosp, 
                  float &v1, float &v2, float &v3);


int grid::ComputeCRTwoMoment(){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  // Indices and field pointers
  int size=1, id, id3[3], id4[4], idl, idr, i,j,k; 
  float *dens, *E_cr, *F1_cr, *F2_cr, *F3_cr, *B[3];

  float *dx = new float[GridRank];

  dx[0] = CellWidth[0][0];
  dx[1] = (GridRank > 1) ? CellWidth[1][0] : 1.0;
  dx[2] = (GridRank > 2) ? CellWidth[2][0] : 1.0;

  for (int dim = 0; dim < GridRank; dim++) 
    size *= GridDimension[dim];

  // Obtain the current (CR) fields
	int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, 
      B1Num, B2Num, B3Num, PhiNum, CRENum, CRF1Num, CRF2Num, CRF3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum,
                                       B1Num, B2Num, B3Num, PhiNum, 
                                       CRENum, CRF1Num, CRF2Num, CRF3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  dens = BaryonField[DensNum];
  E_cr = BaryonField[CRENum];
  F1_cr = BaryonField[CRF1Num];
  F2_cr = BaryonField[CRF2Num];
  F3_cr = BaryonField[CRF3Num];
  B[0] = BaryonField[B1Num];
  B[1] = BaryonField[B2Num];
  B[2] = BaryonField[B3Num];

  // Get system of units
  float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0;
  float VelocityUnits = 1.0, TimeUnits = 1.0;
  double MassUnits = 1.0;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  // big if true
  double diff_units = ((double)LengthUnits)*LengthUnits/((double)TimeUnits);

  // Local cell variables
  float sigma_str[3], sigma_diff[3], sigma_tot[3], tau[3];
  float *v_cr = new float[size*GridRank];

  // for frame transformation:
  float B_mag, Bxy_mag;
  float *B_angles = new float[size*GridRank];

  // for streaming:
  float dP, B_dot_grad_P, inv_sqrt_rho, va[3], va_mag;
  int stream_sign;
  float *v_str;
  if (TRUE) // check streaming
    v_str = new float[size*GridRank];

  // Placeholder values
  float max_opacity = 1/tiny_number;
  float inv_vmax = 1.0/100.0;
  int cr_sound_flag = 1;

  // Set up start and end indexes to cover all of grid except outermost cells.
  int GridStart[] = {0, 0, 0}, GridEnd[] = {0, 0, 0};
  // KEEP? probably yes, because I need left and right cells to be defined
  for (int dim = 0; dim<GridRank; dim++ ) {
    GridStart[dim] = 1;
    GridEnd[dim] = GridDimension[dim]-1;
  }

  for (k = GridStart[2]; k <= GridEnd[2]; k++)
    for (j = GridStart[1]; j <= GridEnd[1]; j++)
      for (i = GridStart[0]; i <= GridEnd[0]; i++) {
        
        // Location in grid arrays
        id = ELT(i,j,k);

        // Indicies for local arrays; e.g., v_cr and B_angles
        id3[0] = id*3;       // x face
        if (GridRank > 1)
          id3[1] = id3[0]+1; // y face
        if (GridRank > 2)
          id3[2] = id3[0]+2; // z face

        id4[0] = id*4;       // b*sin(theta)
        id4[1] = id4[0]+1;   // b*cos(theta)
        id4[2] = id4[0]+2;   // b*sin(phi)
        id4[3] = id4[0]+3;   // b*cos(phi)

        // Magnetic field strength (3D & 2D)
        B_mag = sqrt(B[0][id] * B[0][id] +
                     B[1][id] * B[1][id] +
                     B[2][id] * B[2][id]);
        Bxy_mag = sqrt(B[0][id]*B[0][id]+
                       B[1][id]*B[1][id]);

        /* 1. Set sigma_str and sigma_diff for the cell, 
              in the B-field frame & in units of Vmax */

        // TODO check streaming flag
        if (TRUE) {

          B_dot_grad_P = 0;
          va_mag = 0;
          inv_sqrt_rho = 1.0 / sqrt(dens[id]);

          /* 1a. Calculate the streaming velocity & direction */
          /* TODO: replace grad P with esitmate from reconstruction? */
          for (int dim=0; dim<GridRank; ++dim){ 

            // Get index of left and right cells for this dimension
            if (dim == 0){
              idl = ELT(i-1, j, k);
              idr = ELT(i+1, j, k);
            } else if (dim == 1) {
              idl = ELT(i, j-1, k);
              idr = ELT(i, j+1, k);
            } else if (dim == 2){
              idl = ELT(i, j, k-1);
              idr = ELT(i, j, k+1);
            }

            dP = (E_cr[idr] - E_cr[idl])/3.0;
            B_dot_grad_P += B[dim][id] * dP/dx[dim];

            va[dim] = B[dim][id] / inv_sqrt_rho;
            va_mag += va[dim] * va[dim];
          }

          va_mag = sqrt(va_mag);
          stream_sign = -1.0 * sign(B_dot_grad_P);

          /* 1b. Set the streaming velocity & sigma_str from the Alfven speed 
                 in the B-field frame (B vector is x-axis)*/
          for (int dim=0; dim<GridRank; ++dim){
            v_str[id3[dim]] = stream_sign * va[dim];
          }

          if (va_mag < tiny_number) {
            sigma_str[0] = max_opacity;
          } else {
            sigma_str[0] = fabs(B_dot_grad_P) / 
                          (va_mag * 4.0/3.0*E_cr[id] * inv_vmax);
          }
          sigma_str[1] = max_opacity;
          sigma_str[2] = max_opacity;

        } else { // no streaming
          for (int dim=0; dim<3; ++dim){        
            v_str[dim] = 0;
            sigma_str[dim] = max_opacity;
          }
        }

        for (int dim=0; dim<GridRank; ++dim){ 

          /* 1c. Set the diffusion, sigma_diff. Can be anisotropic,
                or potentially depend on the local plasma properties
                (though this is not yet implemented) */
          // TODO check anisotropy flag
          sigma_diff[dim] = tiny_number; // fucking figure out units later bro

          /* 2. Calculate the transport velocity, v_cr, across each dimension */
        
          /* 2a. Use the cross sections to estimate an optical depth */
          sigma_tot[dim] = 1.0 / (1.0/sigma_str[dim] + 1.0/sigma_diff[dim]);
          tau[dim] = sigma_tot[dim] * dx[dim];
          tau[dim] = tau[dim] * tau[dim]; // /(2*P_xx) ?

          /* 2b. Set the transport velocity from the optical depth */
          if (tau[dim] < 1e-3)
            v_cr[id3[dim]] = sqrt(1.0 - 0.5*tau[dim]);
          else
            v_cr[id3[dim]] = sqrt((1.0 - exp(-tau[dim])) / tau[dim]);
        
        } // End loop over dimensions

        /* 2c. Find angles between B vector and coordinate axes */
        if (B_mag > tiny_number) {
          B_angles[id4[0]] = Bxy_mag/B_mag;
          B_angles[id4[1]] = B[2][id]/B_mag;
        } else {
          B_angles[id4[0]] = 1.0;
          B_angles[id4[1]] = 0.0;
        }

        if (Bxy_mag > tiny_number) {
          B_angles[id4[2]] = B[1][id]/Bxy_mag;
          B_angles[id4[3]] = B[0][id]/Bxy_mag;
        } else {
          B_angles[id4[2]] = 0.0;
          B_angles[id4[3]] = 1.0;
        }

        /* 2d. Rotate v_cr to the simulation frame & finalize */
        InvRotateVec(B_angles[id4[0]], B_angles[id4[1]], B_angles[id4[2]], B_angles[id4[3]], 
                      v_cr[id3[0]], v_cr[id3[1]], v_cr[id3[2]]);

        for (int dim=0; dim<GridRank; ++dim) {
          v_cr[id3[dim]] = fabs(v_cr[id3[dim]]); // each dim is speed across face; no direction needed

          // Add CR sound speed for stability, if specified
          v_cr[id3[dim]] += cr_sound_flag * sqrt( (4.0/3.0) * E_cr[id]/3.0 / dens[id] );
        }

        /* 3. Reconstruct cell faces? */


	} // triple for loop

  delete [] B_angles;
  delete [] v_cr;
  if (TRUE) // check streaming
    delete [] v_str;

  return SUCCESS;  
}

void RotateVec(const float sint, const float cost, 
               const float sinp, const float cosp, 
               float &v1, float &v2, float &v3) {
  /* 
    Rotate from simulation frame to B-field frame,
    with B vector defining the x-axis

    vel1, vel2, vel3 are input
    v1, v2, v3 are output
    The two rotation matrices are
    
    R_1=
    [cos_p  sin_p 0]
    [-sin_p cos_p 0]
    [0       0    1]

    R_2=
    [sin_t  0 cos_t]
    [0      1    0]
    [-cos_t 0 sin_t]
  */

  // First apply R1, then apply R2
  float newv1 =  cosp * v1 + sinp * v2;
  v2 = -sinp * v1 + cosp * v2;

  // now apply R2
  v1 =  sint * newv1 + cost * v3;
  float newv3 = -cost * newv1 + sint * v3;
  v3 = newv3;
}

void InvRotateVec(const float sint, const float cost, 
                  const float sinp, const float cosp, 
                  float &v1, float &v2, float &v3) {
  /* 
    Rotate from B-field frame to simulation frame

    vel1, vel2, vel3 are input
    v1, v2, v3 are output
    The two rotation matrices are

    R_1^-1=
    [cos_p  -sin_p 0]
    [sin_p cos_p 0]
    [0       0    1]

    R_2^-1=
    [sin_t  0 -cos_t]
    [0      1    0]
    [cos_t 0 sin_t]
  */

  // First apply R2^-1, then apply R1^-1
  float newv1 = sint * v1 - cost * v3;
  v3 = cost * v1 + sint * v3;

  // now apply R1^-1
  v1 = cosp * newv1 - sinp * v2;
  float newv2 = sinp * newv1 + cosp * v2;
  v2 = newv2;
}
