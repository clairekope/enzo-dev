/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (SETS OUTFLOW BOUNDARY CONDITIONS FOR GAL SIM)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Munier Salem
/  date:       August, 2013
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <cmath>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "phys_constants.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */
 
int FindField(int f, int farray[], int n);

int GetUnits(float *DensityUnits, float *LengthUnits,
            float *TemperatureUnits, float *TimeUnits,
            float *VelocityUnits, double *MassUnits, FLOAT time);
 
// Set the Left BoundaryValue of the chosen wave direction (set by
// GalaxySimulationRPSWindSpeed) to the appropriate inflow boundary condition.
 
int ExternalBoundary::SetGalaxySimulationBoundary(FLOAT time, class grid *Grid)
{
  if( GalaxySimulationRPSWind == 0 && GalaxySimulationInflow == 0) return SUCCESS;

  /* declarations */

  int i, j, dim, index;
  int NumberOfZones[MAX_DIMENSION], Offset[MAX_DIMENSION];
  float deltime, distance, pos[MAX_DIMENSION];
  const float TwoPi = 6.283185;
 
  /* Compute size of entire mesh. */
 
  int size = 1;
  for (dim = 0; dim < BoundaryRank; dim++)
    size = size*BoundaryDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  /* Determine if we're using metallicity (as a color field) */
  int MetalNum = FindField(Metallicity,BoundaryFieldType,NumberOfBaryonFields);
  int UseMetallicityField = (MetalNum == -1) ? 0 : 1;
  
 
  if (GalaxySimulationRPSWind > 0 ) { 

    /* set the appropriate BoundaryValues on the left side */
 
    for (dim = 0; dim < BoundaryRank; dim++)
      if (BoundaryDimension[dim] != 1) {
  
        /* If the BoundaryValue fields are missing, create them. */
  
        for (int field = 0; field < NumberOfBaryonFields; field++)
          if (BoundaryValue[field][dim][0] == NULL)
            BoundaryValue[field][dim][0] = new float[size/BoundaryDimension[dim]];
  
        /* Compute quantities needed for boundary face loop (below). */
  
        int dim1, dim2;
        dim1 = (dim == 0) ? 1 : 0;
        dim2 = dim1 + 1;
        dim2 = (dim2 == dim) ? dim2+1 : dim2;
        for (i = 0; i < 3; i++) {
          NumberOfZones[i] = max(BoundaryDimension[i] - 2*NumberOfGhostZones,1);
          Offset[i]        = min(NumberOfGhostZones, BoundaryDimension[i]) - 1;
        }
        pos[dim] = 0.0;
  
        /* Loop over the boundary face */

        for (i = 0; i < BoundaryDimension[dim1]; i++)
          for (j = 0; j < BoundaryDimension[dim2]; j++) {
        
            /* Compute the index into the boundary value. */
        
            index = j*BoundaryDimension[dim1] + i;

            // update bndry type (needed for restart runs)
            for( int field = 0 ; field < NumberOfBaryonFields; ++field ){
              BoundaryType[field][dim][0][index] = inflow;  // left bnd
              BoundaryType[field][dim][1][index] = outflow; // right bnd
            }
        
            /* Find the 3D vector from the corner to the current location. */
            /* now supports nonzero left edge */
            pos[dim1] = DomainLeftEdge[dim1]
                      + (float(i-Offset[dim1]))/float(NumberOfZones[dim1]) 
                      * (DomainRightEdge[dim1]-DomainLeftEdge[dim1]);
            pos[dim2] = DomainLeftEdge[dim2]
                      + (float(i-Offset[dim2]))/float(NumberOfZones[dim2]) 
                      * (DomainRightEdge[dim2]-DomainLeftEdge[dim2]);
        
            /* Compute the distance along the wave propogation vector
            *    |d| = |v_wind . x|/|v_wind| */ 		

            float vMag = sqrt(POW(GalaxySimulationRPSWindVelocity[0],2.0) +
                              POW(GalaxySimulationRPSWindVelocity[1],2.0) +
                              POW(GalaxySimulationRPSWindVelocity[2],2.0) );
            distance = 0.0;
            if( vMag > 0.0 )
              distance = fabs(GalaxySimulationRPSWindVelocity[0]*pos[0] +
                              GalaxySimulationRPSWindVelocity[1]*pos[1] +
                              GalaxySimulationRPSWindVelocity[2]*pos[2] )/vMag;

            /* Find the difference between the current time and the time at
              which the wave will reach this point. */
            static int hasNotArrived = 1;
            if( hasNotArrived ) 
              deltime = time - distance/GalaxySimulationRPSWindShockSpeed - GalaxySimulationRPSWindDelay;
            else
              deltime = time - distance/vMag - GalaxySimulationRPSWindDelay; // fluid travels at bulk speed
            if( hasNotArrived && deltime > 0 ) hasNotArrived = 0;

            if( 1 == GalaxySimulationRPSWind ){

              /* Update bounds with simple shock wind */

              if (deltime > 0.0) {  // Shock has arrived, set post-shock values
                BoundaryValue[DensNum][dim][0][index] = GalaxySimulationRPSWindDensity;
                BoundaryValue[TENum][dim][0][index] = GalaxySimulationRPSWindTotalEnergy;
                BoundaryValue[Vel1Num][dim][0][index] = GalaxySimulationRPSWindVelocity[0];
                if (BoundaryRank > 1)
                  BoundaryValue[Vel2Num][dim][0][index] = GalaxySimulationRPSWindVelocity[1];
                if (BoundaryRank > 2)
                  BoundaryValue[Vel3Num][dim][0][index] = GalaxySimulationRPSWindVelocity[2];
              } else { // If not, set pre-shock values
                BoundaryValue[DensNum][dim][0][index] = GalaxySimulationPreWindDensity;
                BoundaryValue[TENum][dim][0]  [index] = GalaxySimulationPreWindTotalEnergy;
                BoundaryValue[Vel1Num][dim][0][index] = GalaxySimulationPreWindVelocity[0];
                if (BoundaryRank > 1)
                  BoundaryValue[Vel2Num][dim][0][index] = GalaxySimulationPreWindVelocity[1];
                if (BoundaryRank > 2)
                  BoundaryValue[Vel3Num][dim][0][index] = GalaxySimulationPreWindVelocity[2];
              }

            } else if( 2 == GalaxySimulationRPSWind ){

              /* Update Bounds w/ table of density and wind velocity components */
          
              static double *ICMDensityTable, *ICMTotalEnergyTable, *ICMVelocityXTable, 
                *ICMVelocityYTable, *ICMVelocityZTable, *ICMTimeTable;	
              static int loadTable = 1,ICMTableSize=0; // only when processor revs up
              if( loadTable ){
          
                /* Find units */
                float DensityUnits,LengthUnits,TemperatureUnits,TimeUnits,VelocityUnits,temperature;
                double MassUnits;
                GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                        &TimeUnits, &VelocityUnits, &MassUnits, time);

                char filename[] = "ICMinflow_data.in"; char line[MAX_LINE_LENGTH]; FILE *fptr;
                if ((fptr = fopen(filename, "r")) == NULL) ENZO_FAIL("Erroring opening ICMinflow_data.in");

                int f_index = 0;
                while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
                  if (line[0] == 'N') {
                    (sscanf(line, "NumberOfSteps = %"ISYM"\n", &ICMTableSize));
                    ICMTimeTable        = new double[ICMTableSize];
                    ICMDensityTable     = new double[ICMTableSize];
                    ICMTotalEnergyTable = new double[ICMTableSize];
                    ICMVelocityXTable   = new double[ICMTableSize];
                    ICMVelocityYTable   = new double[ICMTableSize];
                    ICMVelocityZTable   = new double[ICMTableSize];
                  }
                  if (line[0] != '#' && line[0] != 'N') {
                    
                    // read values
                    sscanf(line,"%lf %lf %lf %lf %lf %lf",&ICMTimeTable[f_index],&ICMDensityTable[f_index], 
                          &temperature,&ICMVelocityXTable[f_index],&ICMVelocityYTable[f_index],&ICMVelocityZTable[f_index]);
                    
                    // convert to code units
                    ICMTimeTable[f_index]        /= TimeUnits;
                    ICMDensityTable[f_index]     /= DensityUnits;
                    ICMVelocityXTable[f_index]   /= (LengthUnits/TimeUnits); 
                    ICMVelocityYTable[f_index]   /= (LengthUnits/TimeUnits); 
                    ICMVelocityZTable[f_index]   /= (LengthUnits/TimeUnits); 
                    ICMTotalEnergyTable[f_index] = temperature/TemperatureUnits/((Gamma-1.0)*0.6);

                    if (HydroMethod != 2) {
                      ICMTotalEnergyTable[f_index] += 0.5*(   pow(ICMVelocityXTable[0],2)
                                                            + pow(ICMVelocityYTable[1],2)
                                                            + pow(ICMVelocityZTable[2],2));
                    }

                    f_index++;
                  } // end data-line if
                } // end file read while

                fclose(fptr); 
                loadTable = 0;
              }// end load table if

              if(deltime < 0.0){

                // use default pre-wind values
                BoundaryValue[DensNum][dim][0][index] = GalaxySimulationPreWindDensity;
                BoundaryValue[TENum  ][dim][0][index] = GalaxySimulationPreWindTotalEnergy;
                BoundaryValue[Vel1Num][dim][0][index] = GalaxySimulationPreWindVelocity[0];
                if (BoundaryRank > 1)
                  BoundaryValue[Vel2Num][dim][0][index] = GalaxySimulationPreWindVelocity[1];
                if (BoundaryRank > 2)
                  BoundaryValue[Vel3Num][dim][0][index] = GalaxySimulationPreWindVelocity[2];

              } else {

                /* interpolate w/ lookup table */

                float t_ratio,v1,v2;
                int i1,i2=-1;
                
                // find times that bracket what we want
                while( ++i2 < ICMTableSize ) if( deltime < ICMTimeTable[i2] ) break;
                i1 = i2-1;
                
                if(i2<ICMTableSize)
                  t_ratio = (deltime - ICMTimeTable[i2])/(ICMTimeTable[i1] - ICMTimeTable[i2]);
                else { // if beyond final time, just use final val
                  i2--;
                  t_ratio = 1.0;
                }

              
                BoundaryValue[DensNum][dim][0][index] = t_ratio*ICMDensityTable[i1]
                                                        + (1.0-t_ratio)*ICMDensityTable[i2];
                BoundaryValue[TENum  ][dim][0][index] = t_ratio*ICMTotalEnergyTable[i1] 
                                                        + (1.0-t_ratio)*ICMTotalEnergyTable[i2];
                BoundaryValue[Vel1Num][dim][0][index] = t_ratio*ICMVelocityXTable[i1]
                                                        + (1.0-t_ratio)*ICMVelocityXTable[i2];
                if (BoundaryRank > 1)
                  BoundaryValue[Vel2Num][dim][0][index] = t_ratio*ICMVelocityYTable[i1]
                                                          + (1.0-t_ratio)*ICMVelocityYTable[i2];
                if (BoundaryRank > 2)
                  BoundaryValue[Vel3Num][dim][0][index] = t_ratio*ICMVelocityZTable[i1]
                                                          + (1.0-t_ratio)*ICMVelocityZTable[i2];

                // update RPS Wind Vector for time delay calc
                if( index == 0.0 ){
                  GalaxySimulationRPSWindVelocity[0] = BoundaryValue[Vel1Num][dim][0][index];
                  GalaxySimulationRPSWindVelocity[1] = BoundaryValue[Vel2Num][dim][0][index];
                  GalaxySimulationRPSWindVelocity[2] = BoundaryValue[Vel3Num][dim][0][index];
                }
          
              }    
            } else {
              ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: GalaxySimulationRPSWind choice invalid");
            }
          } // end loop over boundary slice
      } // end loop over boundary directions
  } // end RPS Wind
  
  if (GalaxySimulationInflow > 0) {

    /* Find units */
    float DensityUnits,LengthUnits,TemperatureUnits,TimeUnits,VelocityUnits,MassUnits;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
             &TimeUnits, &VelocityUnits, &MassUnits, time);

    int target_dim = GalaxySimulationInflowFace%3;
    int target_face = GalaxySimulationInflowFace/3;

    /* Grab info from the abutting grid */
    int grid_ysize, grid_zsize, cell_hwidth;
    float *grid_dens, *grid_tote, *grid_gase, *grid_vel1, *grid_vel2, *grid_vel3;
    grid_ysize = Grid->GetGridDimension(1);
    grid_zsize = Grid->GetGridDimension(2);
    cell_hwidth = 0.5*Grid->GetCellWidth(0, 0); // they're all the same
    grid_dens = Grid->AccessDensity();
    grid_tote = Grid->AccessTotalEnergy();
    grid_gase = Grid->AccessGasEnergy();
    grid_vel1 = Grid->AccessVelocity1();
    grid_vel2 = Grid->AccessVelocity2();
    grid_vel3 = Grid->AccessVelocity3();
    // TODO will need to check for CRs probably

    if (GalaxySimulationInflow == 2) {
      /* TODO Check if inflow is on; update timer state */
    }

    /* Check validity of desired face */
    if (target_dim > BoundaryRank-1)
      ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: GalaxySimulationInflowFace invalid for simulation dimensionality\n");

    if (BoundaryDimension[target_dim] != 1) {

      /* If the BoundaryValue fields are missing, create them. */

      for (int field = 0; field < NumberOfBaryonFields; field++)
        if (BoundaryValue[field][target_dim][target_face] == NULL)
          BoundaryValue[field][target_dim][target_face] = new float[size/BoundaryDimension[target_dim]];

      /* Compute quantities needed for boundary face loop (below). */
      /* Note that we will be computing position differently than with RPS wind */

      int dim1, dim2;
      int grid_index, face_index;

      dim1 = (target_dim == 0) ? 1 : 0;
      dim2 = dim1 + 1;
      dim2 = (dim2 == target_dim) ? dim2+1 : dim2;
      for (i = 0; i < 3; i++) {
        NumberOfZones[i] = max(BoundaryDimension[i] - 2*NumberOfGhostZones,1);
        Offset[i]        = min(NumberOfGhostZones, BoundaryDimension[i]);
      }

      pos[target_dim] = (target_face == 0) ? DomainLeftEdge[target_dim]+cell_hwidth : DomainRightEdge[target_dim]-cell_hwidth;
      face_index = (target_face == 0) ? Offset[target_dim] : Grid->GetGridEndIndex(target_dim); // the abutting edge of nearby active grid

      assert(Offset[target_dim] == Grid->GetGridStartIndex(target_dim));

      /* Loop over the boundary face */

      for (i = 0; i < BoundaryDimension[dim1]; i++)
        for (j = 0; j < BoundaryDimension[dim2]; j++) {
      
          /* Compute the index into the boundary value. */
      
          index = j*BoundaryDimension[dim1] + i;
          
          /* what grid cell abuts this boundary cell? (x*ylength+y)*zlength+z */
          
          grid_index = (dim1 == 0) ? i*grid_ysize : face_index*grid_ysize;
          grid_index += (dim == 0) ? i : 0;
          grid_index += (dim == 1) ? face_index : 0;
          grid_index += (dim == 2) ? j : 0;
          grid_index *= grid_zsize;
          grid_index += (dim2 == 2) ? j : face_index;

          // update bndry type (needed for restart runs)
          for( int field = 0 ; field < NumberOfBaryonFields; ++field ){
            BoundaryType[field][target_dim][target_face][index] = inflow;  // target boundary
            BoundaryType[field][target_dim][target_face^1][index] = outflow; // other boundary
          }

          /* Find current position in box coordinates */
        
          pos[dim1] = DomainLeftEdge[dim1] + cell_hwidth
                    + (float(i-Offset[dim1]))/float(NumberOfZones[dim1]) 
                    * (DomainRightEdge[dim1]-DomainLeftEdge[dim1]);
          pos[dim2] = DomainLeftEdge[dim2] + cell_hwidth
                    + (float(i-Offset[dim2]))/float(NumberOfZones[dim2]) 
                    * (DomainRightEdge[dim2]-DomainLeftEdge[dim2]);

          /* Determine if current cell is within desired circle */
          if (POW(pos[dim1] - GalaxySimulationInflowCenter[dim1], 2) 
            + POW(pos[dim2] - GalaxySimulationInflowCenter[dim2], 2)
            < POW(GalaxySimulationInflowRadius*kpc_cm/LengthUnits, 2)) {
              
            BoundaryValue[DensNum][target_dim][target_face][index] = GalaxySimulationInflowDensity/DensityUnits;
            BoundaryValue[TENum][target_dim][target_face][index]   = GalaxySimulationInflowTemperature/TemperatureUnits 
                                                                   / ((Gamma - 1.0)*0.6);
            if (DualEnergyFormalism)
              BoundaryValue[GENum][target_dim][target_face][index] = GalaxySimulationInflowTemperature/TemperatureUnits 
                                                                   / ((Gamma - 1.0)*0.6);
            BoundaryValue[Vel1Num][target_dim][target_face][index] = 0.0;
            BoundaryValue[Vel2Num][target_dim][target_face][index] = 0.0;
            BoundaryValue[Vel3Num][target_dim][target_face][index] = 0.0;

          } else {
            /* Duplicate outer face of active grid zone */
            BoundaryValue[DensNum][target_dim][target_face][index] = grid_dens[grid_index];
            BoundaryValue[TENum][target_dim][target_face][index] = grid_tote[grid_index];
            if (DualEnergyFormalism)
              BoundaryValue[GENum][target_dim][target_face][index] = grid_gase[grid_index];
            BoundaryValue[Vel1Num][target_dim][target_face][index] = grid_vel1[grid_index];
            BoundaryValue[Vel2Num][target_dim][target_face][index] = grid_vel2[grid_index];
            BoundaryValue[Vel3Num][target_dim][target_face][index] = grid_vel3[grid_index];
          }
        }
    } else {
      ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: desired inflow face has size 1")
    }
  } // end inflow blob

  // update metallicity field
  if( UseMetallicityField )
    BoundaryValue[MetalNum][dim][0][index] = GalaxySimulationGasHaloMetallicity;

  if( BoundaryValue[DensNum][dim][0][index] < 0.0 ) 
    ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: Negative Density");
  if( BoundaryValue[TENum][dim][0][index] < 0.0 ) 
    ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: Negative Total Energy");

  if( BoundaryValue[DensNum][dim][0][index] != BoundaryValue[DensNum][dim][0][index] )  
    ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: Density NaN");
  if( BoundaryValue[TENum][dim][0][index] != BoundaryValue[TENum][dim][0][index] )  
    ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: Total Energy NaN");


  return SUCCESS;
 
}  