/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A STAR PARTICLE TEST)
/
/  written by: Greg Bryan
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
       float *TemperatureUnits, float *TimeUnits,
       float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::TestStarParticleInitializeGrid(float TestStarParticleStarMass,
           float *Initialdt,
           FLOAT TestStarParticleStarVelocity[],
           FLOAT TestStarParticleStarPosition[],
           int NumberOfTestStars, float clusterRadius, 
           char *TestStarInitializationFilename)
{
  /* declarations */

  float CentralMass = 1.0;

  int i, dim, j, k, size, active_size, index, cindex, n;
  float TestInitialdt = *Initialdt;
  int DensNum, TENum, GENum, DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, MetalNum, Vel1Num, Vel2Num, Vel3Num; 
  float *dfield = NULL, *gefield = NULL,
        *v1field = NULL, *v2field = NULL, *v3field = NULL,
        *hifield = NULL, *hiifield = NULL, *heifield = NULL,
        *heiifield = NULL, *heiiifield = NULL, *efield = NULL,
        *h2ifield = NULL, *h2iifield = NULL, *hmfield = NULL, 
        *difield = NULL, *diifield = NULL, *hdifield = NULL,
        *zfield = NULL;
  int gz = NumberOfGhostZones;
  
  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  
  /* Get Units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;
  double MassUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
         &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }


  if (TestStarInitializationFilename){
    /* Initialize a grid based on an input h5 file. Mostly lifted from Grid_PhotonTestInitializeGrid.C - AIW*/
    active_size = 1;
    int ActiveDims[MAX_DIMENSION];
    for (dim = 0; dim < GridRank; dim++) {
      ActiveDims[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
      active_size *= ActiveDims[dim];
    }

    int size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];    

    /* create fields */
    NumberOfBaryonFields = 0;
    FieldType[DensNum = NumberOfBaryonFields++] = Density;
    FieldType[TENum = NumberOfBaryonFields++] = TotalEnergy;
    if (DualEnergyFormalism)
      FieldType[GENum = NumberOfBaryonFields++] = InternalEnergy;
    int ivel = NumberOfBaryonFields;
    FieldType[Vel1Num = NumberOfBaryonFields++] = Velocity1;
    FieldType[Vel2Num = NumberOfBaryonFields++] = Velocity2;
    FieldType[Vel3Num = NumberOfBaryonFields++] = Velocity3;
    if (MultiSpecies) {
      FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
      FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
      FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
      FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
      FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
      FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
      if (MultiSpecies > 1) {
        FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
        FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
        FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
      }
      if (MultiSpecies > 2) {
        FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
        FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
        FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
      }
    }
    if (MetalCooling)
      FieldType[MetalNum    = NumberOfBaryonFields++] = Metallicity; 

    
    this->AllocateGrids();

    // Read Fields
    char *filename;
    hsize_t OutDims[MAX_DIMENSION];
    herr_t h5error;
    hid_t file_id;
    char *delim = "/";
    /* fields we load with this method.  The HDF5 file needs to \
    have one dataset for each field, named as follows. -AIW */

    char *density = "Density";
    char *ge = "GasEnergy";
    char *velocity1 = "x-velocity";
    char *velocity2 = "y-velocity";
    char *velocity3 = "z-velocity";

    char *metal_density = "Metal_Density";

    char *hI_density = "HI_Density";
    char *hII_density = "HII_Density";
    char *heI_density = "HeI_Density";
    char *heII_density = "HeII_Density";
    char *heIII_density = "HeIII_Density";
    char *edensity = "Electron_Density";

    char *h2I_density = "H2I_Density";
    char *h2II_density = "H2II_Density";
    char *hm_density = "HM_Density";

    char *dI_density = "DI_Density";
    char *dII_density = "DII_Density";
    char *hdI_density = "HDI_Density";

    float *dfield = new float [size];
    float *gefield = new float  [size];
    float *v1field = new float  [size];
    float *v2field = new float  [size];
    float *v3field = new float  [size];

    if (MetalCooling) {
      float *zfield = new float  [size];
    }

    if (MultiSpecies) {
      float *hifield = new float  [size];
      float *hiifield = new float  [size];
      float *heifield = new float  [size];
      float *heiifield = new float  [size];
      float *heiiifield = new float  [size];
      float *efield = new float  [size];

      if (MultiSpecies > 1) {
        float *hmfield = new float  [size];
        float *h2ifield = new float  [size];
        float *h2iifield = new float  [size];
        
      }
      if (MultiSpecies > 2){
        float *difield = new float  [size];
        float *diifield = new float  [size];
        float *hdifield = new float  [size];
      }
    }

    /////////////////////////////////////////////////////////////////
    printf("TestStarParticleInitialize\n");
    printf("Active size = %d, OutDim = %d, gz = %d, Grid dim = %ld\n", ActiveDims[0], OutDims[0], gz, 
                  (GridDimension[0]));
    fflush(stdout);
    ////////////////////////////////////////////////////////////////

    filename = strtok(TestStarInitializationFilename, delim);
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      OutDims[GridRank-dim-1] = GridEndIndex[dim] - GridStartIndex[dim] + 1;

    if (file_id == -1) ENZO_FAIL("Error opening field file.");
    
    this->read_dataset(GridRank, OutDims, density, file_id,
           HDF5_REAL, dfield, FALSE, NULL, NULL);

    this->read_dataset(GridRank, OutDims, ge, file_id,
                      HDF5_REAL, gefield, FALSE, NULL, NULL);

    this->read_dataset(GridRank, OutDims, velocity1, file_id,
                      HDF5_REAL, v1field, FALSE, NULL, NULL);

    this->read_dataset(GridRank, OutDims, velocity2, file_id,
                      HDF5_REAL, v2field, FALSE, NULL, NULL);

    this->read_dataset(GridRank, OutDims, velocity3, file_id,
                      HDF5_REAL, v3field, FALSE, NULL, NULL);

    if (MultiSpecies) {
      this->read_dataset(GridRank, OutDims, hI_density, file_id,
                        HDF5_REAL, hifield, FALSE, NULL, NULL);

      this->read_dataset(GridRank, OutDims, hII_density, file_id,
                        HDF5_REAL, hiifield, FALSE, NULL, NULL);

      this->read_dataset(GridRank, OutDims, heI_density, file_id,
                        HDF5_REAL, heifield, FALSE, NULL, NULL);

      this->read_dataset(GridRank, OutDims, heII_density, file_id,
                        HDF5_REAL, heiifield, FALSE, NULL, NULL);
      
      this->read_dataset(GridRank, OutDims, heIII_density, file_id,
                        HDF5_REAL, heiiifield, FALSE, NULL, NULL);

      this->read_dataset(GridRank, OutDims, edensity, file_id,
                        HDF5_REAL, efield, FALSE, NULL, NULL);

      if (MultiSpecies > 1){
        this->read_dataset(GridRank, OutDims, h2I_density, file_id,
                          HDF5_REAL, h2ifield, FALSE, NULL, NULL);

        this->read_dataset(GridRank, OutDims, h2II_density, file_id,
                          HDF5_REAL, h2iifield, FALSE, NULL, NULL);

        this->read_dataset(GridRank, OutDims, hm_density, file_id,
                          HDF5_REAL, hmfield, FALSE, NULL, NULL);
      }
      
      if (MultiSpecies > 2){
        this->read_dataset(GridRank, OutDims, dI_density, file_id,
                          HDF5_REAL, difield, FALSE, NULL, NULL);

        this->read_dataset(GridRank, OutDims, dII_density, file_id,
                          HDF5_REAL, diifield, FALSE, NULL, NULL);

        this->read_dataset(GridRank, OutDims, hdI_density, file_id,
                          HDF5_REAL, hdifield, FALSE, NULL, NULL);
      }
    }

    if (MetalCooling)
      this->read_dataset(GridRank, OutDims, metal_density, file_id,
                        HDF5_REAL, zfield, FALSE, NULL, NULL);

    h5error = H5Fclose(file_id);
    if (h5error == -1) ENZO_FAIL("Error closing initial conditions file.");

    for (n=0; n < size; n++){
      BaryonField[DensNum][n] = dfield[n];
      BaryonField[GENum][n] = gefield[n];
      BaryonField[TENum][n] = BaryonField[GENum][n]; //0.5 * density_field[cindex] 
      BaryonField[Vel1Num][n] = v1field[n];
      BaryonField[Vel2Num][n] = v2field[n];
      BaryonField[Vel3Num][n] = v3field[n];
    }

    if (MetalCooling) {
      for (n=0; n < size; n++)
        BaryonField[MetalNum][n] = zfield[n];
    }

    if (MultiSpecies) {
      for (n=0; n < size; n++){
        BaryonField[HINum][n] = hifield[n];
        BaryonField[HIINum][n] = hiifield[n];
        BaryonField[HeINum][n] = heifield[n];
        BaryonField[HeIINum][n] = heiifield[n];
        BaryonField[HeIIINum][n] = heiiifield[n];
        BaryonField[DeNum][n] = efield[n];
      }
      if (MultiSpecies > 1){
        for (n=0; n < size; n++){
          BaryonField[H2INum][n] = h2ifield[n];
          BaryonField[H2IINum][n] = h2iifield[n];
          BaryonField[HMNum][n] = hmfield[n];
        }
      }
      if (MultiSpecies > 2){
        for (n=0; n < size; n++){
          BaryonField[DINum][n] = difield[n];
          BaryonField[DIINum][n] = diifield[n];
          BaryonField[HDINum][n] = hdifield[n];
        }
      }
    }
    delete [] dfield;
    delete [] gefield;
    delete [] v1field;
    delete [] v2field;
    delete [] v3field;

    if (MetalCooling)
      delete [] zfield;

    if (MultiSpecies){
      delete [] hifield;
      delete [] hiifield;
      delete [] heifield;
      delete [] heiifield;
      delete [] heiiifield;
      delete [] efield;
      if (MultiSpecies > 1) {
        delete [] hmfield;
        delete [] h2ifield;
        delete [] h2iifield;
      }
      if (MultiSpecies > 2) {
        delete [] difield;
        delete [] diifield;
        delete [] hdifield;
      }
    }
  }
  printf("Central Mass: %f \n",CentralMass);

  /* Get Units. */


  /* Set Central Mass in simulation units */

  CentralMass = TestStarParticleStarMass*1.99e33* pow(LengthUnits*CellWidth[0][0],-3.0)/DensityUnits;


  /* Set number of particles for this grid and allocate space. */

  NumberOfParticles = NumberOfTestStars;
  NumberOfParticleAttributes = 4;
  this->AllocateNewParticles(NumberOfParticles);
  printf("Allocated %d particles\n", NumberOfParticles);


  /* Set particle IDs and types */

  for (i = 0; i < NumberOfParticles; i++) {
    ParticleNumber[i] = i;
    if (STARFEED_METHOD(POP3_STAR))
      ParticleType[i] = -1*PARTICLE_TYPE_SINGLE_STAR;
    else
      ParticleType[i] = PARTICLE_TYPE_STAR;
  }
  float p1;

  /* Set central particle. */
  for (i = 0; i <= NumberOfParticles; i++){
    for (dim = 0; dim < GridRank; dim++) {
      if (NumberOfParticles == 1){
        p1 = TestStarParticleStarPosition[dim];
      }else{
        int rng = clusterRadius*200;
        p1 = float(rand() % rng)/100.0+(0.5-clusterRadius);
      }
        ParticlePosition[dim][i] = p1 *
        (DomainLeftEdge[dim]+DomainRightEdge[dim]) + 0.5*CellWidth[0][0];
      ParticleVelocity[dim][i] = TestStarParticleStarVelocity[dim]*1e5*TimeUnits/LengthUnits;
    }
    ParticleMass[i] = CentralMass;
    ParticleAttribute[0][i] = Time+1e-7; //creation time:make sure it is non-zero
    if (STARFEED_METHOD(UNIGRID_STAR)) ParticleAttribute[1][i] = 10.0 * Myr_s/TimeUnits;
    if (STARFEED_METHOD(MOM_STAR))
      if(StarMakerExplosionDelayTime >= 0.0)
        ParticleAttribute[1][i] = 1.0;
      else
        ParticleAttribute[1][i] =10.0 * Myr_s/TimeUnits;
    if (STARFEED_METHOD(MECHANICAL)) {
      if (StarParticleRadiativeFeedback){
        ParticleAttribute[1][i] = 25 * Myr_s/TimeUnits; // radiate for 25 Myr
      }
      else{
       ParticleAttribute[1][i] = 3.953;
      }
    }
  ParticleAttribute[2][i] = 0.0;  // Metal fraction
  ParticleAttribute[3][i] = 0.0;  // metalfSNIa
  }
  return SUCCESS;
}

