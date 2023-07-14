 /*
    Routine actually checks to see whether the input grid 
    is capable of star formation

    07/2019: Azton Wells
 */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#include "phys_constants.h"
    
    int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, float Time);
#define PASS 1;
int FindField(int field, int farray[], int numfields);

int CheckCreationCriteria(float* Density, float* Metals,
                        float* Temperature,float* DMField,
                        float* Vel1, float* Vel2, float* Vel3, float* TotE,
                        float* H2II, float* H2, float* Mu,
                        float* CoolingTime, int* GridDim,
                        float* shieldedFraction, float* freeFallTime, 
                        float* dynamicalTime, int i, int j, int k,  
                        float Time, float* RefinementField, float CellWidth,
                        bool* gridShouldFormStars, bool* notEnoughMetals, 
                        int continuingFormation, int* seedIndex)
{  
    bool use_F2 = true; // use FIRE-2 methods of self-shielded fractionn and virial parameter
    float Zsolar = SolarMetalFractionByMass;
    float maxZ = 0.0;
    bool internal_debug = false;
    int status = PASS;
    float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
                TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
        fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;    
    } 

    float EnergyUnits = DensityUnits * pow(LengthUnits * CellWidth, 3) * pow(VelocityUnits, 2.0); //[g cm^2/s^2] -> code_energy    
    MassUnits = DensityUnits*pow(LengthUnits*CellWidth, 3);
    int index = i+j*GridDim[0]+k*GridDim[0]*GridDim[1];
    int iminus = index-1;
    int iplus = index+1;
    int jminus =i+(j-1)*GridDim[0]+k*GridDim[0]*GridDim[1];
    int jplus = i+(j+1)*GridDim[0]+k*GridDim[0]*GridDim[1];
    int kminus = i+j*GridDim[0]+(k-1)*GridDim[0]*GridDim[1];
    int kplus = i+j*GridDim[0]+(k+1)*GridDim[0]*GridDim[1];
   
    /*
    Checking creation criteria!
    */

    // if this isnt finest grid in this space, continue
    if (RefinementField[index] != 0) status = FAIL;

    /* Baryon overdensity. */
    if (StarMakerOverDensityThreshold > 0){
        
        if (Density[index] < StarMakerOverDensityThreshold) 
        {
            status = FAIL;
        }
    }
    else if (StarMakerOverDensityThreshold < 0) 
    { // checking number density
        float nb = Density[index]*DensityUnits/mh/Mu[index];
        if (nb < -1*StarMakerOverDensityThreshold)
            status = FAIL;
            
    }

    // check that metals exist in enough quantity
    if (MechStarsCriticalMetallicity > -1) {
        if (Metals[index]/Density[index] < MechStarsCriticalMetallicity) // metallicity in absolute
            status = FAIL;
    }

    /* in addition to the converging flow check, we check
        the virial parameter of the gas to see if it is 
        locally gravitationally bound*/

    float div = 0.0; //divergence
    float alpha = 0.0; //virial parameter
    float vfactor= 0.0; //norm of velocity gradient tensor
    float cSound = 0.0; //sound speed

    float dxvx, dxvy, dxvz, dyvx, dyvy, dyvz, dzvx, dzvy, dzvz;
    
    if (StarMakerVelDivCrit) {

        dxvx = (Vel1[iplus] - Vel1[iminus])/2.0;
        dxvy = (Vel2[iplus] - Vel2[iminus])/2.0;
        dxvz = (Vel3[iplus] - Vel3[iminus])/2.0;
        
        dyvx = (Vel1[jplus] - Vel1[jminus])/2.0;
        dyvy = (Vel2[jplus] - Vel2[jminus])/2.0;
        dyvz = (Vel3[jplus] - Vel3[jminus])/2.0;
        
        dzvx = (Vel1[kplus] - Vel1[kminus])/2.0;
        dzvy = (Vel2[kplus] - Vel2[kminus])/2.0;
        dzvz = (Vel3[kplus] - Vel3[kminus])/2.0;

        /* Chck for converging flow */
        div = dxvx+dyvy+dzvz;
        if (div > 0.0) status = FAIL;

    }

    float totalDensity = (Density[index]+DMField[index])*DensityUnits;

    if (StarMakerSelfBoundCrit){

        /* check for virial parameter */
        vfactor = (dxvx*dxvx + dxvy*dxvy + dxvz*dxvz 
                 + dyvx*dyvx + dyvy*dyvy + dyvz*dyvz
                 + dzvx*dzvx + dzvy*dzvy + dzvz*dzvz) / CellWidth / CellWidth;
        
        /* approximate taking gas as monatomic and mu = 0.6*/
        /* Gravitational constant [cm3g-1s-2]*/
        float Gcode = GravConst*DensityUnits*pow(TimeUnits,2);
        float KBcode = kboltz*MassUnits/(LengthUnits*CellWidth)/pow(TimeUnits,2);
        cSound = sqrt(5/3*kboltz*Temperature[index]/mh/0.6)/VelocityUnits;
        alpha = ((vfactor) + pow(cSound/(CellWidth), 2.0))
                / (8.0 * M_PI* Gcode * Density[index]);

        // alternative formula for virial parameter
        float TE = TotE[index] * Density[index] * EnergyUnits; // total energy of cell
        float PE =  GravConst * pow(totalDensity,2) * pow(CellWidth*LengthUnits, 5); // Approx grav PE of cell, taking r = dx
        float AltAlpha = TE / PE; // canonically, 2 KE / PE, but we explicitly include thermal+internal energy in TE
    
        if (alpha < 20 && internal_debug) 
            fprintf(stdout, "STARSS_CHCR: Compare alphas: F3 = %f; Energy method = %f (G, Gcode, rho, mcell, TE, PE = %e %e %f %e %e %e\n", 
                            alpha, AltAlpha, GravConst, Gcode, Density[index], Density[index]*MassUnits/SolarMass, TE, PE);

        if (MechStarsUseSimpleVirialCriterion) {
            if (AltAlpha > 1.0) status = FAIL;
        } else {
            if (alpha > 1.0) status = FAIL;
        }
    }

    /* Is cooling time < dynamical time or temperature < 1e4 */
    if (StarMakerThermalCrit) {
        *dynamicalTime = pow(3.0*pi/32.0/GravConst/totalDensity, 0.5); //seconds
        if (Temperature[index] > StarMakerTemperatureThreshold)
        {
        if (MultiSpecies > 0) 
                return FAIL; //no hot gas forming stars!
        else if (*dynamicalTime/TimeUnits < CoolingTime[index]) 
                return FAIL;   
        }
    }

    /* is gas mass > critical jeans mass? */
    if (StarMakerUseJeansMass) {
        float baryonMass = Density[index]*DensityUnits
                *LengthUnits*LengthUnits*LengthUnits
                *CellWidth*CellWidth*CellWidth
                /SolarMass;
        float IsoSndSpeed = 1.3095e8 * Temperature[index];
        float jeansMass = pi/(6.0*pow(Density[index]*DensityUnits, 0.5))
                *pow(pi*IsoSndSpeed/GravConst, 1.5)/SolarMass;
        
        float altJeans = 2 * pow(cSound*VelocityUnits/1e5/0.2, 3) * pow((Density[index]*DensityUnits/(mh/0.6)/1e3), -0.5);
        
        if (jeansMass > max(baryonMass, 1e3)) status = FAIL;
    }

    if (MechStarsUseAnalyticShieldedFraction && MechStarsUseMeasuredShieldedFraction)
        ENZO_FAIL("MechStars_CheckCreationCriteria: cannot use both MechStarsUseAnalyticShieldedFraction and MechStarsUseMeasuredShieldedFraction");

    if (MechStarsUseAnalyticShieldedFraction == 1) // allows for other analytic formulae
    {

        /* Is self Shielded fraction > 0.0 by Krumholz & Gnedin */
        float gradRho = (Density[index+1]-Density[index-1])
                        *(Density[index+1]-Density[index-1]);
        gradRho += (Density[jplus]-Density[jminus])
                    *(Density[jplus]-Density[jminus]);
        gradRho +=  (Density[kplus]-Density[kminus])
                    *(Density[kplus]-Density[kminus]);
        gradRho = pow(gradRho, 0.5);

        // factors were given in physical units
        float TauFactor = 434.8/*cm**2/g*/ * MassUnits/pow(LengthUnits*CellWidth, 2); // cm**2/g
        float Tau = TauFactor * Density[index] *(CellWidth+Density[index]/gradRho);

        float Phi = 0.756*pow(1+3.1*Metals[index]/Density[index]/Zsolar, 0.365);

        float Psi = 0.6*Tau*(0.01+Metals[index]/Density[index]/Zsolar)/
                    log(1+0.6*Phi+0.01*Phi*Phi);
        *shieldedFraction = 1.0 - 3.0/(1.0+4.0*Psi);

        if (*shieldedFraction < 0) status = FAIL;

    } else if (MechStarsUseMeasuredShieldedFraction) {
        if (MultiSpecies < 2){
            ENZO_FAIL("MechStars_CheckCreationCriteria: MechStarsUseMeasuredShieldedFraction = 1 requires MultiSpecies > 1\n");
        }

        *shieldedFraction = H2II[index] / (H2II[index] + H2[index]);

    } else {
        *shieldedFraction = 1.0; // use full mass if not shielding
    }
    
    *freeFallTime = pow(3*(pi/(32*GravConst*Density[index]*DensityUnits)), 0.5)/TimeUnits; // that theres code-time

    return status;

}
