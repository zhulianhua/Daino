
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_UserCriteria
// Description :  Check if the element (i,j,k) of the input data satisfies the user-defined flag criteria
//
// Note        :  Users can put their favorite flag criteria in this function
//
// Parameter   :  i,j,k       : Indices of the targeted element in the patch ptr[ patch->FluSg[lv] ][lv][PID]
//                lv          : Refinement level of the targeted patch
//                PID         : ID of the targeted patch
//                Threshold   : User-provided threshold for the flag operation, which is loaded from the 
//                              file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_UserCriteria( const int i, const int j, const int k, const int lv, const int PID, const real Threshold )
{

   /*
   const double dh     = patch->dh[lv];                                              // grid size
   const int    scale  = patch->scale[lv];                                           // grid scale
   const real   Pos[3] = { ( patch->ptr[0][lv][PID]->corner[0]/scale + i + 0.5 )*dh, // x,y,z position
                           ( patch->ptr[0][lv][PID]->corner[1]/scale + j + 0.5 )*dh,
                           ( patch->ptr[0][lv][PID]->corner[2]/scale + k + 0.5 )*dh  };  

   const real (*Rho )[PS1][PS1] = patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[DENS];   // density
   const real (*MomX)[PS1][PS1] = patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[MOMX];   // momentum x
   const real (*MomY)[PS1][PS1] = patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[MOMY];   // momentum y
   const real (*MomZ)[PS1][PS1] = patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[MOMZ];   // momentum z
   const real (*Egy )[PS1][PS1] = patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[ENGY];   // total energy
#  ifdef GRAVITY
   const real (*Pot )[PS1][PS1] = patch->ptr[ patch->PotSg[lv] ][lv][PID]->pot;           // potential
#  endif
   */

   bool Flag = false;


// put your favorite flag criteria here
// ##########################################################################################################

// Example 1 : flag if the velocity exceeds the given threshold
/*
   const real Vel = sqrt(  ( MomX[k][j][i]*MomX[k][j][i] + MomY[k][j][i]*MomY[k][j][i] + 
                             MomZ[k][j][i]*MomZ[k][j][i] )  ) / Rho[k][j][i];
   Flag = Vel > Threshold;
*/


// Example 2 : flag if the grid is within a sphere with the radius eqaul to the input "Threshold" and the origin 
//             in the center of the simulation box
/*
   const real Center[3] = { 0.5*patch->BoxSize[0], 0.5*patch->BoxSize[1], 0.5*patch->BoxSize[2] };
   const real dr[3]     = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
   const real Radius    = sqrt( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );

   Flag = Radius < Threshold;
*/

// ##########################################################################################################


   return Flag;

} // FUNCTION : Flag_UserCriteria

