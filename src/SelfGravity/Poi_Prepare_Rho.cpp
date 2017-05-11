
#include "DAINO.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_Prepare_Rho
// Description :  Prepare the input array "h_Rho_Array_P" for the Poisson solver 
//
// Note        :  Invoke the function "Prepare_PatchGroupData"
//
// Parameter   :  lv             : Targeted refinement level
//                PrepTime       : Targeted physical time to prepare the coarse-grid data
//                h_Rho_Array_P  : Host array to store the prepared data
//                NPG            : Number of patch groups to be prepared at a time
//                PID0_List      : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Poi_Prepare_Rho( const int lv, const double PrepTime, real h_Rho_Array_P[][RHO_NXT][RHO_NXT][RHO_NXT], 
                      const int NPG, const int *PID0_List )
{

   const bool IntPhase_No = false;

   Prepare_PatchGroupData( lv, PrepTime, &h_Rho_Array_P[0][0][0][0], RHO_GHOST_SIZE, NPG, PID0_List, _DENS,
                           OPT__RHO_INT_SCHEME, UNIT_PATCH, NSIDE_26, IntPhase_No );

// subtract the background density to be consistent with the base-level FFT solver
#  pragma omp parallel for
   for (int TID=0; TID<8*NPG; TID++)
   for (int k=0; k<RHO_NXT; k++)
   for (int j=0; j<RHO_NXT; j++)
   for (int i=0; i<RHO_NXT; i++)
      h_Rho_Array_P[TID][k][j][i] -= AveDensity;

} // FUNCTION : Poi_Prepare_Rho



#endif // #ifdef GRAVITY
