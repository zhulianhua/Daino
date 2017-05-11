
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_Prepare
// Description :  Prepare the input array "Flu_Array_F_In" for the fluid solver 
//
// Note        :  Invoke the function "Prepare_PatchGroupData"
//
// Parameter   :  lv                : Targeted refinement level
//                PrepTime          : Targeted physical time to prepare the coarse-grid data
//                h_Flu_Array_F_In  : Host array to store the prepared data
//                NPG               : Number of patch groups to be prepared at a time
//                PID0_List         : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Flu_Prepare( const int lv, const double PrepTime, real h_Flu_Array_F_In[][FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT], 
                  const int NPG, const int *PID0_List )
{

   const bool IntPhase_No = false;

#  if ( MODEL == ELBDM )
   Prepare_PatchGroupData( lv, PrepTime, &h_Flu_Array_F_In[0][0][0], FLU_GHOST_SIZE, NPG, PID0_List, 
                           _REAL | _IMAG, OPT__FLU_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, OPT__INT_PHASE );
#  else
   Prepare_PatchGroupData( lv, PrepTime, &h_Flu_Array_F_In[0][0][0], FLU_GHOST_SIZE, NPG, PID0_List, 
                           _FLU,          OPT__FLU_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No );
#  endif

} // FUNCTION : Flu_Prepare
