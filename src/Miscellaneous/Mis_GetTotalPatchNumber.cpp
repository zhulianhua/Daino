
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTotalPatchNumber
// Description :  Get the total number of real patches among all DAINO ranks 
//
// Note        :  For the out-of-core computing, we must record the number of patches of each OOC rank
//                in "ooc.NPatch" in advance
//
// Parameter   :  lv : Targeted refinement level
//-------------------------------------------------------------------------------------------------------
void Mis_GetTotalPatchNumber( const int lv )
{

   if ( lv < 0  ||  lv > NLEVEL-1 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );


#  ifdef OOC
   int NPatch_local = 0;
   for (int r=0; r<ooc.NRank; r++)     NPatch_local += ooc.NRealPatch[r][lv];

   MPI_Allreduce( &NPatch_local, &NPatchTotal[lv], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

#  else

   MPI_Allreduce( &patch->NPatchComma[lv][1], &NPatchTotal[lv], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
#  endif

} // FUNCTION : Mis_GetTotalPatchNumber
