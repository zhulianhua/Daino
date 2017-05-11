
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_AllocateFluxArray
// Description :  Allocate flux arrays for the coarse-grid patches (at level lv ) adjacent to the 
//                coarse-fine boundaries (including the buffer patches)
//
// Parameter   :  lv : Coarse-grid level
//-------------------------------------------------------------------------------------------------------
void Flu_AllocateFluxArray( const int lv )
{

// check
   if ( !patch->WithFlux )
      Aux_Message( stderr, "WARNING : why invoking %s when patch->WithFlux is off ??\n", __FUNCTION__ );


// deallocate the flux arrays allocated previously
#  pragma omp parallel for
   for (int PID=0; PID<patch->NPatchComma[lv][7]; PID++)     patch->ptr[0][lv][PID]->fdelete();


// allocate flux arrays for the real patches
   if ( patch->NPatchComma[lv+1][7] != 0 )
   {
#     pragma omp parallel for
      for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
      {
         if ( patch->ptr[0][lv][PID]->son == -1 )
         {
            for (int s=0; s<6; s++)
            {
               if ( patch->ptr[0][lv][PID]->sibling[s] != -1 )
               {
                  if ( patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[s] ]->son != -1 )
                     patch->ptr[0][lv][PID]->fnew( s );
               }
            }
         }
      }
   }


// allocate flux arrays for the buffer patches
   if ( patch->NPatchComma[lv+1][7] != 0 )   Flu_AllocateFluxArray_Buffer( lv );

   
// get the PIDs for sending/receiving fluxes to/from neighboring ranks
   Buf_RecordExchangeFluxPatchID( lv );

} // Flu_AllocateFluxArray
