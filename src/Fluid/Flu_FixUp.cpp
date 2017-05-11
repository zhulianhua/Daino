
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_FixUp
// Description :  1. Use the corrected coarse-fine boundary fluxes to fix the data at level "lv"
//                2. Use the average data at level "lv+1" to replace the data at level "lv"
//
// Note        :  1. Also include the fluxes from neighbor ranks
//                2. The boundary fluxes must be received in advance by invoking the function "Buf_GetBufferData"
//
// Parameter   :  lv : Targeted refinement level
//                dt : Time interval to advance solution
//-------------------------------------------------------------------------------------------------------
void Flu_FixUp( const int lv, const double dt )
{

   const real Const = dt / patch->dh[lv];

   real (*FluxPtr)[PATCH_SIZE][PATCH_SIZE] = NULL;


// a. correct the coarse-fine boundary fluxes
   if ( OPT__FIXUP_FLUX )
   {
//    check
      if ( !patch->WithFlux )    
         Aux_Error( ERROR_INFO, "patch->WithFlux is off -> no flux array is allocated for OPT__FIXUP_FLUX !!\n" );

#     pragma omp parallel for private( FluxPtr )
      for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
      {
//       a1. sum up the coarse-grid and fine-grid fluxes for the debug mode
#        ifdef DAINO_DEBUG
         for (int s=0; s<6; s++)
         {
            FluxPtr = patch->ptr[0][lv][PID]->flux[s];

            if ( FluxPtr != NULL )
            {
               for (int v=0; v<NCOMP; v++)
               for (int m=0; m<PS1; m++)
               for (int n=0; n<PS1; n++)
                  FluxPtr[v][m][n] += patch->ptr[0][lv][PID]->flux_debug[s][v][m][n];
            }
         }
#        endif


//       a2. correct fluid variables by the difference between coarse-grid and fine-grid fluxes
         if ( NULL != (FluxPtr = patch->ptr[0][lv][PID]->flux[0]) )
         {
            for (int v=0; v<NCOMP; v++)
            for (int k=0; k<PATCH_SIZE; k++)
            for (int j=0; j<PATCH_SIZE; j++)
               patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[v][k][j][           0] -= FluxPtr[v][k][j] * Const;
         }

         if ( NULL != (FluxPtr = patch->ptr[0][lv][PID]->flux[1]) )
         {
            for (int v=0; v<NCOMP; v++)
            for (int k=0; k<PATCH_SIZE; k++)
            for (int j=0; j<PATCH_SIZE; j++)
               patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[v][k][j][PATCH_SIZE-1] += FluxPtr[v][k][j] * Const;
         }

         if ( NULL != (FluxPtr = patch->ptr[0][lv][PID]->flux[2]) )
         {
            for (int v=0; v<NCOMP; v++)
            for (int k=0; k<PATCH_SIZE; k++)
            for (int i=0; i<PATCH_SIZE; i++)
               patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[v][k][           0][i] -= FluxPtr[v][k][i] * Const;
         }

         if ( NULL != (FluxPtr = patch->ptr[0][lv][PID]->flux[3]) )
         {
            for (int v=0; v<NCOMP; v++)
            for (int k=0; k<PATCH_SIZE; k++)
            for (int i=0; i<PATCH_SIZE; i++)
               patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[v][k][PATCH_SIZE-1][i] += FluxPtr[v][k][i] * Const;
         }

         if ( NULL != (FluxPtr = patch->ptr[0][lv][PID]->flux[4]) )
         {
            for (int v=0; v<NCOMP; v++)
            for (int j=0; j<PATCH_SIZE; j++)
            for (int i=0; i<PATCH_SIZE; i++)
               patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[v][           0][j][i] -= FluxPtr[v][j][i] * Const;
         }

         if ( NULL != (FluxPtr = patch->ptr[0][lv][PID]->flux[5]) )
         {
            for (int v=0; v<NCOMP; v++)
            for (int j=0; j<PATCH_SIZE; j++)
            for (int i=0; i<PATCH_SIZE; i++)
               patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[v][PATCH_SIZE-1][j][i] += FluxPtr[v][j][i] * Const;
         }

      } // for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)


//    a3. reset all flux arrays (in both real and buffer patches) to zero for the debug mode
#     ifdef DAINO_DEBUG
#     pragma omp parallel for private( FluxPtr )
      for (int PID=0; PID<patch->NPatchComma[lv][27]; PID++)
      {
         for (int s=0; s<6; s++)
         {
            FluxPtr = patch->ptr[0][lv][PID]->flux[s];
            if ( FluxPtr != NULL )
            {
               for (int v=0; v<NCOMP; v++)
               for (int m=0; m<PS1; m++)
               for (int n=0; n<PS1; n++)
                  FluxPtr[v][m][n] = 0.0;
            }

            FluxPtr = patch->ptr[0][lv][PID]->flux_debug[s];
            if ( FluxPtr != NULL )
            {
               for (int v=0; v<NCOMP; v++)
               for (int m=0; m<PS1; m++)
               for (int n=0; n<PS1; n++)
                  FluxPtr[v][m][n] = 0.0;
            }
         }
      }
#     endif

   } // if ( OPT__FIXUP_FLUX )


// b. average over the data at level "lv+1" to correct the data at level "lv"
   if ( OPT__FIXUP_RESTRICT )    
   {
      Flu_Restrict( lv, patch->FluSg[lv+1], patch->FluSg[lv], NULL_INT, NULL_INT, _FLU );

#     ifdef LOAD_BALANCE 
      LB_GetBufferData( lv, patch->FluSg[lv], NULL_INT, DATA_RESTRICT, _FLU, NULL_INT );
#     endif
   }

} // FUNCTION : Flu_FixUp
