
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Refine
// Description :  Perform the refine operation during initialization (including the buffer patches)
//
// Note        :  Allocate patches at level "lv+1"
//
// Parameter   :  lv : Targeted level to be refined (lv --> lv+1)
//-------------------------------------------------------------------------------------------------------
void Init_Refine( const int lv )
{

   if ( lv == NLEVEL-1 )   Aux_Error( ERROR_INFO, "refine the maximum level !!\n" );


   const int Width = PATCH_SIZE*patch->scale[lv+1];
   bool AllocData[8];         // allocate data or not
   int *Cr;    

   for (int m=0; m<27; m++)   
   {

//    all real patches must store physical data 
      if ( m == 0 )
         for (int LocalID=0; LocalID<8; LocalID++)    AllocData[LocalID] = true;

//    only the outer buffer patches do NOT need to store physical data 
      else
      {
         for (int LocalID=0; LocalID<8; LocalID++)
         {
#           ifdef OOC
               AllocData[LocalID] = true;
#           else
            if (  TABLE_01( m-1, 'x', -1, 0, 1 ) == TABLE_02( LocalID, 'x', -1, 1 )  ||
                  TABLE_01( m-1, 'y', -1, 0, 1 ) == TABLE_02( LocalID, 'y', -1, 1 )  ||
                  TABLE_01( m-1, 'z', -1, 0, 1 ) == TABLE_02( LocalID, 'z', -1, 1 )      )  
               AllocData[LocalID] = false;

            else
               AllocData[LocalID] = true;
#           endif
         }
      }


      for (int PID=patch->NPatchComma[lv][m]; PID<patch->NPatchComma[lv][m+1]; PID++)
      {
         if ( patch->ptr[0][lv][PID]->flag )
         {

//          construct relation : father -> child   
            patch->ptr[0][lv][PID]->son = patch->num[lv+1];


//          allocate child patches and construct relation : child -> father
            Cr = patch->ptr[0][lv][PID]->corner;

            patch->pnew( lv+1, Cr[0],       Cr[1],       Cr[2],       PID, AllocData[0], AllocData[0] );
            patch->pnew( lv+1, Cr[0]+Width, Cr[1],       Cr[2],       PID, AllocData[1], AllocData[1] );
            patch->pnew( lv+1, Cr[0],       Cr[1]+Width, Cr[2],       PID, AllocData[2], AllocData[2] );
            patch->pnew( lv+1, Cr[0],       Cr[1],       Cr[2]+Width, PID, AllocData[3], AllocData[3] );
            patch->pnew( lv+1, Cr[0]+Width, Cr[1]+Width, Cr[2],       PID, AllocData[4], AllocData[4] );
            patch->pnew( lv+1, Cr[0],       Cr[1]+Width, Cr[2]+Width, PID, AllocData[5], AllocData[5] );
            patch->pnew( lv+1, Cr[0]+Width, Cr[1],       Cr[2]+Width, PID, AllocData[6], AllocData[6] );
            patch->pnew( lv+1, Cr[0]+Width, Cr[1]+Width, Cr[2]+Width, PID, AllocData[7], AllocData[7] );

            
//          record the number of buffer patches in each sibling direction
            patch->NPatchComma[lv+1][m+1] += 8;

         } // if ( patch->ptr[0][lv][PID]->flag )
      } // for (int PID=patch->NPatchComma[lv][s+1]; PID<patch->NPatchComma[lv][s+2]; PID++)

      for (int n=m+2; n<28; n++)    patch->NPatchComma[lv+1][n] = patch->num[lv+1];

   } // for (int m=0; m<27; m++)


// record the numbers of real and data patches for the out-of-core computing
#  ifdef OOC
   ooc.NDataPatch[ooc.Rank][lv  ] = 0;
   ooc.NDataPatch[ooc.Rank][lv+1] = 0;
   ooc.NRealPatch[ooc.Rank][lv+1] = patch->NPatchComma[lv+1][1];

   for (int PID=0; PID<patch->NPatchComma[lv  ][1]; PID++)
      if ( patch->ptr[0][lv  ][PID]->son == -1 )      ooc.NDataPatch[ooc.Rank][lv  ]++; 

   for (int PID=0; PID<patch->NPatchComma[lv+1][1]; PID++)
      if ( patch->ptr[0][lv+1][PID]->son == -1 )      ooc.NDataPatch[ooc.Rank][lv+1]++; 
#  endif

// set up the BounP_IDMap for the level just created
   Buf_RecordBoundaryPatch( lv+1 );

// construct the sibling relation for the level just created (including the buffer patches)
   SiblingSearch( lv+1 );

// get the patch IDs for sending and receiving data between neighboring ranks
   Buf_RecordExchangeDataPatchID( lv+1 );

// allocate flux arrays at the level "lv"
   if ( patch->WithFlux )
   Flu_AllocateFluxArray( lv );

} // FUNCTION : Init_Refine
