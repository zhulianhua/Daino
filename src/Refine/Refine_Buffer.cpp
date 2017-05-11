
#include "DAINO.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Refine_Buffer
// Description :  Construct buffer patches at level "lv+1" according to the flagging result of buffer patches
//                at level "lv" 
//
// Parameter   :  lv          : Targeted level to be refined
//                SonTable    : Table recording the linking index of each buffer father patch to GrandTable
//                GrandTable  : Table recording the patch IDs of grandson buffer patches
//-------------------------------------------------------------------------------------------------------
void Refine_Buffer( const int lv, const int *SonTable, const int *GrandTable )
{

   const int Width = PATCH_SIZE * patch->scale[lv+1]; // scale of a single patch at level "lv+1"
   bool AllocData[8];                                 // allocate data or not
   int *Cr;                                           // corner coordinates


// a. allocate buffer patches
// ------------------------------------------------------------------------------------------------
   for (int s=0; s<26; s++)   
   {

//    determine the buffer patches that must store physical data
      for (int LocalID=0; LocalID<8; LocalID++)
      {
#        ifdef OOC
            AllocData[LocalID] = true;
#        else
         if (  TABLE_01( s, 'x', -1, 0, 1 ) == TABLE_02( LocalID, 'x', -1, 1 )  ||
               TABLE_01( s, 'y', -1, 0, 1 ) == TABLE_02( LocalID, 'y', -1, 1 )  ||
               TABLE_01( s, 'z', -1, 0, 1 ) == TABLE_02( LocalID, 'z', -1, 1 )      )  
            AllocData[LocalID] = false;

         else
            AllocData[LocalID] = true;
#        endif
      }


      for (int PID=patch->NPatchComma[lv][s+1]; PID<patch->NPatchComma[lv][s+2]; PID++)
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
            patch->NPatchComma[lv+1][s+2] += 8;

         } // if ( patch->ptr[0][lv][PID]->flag )
      } // for (int PID=patch->NPatchComma[lv][s+1]; PID<patch->NPatchComma[lv][s+2]; PID++)

      for (int n=s+3; n<28; n++)    patch->NPatchComma[lv+1][n] = patch->num[lv+1];

   } // for (int s=0; s<26; s++)


// check the patch->NPatchComma recording
   if ( patch->NPatchComma[lv+1][27] != patch->num[lv+1] )
      Aux_Error( ERROR_INFO, "patch->NPatchComma[%d][27] (%d) != patch->num[%d] (%d) !!\n",
                 lv+1, patch->NPatchComma[lv+1][27], lv+1, patch->num[lv+1] );


// b. construct relations : son <-> grandson
// ------------------------------------------------------------------------------------------------
   if ( lv < NLEVEL-2 )
   {
      /*
      for (int s=0; s<26; s++)   
      {
         if ( patch->NPatchComma[lv+2][s+1] == patch->NPatchComma[lv+2][s+2] )     continue;

         for (int SonPID0=patch->NPatchComma[lv+1][s+1]; SonPID0<patch->NPatchComma[lv+1][s+2]; SonPID0+=8)
         */
         for (int SonPID0=patch->NPatchComma[lv+1][1]; SonPID0<patch->NPatchComma[lv+1][27]; SonPID0+=8)
         {
            const int FaID      = patch->ptr[0][lv+1][SonPID0]->father - patch->NPatchComma[lv][1];
            const int TempSonID = SonTable[FaID];

            if ( TempSonID != -1 )
            {
               for (int t=0; t<8; t++)
               {
                  const int SonPID    = SonPID0 + t;
                  const int GrandPID0 = GrandTable[ TempSonID + t ];
                 
                  if ( GrandPID0 != -1 )
                  {
                     patch->ptr[0][lv+1][SonPID]->son = GrandPID0;
                 
                     for (int GrandPID=GrandPID0; GrandPID<GrandPID0+8; GrandPID++) 
                        patch->ptr[0][lv+2][GrandPID]->father = SonPID;
                  }
               }
            }

         } // for (int SonPID0=patch->NPatchComma[lv+1][s+1]; SonPID0<patch->NPatchComma[lv+1][s+2]; SonPID0+=8)
//      } // for (int s=0; s<26; s++)
   } // if ( lv < NLEVEL-2 )

} // FUNCTION : Refine_Buffer



#endif // #ifndef SERIAL
