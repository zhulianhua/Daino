#include "GetCube.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Buf_AllocateBufferPatch
// Description :  Allocate the buffer patches at level "lv"
//
// Note        :  Currently it only works for the function "Init_Reload"
//
// Parameter   :  lv : The targeted refinement lv to allocate buffer patches
//-------------------------------------------------------------------------------------------------------
void Buf_AllocateBufferPatch( const int lv )
{

// invoke the function "Buf_AllocateBufferPatch_Base" for the base level
   if ( lv == 0 )
   {
      Buf_AllocateBufferPatch_Base();

      return;
   }


   int BounPID, BuffPID, SonPID, RefinePos, TargetID, NSend[26], NRecv[26], NBuffer[26];
   int *Send_PosList[26], *Recv_PosList[26];
   bool AllocateData[8];   // allocate data or not


   bool AllocateData2[8];     // check whether of not the patch is within the candidate box


// set up the Send_PosList by scanning over the BounP_IDList at level "lv-1"
   for (int s=0; s<26; s++)
   {
//    initialize counter as zero
      NSend[s] = 0;
      
//    allocate the maximum necessary memory
      Send_PosList[s] = new int [ ParaVar.BounP_NList[lv-1][s] ];

      for (int ID=0; ID<ParaVar.BounP_NList[lv-1][s]; ID++) 
      {
         RefinePos = ParaVar.BounP_PosList[lv-1][s][ID];
         BounPID   = ParaVar.BounP_IDList [lv-1][s][ID];
         SonPID    = patch.ptr[lv-1][BounPID]->son;
      
         if ( SonPID != -1 )
         {
            Send_PosList[s][ NSend[s] ] = RefinePos;
            NSend[s] ++;
         }
      }

   } // for (int s=0; s<26; s+=2)


// get the position of buffer patches to be refined at level "lv"
   MPI_ExchangeBufferPosition( NSend, NRecv, Send_PosList, Recv_PosList );


// allocate the buffer patches at level "lv"
   const int scale = patch.scale[lv];
   const int Disp  = PATCH_SIZE*scale;
   int *Corner;

   for (int s=0; s<26; s++)
   {

//    determine the buffer patches that must store physical data
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         if (  TABLE_01( s, 'x', -1, 0, 1 ) == TABLE_02( LocalID, 'x', -1, 1 )  ||
               TABLE_01( s, 'y', -1, 0, 1 ) == TABLE_02( LocalID, 'y', -1, 1 )  ||
               TABLE_01( s, 'z', -1, 0, 1 ) == TABLE_02( LocalID, 'z', -1, 1 )      )  

            AllocateData[LocalID] = false;

         else

            AllocateData[LocalID] = true;
      }


//    find the targeted buffer patches to be refined
      NBuffer[s] = 0;
      TargetID   = 0;
      
      for (int ID=0; ID<NRecv[s]; ID++)
      {
         RefinePos = Recv_PosList[s][ID];

         while ( ParaVar.BounP_PosList[lv-1][s][TargetID] != RefinePos )  
         {
            TargetID++;

//          the BounPID must exist due to the proper-nesting condition
            if ( TargetID >= ParaVar.BounP_NList[lv-1][s] )
            {
               fprintf( stderr, "ERROR : \"incorrect parameter %s = %d\" !!\n", "TargetID", TargetID );
               fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", 
                        __FILE__, __LINE__,  __FUNCTION__  );
               MPI_Exit();
            }
         }

         BounPID = ParaVar.BounP_IDList[lv-1][s][TargetID];
         BuffPID = patch.ptr[lv-1][BounPID]->sibling[s];
      

//       the BuffPID must exist since that the flagging status should be the same over all processes
         if ( BuffPID == -1 )
         {
            fprintf( stderr, "ERROR : \"incorrect parameter %s = %d\" !!\n", "BuffPID", BuffPID );
            fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", 
                     __FILE__, __LINE__,  __FUNCTION__  );
            MPI_Exit();
         }


//       create buffer patches
         Corner = patch.ptr[lv-1][BuffPID]->corner;
      
         patch.ptr[lv-1][BuffPID]->son = patch.num[lv];



         for (int LocalID=0; LocalID<8; LocalID++)
         {
            const int PS = PATCH_SIZE*patch.scale[lv];
            int Cr[3];

            for (int d=0; d<3; d++)    Cr[d] = Corner[d] + TABLE_02( LocalID, 'x'+d, 0, Disp );

            AllocateData2[LocalID] = ( AllocateData[LocalID] && WithinCandidateBox(Cr,PS,CanBuf) );
         }


      
         patch.pnew( lv, Corner[0],      Corner[1],      Corner[2],      BuffPID, AllocateData2[0] );
         patch.pnew( lv, Corner[0]+Disp, Corner[1],      Corner[2],      BuffPID, AllocateData2[1] );
         patch.pnew( lv, Corner[0],      Corner[1]+Disp, Corner[2],      BuffPID, AllocateData2[2] );
         patch.pnew( lv, Corner[0],      Corner[1],      Corner[2]+Disp, BuffPID, AllocateData2[3] );
         patch.pnew( lv, Corner[0]+Disp, Corner[1]+Disp, Corner[2],      BuffPID, AllocateData2[4] );
         patch.pnew( lv, Corner[0],      Corner[1]+Disp, Corner[2]+Disp, BuffPID, AllocateData2[5] );
         patch.pnew( lv, Corner[0]+Disp, Corner[1],      Corner[2]+Disp, BuffPID, AllocateData2[6] );
         patch.pnew( lv, Corner[0]+Disp, Corner[1]+Disp, Corner[2]+Disp, BuffPID, AllocateData2[7] );
      
         NBuffer[s] += 8;

       } // for (int ID=0; ID<NRecv[s]; ID++)
      
       delete [] Send_PosList[s];
       delete [] Recv_PosList[s];

   } // for (int s=0; s<26; s++)


// set up the NPatchComma list
   for (int s=0; s<26; s++)   NPatchComma[lv][s+2] = NPatchComma[lv][s+1] + NBuffer[s];

   if ( NPatchComma[lv][27] != patch.num[lv] )
   {
      fprintf( stderr, "ERROR : \"NPatchComma[%d][27] (%d) != patch.num[%d] (%d)\" !!\n",
               lv, NPatchComma[lv][27], lv, patch.num[lv] ); 
      fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
      MPI_Exit();
   }

}
