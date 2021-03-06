
#include "DAINO.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Buffer
// Description :  Flag the buffer patches (also the boundary patches)
//
// Note        :  a. All flags should be initialized as "false" by calling the function "Flag_Real" in advance
//                b. The arrays "ParaVar->BuffFlag_NList and ParaVar->BuffFlag_PosList" must be prepared 
//                   in advance by calling the functions "Buf_RecordBoundaryFlag and MPI_ExchangeBoundaryFlag" 
//                c. No OpenMP directives are applied in this function
//
// Parameter   :  lv : Targeted refinement level to be flagged
//-------------------------------------------------------------------------------------------------------
void Flag_Buffer( const int lv )
{

// MirrorSib : the mirror-symmetric sibling index
   const int MirrorSib[26] = { 1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18 };
   int BounPID, BuffPID, FlagPos, SibPID, TargetID, FlagLayer, Sib;

   for (int s=0; s<26; s++)
   {
      FlagLayer = TABLE_05( s );
      TargetID  = 0;

      for (int ID=0; ID<patch->ParaVar->BuffFlag_NList[lv][s]; ID+=FlagLayer)
      {
         FlagPos = patch->ParaVar->BuffFlag_PosList[lv][s][ID];

         while ( patch->ParaVar->BounP_PosList[lv][s][TargetID] != FlagPos )  
         {
            TargetID++;

//          the BounPID must exist due to the proper-nesting condition
#           ifdef DAINO_DEBUG
            if ( TargetID >= patch->ParaVar->BounP_NList[lv][s] )
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TargetID", TargetID );
#           endif
         }

         BounPID = patch->ParaVar->BounP_IDList[lv][s][TargetID];
         BuffPID = patch->ptr[0][lv][BounPID]->sibling[s];


//       the BuffPID must exist since that the flagging status should be the same over all processes
#        ifdef DAINO_DEBUG
         if ( BuffPID == -1 )    Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "BuffPID", BuffPID );
#        endif


//       flag the patch with layer = 0
         patch->ptr[0][lv][BuffPID]->flag = true;


//       flag the patch with layer > 0
         for (int n=1; n<FlagLayer; n++)
         {
            if ( patch->ParaVar->BuffFlag_PosList[lv][s][ID+n] == -999 )
            {
               Sib    = TABLE_06( MirrorSib[s], n );
               SibPID = patch->ptr[0][lv][BuffPID]->sibling[Sib];

               patch->ptr[0][lv][SibPID]->flag = true;
            }
         }

      } // for (int ID=0; ID<BuffFlag_NList[lv][s]; ID+=FlagLayer)


//    deallocate the memory previously allocated by "MPI_ExchangeBoundaryFlag or OOC_ExchangeBoundaryFlag"
      if ( patch->ParaVar->BuffFlag_PosList[lv][s] != NULL )
      {
         delete [] patch->ParaVar->BuffFlag_PosList[lv][s];
         patch->ParaVar->BuffFlag_PosList[lv][s] = NULL;
      }

      patch->ParaVar->BuffFlag_NList[lv][s] = 0;

   } // for (int s=0; s<26; s++)

} // FUNCTION : Flag_Buffer



#endif // #ifndef SERIAL
