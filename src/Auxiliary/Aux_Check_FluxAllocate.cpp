
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_FluxAllocate
// Description :  Verify that the flux arrays are properly allocated at all coarse-fine boundaries
//
// Note        :  This function will not check anything if "patch->WithFlux == false", in which no flux
//                array is allocated
//
// Parameter   :  lv       : Targeted refinement level
//                comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Aux_Check_FluxAllocate( const int lv, const char *comment )
{

// check
   if ( lv == NLEVEL-1 )
   {
      Aux_Message( stderr, "WARNING : function \"%s\" should NOT be applied to the highest level !!\n",
                   __FUNCTION__ );
      return;
   }

   if ( !patch->WithFlux )
   {
      Aux_Message( stderr, "WARNING : function \"%s\" is useless since no flux is required !!\n", __FUNCTION__ );
      OPT__CK_FLUX_ALLOCATE = false;
      return;
   }


   int Pass = true;
   int SonPID, SibPID, SibSonPID;
   real (*FluxPtr)[PATCH_SIZE][PATCH_SIZE] = NULL;

   for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   {
      if ( MPI_Rank == TargetRank )
      {
//       check all "real" patches
         for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
         {
            SonPID = patch->ptr[0][lv][PID]->son;

            for (int sib=0; sib<6; sib++)    
            {
               FluxPtr = patch->ptr[0][lv][PID]->flux[sib];

               if ( SonPID == -1 )
               {
                  SibPID = patch->ptr[0][lv][PID]->sibling[sib];

                  if ( SibPID != -1 )
                  {
                     SibSonPID = patch->ptr[0][lv][SibPID]->son;

                     if ( SibSonPID == -1  &&  FluxPtr != NULL )
                     {
                        if ( Pass )
                        {
                           Aux_Message( stderr, "\"%s\" : <%s> FAILED at lv %2d, Time = %13.7e, Step = %ld !!\n",
                                        comment, __FUNCTION__, lv, Time[lv], Step );
                           Aux_Message( stderr, "%4s\t%7s\t%7s\t\t%19s\t%7s\t\t%7s\n", 
                                        "Rank", "PatchID", "SibID", "Patch Corner", "Surplus", "Missing" );

                           Pass = false;
                        }

                        Aux_Message( stderr, "%4d\t%7d\t%7d\t\t(%5d,%5d,%5d)\t%7s\t\t%7s\n",
                                     MPI_Rank, PID, sib, patch->ptr[0][lv][PID]->corner[0],
                                                         patch->ptr[0][lv][PID]->corner[1],
                                                         patch->ptr[0][lv][PID]->corner[2], "O", "X" );
                     }

                     if ( SibSonPID != -1  &&  FluxPtr == NULL )
                     {
                        if ( Pass )
                        {
                           Aux_Message( stderr, "\"%s\" : <%s> FAILED at lv %2d, Time = %13.7e, Step = %ld !!\n",
                                        comment, __FUNCTION__, lv, Time[lv], Step );
                           Aux_Message( stderr, "%4s\t%7s\t%7s\t\t%19s\t%7s\t\t%7s\n", 
                                        "Rank", "PatchID", "SibID", "Patch Corner", "Surplus", "Missing" );

                           Pass = false;
                        }

                        Aux_Message( stderr, "%4d\t%7d\t%7d\t\t(%5d,%5d,%5d)\t%7s\t\t%7s\n",
                                     MPI_Rank, PID, sib, patch->ptr[0][lv][PID]->corner[0],
                                                         patch->ptr[0][lv][PID]->corner[1],
                                                         patch->ptr[0][lv][PID]->corner[2], "X", "O" );
                     }
                  } // if ( SibPID != -1 )
               } // if ( SonPID != -1 )

               else // if ( SonPID != -1 ) ... else ...
               {
                  if ( FluxPtr != NULL )
                  {
                     if ( Pass )
                     {
                        Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                     comment, __FUNCTION__, lv, Time[lv], Step );
                        Aux_Message( stderr, "%4s\t%7s\t%7s\t\t%19s\t%7s\t\t%7s\n", 
                                     "Rank", "PatchID", "SibID", "Patch Corner", "Surplus", "Missing" );

                        Pass = false;
                     }

                     Aux_Message( stderr, "%4d\t%7d\t%7d\t\t(%5d,%5d,%5d)\t%7s\t\t%7s\n",
                                  MPI_Rank, PID, sib, patch->ptr[0][lv][PID]->corner[0],
                                                      patch->ptr[0][lv][PID]->corner[1],
                                                      patch->ptr[0][lv][PID]->corner[2], "O", "X" );
                  }
               } // if ( SonPID != -1 ) ... else ... 
            } // for (int sib=0; sib<6; sib++)
         } // for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)


//       check "buffer" patches (sib = 0 ~ 5) (not supported in LOAD_BALANCE)
#        ifndef LOAD_BALANCE
         const int MirrorSib[6] = { 1, 0, 3, 2, 5, 4 };
         int PID, Table[4];

         for (int sib=0; sib<6; sib++)
         {
            for (int t=0; t<4; t++)    Table[t] = TABLE_03(sib, t);

            for (int PID0=patch->NPatchComma[lv][1+sib]; PID0<patch->NPatchComma[lv][2+sib]; PID0+=8)
            for (int t=0; t<4; t++)
            {
               PID    = PID0 + Table[t];
               SonPID = patch->ptr[0][lv][PID]->son;

               FluxPtr = patch->ptr[0][lv][PID]->flux[ MirrorSib[sib] ];

               if ( SonPID == -1 )
               {
                  SibPID = patch->ptr[0][lv][PID]->sibling[ MirrorSib[sib] ];

                  if ( SibPID != -1 )
                  {
                     SibSonPID = patch->ptr[0][lv][SibPID]->son;

                     if ( SibSonPID == -1  &&  FluxPtr != NULL )
                     {
                        if ( Pass )
                        {
                           Aux_Message( stderr, "\"%s\" : <%s> FAILED at lv %2d, Time = %13.7e, Step = %ld !!\n",
                                        comment, __FUNCTION__, lv, Time[lv], Step );
                           Aux_Message( stderr, "%4s\t%7s\t%7s\t\t%19s\t%7s\t\t%7s\n", 
                                        "Rank", "PatchID", "SibID", "Patch Corner", "Surplus", "Missing" );

                           Pass = false;
                        }

                        Aux_Message( stderr, "%4d\t%7d\t%7d\t\t(%5d,%5d,%5d)\t%7s\t\t%7s\n",
                                     MPI_Rank, PID, sib, patch->ptr[0][lv][PID]->corner[0],
                                                         patch->ptr[0][lv][PID]->corner[1],
                                                         patch->ptr[0][lv][PID]->corner[2], "O", "X" );
                     }

                     if ( SibSonPID != -1  &&  FluxPtr == NULL )
                     {
                        if ( Pass )
                        {
                           Aux_Message( stderr, "\"%s\" : <%s> FAILED at lv %2d, Time = %13.7e, Step = %ld !!\n",
                                        comment, __FUNCTION__, lv, Time[lv], Step );
                           Aux_Message( stderr, "%4s\t%7s\t%7s\t\t%19s\t%7s\t\t%7s\n", 
                                        "Rank", "PatchID", "SibID", "Patch Corner", "Surplus", "Missing" );

                           Pass = false;
                        }

                        Aux_Message( stderr, "%4d\t%7d\t%7d\t\t(%5d,%5d,%5d)\t%7s\t\t%7s\n",
                                     MPI_Rank, PID, sib, patch->ptr[0][lv][PID]->corner[0],
                                                         patch->ptr[0][lv][PID]->corner[1],
                                                         patch->ptr[0][lv][PID]->corner[2], "X", "O" );
                     }

                  } // if ( SibPID != -1 )
               } // if ( SonPID != -1 )

               else
               {
                  if ( FluxPtr != NULL )
                  {
                     if ( Pass )
                     {
                        Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                     comment, __FUNCTION__, lv, Time[lv], Step );
                        Aux_Message( stderr, "%4s\t%7s\t%7s\t\t%19s\t%7s\t\t%7s\n", 
                                     "Rank", "PatchID", "SibID", "Patch Corner", "Surplus", "Missing" );

                        Pass = false;
                     }

                     Aux_Message( stderr, "%4d\t%7d\t%7d\t\t(%5d,%5d,%5d)\t%7s\t\t%7s\n",
                                  MPI_Rank, PID, sib, patch->ptr[0][lv][PID]->corner[0],
                                                      patch->ptr[0][lv][PID]->corner[1],
                                                      patch->ptr[0][lv][PID]->corner[2], "O", "X" );
                  }
               } // if ( SonPID != -1 ) ... else

            } // for (int PID0=patch->NPatchComma[lv][1+sib]; PID0<patch->NPatchComma[lv][2+sib]; PID0+=8) 
              // for (int t=0; t<4; t++)
         } // for (int sib=0; sib<6; sib++)
#        endif // #ifndef LOAD_BALANCE

      } // if ( MPI_Rank == TargetRank )

      MPI_Bcast( &Pass, 1, MPI_INT, TargetRank, MPI_COMM_WORLD );

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   

   if ( Pass )
   {
      if ( MPI_Rank == 0 )   
         Aux_Message( stdout, "\"%s\" : <%s> PASSED at level %2d, Time = %13.7e, Step = %ld\n", 
                      comment, __FUNCTION__, lv, Time[lv], Step );
   }

   else
      MPI_Exit();

} // FUNCTION : Aux_Check_FluxAllocate
