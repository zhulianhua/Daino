
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_ProperNesting
// Description :  Verify the proper-nesting condition at the input level
//                --> Ensure that all patches at level "lv" are properly enclosed by patches at level "lv-1"
//
// Parameter   :  lv       : Targeted refinement level
//                comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Aux_Check_ProperNesting( const int lv, const char *comment )
{

// check
   if ( lv == 0 )
   {
      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : the function \"%s\" should NOT be applied to the base level !!\n",
                      __FUNCTION__ );
      return;
   }


   int Pass = true;
   int FaPID, FaSibPID;

   for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   {
      if ( MPI_Rank == TargetRank )
      {
         for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
         {
            FaPID = patch->ptr[0][lv][PID]->father;

            if ( FaPID == -1 )
            {
               if ( Pass )
               {
                  Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                               comment, __FUNCTION__, lv, Time[lv], Step );
                  Aux_Message( stderr, "%4s\t%7s\t%7s\t\t%19s\n", 
                               "Rank", "PID", "FaPID", "Coordinate" );

                  Pass = false;
               }

               Aux_Message( stderr, "%4d\t%7d\t%7d\t\t(%5d,%5d,%5d)\n",
                            MPI_Rank, PID, FaPID, patch->ptr[0][lv][PID]->corner[0],  
                                                  patch->ptr[0][lv][PID]->corner[1],  
                                                  patch->ptr[0][lv][PID]->corner[2]  ); 
            } // if ( FaPID == -1 )

            for (int sib=0; sib<26; sib++)
            {
               FaSibPID = patch->ptr[0][lv-1][FaPID]->sibling[sib];

               if ( FaSibPID == -1 )
               {
                  if ( Pass )
                  {
                     Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                  comment, __FUNCTION__, lv, Time[lv], Step );
                     Aux_Message( stderr, "%4s\t%7s\t%7s\t%7s\t\t%19s\n", 
                                  "Rank", "PID", "FaPID", "Sib", "Coordinate" );

                     Pass = false;
                  }

                  Aux_Message( stderr, "%4d\t%7d\t%7d\t%7d\t\t(%5d,%5d,%5d)\n",
                               MPI_Rank, PID, FaPID, sib, patch->ptr[0][lv][PID]->corner[0],  
                                                          patch->ptr[0][lv][PID]->corner[1],  
                                                          patch->ptr[0][lv][PID]->corner[2]  ); 
               } // if ( FaSibPID == -1 )

            } // for (int sib=0; sib<26; sib++)
         } // for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
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

} // FUNCTION : Aux_Check_ProperNesting
