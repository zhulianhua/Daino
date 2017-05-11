
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_PatchCorner
// Description :  Output the corner coordinates of all real and buffer patches, which can be visualized
//                by "splot" in gnuplot
//
// Note        :  An empty line is inserted between the data of real and buffer patches
//                --> For gnuplot plot, one can use (1) "every :::0::0" to plot the real   patches
//                                                  (2) "every :::1::1" to plot the buffer patches
//
// Parameter   :  lv       : Targeted refinement level 
//                comment  : String to attach to the end of the file name
//-------------------------------------------------------------------------------------------------------
void Output_PatchCorner( const int lv, const char *comment )
{

   char FileName[100];
   sprintf( FileName, "PatchCorner_%d%d%d%d_%d%d", MPI_Rank/1000, MPI_Rank%1000/100, MPI_Rank%100/10, MPI_Rank%10,
                                                   lv/10, lv%10 );
   if ( comment != NULL )       
   {
      strcat( FileName, "_" );
      strcat( FileName, comment );
   }

   FILE *File_Check = fopen( FileName, "r" );
   if ( File_Check != NULL )
   {
      Aux_Message( stderr, "WARNING : the file \"%s\" already exists and will be overwritten !!\n", FileName );
      fclose( File_Check );
   }


   const int NReal = patch->NPatchComma[lv][1];
   const int NBuff = patch->NPatchComma[lv][27] - patch->NPatchComma[lv][1];

   FILE *File = fopen( FileName, "w" );

   fprintf( File, "Time %13.7e  Step %ld  Counter %u  Rank %d  Level %d  NRealPatch %d  NBufferPatch %d\n", 
            Time[lv], Step, AdvanceCounter[lv], MPI_Rank, lv, NReal, NBuff );
   fprintf( File, "=========================================================================================\n" );
   fprintf( File, "%8s   %11s   %11s   %11s\n", "PID", "Corner[x]", "Corner[y]", "Corner[z]" );

// output real patches
   for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
      fprintf( File, "%8d   %11d   %11d   %11d\n", PID, patch->ptr[0][lv][PID]->corner[0], 
                                                        patch->ptr[0][lv][PID]->corner[1], 
                                                        patch->ptr[0][lv][PID]->corner[2] );
   fprintf( File, "\n" );

// output buffer patches
   for (int PID=patch->NPatchComma[lv][1]; PID<patch->NPatchComma[lv][27]; PID++)
      fprintf( File, "%8d   %11d   %11d   %11d\n", PID, patch->ptr[0][lv][PID]->corner[0], 
                                                        patch->ptr[0][lv][PID]->corner[1], 
                                                        patch->ptr[0][lv][PID]->corner[2] );

   fclose( File );

} // FUNCTION : Output_PatchCorner


