
#include "DAINO.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_BoundaryFlagList
// Description :  Output BounP_PosList[lv] or BuffP_PosList[lv]
//
// option      :  0 -> BuffP_PosList[lv]
//             :  1 -> BounP_PosList[lv]
//-------------------------------------------------------------------------------------------------------
void Output_BoundaryFlagList( const int option, const int lv, const char *comment )
{

// check
   if ( option != 0  &&  option != 1 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "option", option );


   char FileName[100];
   if ( option )  sprintf( FileName, "BoundaryFlagList_%d_%d", MPI_Rank, lv );
   else           sprintf( FileName, "BufferFlagList_%d_%d", MPI_Rank, lv );

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


   FILE *File = fopen( FileName, "w" );

   fprintf( File, "Time = %13.7e  Step = %ld  Rank = %d  Level = %d\n\n", Time[0], Step, MPI_Rank, lv );


   for (int s=0; s<26; s++)
   {
      const int NP    = ( option ) ? patch->ParaVar->BounFlag_NList  [lv][s] : 
                                     patch->ParaVar->BuffFlag_NList  [lv][s];
      const int *List = ( option ) ? patch->ParaVar->BounFlag_PosList[lv][s] : 
                                     patch->ParaVar->BuffFlag_PosList[lv][s];
      const int FlagLayer = TABLE_05( s );

      fprintf( File, "Face = %d     Length = %d\n", s, NP );

      for (int P=0; P<NP; P++)  
      {
         fprintf( File, "%5d ", List[P] );

         if ( (P+1)%FlagLayer == 0 )   fprintf( File, "  ||  " );
      }

      fprintf( File, "\n\n" );
   }

   fclose( File );

} // FUNCTION : Output_BoundaryFlagList



#endif // #ifndef SERIAL
