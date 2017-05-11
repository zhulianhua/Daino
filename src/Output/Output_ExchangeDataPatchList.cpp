
#include "DAINO.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_ExchangeDataPatchList
// Description :  Output SendP_IDList[lv] or RecvP_IDList[lv]
//
// Option      :  0 -> RecvP_IDList[lv]
//             :  1 -> SendP_IDList[lv]
//-------------------------------------------------------------------------------------------------------
void Output_ExchangeDataPatchList( const int option, const int lv, const char *comment )
{

// check
   if ( option != 0  &&  option != 1 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "option", option );


   char FileName[100];
   if ( option )  sprintf( FileName, "SendDataPatchList_%d_%d", DAINO_RANK, lv );
   else           sprintf( FileName, "RecvDataPatchList_%d_%d", DAINO_RANK, lv );

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

   fprintf( File, "Time = %13.7e  Step = %ld  Rank = %d  Level = %d\n\n", Time[0], Step, DAINO_RANK, lv );


   for (int s=0; s<26; s++)
   {
      const int NP    = ( option )  ?  patch->ParaVar->SendP_NList [lv][s] : patch->ParaVar->RecvP_NList [lv][s];
      const int *List = ( option )  ?  patch->ParaVar->SendP_IDList[lv][s] : patch->ParaVar->RecvP_IDList[lv][s];

      fprintf( File, "Face = %d     Length = %d\n", s, NP );

      for (int P=0; P<NP; P++)      fprintf( File, "%5d ", List[P] );

      fprintf( File, "\n\n" );
   }

   fclose( File );

} // FUNCTION : Output_ExchangeDataPatchList



#endif // #ifndef SERIAL
