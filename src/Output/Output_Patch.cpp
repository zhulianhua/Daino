
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_Patch
// Description :  Output the data of a single patch
//
// Example     :  char comment[10];
//                sprintf( comment, "Step%d", AdvanceCounter[6] );
//                Output_Patch( 6, 5560, comment );
//
// Parameter   :  lv       : Targeted refinement level 
//                PID      : Targeted patch index
//                FluSg    : Sandglass of the fluid data
//                PotSg    : Sandglass of the potential data
//                comment  : String to attach to the end of the file name
//-------------------------------------------------------------------------------------------------------
void Output_Patch( const int lv, const int PID, const int FluSg, const int PotSg, const char *comment )
{

// check
   if ( lv < 0  ||  lv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

   if ( PID < 0  ||  PID >= MAX_PATCH )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d (MAX_PATCH = %d) !!\n", "PID", PID, MAX_PATCH );

   if ( FluSg < 0  ||  FluSg >= 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FluSg", FluSg );

#  ifdef GRAVITY
   if ( PotSg < 0  ||  PotSg >= 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "PotSg", PotSg );
#  endif

   if ( patch->ptr[0][lv][PID] == NULL )
   {
      Aux_Message( stderr, "WARNING : level = %d, PID = %d does NOT exist !!\n", lv, PID );
      return;
   }


   patch_t *Relation = patch->ptr[    0][lv][PID];
   patch_t *FluData  = patch->ptr[FluSg][lv][PID];
#  ifdef GRAVITY
   patch_t *PotData  = patch->ptr[PotSg][lv][PID];
#  endif

   char FileName[100];
   sprintf( FileName, "Patch_r%d_lv%d_p%d", DAINO_RANK, lv, PID );
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


// output patch information
   FILE *File = fopen( FileName, "w" );

   fprintf( File, "Rank %d  Lv %d  PID %d  Local ID %d  FluSg %d  PotSg %d  Time %13.7e  Step %ld  Counter %u\n", 
            DAINO_RANK, lv, PID, PID%8, FluSg, PotSg, Time[lv], Step, AdvanceCounter[lv] );

   fprintf( File, "Father %d  Son %d  Corner (%4d,%4d,%4d)", Relation->father, Relation->son, 
            Relation->corner[0], Relation->corner[1], Relation->corner[2] );
#  ifdef LOAD_BALANCE
   fprintf( File, "  LB_Idx %ld  PaddedCr1D %ld", Relation->LB_Idx, Relation->PaddedCr1D );
#  endif
   fprintf( File, "\n\n" );

   fprintf( File, "Sibling, Sibling->Son, and Father->Sibling Lists :\n" );

   int Sib, FaSib, SibSon, Fa;
   for (int S=0; S<26; S++)   
   {
      Fa     = Relation->father;
      Sib    = Relation->sibling[S];
      FaSib  = (  Fa == -1 ) ? -1 : ( patch->ptr[0][lv-1][Fa] != NULL ) ? 
                                      patch->ptr[0][lv-1][Fa]->sibling[S] : -999;
      SibSon = ( Sib == -1 ) ? -1 : patch->ptr[0][lv][Sib]->son; 

      fprintf( File, "Sib[%2d] = %6d     Sib_Son = %6d     Fa_Sib[%2d] = %6d\n", 
               S, Sib, SibSon, S, FaSib );
   }
   fprintf( File, "\n" );



// check whether or not the targeted patch stores physical data
   if ( FluData->fluid == NULL )
      Aux_Message( stderr, "WARNING : Lv = %d, PID = %d does NOT store fluid data !!\n", lv, PID );
#  ifdef GRAVITY
   if ( PotData->pot == NULL )
      Aux_Message( stderr, "WARNING : Lv = %d, PID = %d does NOT store potential data !!\n", lv, PID );
#  endif


// output header
#  if   ( MODEL == HYDRO )
#  ifdef GRAVITY
   fprintf( File, "(%2s,%2s,%2s)%14s%14s%14s%14s%14s%14s%14s\n", "i", "j", "k", "Density", "Px", "Py", "Pz",
                                                                 "Energy", "Pressure", "Potential" );
#  else
   fprintf( File, "(%2s,%2s,%2s)%14s%14s%14s%14s%14s%14s\n",     "i", "j", "k", "Density", "Px", "Py", "Pz",
                                                                 "Energy", "Pressure" );
#  endif

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
#  ifdef GRAVITY
   fprintf( File, "(%2s,%2s,%2s)%14s%14s%14s%14s\n", "i", "j", "k", "Density", "Real", "Imag", "Potential" );
#  else
   fprintf( File, "(%2s,%2s,%2s)%14s%14s%14s\n",     "i", "j", "k", "Density", "Real", "Imag" );
#  endif

#  else
#  warning : WARNING : DO YOU WANT TO ADD the FILE HEADER HERE FOR THE NEW MODEL ??
#  endif // MODEL


// output data
   real u[NCOMP]; 
#  ifdef GRAVITY
   real pot=NULL_REAL;
#  endif

   for (int k=0; k<PATCH_SIZE; k++)
   for (int j=0; j<PATCH_SIZE; j++)
   for (int i=0; i<PATCH_SIZE; i++)
   {
      if ( FluData->fluid != NULL )    for (int v=0; v<NCOMP; v++)   u[v] = FluData->fluid[v][k][j][i];

#     ifdef GRAVITY
      if ( PotData->pot   != NULL )    pot = PotData->pot[k][j][i];
#     endif

//    output cell indices      
      fprintf( File, "(%2d,%2d,%2d)", i, j, k );

      if ( FluData->fluid != NULL )
      {
//       output all variables in the fluid array
         for (int v=0; v<NCOMP; v++)   fprintf( File, " %13.6e", u[v] );

//       output pressure in HYDRO
#        if   ( MODEL == HYDRO )
         fprintf( File, " %13.6e", ( u[ENGY]-0.5*(u[MOMX]*u[MOMX]+u[MOMY]*u[MOMY]+u[MOMZ]*u[MOMZ])/u[DENS] )*
                                   (GAMMA-1.0) );
#        elif ( MODEL == MHD )
#        warning : WAIT MHD !!!
#        endif // MODEL
      }

      else
      {
//       output empty strings if the fluid array is not allocated
         for (int v=0; v<NCOMP; v++)   fprintf( File, " %13s", "" );

#        if   ( MODEL == HYDRO )
         fprintf( File, " %13s", "" );
#        elif ( MODEL == MHD )
#        warning : WAIT MHD !!!
#        endif // MODEL
      }

//    output potential
#     ifdef GRAVITY
      if ( PotData->pot != NULL )
      fprintf( File, " %13.6e", pot );
#     endif

      fprintf( File, "\n" );
   } // i,j,k

   fclose( File );

} // FUNCTION : Output_Patch
