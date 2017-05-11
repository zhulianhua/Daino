
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_PatchCount
// Description :  Count the total number of patches at each level
//
// Node        :  OPT__PATCH_COUNT = (1,3/2,4) --> count the number of patches per OOC/MPI rank
//-------------------------------------------------------------------------------------------------------
void Aux_PatchCount()
{

// check 
   if ( OPT__PATCH_COUNT == 0  ||  OPT__PATCH_COUNT > 4 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "OPT__PATCH_COUNT", OPT__PATCH_COUNT );


   const char FileName[] = "Record__PatchCount";
   static bool FirstTime = true;

   if ( MPI_Rank == 0  &&  FirstTime )
   {
      FILE *File_Check = fopen( FileName, "r" );
      if ( File_Check != NULL )  
      {
         Aux_Message( stderr, "WARNING : the file \"%s\" already exists !!\n", FileName );
         fclose( File_Check );
      }
      FirstTime = false;
   }


#  ifdef OOC
   const int OOC_NRank = ooc.NRank;
#  else
   const int OOC_NRank = 1;
#  endif
   long   MaxNPatch;
   int    NRank;
   int    (*NPatch_Local)[NLEVEL], (*NPatch_Gather)[NLEVEL];
   double (*Coverage_Gather)[NLEVEL], Coverage_Total[NLEVEL];


// a. gather information from all DAINO ranks
   if ( OPT__PATCH_COUNT == 2  ||  OPT__PATCH_COUNT == 4 )     
   {
      NRank        = MPI_NRank;
      NPatch_Local = new int [1][NLEVEL];

      for (int lv=0; lv<NLEVEL; lv++)
      {
         NPatch_Local[0][lv] = 0;

         for (int r=0; r<OOC_NRank; r++)     
         {
#           ifdef OOC
            NPatch_Local[0][lv] += ooc.NRealPatch[r][lv];
#           else
            NPatch_Local[0][lv] = patch->NPatchComma[lv][1];
#           endif
         }
      }
   }

   else  // OPT__PATCH_COUNT == 1 or 3                                            
   {
      NRank        = DAINO_NRANK;
      NPatch_Local = new int [OOC_NRank][NLEVEL];

      for (int lv=0; lv<NLEVEL; lv++)
      {
         for (int r=0; r<OOC_NRank; r++)     
         {
#           ifdef OOC
            NPatch_Local[r][lv] = ooc.NRealPatch[r][lv];
#           else
            NPatch_Local[0][lv] = patch->NPatchComma[lv][1];
#           endif
         }
      }
   }

   NPatch_Gather   = new int    [NRank][NLEVEL];
   Coverage_Gather = new double [NRank][NLEVEL];

#  ifdef SERIAL
   for (int r=0; r<NRank/MPI_NRank; r++)
   for (int lv=0; lv<NLEVEL; lv++)     NPatch_Gather[r][lv] = NPatch_Local[r][lv];
#  else
   MPI_Gather( NPatch_Local, NRank/MPI_NRank*NLEVEL, MPI_INT, NPatch_Gather, NRank/MPI_NRank*NLEVEL, 
               MPI_INT, 0, MPI_COMM_WORLD ); 
#  endif


   if ( MPI_Rank == 0 )
   {
//    b. evaluate the coverage ratio of each OOC/MPI rank
      for (int lv=0; lv<NLEVEL; lv++)
      {
         MaxNPatch = (NX0_TOT[0]/PATCH_SIZE)*(NX0_TOT[1]/PATCH_SIZE)*(NX0_TOT[2]/PATCH_SIZE)*(long)(1L<<3*lv);

         if ( OPT__PATCH_COUNT == 2  ||  OPT__PATCH_COUNT == 4 )  MaxNPatch *= OOC_NRank;

         for (int r=0; r<NRank; r++)
            Coverage_Gather[r][lv] = 100.0*NPatch_Gather[r][lv]/MaxNPatch*MPI_NRank;
      }


//    c. evaluate the coverage ratio of the entire simulation domain
      for (int lv=0; lv<NLEVEL; lv++)
      {
         MaxNPatch = (NX0_TOT[0]/PATCH_SIZE)*(NX0_TOT[1]/PATCH_SIZE)*(NX0_TOT[2]/PATCH_SIZE)*(long)(1L<<3*lv);

         Coverage_Total[lv] = 100.0*NPatchTotal[lv]/MaxNPatch;
      }


//    d. write to the file "Record__PatchCount"
      FILE *File = fopen( FileName, "a" );

      fprintf( File, "Time = %13.7e,  Step = %7ld\n\n", Time[0], Step );

      fprintf( File, "%4s", "Rank" );
      for (int lv=0; lv<NLEVEL; lv++)     fprintf( File, "%12s%2d ", "Level", lv );
      fprintf( File, "\n" );

      for (int r=0; r<NRank; r++)
      {
         fprintf( File, "%4d", r );
         for (int lv=0; lv<NLEVEL; lv++)  fprintf( File, "%6d(%6.2lf%%)", NPatch_Gather[r][lv], 
                                                                          Coverage_Gather[r][lv] );
         fprintf( File, "\n" );
      }

      fprintf( File, "-------------------------------------------------------------------------------------" );
      fprintf( File, "----------------------------\n" );
      fprintf( File, "%4s", "Sum:" );
      for (int lv=0; lv<NLEVEL; lv++)     fprintf( File, "%6d(%6.2lf%%)", NPatchTotal[lv], Coverage_Total[lv] );
      fprintf( File, "\n" );


//    e. get the load and load-imbalance factor at each level
      int  Load_Max[NLEVEL];
      real Load_Ave[NLEVEL], Load_Imb[NLEVEL];  // Imb : Imbalance
      real WLoad_Ave=0.0, WLoad_Max=0.0, WLI;   // WLoad : weighted load

#     ifdef LOAD_BALANCE
      const real *Weighting = patch->LB->WLI_Weighting;
#     else
      real Weighting[NLEVEL];
      for (int lv=0; lv<NLEVEL; lv++)  Weighting[lv] = real( 1UL << lv );  // 2^lv
#     endif

//    e1. estimate the load-imbalance factor at each level
      for (int lv=0; lv<NLEVEL; lv++)
      {
         Load_Max[lv] = -1;
         Load_Ave[lv] = (real)NPatchTotal[lv] / MPI_NRank;

         for (int r=0; r<MPI_NRank; r++)
            if ( NPatch_Gather[r][lv] > Load_Max[lv] )   Load_Max[lv] = NPatch_Gather[r][lv];

         Load_Imb[lv] = (Load_Max[lv]==0) ? 0.0 : ( Load_Max[lv] - Load_Ave[lv] ) / Load_Ave[lv];
      }

//    e2. estimate the weighted load-imbalance factor of patches at all levels
      for (int lv=0; lv<NLEVEL; lv++)
      {
         WLoad_Max += Weighting[lv]*Load_Max[lv];
         WLoad_Ave += Weighting[lv]*Load_Ave[lv];
      }

      WLI = ( WLoad_Max - WLoad_Ave ) / WLoad_Ave;

//    e3. record the load-imbalance factors
      fprintf( File, "%4s", "LIM:" );
      for (int lv=0; lv<NLEVEL; lv++)  fprintf( File, "%6d(%6.2f%%)", Load_Max[lv], 100.0*Load_Imb[lv] );
      fprintf( File, "\n" );
      fprintf( File, "Weighted load-imbalance factor = %6.2f%%\n", 100.0*WLI );

//    e4. record WLI for LOAD_BALANCE
#     ifdef LOAD_BALANCE
      patch->LB->WLI = WLI;
#     endif // #ifdef LOAD_BALANCE


      fprintf( File, "-------------------------------------------------------------------------------------" );
      fprintf( File, "----------------------------\n" );
      fprintf( File, "\n\n" );

      fclose( File );

   } // if ( MPI_Rank == 0 )


// f. broadcast WLI to all ranks
#  ifdef LOAD_BALANCE
#  ifdef FLOAT8
   MPI_Bcast( &patch->LB->WLI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#  else
   MPI_Bcast( &patch->LB->WLI, 1, MPI_FLOAT , 0, MPI_COMM_WORLD );
#  endif
#  endif


   delete [] NPatch_Local;
   delete [] NPatch_Gather;
   delete [] Coverage_Gather;

} // FUNCTION : Aux_PatchCount
