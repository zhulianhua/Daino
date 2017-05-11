
#include "DAINO.h"

#if ( MODEL == MHD )
#warning : WAIT MHD !!!
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_Conservation
// Description :  Verify the conservation laws
//                --> HYDRO : check mass, momenum, and energy
//                    MHD   : check mass, momenum, and energy
//                    ELBDM : check mass only
//
// Note        :  1. This check only works with the models HYDRO, MHD, and ELBDM
//                2. The values measured during the first time this function is invoked will be taken as the
//                   reference values to estimate errors
//                3. The following gnuplot command can be used to plot "error vs. time", assuming 
//                   "TVar = [0 ... NVar-1]"
//
//                   plot 'Record__Conservation' u 1:7 every NCOMP+1::(2+TVar) w lp ps 4
//
// Parameter   :  Output2File : true --> Output results to file instead of showing on the screen
//                comment     : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Aux_Check_Conservation( const bool Output2File, const char *comment )
{

   static bool FirstTime = true;
   const char *FileName  = "Record__Conservation";


// check
#  if ( MODEL != HYDRO  &&  MODEL != MHD  &&  MODEL != ELBDM )
   Aux_Message( stderr, "Warning : function \"%s\" is supported only in the models HYDRO, MHD, and ELBDM !!\n", 
                __FUNCTION__ );
   OPT__CK_CONSERVATION = false;
   return;
#  endif

   if ( FirstTime  &&  MPI_Rank == 0  &&  Output2File )
   {
      FILE *File_Check = fopen( FileName, "r" );

      if ( File_Check != NULL )  
      {
         Aux_Message( stderr, "WARNING : the file \"%s\" already exists !!\n", FileName );
         fclose( File_Check );
      }
   }


#  if   ( MODEL == HYDRO )
   const int NVar = NCOMP;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   const int NVar = 1;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif

   double dV, Total_local[NVar], Total_sum[NVar], Total_lv[NVar]; // dV : cell volume at each level
   int    Sg;
   FILE  *File = NULL;


// output message if Output2File is off
   if ( MPI_Rank == 0  &&  !Output2File )
   {
      if ( FirstTime )  
         Aux_Message( stdout, "\"%s\" : <%s> referencing at Time = %13.7e, Step = %7ld\n", 
                      comment, __FUNCTION__, Time[0], Step );
      else        
         Aux_Message( stdout, "\"%s\" : <%s> checking at Time = %13.7e, Step = %7ld\n", 
                      comment, __FUNCTION__, Time[0], Step );
   }


// measure the total amount of the targeted variables
   for (int v=0; v<NVar; v++)
   {
      Total_local[v] = 0.0;
      Total_sum  [v] = 0.0;
   }

   for (int lv=0; lv<NLEVEL; lv++)
   {  
      for (int v=0; v<NVar; v++)    Total_lv[v] = 0.0;

      dV = patch->dh[lv] * patch->dh[lv] * patch->dh[lv];
      Sg = patch->FluSg[lv];

      for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
      {  
         if ( patch->ptr[0][lv][PID]->son == -1 )
         {
#           if   ( MODEL == HYDRO )
            for (int v=0; v<NVar; v++)
            for (int k=0; k<PATCH_SIZE; k++)
            for (int j=0; j<PATCH_SIZE; j++)
            for (int i=0; i<PATCH_SIZE; i++)
               Total_lv[v] += (double)patch->ptr[Sg][lv][PID]->fluid[v][k][j][i];

#           elif ( MODEL == MHD )
#           warning : WAIT MHD !!!

#           elif ( MODEL == ELBDM )
            for (int k=0; k<PATCH_SIZE; k++)
            for (int j=0; j<PATCH_SIZE; j++)
            for (int i=0; i<PATCH_SIZE; i++)
               Total_lv[0] += (double)patch->ptr[Sg][lv][PID]->fluid[DENS][k][j][i];
#           endif // MODEL
         }
      } // for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)

//    sum over NLEVEL levels
      for (int v=0; v<NVar; v++)
      {
         Total_lv   [v] *= dV;
         Total_local[v] += Total_lv[v];
      }
   } // for (int lv=0; lv<NLEVEL; lv++)


// sum over all ranks
   MPI_Reduce( Total_local, Total_sum, NVar, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );


// output
   if ( MPI_Rank == 0 )
   {
      static double RefTotal[NVar];
      double AbsErr[NVar], RelErr[NVar];
   
//    record the reference values
      if ( FirstTime )
      {  
         for (int v=0; v<NVar; v++)    RefTotal[v] = Total_sum[v];

         if ( Output2File )
         {
            File = fopen( FileName, "a" );
            Aux_Message( File, "%16s%10s%13s%18s%18s%18s%18s\n", 
                         "Time", "Step", "Attribute", "Evaluate", "Reference", "Absolute Error", "Error" );
            Aux_Message( File, "-----------------------------------------------------------------------------" );
            Aux_Message( File, "----------------------------------\n" );
            fclose( File );
         }
      }

//    calculate errors
      else
      {
         for (int v=0; v<NVar; v++)
         {
            AbsErr[v] = Total_sum[v] - RefTotal[v];
            RelErr[v] = AbsErr[v] / RefTotal[v]; 
         }
   
         if ( Output2File )   
         {
            File = fopen( FileName, "a" );

            for (int v=0; v<NVar; v++)
               Aux_Message( File, "%16.7e%10ld%13d%18.7e%18.7e%18.7e%18.7e\n", 
                            Time[0], Step, v, Total_sum[v], RefTotal[v], AbsErr[v], RelErr[v] );

            Aux_Message( File, "-----------------------------------------------------------------------------" );
            Aux_Message( File, "----------------------------------\n" );
            fclose( File ); 
         }

         else
         {
            Aux_Message( stdout, "%13s%20s%20s%20s%20s\n", 
                         "Attribute", "Evaluate", "Reference", "Absolute Error", "Error" );

            for (int v=0; v<NVar; v++)
               Aux_Message( stdout, "%13d%20.7e%20.7e%20.7e%20.7e\n", 
                            v, Total_sum[v], RefTotal[v], AbsErr[v], RelErr[v] );
         }
      } // if ( FirstTime ) ... else ...
   } // if ( MPI_Rank == 0 )


   if ( FirstTime )  FirstTime = false;

} // FUNCTION : Aux_Check_Conservation
