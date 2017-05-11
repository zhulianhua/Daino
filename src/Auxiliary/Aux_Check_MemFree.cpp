
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_MemFree
// Description :  Check the total free memory and terminate the program if it is below a given threshold 
//
// Note        :  1. The minimum free memory is specified by the variable "OPT__CK_MEMFREE"
//                2. This check will be performed every "global step"
//                   --> included in the function "Aux_Check"
//                3. The total free memory is estimated as the sum of the free, buffer, and cached memories
// 
// Parameter   :  MinMemFree_Total  : Minimum total free memory (in GB)
//                comment           : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Aux_Check_MemFree( const real MinMemFree_Total, const char *comment )
{

   const int  StrSize               = 128;
   const char FileName_Mem[StrSize] = "/proc/meminfo";

   char   Useless[2][StrSize], *line=NULL;
   bool   GetMemTotal=false, GetMemFree=false, GetBuffers=false, GetCached=false;
   char   MemTotal[StrSize], MemFree[StrSize], Buffers[StrSize], Cached[StrSize];
   real   MemTotal_float, MemFree_float, Buffers_float, Cached_float;
   real   MemFree_Total;
   size_t len=0;


// 1. read the memory information
   FILE *MemFile = fopen( FileName_Mem, "r" );
   if ( MemFile == NULL )
   {
      Aux_Message( stderr, "WARNING : memory information file \"%s\" does not exist (Rank %d) !!\n", 
                   FileName_Mem, MPI_Rank );
      return;
   }

   while (  !GetMemTotal  ||  !GetMemFree  ||  !GetBuffers  ||  !GetCached  )
   {
      if ( getline( &line, &len, MemFile ) == -1 )
      {
         Aux_Message( stderr, "WARNING : some memory information is not found at Rank %d ", MPI_Rank );
         Aux_Message( stderr, "(MemTotal: %s, MemFree: %s, Buffers: %s, Cached: %s)\n", 
                      GetMemTotal? "OK":"NO", GetMemFree? "OK":"NO", GetBuffers? "OK":"NO", GetCached? "OK":"NO" );
         break;
      }

      if      ( strncmp( line, "MemTotal:", 9 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], MemTotal, Useless[1] );
         GetMemTotal = true;
      }

      else if ( strncmp( line, "MemFree:", 8 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], MemFree, Useless[1] );
         GetMemFree = true;
      }

      else if ( strncmp( line, "Buffers:", 8 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], Buffers, Useless[1] );
         GetBuffers = true;
      }

      else if ( strncmp( line, "Cached:", 7 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], Cached, Useless[1] );
         GetCached = true;
      }
   } // while (  !GetMemTotal  ||  !GetMemFree  ||  !GetBuffers  ||  !GetCached  )

   fclose( MemFile );

   if ( line != NULL )  free( line );


// 2. estimate the total free memory
   MemTotal_float = (real)atof( MemTotal )*1.e-6;  // char (kB) --> real (GB)
   MemFree_float  = (real)atof( MemFree  )*1.e-6;
   Buffers_float  = (real)atof( Buffers  )*1.e-6;
   Cached_float   = (real)atof( Cached   )*1.e-6;

   MemFree_Total = MemFree_float + Buffers_float + Cached_float;


// 3. terminate the program if the total free memory is below the threshold
   int Terminate_global, Terminate_local;

   if ( MemFree_Total < MinMemFree_Total )   
   {
      Terminate_local = true;

      Aux_Message( stderr, "\"%s\" : <%s> FAILED at Time = %13.7e, Step = %ld !!\n", comment, __FUNCTION__, Time[0], Step );
      Aux_Message( stderr, "   Total free memory (%5.2f GB) is below the threshold (%5.2f GB) at Rank %d !!\n",
                   (float)MemFree_Total, (float)MinMemFree_Total, MPI_Rank );
      Aux_Message( stderr, "   MemTotal %5.2f GB, MemFree %5.2f GB, Buffers %5.2f GB, Cached %5.2f GB\n",
                   (float)MemTotal_float, (float)MemFree_float, (float)Buffers_float, (float)Cached_float );
      Aux_Message( stderr, "\n   The program is going to be terminated automatically ...\n\n" );
   }
   else
      Terminate_local = false;

// the program will be terminated as long as ONE process has detected the low free memory
   MPI_Allreduce( &Terminate_local, &Terminate_global, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD );

   if ( Terminate_global )
   {
//    output data
      Output_DumpData( 2 );

      End_DAINO();
   }

} // FUNCTION : Aux_Check_MemFree


