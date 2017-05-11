
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree
// Description :  Free memory allocated by the function "Init_MemoryAllocate"
//-------------------------------------------------------------------------------------------------------
void End_MemFree()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );


// a. deallocate load-balance variables
#  ifdef LOAD_BALANCE
   if ( patch->LB != NULL )
   {
      delete patch->LB;
      patch->LB = NULL;
   }
#  endif


// b. deallocate the AMR structure
   if ( patch != NULL )   
   {
      delete patch;
      patch = NULL;
   }
#  ifdef OOC
   if ( ooc.NRank > 1  &&  OOC_patch != NULL )   
   {
      delete OOC_patch;
      OOC_patch = NULL;
   }
#  endif


// c. deallocate the BaseP
   if ( BaseP != NULL )   
   {
      delete [] BaseP;
      BaseP = NULL;
   }


// d. deallocate arrays for GPU (or CPU) solvers
#  ifdef GPU
      CUAPI_MemFree_Fluid( GPU_NSTREAM );
#  else
      End_MemFree_Fluid();
#  endif

#  ifdef GRAVITY
#     ifdef GPU
         CUAPI_MemFree_PoissonGravity();
#     else
         End_MemFree_PoissonGravity();
#     endif
#  endif


// e. deallocate the dump table
   if ( DumpTable != NULL )   
   {
      delete [] DumpTable;
      DumpTable = NULL;
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   
} // FUNCTION : End_MemFree
