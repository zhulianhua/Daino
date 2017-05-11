
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate
// Description :  Allocate memory for several global arrays 
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_MemAllocate ... \n" );


// a. allocate the BaseP
   const int NPatch1D[3] = { NX0[0]/PATCH_SIZE+4, NX0[1]/PATCH_SIZE+4, NX0[2]/PATCH_SIZE+4 };
   BaseP = new int [ NPatch1D[0]*NPatch1D[1]*NPatch1D[2] ];


// b. allocate memory for all GPU (or CPU) solvers (including the global memory in GPU)
#  ifdef GPU
      CUAPI_MemAllocate_Fluid( FLU_GPU_NPGROUP, GPU_NSTREAM );
#  else
      Init_MemAllocate_Fluid ( FLU_GPU_NPGROUP );
#  endif

#  ifdef GRAVITY
#     ifdef GPU
         CUAPI_MemAllocate_PoissonGravity( POT_GPU_NPGROUP );
#     else
         Init_MemAllocate_PoissonGravity ( POT_GPU_NPGROUP );
#     endif
#  endif


// c. allocate load-balance variables
#  ifdef LOAD_BALANCE
   patch->LB = new LB_t( MPI_NRank, NX0_TOT, LB_INPUT__WLI_MAX );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_MemAllocate ... done\n" );

} // FUNCTION : Init_MemAllocate
