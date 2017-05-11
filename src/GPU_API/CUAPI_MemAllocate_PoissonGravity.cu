
#include "DAINO.h"

#if ( defined GPU  &&  defined GRAVITY )



extern real (*d_Rho_Array_P    )[ RHO_NXT*RHO_NXT*RHO_NXT ];
extern real (*d_Pot_Array_P_In )[ POT_NXT*POT_NXT*POT_NXT ];
extern real (*d_Pot_Array_P_Out)[ GRA_NXT*GRA_NXT*GRA_NXT ];
extern real (*d_Flu_Array_G    )[GRA_NIN][ PATCH_SIZE*PATCH_SIZE*PATCH_SIZE ];




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemAllocate_PoissonGravity
// Description :  Allocate device and host memory for the Poisson and Gravity solvers
//
// Parameter   :  Pot_NPatchGroup   : Number of patch groups evaluated simultaneously by GPU 
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemAllocate_PoissonGravity( const int Pot_NPatchGroup )
{
   
   const long Pot_NPatch        = 8*Pot_NPatchGroup;
   const long Rho_MemSize_P     = sizeof(real)*Pot_NPatch*RHO_NXT   *RHO_NXT   *RHO_NXT;
   const long Pot_MemSize_P_In  = sizeof(real)*Pot_NPatch*POT_NXT   *POT_NXT   *POT_NXT;
   const long Pot_MemSize_P_Out = sizeof(real)*Pot_NPatch*GRA_NXT   *GRA_NXT   *GRA_NXT;
   const long Flu_MemSize_G     = sizeof(real)*Pot_NPatch*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*GRA_NIN;


// output the total memory requirement
   long TotalSize = Rho_MemSize_P + Pot_MemSize_P_In + Pot_MemSize_P_Out + Flu_MemSize_G;

   if ( MPI_Rank == 0 )
      Aux_Message( stdout, "NOTE : total memory requirement in GPU Poisson and gravity solver = %ld MB\n", 
                   TotalSize/(1<<20) ); 


// allocate the device memory
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Rho_Array_P,     Rho_MemSize_P     )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Pot_Array_P_In,  Pot_MemSize_P_In  )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Pot_Array_P_Out, Pot_MemSize_P_Out )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Flu_Array_G,     Flu_MemSize_G     )  );


// allocate the host memory by CUDA
   for (int t=0; t<2; t++)
   {
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Rho_Array_P    [t], Rho_MemSize_P     )  );
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Pot_Array_P_In [t], Pot_MemSize_P_In  )  );
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Pot_Array_P_Out[t], Pot_MemSize_P_Out )  );
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Flu_Array_G    [t], Flu_MemSize_G     )  );
   }

} // FUNCTION : CUAPI_MemAllocate_PoissonGravity



#endif // #if ( defined GPU  &&  defined GRAVITY )
