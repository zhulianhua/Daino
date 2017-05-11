
#ifndef GPU

#include "DAINO.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_Fluid
// Description :  Allocate memory for the fluid solver
//
// Note        :  Work when using CPUs only
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_Fluid( const int Flu_NPatchGroup )
{

   for (int t=0; t<2; t++)
   {
      h_Flu_Array_F_In       [t] = new real [Flu_NPatchGroup][FLU_NIN ][  FLU_NXT   *FLU_NXT   *FLU_NXT   ];
      h_Flu_Array_F_Out      [t] = new real [Flu_NPatchGroup][FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE];

      if ( patch->WithFlux )
      h_Flux_Array           [t] = new real [Flu_NPatchGroup][9][NCOMP][4*PATCH_SIZE*PATCH_SIZE];

      if ( OPT__ADAPTIVE_DT )
      h_MinDtInfo_Fluid_Array[t] = new real [Flu_NPatchGroup];
   }

} // FUNCTION : Init_MemAllocate_Fluid



#endif // #ifndef GPU
