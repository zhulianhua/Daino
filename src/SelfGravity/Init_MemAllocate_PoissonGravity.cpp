
#include "DAINO.h"

#if ( !defined GPU  &&  defined GRAVITY )




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_PoissonGravity
// Description :  Allocate memory for the Poisson and Gravity solvers
//
// Note        :  Only work when using CPUs only
//
// Parameter   :  Pot_NPatchGroup   : Number of patch groups to be evaluated at a time 
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_PoissonGravity( const int Pot_NPatchGroup )
{
   
   const int Pot_NPatch = 8*Pot_NPatchGroup;

   for (int t=0; t<2; t++)
   {
      h_Rho_Array_P    [t] = new real [Pot_NPatch][RHO_NXT][RHO_NXT][RHO_NXT]; 
      h_Pot_Array_P_In [t] = new real [Pot_NPatch][POT_NXT][POT_NXT][POT_NXT];
      h_Pot_Array_P_Out[t] = new real [Pot_NPatch][GRA_NXT][GRA_NXT][GRA_NXT];
      h_Flu_Array_G    [t] = new real [Pot_NPatch][GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
   }

} // FUNCTION : Init_MemAllocate_PoissonGravity



#endif // #if ( !defined GPU  &&  defined GRAVITY )
