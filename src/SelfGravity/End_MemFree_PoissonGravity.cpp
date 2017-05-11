
#include "DAINO.h"

#if ( !defined GPU  &&  defined GRAVITY )




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree_PoissonGravity
// Description :  Free memory previously allocated by the function "Init_MemAllocate_PoissonGravity"
//
// Parameter   :  None 
//-------------------------------------------------------------------------------------------------------
void End_MemFree_PoissonGravity()
{

   for (int t=0; t<2; t++)
   {
      if ( h_Rho_Array_P    [t] != NULL )    delete [] h_Rho_Array_P    [t];
      if ( h_Pot_Array_P_In [t] != NULL )    delete [] h_Pot_Array_P_In [t];
      if ( h_Pot_Array_P_Out[t] != NULL )    delete [] h_Pot_Array_P_Out[t];
      if ( h_Flu_Array_G    [t] != NULL )    delete [] h_Flu_Array_G    [t];

      h_Rho_Array_P    [t] = NULL; 
      h_Pot_Array_P_In [t] = NULL;
      h_Pot_Array_P_Out[t] = NULL;
      h_Flu_Array_G    [t] = NULL;
   }

} // FUNCTION : End_MemFree_PoissonGravity



#endif // #if ( !defined GPU  &&  defined GRAVITY )
