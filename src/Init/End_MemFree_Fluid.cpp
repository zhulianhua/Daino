
#ifndef GPU

#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree_Fluid
// Description :  Free memory previously allocated by the function "Init_MemAllocate_Fluid"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_MemFree_Fluid()
{

   for (int t=0; t<2; t++)
   {
      if ( h_Flu_Array_F_In       [t] != NULL )    delete [] h_Flu_Array_F_In       [t];
      if ( h_Flu_Array_F_Out      [t] != NULL )    delete [] h_Flu_Array_F_Out      [t];
      if ( h_Flux_Array           [t] != NULL )    delete [] h_Flux_Array           [t];
      if ( h_MinDtInfo_Fluid_Array[t] != NULL )    delete [] h_MinDtInfo_Fluid_Array[t];

      h_Flu_Array_F_In       [t] = NULL; 
      h_Flu_Array_F_Out      [t] = NULL;
      h_Flux_Array           [t] = NULL;
      h_MinDtInfo_Fluid_Array[t] = NULL;
   }

} // FUNCTION : End_MemFree_Fluid



#endif // #ifndef GPU
