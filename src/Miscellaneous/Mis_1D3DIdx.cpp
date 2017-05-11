
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Idx1D2Idx3D
// Description :  Transform a 1D index to 3D indices 
//
// Parameter   :  Size  : Size in all three spatial directions
//                Idx1D : Input 1D index 
//                Idx3D : Output 3D index
//-------------------------------------------------------------------------------------------------------
void Mis_Idx1D2Idx3D( const int Size[], const int Idx1D, int Idx3D[] )
{

   const int SizeXY = Size[0]*Size[1];

   Idx3D[0] = Idx1D % Size[0];
   Idx3D[1] = Idx1D % SizeXY / Size[0];
   Idx3D[2] = Idx1D / SizeXY;

} // FUNCTION : Mis_Idx1D2Idx3D



//-------------------------------------------------------------------------------------------------------
// overloaded function for "long" Idx1D
//-------------------------------------------------------------------------------------------------------
void Mis_Idx1D2Idx3D( const int Size[], const long Idx1D, int Idx3D[] )
{

   const long SizeXY = (long)Size[0]*Size[1];

   Idx3D[0] = Idx1D % Size[0];
   Idx3D[1] = Idx1D % SizeXY / Size[0];
   Idx3D[2] = Idx1D / SizeXY;

} // FUNCTION : Mis_Idx1D2Idx3D



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Idx3D2Idx1D
// Description :  Transform 3D indices to a 1D index
//
// Note        :  The macro "IDX321" does the same job
//                
// Parameter   :  Size  : Size in all three spatial directions
//                Idx3D : Output 3D index
//
// Return      :  1D index 
//-------------------------------------------------------------------------------------------------------
long Mis_Idx3D2Idx1D( const int Size[], const int Idx3D[] )
{

   return   ( (long)Idx3D[2]*Size[1] + Idx3D[1] )*Size[0] + Idx3D[0];

} // FUNCTION : Mis_Idx3D2Idx1D
