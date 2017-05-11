
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_BinarySearch
// Description :  Use binary search to locate the position of an input number in a sorted array 
//
// Note        :  1. "Array" must be sorted in advance in ascending numerical order
//                2. If there are multiple elements matching Key, the return index can be any of them
//                3. An overloaded function for the "long" Array and Key is also created
//
// Return      :     match --> array index
//                no match --> -1
//
// Parameter   :  Array : Sorted look-up integer array (in ascending numerical order)
//                Min   : Minimum array index for searching
//                Max   : Maximum array index for searching
//                Key   : Integer number to search for
//-------------------------------------------------------------------------------------------------------
int Mis_BinarySearch( const int Array[], int Min, int Max, const int Key )
{

   int Mid = -1;

   while ( Min <= Max )
   {
      Mid = ( Min + Max ) / 2;

      if      ( Array[Mid] > Key )  Max = Mid-1;
      else if ( Array[Mid] < Key )  Min = Mid+1;
      else                          return Mid;
   }

   return -1;

} // FUNCTION : Mis_BinarySearch



//-------------------------------------------------------------------------------------------------------
// overloaded function for "long" Array and Key
//-------------------------------------------------------------------------------------------------------
int Mis_BinarySearch( const long Array[], int Min, int Max, const long Key )
{

   int Mid = -1;

   while ( Min <= Max )
   {
      Mid = ( Min + Max ) / 2;

      if      ( Array[Mid] > Key )  Max = Mid-1;
      else if ( Array[Mid] < Key )  Min = Mid+1;
      else                          return Mid;
   }

   return -1;

} // FUNCTION : Mis_BinarySearch
