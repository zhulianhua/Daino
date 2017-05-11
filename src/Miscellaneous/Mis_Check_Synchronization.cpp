
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Check_Synchronization
// Description :  Verify the time synchronization between two input times
//
// Parameter   :  Time1    : The first  input time
//                Time2    : The second input time
//                comment  : You can put the location where this function is invoked in this string
//                Verbose  : Output the warning message if the check fails
//-------------------------------------------------------------------------------------------------------
bool Mis_Check_Synchronization( const double Time1, const double Time2, const char *comment, const bool Verbose )
{

   const double TolErr = 1.e-12;

   double RelErr;

   if ( Time1 == 0.0 )     RelErr = fabs(  Time1 - Time2 );
   else                    RelErr = fabs( (Time1 - Time2)/Time1 );

   if ( RelErr > TolErr )
   {
      if ( Verbose )
      {
         Aux_Message( stderr, "WARNING : \"%s\" : <%s> FAILED at rank %d !!\n", 
                      comment, __FUNCTION__, DAINO_RANK );
         Aux_Message( stderr, "          Time1 = %20.14e vs. Time2 = %20.14e --> RelErr = %20.14e !!\n",
                      Time1, Time2, RelErr );
      }

      return false;
   }

   else
      return true;

}
