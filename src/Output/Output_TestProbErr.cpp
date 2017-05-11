
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_TestProbErr
// Description :  Compare and output the numerical and analytical solutions for the chosen test problem 
//
// Note        :  Most test problems supported in DAINO should contain this file in the corresponding 
//                test problem directories 
//                (e.g., DAINO/test_problem/Model_ELBDM/Jeans_Instability/Output_TestProbErr.cpp). 
//                Please copy the file to here directly.
//
// Parameter   :  BaseOnly :  Only output the base-level data
//-------------------------------------------------------------------------------------------------------
void Output_TestProbErr( const bool BaseOnly )
{  

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );

// ===========================================================================================
// please replace this file by the one with the same file name in the test problem directories
// ===========================================================================================

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_TestProbErr
